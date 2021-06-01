import math
import os
import sys
from collections import defaultdict
from operator import itemgetter
from sys import stderr
# try:
#     import cPickle as pickle
# except:
#     import pickle

from igraph import Graph
from visualize import visualize_separate, visualize_grid

# from .visualize import visualize_separate, visualize_grid


def read(path):
    with open(path, 'r') as f:
        lines = f.read().splitlines()
    return lines


class Compressor:
    """ Parameters and state of the GraphZip model """

    def __init__(self, batch_size=10, dict_size=math.inf, directed=False):
        """ Initialize state and set parameters """
        # Should remain constant
        _directed = directed

        # For keeping track of internal stats
        _compress_count = 0
        _lines_read = 0
        _dict_trimmed = 0

        # Max dictionary size. If we go over size, keep only the top patterns
        dict_size = dict_size

        # The window size for which we process the graph stream
        batch_size = batch_size

        # True if we don't just want subisomorphisms, but strict subgraphs
        match_strict = True

        # True if we want to allow edges to be added without the vertices
        # being declared in the same batch
        # However, in order for the edge to be added, we still need the
        # label info of the vertex we're adding (in the vid_to_label dict)
        add_implicit_vertices = True

        # Mapping from vid's to label
        # We need this because a batch (from a .graph file) may just contain
        # e label vid1 vid2, and if we want to add a single edge vid1 vid2 we
        # don't want to loose label information
        vid_to_label = dict()

        # If true, then we only keep vertex->label mapping information
        # (from 'v x x' lines) for the duration of a file
        # This means that any vertex mentioned in an edge must have been
        # previously declared in that same .graph file
        # If false, then we keep the vertex->label mapping forever
        # i.e. if the vertex is declared in ANY previous file, it can be
        # mentioned in an edge at any later time
        self.label_history_per_file = False  # XXX

        # Initialize patterns list P of (graph,count,score) tuples
        self.P = []


    def safe_add_edge(self, graph, source, target, **kwds):
        """ Makes sure vertices are added to a graph before adding an edge

        iGraph wrapper function around graph.add_edge()
        Note: we still need the id->label mapping of the vertex to proceed

        """
        # print(f"try to add from {source} to {target}")
        try:
            graph.vs.find(name=source)
        except ValueError:
            graph.add_vertex(source)

        try:
            graph.vs.find(name=target)
        except ValueError:
            graph.add_vertex(target)

        # don't add duplicate edges
        if not graph.are_connected(source, target):
            graph.add_edge(source, target, **kwds)

    def get_score(self, graph, count):
        """ Calculate a pattern's compression score
        Compression scores are used to rank the pattern dictionary
        """
        size = len(graph.es)
        return (size-1) * (count-1)

    def trim_dictionary(self, threshold_multiplier=2):
        """ Trim pattern dictionary to theta if size > multiplier*theta
        By default, the multiplier is set to 2, so the dictionary is trimmed to
        size theta when it exceeds size 2*theta
        """
        if len(self.P) > threshold_multiplier * self.dict_size:
            self._dict_trimmed += 1
            self.P = sorted(self.P, key=itemgetter(2), reverse=True)
            # Should be faster than `self.P = self.P[:self.dict_size]'
            del self.P[self.dict_size:]

    def update_dictionary(self, pattern):
        """ Update the pattern dictionary with a new graph
        Iterate through the pattern dictionary:
            If the pattern already exists, update its counter
            Otherwise, add the new pattern to the dictionary
        """
        # c2 = pattern.vs['label']
        c2_edge = pattern.es['label']

        for i, (graph, count, score) in enumerate(self.P):
            # c1 = graph.vs['label']
            c1_edge = graph.es['label']

            if len(c1_edge) != len(c2_edge):
                continue

            if graph.isomorphic_vf2(pattern):
                # match found, update counter and score
                new_count = count+1
                new_score = self.get_score(graph, new_count)
                self.P[i] = (graph, new_count, new_score)
                return

        self.trim_dictionary()

        # If pattern is not in dictionary, add it
        count = 1
        score = self.get_score(pattern, count)
        self.P.append((pattern, count, score))

    def iterate_batch(self, G_batch):
        """ `Compress` a single graph stream object G_batch """

        if len(G_batch.es) == 0:
            return
        taken = defaultdict(lambda: False)  # XXX move to under vmap in maps?

        # Keep track of all the extended patterns
        # then update dictionary at the end
        new_patterns = []

        # For each pattern-graph p in P
        for p, c, s in self.P:

            print(p)
            print('\n------------------------\n')
            # Get all subgraphs matching pattern p in the batch's graph
            if self.match_strict:
                if len(p.es) >= len(G_batch.es):
                    maps = []
                else:
                    maps = G_batch.get_subisomorphisms_vf2(p,
                                                           edge_color1=G_batch.es['label'],
                                                           edge_color2=p.es['label'])
            else:
                print("Getting loose embeddings (no label match)", file=stderr)
                maps = G_batch.get_subisomorphisms_vf2(p)

            maps
            # For each instance of i of p in B
            # 'vmap' is a mapping of p's vertex indices to G's v. indices
            counter = 0
            # 그래프에서 찾은 subgraph set에 속한 각 subgraph에 대해 반복 수행
            for v_map in maps:
                counter += 1

                # Extend the pattern by "one degree/layer" as p_new
                p_new = None

                # respective vertex indices for the mappings
                # e.g.
                # p_v: 0, 1, 2, 3, 4, 5, ...
                # G_v: 4, 3, 2, 5, 1, 0, ...
                Gv_to_pv = {G_v: p_v for p_v, G_v in enumerate(v_map)}
                #각 subgraph를 반복 수행하면서 버텍스를 인덱스로 매칭시켜준다

                # Iterate through the mapped vertices in a single embedding
                # Equivalent to `for G_v,p_v in zip(v_map,range(len(v_map)))'
                #각 매칭된 버텍스들을 반복하면서 체크해준다
                for p_v, G_v in enumerate(v_map):

                    # Check whether the equivalent node on the big graph has
                    # extra incident edges not on the dictionary pattern
                    # XXX Incident takes keyword 'mode' when directed,
                    #     defaults to only OUT edges (not ALL)
                    G_edges = G_batch.es[G_batch.incident(G_v)] # 그래프의 엣지를 인덱스로 비교
                    p_edges = p.es[p.incident(p_v)] # 사전의 엣지를 인덱스로 비교

                    # (p MUST be subgraph of G) P는 반드시 그래프 G의 서브그래프
                    if len(G_edges) <= len(p_edges):
                        continue

                    # Find which extra edges we need to add
                    for Ge in G_edges: #G_batch 그래프의 v의 인접한 vertex들에 대하여

                        # If we don't want to re-use edges in pattern building
                        #  if taken[Ge]:
                        #     continue

                        # One of these should be the same as G_v
                        Gv_source_index = Ge.source
                        Gv_target_index = Ge.target
                        Gv_source = G_batch.vs[Gv_source_index]
                        Gv_target = G_batch.vs[Gv_target_index]

                        # First possibility is the extending edge leads to a
                        # vertex not on the pattern, in which case it does not
                        # exist in the mapping
                        # Gv_to_pv : name과 인덱스번호를 매핑시켜줌
                        if Gv_source_index not in Gv_to_pv:  # not in map
                            # In which case we want to add the vertex, then add
                            # the edge extending to that vertex
                            # p_new:새로운 패턴을 확장해서 넣어주는부분
                            if p_new is None: # source vertex에대한 첫루프는 여기서 시작
                                    p_new = p.copy()
                            # p_new.add_vertex(label=Gv_source['label'])
                            pv_target_index = Gv_to_pv[Gv_target_index]
                            pv_source_index = p_new.vcount()-1
                            p_new.add_edge(pv_source_index,
                                           pv_target_index,
                                           label=Ge['label'])
                        elif Gv_target_index not in Gv_to_pv:  # not in map
                            if p_new is None: # target vertex에 대한 첫루프는 여기서 시작
                                p_new = p.copy()
                            # p_new.add_vertex(label=Gv_target['label'])
                            pv_source_index = Gv_to_pv[Gv_source_index]
                            pv_target_index = p_new.vcount()-1
                            p_new.add_edge(pv_source_index,
                                           pv_target_index,
                                           label=Ge['label'])

                        # Second possibility is that both vertices exist, but
                        # the edge only exists in the larger (batch) graph
                        # e.g., closing a cycle
                        else:  # source index나 target index가 있는경우(=확장가능한경우)
                            pv_source_index = Gv_to_pv[Gv_source_index]
                            pv_target_index = Gv_to_pv[Gv_target_index]
                            # 패턴 p(subgraph)가 사전에 존재하지 않으면,
                            if not p.are_connected(pv_source_index,
                                                   pv_target_index):
                                if p_new is None:
                                    p_new = p.copy()
                                self.safe_add_edge(p_new,
                                                   pv_source_index,
                                                   pv_target_index,
                                                   label=Ge['label'])

                        # Third possibility is that the edge exists in both the
                        # pattern and the larger (batch) graph (do nothing)

                        # Mark the subgraph and the extra edge as taken
                        taken[Ge] = True

                # Add the new pattern to the dictionary
                if p_new is not None:
                    new_patterns.append(p_new)

        for g in new_patterns:
            self.update_dictionary(g)

        # Add remaining edges in B as single-edge patterns in P
        for e in G_batch.es:
            if e not in taken:
                source = G_batch.vs[e.source]
                target = G_batch.vs[e.target]
                single_edge = Graph(directed=self._directed)
                # Don't need the safe methods here since it's a fresh graph
                single_edge.add_vertex(source)
                single_edge.add_vertex(target)
                single_edge.add_edge(0, 1, label=e['label'])
                # self.update_dictioonary(single_edge)
                # single_edge.add_edge(0, 1)
                self.update_dictionary(single_edge)

    def iterate_batch_hb(self, G_batch):
        """ `Compress` a single graph stream object G_batch """

        if len(G_batch.es) == 0:
            return
        taken = defaultdict(lambda: False)  # XXX move to under vmap in maps?
        new_patterns = []

        for p, c, s in self.P:
            print(p)
            # Get all subgraphs matching pattern p in the batch's graph
            if self.match_strict:
                if len(p.es) >= len(G_batch.es):
                    maps = []
                else:
                    maps = G_batch.get_subisomorphisms_vf2(p,
                                            color1=G_batch.vs['label'],
                                            color2=p.vs['label'],
                                            edge_color1=G_batch.es['label'],
                                            edge_color2=p.es['label'])
            else:
                print("Getting loose embeddings (no label match)", file=stderr)
                maps = G_batch.get_subisomorphisms_vf2(p)

            # For each instance of i of p in B
            # 'vmap' is a mapping of p's vertex indices to G's v. indices
            counter = 0
            for v_map in maps:
                counter += 1

                # Extend the pattern by "one degree/layer" as p_new
                p_new = None

                # respective vertex indices for the mappings
                # e.g.
                # p_v: 0, 1, 2, 3, 4, 5, ...
                # G_v: 4, 3, 2, 5, 1, 0, ...
                Gv_to_pv = {G_v: p_v for p_v, G_v in enumerate(v_map)}

                # Iterate through the mapped vertices in a single embedding
                # Equivalent to `for G_v,p_v in zip(v_map,range(len(v_map)))'
                for p_v, G_v in enumerate(v_map):

                    # Check whether the equivalent node on the big graph has
                    # extra incident edges not on the dictionary pattern
                    # XXX Incident takes keyword 'mode' when directed,
                    #     defaults to only OUT edges (not ALL)
                    G_edges = G_batch.es[G_batch.incident(G_v)]
                    p_edges = p.es[p.incident(p_v)]

                    # (p MUST be subgraph of G)
                    if len(G_edges) <= len(p_edges):
                        continue

                    # Find which extra edges we need to add
                    for Ge in G_edges:

                        # If we don't want to re-use edges in pattern building
                        #  if taken[Ge]:
                        #     continue

                        # One of these should be the same as G_v
                        Gv_source_index = Ge.source
                        Gv_target_index = Ge.target
                        Gv_source = G_batch.vs[Gv_source_index]
                        Gv_target = G_batch.vs[Gv_target_index]

                        # First possibility is the extending edge leads to a
                        # vertex not on the pattern, in which case it does not
                        # exist in the mapping
                        if Gv_source_index not in Gv_to_pv:  # not in map
                            # In which case we want to add the vertex, then add
                            # the edge extending to that vertex
                            if p_new is None:
                                    p_new = p.copy()
                            p_new.add_vertex(label=Gv_source['label'])
                            pv_target_index = Gv_to_pv[Gv_target_index]
                            pv_source_index = p_new.vcount()-1
                            p_new.add_edge(pv_source_index,
                                           pv_target_index,
                                           label=Ge['label'])
                        elif Gv_target_index not in Gv_to_pv:  # not in map
                            if p_new is None:
                                p_new = p.copy()
                            p_new.add_vertex(label=Gv_target['label'])
                            pv_source_index = Gv_to_pv[Gv_source_index]
                            pv_target_index = p_new.vcount()-1
                            p_new.add_edge(pv_source_index,
                                           pv_target_index,
                                           label=Ge['label'])

                        # Second possibility is that both vertices exist, but
                        # the edge only exists in the larger (batch) graph
                        # e.g., closing a cycle
                        else:
                            pv_source_index = Gv_to_pv[Gv_source_index]
                            pv_target_index = Gv_to_pv[Gv_target_index]

                            if not p.are_connected(pv_source_index,
                                                   pv_target_index):
                                if p_new is None:
                                    p_new = p.copy()
                                self.safe_add_edge(p_new,
                                                   pv_source_index,
                                                   pv_target_index,
                                                   label=Ge['label'])

                        # Third possibility is that the edge exists in both the
                        # pattern and the larger (batch) graph (do nothing)

                        # Mark the subgraph and the extra edge as taken
                        taken[Ge] = True

                # Add the new pattern to the dictionary
                if p_new is not None:
                    new_patterns.append(p_new)

        for g in new_patterns:
            self.update_dictionary(g)

        # Add remaining edges in B as single-edge patterns in P
        for e in G_batch.es:
            if e not in taken:
                source = G_batch.vs[e.source]
                target = G_batch.vs[e.target]
                single_edge = Graph(directed=self._directed)
                # Don't need the safe methods here since it's a fresh graph
                single_edge.add_vertex(label=source['label'])
                single_edge.add_vertex(label=target['label'])
                single_edge.add_edge(0, 1, label=e['label'])
                self.update_dictionary(single_edge)

    def parse_line_hb(self, line_str, G_batch):
        if line_str == '':
            return
        raw = line_str.strip().split()
        # 지금은.. time 까지만 고려
        # 추가로 다른 레이블 고려
        e_source_id, e_dest_id, e_label = raw[0], raw[1], int(raw[2])
        try:
            if not G_batch.are_connected(e_source_id, e_dest_id):
                # We can only add an edge if the vertices already exist
                if self.add_implicit_vertices:
                    self.safe_add_edge(G_batch,
                                       e_source_id,
                                       e_dest_id,
                                       label=e_label)
                else:
                    G_batch.add_edge(e_source_id,
                                     e_dest_id,
                                     label=e_label)
        except ValueError:
            # We can only add an edge if the vertices already exist
            if self.add_implicit_vertices:
                self.safe_add_edge(G_batch,
                                   e_source_id, e_dest_id, label=e_label)
            else:
                print("Error: vertex in line DNE:\n%s" % line_str,
                      file=stderr)
                raise

    def visualize_dictionary(self, fout, top=True, n=None):
        """ Visualize the dictionary as a grid of graphs (in an SVG file) """
        visualize_grid(fout, self.P, top, n)

    def visualize_dictionary_separate(self, fout, n=None):
        """ Save the top-N dictionary patterns as separate SVG files """
        visualize_separate(fout, self.P, n)

    def compress_file_hb(self, fin, fout=None):
        self._compress_count += 1
        with open(fin, 'r') as f:
            G_batch = Graph(directed=self._directed)
            line_count = 0
            edge_count = 0

            for line in f:
                line_count += 1
                edge_count += 1

                if (line_count % 1000 == 0):
                    print("Read %d lines (%d edges) from %s" %(line_count, edge_count, fin), file=stderr, end='\r')

                self.parse_line_hb(line, G_batch) # g_batch에 그래프 정보 입력
                # Only process our "batch" once we've reached a certain size
                if (edge_count == 0 or (edge_count % self.batch_size) != 0):
                    continue

                # Processed the batch, then create a fresh stream object/graph
                self.iterate_batch(G_batch)
                G_batch = Graph(directed=self._directed)

            # Process the leftovers (if any)
            if len(G_batch.es) > 0:
                self.iterate_batch(G_batch)

        print("Read %d lines (%d edges) from %s" %  # final count
              (line_count, edge_count, fin), file=stderr, end='\r')

        if self.label_history_per_file:
            # Wipe vid->label mapping
            self.vid_to_label = dict()


file1 = 'data/50.txt'
file2 = 'data/1k.txt'

comp4 = Compressor(batch_size=20, dict_size=5)
comp4.compress_file_hb(file1)

comp4.visualize_dictionary("t1.out")
# comp4.visualize_dictionary_separate("t2.out")

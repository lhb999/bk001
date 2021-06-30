""" Python implementation of the GraphZip algorithm """
import pyfpgrowth
import logging
import traceback
logging.getLogger().setLevel(logging.DEBUG)
import math
from collections import defaultdict
from operator import itemgetter
from sys import stderr
try:
    import cPickle as pickle
except:
    import pickle
from igraph import Graph
import pyfpgrowth
from itertools import chain
# from .visualize import visualize_separate, visualize_grid

class Compressor2:
    """ Parameters and state of the GraphZip model """

    def __init__(self, batch_size=10, dict_size=math.inf, directed=False):
        """ Initialize state and set parameters """
        # Should remain constant
        self._directed = directed

        # For keeping track of internal stats
        self._compress_count = 0
        self._lines_read = 0
        self._dict_trimmed = 0
        self._index_to_path = {}
        self._expanded_pattern_list = []
        self._time_stamp = 0
        self._pattern_from = {}
        self._pattern_id_inc = 0

        # Max dictionary size. If we go over size, keep only the top patterns
        self.dict_size = dict_size

        # The window size for which we process the graph stream
        self.batch_size = batch_size

        # True if we don't just want subisomorphisms, but strict subgraphs
        self.match_strict = True

        # True if we want to allow edges to be added without the vertices
        # being declared in the same batch
        # However, in order for the edge to be added, we still need the
        # label info of the vertex we're adding (in the vid_to_label dict)
        self.add_implicit_vertices = True

        # Mapping from vid's to label
        # We need this because a batch (from a .graph file) may just contain
        # e label vid1 vid2, and if we want to add a single edge vid1 vid2 we
        # don't want to loose label information
        self.vid_to_label = dict()

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

    def append_expanded_pattern(self, j, k):
        if(j != k):
            tmp_p = [self._time_stamp, j, k]
            self._expanded_pattern_list.append(tmp_p)

    def pattern_score_regularization(self):
        tmp_p = self.P
        try:
            max_score = max(list(filter(lambda x: x[3] > 0, tmp_p)))[0]
        except Exception as e:
            max_score = 1
        # print("MAX -> ", max_score)
        for i, (pid, p, c, s) in enumerate(tmp_p):
            new_score = round(s / max_score) * 100
            tmp_p[i] = (pid, p, c, new_score)

    def pattern_score_update_with_FPtree(self, filtered_pattern):
        for i, (pid, p, c, old_score) in enumerate(filtered_pattern):
            new_score = old_score + 100
            filtered_pattern[i] = (pid, p, c, new_score)

    def pattern_score_update_with_timeStamp(self):
        for (time_stamp, pid1, pid2) in self._expanded_pattern_list:
            g1 = list(filter(lambda x: (x[1][0] == pid1) or (x[1][0] == pid2), enumerate(self.P)))
            if g1 is not None:
                for idx, (pid, g, c, s) in g1:
                    new_score = s + time_stamp
                    self.P[idx] = (pid, g, c, new_score)

    def pattern_pruner(self, fp_frequent_factor=2, fp_association_factor=0.5):
        # SCORE 정규화
        self.pattern_score_regularization()

        # FP 트리 기반 SCORING
        check1 = self.P.__len__()
        transactions = self._expanded_pattern_list[1:]
        patterns = pyfpgrowth.find_frequent_patterns(transactions, fp_frequent_factor)
        rules = pyfpgrowth.generate_association_rules(patterns, fp_association_factor)
        unnested = list(chain.from_iterable(rules))
        parsed_list = list(map(lambda x: x[:1][0], self.P))
        inter_section = list(set(unnested) & set(parsed_list))
        filtered = list(filter(lambda x: inter_section.__contains__(x[:1][0]), self.P))
        check2 = filtered.__len__()
        self.pattern_score_update_with_FPtree(filtered)
        # 여기까지 FP ~
        # Time Stamp 기반 Scoring
        self.pattern_score_update_with_timeStamp()
        sorted_p = sorted(self.P, key=lambda p: p[3], reverse=True)[0:self.dict_size]
        # print(f"Pattern Pruned => {check1} ->  {sorted_p.__len__()}")
        # sorted(students, key=lambda student: student[2])
        return sorted_p

    def save_state(self, fout):
        """ Save the compressor state as a pickle file
        Object format is:
          (_compress_count,_lines_read,_dict_trimmed,P)
        """

        with open(fout, 'wb') as pickle_file:
            pickle.dump((self._compress_count,
                         self._lines_read,
                         self._dict_trimmed,
                         self.P),
                        pickle_file)

    def import_state(self, fin):
        """ Import a saved compressor state generated using save_state()
        Expects object format:
          (_compress_count,_lines_read,_dict_trimmed,P)
        """
        with open(fin, 'rb') as pickle_file:
            (self._compress_count,
             self._lines_read,
             self._dict_trimmed,
             self.P) = pickle.load(pickle_file)

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

    def update_dictionary_idx(self, pattern, pi):
        """ Update the pattern dictionary with a new graph
        Iterate through the pattern dictionary:
            If the pattern already exists, update its counter
            Otherwise, add the new pattern to the dictionary
        """
        # j = 이전
        c2 = pattern.vs['label']
        c2_edge = pattern.es['label']

        for i, (pid, graph, count, score) in enumerate(self.P):
            # 이미 존재하는 패턴일 때
            c1 = graph.vs['label']
            c1_edge = graph.es['label']

            if len(c1) != len(c2) or len(c1_edge) != len(c2_edge):
                continue

            if graph.isomorphic_vf2(pattern,
                                    color1=c1, color2=c2,
                                    edge_color1=c1_edge, edge_color2=c2_edge):
                # match found, update counter and score
                new_count = count+1
                # self.append_expanded_pattern(pid, pi)
                new_score = self.get_score(graph, new_count)
                pattern_from = self._pattern_from[pid]
                tmp_path = [self._time_stamp,pattern_from, pid]
                self._expanded_pattern_list.append(tmp_path)
                # print(f"{pattern_from} -> {pid}")
                self.P[i] = (pid, graph, new_count, new_score)
                return

        # self.trim_dictionary()

        # If pattern is not in dictionary, add it
        count = 1
        score = self.get_score(pattern, count)

        current_pattern_id = self._pattern_id_inc

        self._pattern_from[current_pattern_id] = pi
        # print(f"{pi} -> {current_pattern_id}")
        self.append_expanded_pattern(pi, current_pattern_id)
        self.P.append((current_pattern_id, pattern, count, score))
        self._pattern_id_inc += 1
        # self.P.append((pattern, count, score))


    def update_dictionary(self, pattern):
        """ Update the pattern dictionary with a new graph
        Iterate through the pattern dictionary:
            If the pattern already exists, update its counter
            Otherwise, add the new pattern to the dictionary
        """
        # j = 이전
        c2 = pattern.vs['label']
        c2_edge = pattern.es['label']

        for i, (pid, graph, count, score) in enumerate(self.P):
            # 이미 존재하는 패턴일 때
            c1 = graph.vs['label']
            c1_edge = graph.es['label']

            if len(c1) != len(c2) or len(c1_edge) != len(c2_edge):
                continue

            if graph.isomorphic_vf2(pattern,
                                    color1=c1, color2=c2,
                                    edge_color1=c1_edge, edge_color2=c2_edge):
                # match found, update counter and score
                new_count = count+1
                # self.append_expanded_pattern(pid, pi)
                new_score = self.get_score(graph, new_count)
                self.P[i] = (pid, graph, new_count, new_score)
                return

        # self.trim_dictionary()

        # If pattern is not in dictionary, add it
        count = 1
        score = self.get_score(pattern, count)

        current_pattern_id = self._pattern_id_inc
        self.P.append((current_pattern_id, pattern, count, score))
        self._pattern_id_inc += 1
        # self.P.append((pattern, count, score))


    def iterate_batch(self, G_batch):
        """ `Compress` a single graph stream object G_batch """

        if len(G_batch.es) == 0:
            return
        taken = defaultdict(lambda: False)  # XXX move to under vmap in maps?

        # Keep track of all the extended patterns
        # then update dictionary at the end
        new_patterns = []

        # For each pattern-graph p in P
        for i, (pid, p, c, s) in enumerate(self.P):

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
                    new_patterns.append((p_new, pid))
                    # new_patterns.append(p_new)

        for (g, ii) in new_patterns:
            # self.update_dictionary(g)
            self.update_dictionary_idx(g, ii)

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

    def compress_file(self, fin, fout=None):
        """ Run GraphZip on a graph specified by an input (.graph) file
        For larger files we process input stream in chunks specified by
        parameter alpha ('batch_size')
        """
        self._compress_count += 1

        with open(fin, 'r') as f:
            G_batch = Graph(directed=self._directed)
            line_count = 0
            edge_count = 0

            for line in f:
                line_count += 1
                # if (line_count % 1000 == 0):
                #     # print("Read %d lines (%d edges) from %s" %
                #     #       (line_count, edge_count, fin), file=stderr, end='\n')
                if line[0] == 'e':
                    edge_count += 1

                # Add the vertex/edge to our graph stream object (G_batch)
                self.parse_line(line, G_batch)

                # Only process our "batch" once we've reached a certain size
                if (edge_count == 0 or (edge_count % self.batch_size) != 0):
                    continue

                # print(f"Continue Process => line: {line_count}, edge count: {edge_count} time stamp: {self._time_stamp}")
                if(self.P.__len__() >= self.dict_size):
                    self.P = self.pattern_pruner()
                # Processed the batch, then create a fresh stream object/graph
                self.iterate_batch(G_batch)
                self._time_stamp += 1
                G_batch = Graph(directed=self._directed)

            # Process the leftovers (if any)
            if len(G_batch.es) > 0:
                self.iterate_batch(G_batch)

        # print("Read %d lines (%d edges) from %s" %  # final count
        #       (line_count, edge_count, fin), file=stderr, end='\r')

        if self.label_history_per_file:
            # Wipe vid->label mapping
            self.vid_to_label = dict()

    def safe_add_edge(self, graph, source, target, **kwds):
        """ Makes sure vertices are added to a graph before adding an edge
        iGraph wrapper function around graph.add_edge()
        Note: we still need the id->label mapping of the vertex to proceed
        """
        try:
            graph.vs.find(name=source)
        except ValueError:
            graph.add_vertex(source, label=self.vid_to_label[source])

        try:
            graph.vs.find(name=target)
        except ValueError:
            graph.add_vertex(target, label=self.vid_to_label[target])

        # don't add duplicate edges
        if not graph.are_connected(source, target):
            graph.add_edge(source, target, **kwds)

    def parse_line(self, line_str, G_batch):
        """ Parse a line from .graph input and update the graph accordingly
        Line format is either:
            a) 'v id label'
            b) 'e v1 v2 label'
            c) '%  comment to ignore'
        Where label must be an int
        Generally, we assume that the .graph file is formatted correctly
        e.g., we don't cover the case where a vertex is defined more than once
        We index the nodes on IDs, since labels aren't guaranteed to be unique
        """
        if line_str == '':
            return

        raw = line_str.strip().split()

        # lines that begin with % are comments
        if raw[0] == '%':
            return

        # Vertex case
        if (raw[0] == 'v'):
            v_id, v_label = raw[1], int(raw[2])
            self.vid_to_label[v_id] = v_label

        # Edge case
        elif (raw[0] == 'e' or raw[0] == 'u' or raw[0] == 'd'):
            e_type = raw[0]  # d=directed, u=undirected, e=directed unless flag
            e_source_id, e_dest_id, e_label = raw[1], raw[2], int(raw[3])

            # Don't add an edge if it already exists
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

            # Means that one of the vertices DNE (also not connected)
            except ValueError:
                # We can only add an edge if the vertices already exist
                if self.add_implicit_vertices:
                    self.safe_add_edge(G_batch,
                                       e_source_id, e_dest_id, label=e_label)
                else:
                    print("Error: vertex in line DNE:\n%s" % line_str,
                          file=stderr)
                    raise

        else:
            # Generic parsing problem
            raise ValueError("Error: could not parse line:\n%s" % raw)

    # XXX Convenience methods
    # XXX Should eventually remove and just use vis. module directly (SRP)
    # def visualize_dictionary(self, fout, top=True, n=None):
    #     """ Visualize the dictionary as a grid of graphs (in an SVG file) """
    #     visualize_grid(fout, self.P, top, n)
    #
    # def visualize_dictionary_separate(self, fout, n=None):
    #     """ Save the top-N dictionary patterns as separate SVG files """
    #     visualize_separate(fout, self.P, n)

file7 = 'data/SUBGEN/4PATH/4PATH_1_5_50cx.graph'
file8 = 'data/SUBGEN/4PATH/100K.graph'
file9 = 'data/SUBGEN/4PATH/10K.graph'
file0 = 'data/SUBGEN/4PATH/10.graph'




def write_dictionary(model, fout=None):
    """ Print the pattern dictionary in readable .graph format
    If fout == None, print to stdout
    Otherwise, write to fout/file
    """
    model.P = sorted(model.P, key=itemgetter(2), reverse=True)
    for i, (pid, g, c, s) in enumerate(model.P):
        if fout is None:
            print("%% Pattern %d" % (i + 1))
            print("%% Score:  %d" % s)
            print("%% Count:  %d" % c)
            for i, v in enumerate(g.vs):
                print("v %d %d" % (i, v['label']))
            for e in g.es:
                print("e %d %d %d" % (e.source, e.target, e['label']))
        else:
            fout.write("%% Pattern %d\n" % i)
            fout.write("%% Score:  %d\n" % s)
            fout.write("%% Count:  %d\n" % c)
            for i, v in enumerate(g.vs):
                fout.write("v %d %d\n" % (i, v['label']))
            for e in g.es:
                fout.write("e %d %d %d\n" %
                           (e.source, e.target, e['label']))

import time
import os

sizes = [5]

# print
# n, "Bytes"  # 바이트 단
for size in sizes:
    start = time.time()  # 시작 시간 저장
    for test_case in range(5):
        fname = f"test_a_{size}_{test_case}.out"
        comp4 = Compressor2(batch_size=200, dict_size=size)
        comp4.compress_file(file9)
        with open(fname, 'w') as fout:
            write_dictionary(comp4, fout)
    n = os.path.getsize(fname)
    ttime = round(time.time() - start, 2)
    # print(f"pattern size: {size}")
    print(ttime, n)



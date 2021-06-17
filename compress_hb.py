import math
import os
import sys
from collections import defaultdict
from operator import itemgetter
from sys import stderr
from pandas.util import hash_pandas_object
import pyfpgrowth
import logging
import traceback
logging.getLogger().setLevel(logging.DEBUG)

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
        self._directed = directed

        # For keeping track of internal stats
        self._compress_count = 0
        self._lines_read = 0
        self._dict_trimmed = 0
        self._batch_count = 0
        self._max_pattern_size = 100
        self._update_call_counter = 0
        self._index_to_path = {}
        self._expanded_pattern_list = []
        self._time_stamp = 0
        self._pattern_from = {}

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
        self.label_to_id = dict()

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

    def get_label_idx(self, label):
        if label not in self.label_to_id:
            self.label_to_id[label] = len(self.label_to_id)
        # print(f"label: {label} ==> {self.label_to_id[label]}") # FOR DEBUG
        return self.label_to_id[label]

    def get_or_create_path(self, sp):
        pf = self._pattern_from

        k = sp
        if(pf.__contains__(k)):
            return pf[k]
            ...
        else:
            # print(f"search with : {k}")
            return k
            ...


    def safe_add_edge(self, graph, source, target, **kwds):
        """
        Makes sure vertices are added to a graph before adding an edge
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
            print("TRIM!")
            self._dict_trimmed += 1
            self.P = sorted(self.P, key=itemgetter(2), reverse=True)
            # Should be faster than `self.P = self.P[:self.dict_size]'
            del self.P[self.dict_size:]

    def pattern_pruner(self):
        transactions = self._expanded_pattern_list[1:]
        patterns = pyfpgrowth.find_frequent_patterns(transactions, 3)
        rules = pyfpgrowth.generate_association_rules(patterns, 0.5)

        print("TEST")

    def update_dictionary(self, pattern):
        self._update_call_counter += 1
        # print(f"update function call counter ==> {self._update_call_counter}")
        """ Update the pattern dictionary with a new graph
        Iterate through the pattern dictionary:
            If the pattern already exists, update its counter
            Otherwise, add the new pattern to the dictionary
        """
        # pattern -> 그래프
        # 해당 함수에서는 패턴을 넣는것이 아니라 패턴을 추출한 뒤에
        # 서브그래프 매칭하여 패턴 일치(검색) 수, 점수(그래프 크기 고려) 측정

        c2_edge = pattern.es['label']

        for i, (graph, count, score) in enumerate(self.P):
            # c1 = graph.vs['label']
            c1_edge = graph.es['label']

            if len(c1_edge) != len(c2_edge):
                continue

            if graph.isomorphic_vf2(pattern):
                # match found, update counter and score
                # tmp_p.append(i)
                new_count = count+1
                new_score = self.get_score(graph, new_count)
                self.P[i] = (graph, new_count, new_score)
                return

        self.trim_dictionary()

        # If pattern is not in dictionary, add it
        count = 1
        score = self.get_score(pattern, count)
        self.P.append((pattern, count, score))

    def update_dictionary_idx(self, pattern, j, k):
        self._update_call_counter += 1
        # print(f"update function call counter ==> {self._update_call_counter}")
        """ Update the pattern dictionary with a new graph
        Iterate through the pattern dictionary:
            If the pattern already exists, update its counter
            Otherwise, add the new pattern to the dictionary
        """
        # pattern -> 그래프
        # 해당 함수에서는 패턴을 넣는것이 아니라 패턴을 추출한 뒤에
        # 서브그래프 매칭하여 패턴 일치(검색) 수, 점수(그래프 크기 고려) 측정

        c2_edge = pattern.es['label']

        for i, (graph, count, score) in enumerate(self.P):
            # c1 = graph.vs['label']
            c1_edge = graph.es['label']

            if len(c1_edge) != len(c2_edge):
                continue

            if graph.isomorphic_vf2(pattern):
                # match found, update counter and score
                # tmp_p.append(i)
                new_count = count+1
                self.append_expanded_pattern(j, i)
                new_score = self.get_score(graph, new_count)
                self.P[i] = (graph, new_count, new_score)
                return

        self.trim_dictionary()

        # If pattern is not in dictionary, add it
        count = 1
        score = self.get_score(pattern, count)

        self._pattern_from[k] = j
        self.append_expanded_pattern(j, k)

        self.P.append((pattern, count, score))

    # def add_pattern_fp_growth(self, old_pp, new_pp):
    #     # 현재는 그래프의 ID 를 가지고 키로 사용함 문제 발생 가능성 있음
    #     old_p = id(old_pp)
    #     new_p = id(new_pp)
    #     if old_p not in self._pattern_to_path:
    #         tmp_idx = len(self._pattern_to_path)
    #         self._pattern_to_index[old_p] = tmp_idx
    #         self._pattern_to_path[old_p] = [tmp_idx]
    #     else:
    #         old_idx = self._pattern_to_index[old_p]
    #         new_idx = old_idx + 1
    #         old_path = self._pattern_to_path[old_p]
    #         new_path = old_path.append(new_idx)
    #         self._pattern_to_path[new_p] = new_path

    def print_graph_edges(self, Es):
        for E in Es:
            print(f"{E.tuple} {E.attributes()}")

    def print_vertexes(self, vs):
        for i in range(len(vs)):
            print(f"{i} -> {vs[i]}")

    def print_graph_patterns(self):
        object_graph = self.P
        for i in range(len(object_graph)):
            p, c, s = object_graph[i]
            print(f"Pattern Graph #[{i}] Count: [{c}] Score: [{s}]")
            print("------------graph detail------------")
            for E in p.es:
                print(f"Edges: {E.tuple} Attr: {E.attributes()}")
            print("------------------------------------")

    def iterate_batch(self, G_batch):
        """ `Compress` a single graph stream object G_batch """
        # 기존 알고리즘과 다르게 한번의 배치처리로 패턴의 MAX SIZE 미만의 반복 수행
        if len(G_batch.es) == 0:
            return
        taken = defaultdict(lambda: False)
        new_patterns = []
        # print(f"batch loop ... {self._batch_count} len P : {len(self.P)}")
        self._batch_count += 1

        # 아래는 2회 이상일 때
        for i, (p, c, s) in enumerate(self.P):
            # 패턴 그래프의 엣지 수를 세어 일정크기 이상이면 루프 돌지 않음
            if len(p.es) >= self._max_pattern_size:
                print(f"Pattern P is greater than _max_pattern_size ...")
                continue

            if len(p.es) >= len(G_batch.es):
                # 첫번째 루프에는 정해진 패턴이 없기때문에 맵이 빈 상태로 시작
                # 저장된 패턴 변수 P는 update dictionary 함수에서 갱신
                maps = []
            else:
                # 두번째 이후의 루프에서는 패턴이 저장되어있고 이를
                maps = G_batch.get_subisomorphisms_vf2(p)
            # For each instance of i of p in B
            # 'vmap' is a mapping of p's vertex indices to G's v. indices
            counter = 0
            # 그래프에서 찾은 subgraph set에 속한 각 subgraph에 대해 반복 수행
            # 첫 루프에서는 돌지 않음
            # print(f"subgraph match result => {len(maps)}")
            for v_map in maps: # 이미 존재하는 패턴에 대해서 도는거
                counter += 1
                # Extend the pattern by "one degree/layer" as p_new
                p_new = None
                # respective vertex indices for the mappings
                # e.g.
                # p_v: 0, 1, 2, 3, 4, 5, ...
                # G_v: 4, 3, 2, 5, 1, 0, ...
                Gv_to_pv = {G_v: p_v for p_v, G_v in enumerate(v_map)}
                # print(f"gv to pv info => {Gv_to_pv}")
                #각 subgraph를 반복 수행하면서 버텍스를 인덱스로 매칭시켜준다

                # Iterate through the mapped vertices in a single embedding
                # Equivalent to `for G_v,p_v in zip(v_map,range(len(v_map)))'
                #각 매칭된 버텍스들을 반복하면서 체크해준다
                for p_v, G_v in enumerate(v_map):
                    G_edges = G_batch.es[G_batch.incident(G_v)] # 그래프의 엣지를 인덱스로 비교
                    p_edges = p.es[p.incident(p_v)] # 사전의 엣지를 인덱스로 비교
                    # print(f"G_v ==> {G_v}")
                    # self.print_graph_edges(G_edges)

                    # (p MUST be subgraph of G) P는 반드시 그래프 G의 서브그래프
                    # print(len(G_edges) <= len(p_edges))
                    if len(G_edges) <= len(p_edges):
                        continue

                    # Find which extra edges we need to add
                    for Ge in G_edges: # G_batch 그래프의 v의 인접한 edge 들에 대하여
                        # One of these should be the same as G_v
                        # Ge 는 정점 V의 인접인 엣지 정보
                        Gv_source_index = Ge.source
                        Gv_target_index = Ge.target
                        # print(f"ori -> {Ge.tuple}, and -> ({Gv_source_index},{Gv_target_index})")
                        # Gv_source = G_batch.vs[Gv_source_index]
                        # Gv_target = G_batch.vs[Gv_target_index]

                        if Gv_source_index not in Gv_to_pv:  # not in map
                            # p_new:새로운 패턴을 확장해서 넣어주는부분
                            # 소스 노드 인덱스가 없을때 만들어서 넣어주는부분
                            if p_new is None:
                                p_new = p.copy()
                            # p_new.add_vertex(pv_source_index)
                            pv_target_index = Gv_to_pv[Gv_target_index]
                            pv_source_index = p_new.vcount()
                            Gv_to_pv[Gv_source_index] = pv_source_index
                            # print(f"new edge info 1-> {pv_source_index, pv_target_index}")

                            # p_new.add_edge(pv_source_index, pv_target_index, label=Ge['label'])
                            self.safe_add_edge(p_new,
                                               pv_source_index,
                                               pv_target_index,
                                               label=Ge['label'])
                        elif Gv_target_index not in Gv_to_pv:  # not in map
                            # 타겟 노드 인덱스가 없을때 만들어서 넣어주는부분
                            if p_new is None:
                                p_new = p.copy()
                            # p_new.add_vertex(pv_target_index)
                            pv_source_index = Gv_to_pv[Gv_source_index]
                            pv_target_index = p_new.vcount()
                            Gv_to_pv[Gv_target_index] = pv_target_index
                            # print(f"new edge info 2-> {pv_source_index, pv_target_index}")

                            # p_new.add_edge(pv_source_index, pv_target_index, label=Ge['label'])
                            self.safe_add_edge(p_new,
                                               pv_source_index,
                                               pv_target_index,
                                               label=Ge['label'])

                        else:  # source index나 target index가 있는경우(=확장가능한경우)
                            pv_source_index = Gv_to_pv[Gv_source_index]
                            pv_target_index = Gv_to_pv[Gv_target_index]
                            # print(f"1 >>> {Gv_source_index} => {pv_source_index}")
                            # print(f"2 >>> {Gv_target_index} => {pv_target_index}")

                            # print('DEBUG IS CON? -> ', p.are_connected(pv_source_index, pv_target_index))
                            # 패턴 p(subgraph)가 사전에 존재하지 않으면,
                            # 소문자 p는 이미 저장된 패턴

                            # self.print_graph_edges(p.es)
                            # print("--")
                            # print(f"\n{pv_source_index} -> {pv_target_index}")

                            # for key in Gv_to_pv.keys():
                            #     print("[", key, "]:[", Gv_to_pv[key], "]")

                            # print(p)
                            # self.print_vertexes(p.vs)
                            try:
                                if not p.are_connected(pv_source_index, pv_target_index):
                                    # print("not connected")
                                    if p_new is None:
                                        p_new = p.copy()
                                    self.safe_add_edge(p_new,
                                                       pv_source_index,
                                                       pv_target_index,
                                                       label=Ge['label'])
                            except Exception as e:
                                print(repr(e))
                                traceback.print_exc()
                                print(Gv_to_pv)
                                print(f"1 >>> {Gv_source_index} => {pv_source_index}")
                                print(f"2 >>> {Gv_target_index} => {pv_target_index}")

                        # Third possibility is that the edge exists in both the
                        # pattern and the larger (batch) graph (do nothing)

                        # Mark the subgraph and the extra edge as taken
                        taken[Ge] = True

                # Add the new pattern to the dictionary
                if p_new is not None:
                    # print(f"new pattern size -> {len(p_new.vs)}")
                    # 여기에 Src Idx to Tgt Idx 추가
                    target_idx = len(self.P)
                    # tmp_p = [self._time_stamp, i, target_idx]
                    # self._pattern_list.append(tmp_p)
                    # print(f">>>  {i} -> {target_idx}")
                    # self.add_pattern_fp_growth(p, p_new)
                    new_patterns.append((p_new, i, target_idx))

                print(f"loop count : {counter}")
        for (g, i, j) in new_patterns:
            # self.update_dictionary(g)
            self.update_dictionary_idx(g, i, j)

        # Add remaining edges in B as single-edge patterns in P
        # LOOP 1 => 첫 패턴에서 모든 단일 엣지가 아래에서 입력

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

    def parse_line_hb(self, line_str, G_batch):
        if line_str == '':
            return
        raw = line_str.strip().split()
        # src, tgt, lbl(MT, RT, RE)
        # 소스, 타겟, 레이블(멘션, 리트윗, 리플)
        # 리플은 문자열을 숫자로 치환하여 사용함 (딕셔너리 사용)
        # timestamp 고려하지않은 상태로 진행
        e_source_id, e_dest_id = raw[0], raw[1]
        e_label = self.get_label_idx(raw[2])
        try:
            if not G_batch.are_connected(e_source_id, e_dest_id):
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
            # print("no vertex info ... add vertexes and edges")
            # We can only add an edge if the vertices already exist
            if self.add_implicit_vertices:
                self.safe_add_edge(G_batch,
                                   e_source_id, e_dest_id, label=e_label)
            else:
                print("Error: vertex in line DNE:\n%s" % line_str,
                      file=stderr)
                raise

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
                if (edge_count % self.batch_size) != 0:
                    continue

                # 이 부분에 패턴 업데이트 하는 알고리즘

                # -------------------------------------

                # Processed the batch, then create a fresh stream object/graph
                self.iterate_batch(G_batch)
                G_batch = Graph(directed=self._directed)
                self._time_stamp += 1

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
file3 = 'data/4rec.txt'
file4 = 'data/10.txt'
file5 = 'rand.txt'
file6 = 'data/144.txt'

comp4 = Compressor(batch_size=10, dict_size=10)
comp4.compress_file_hb(file6)

# 첫번째 원소는 timestamp 의미함
transactions = comp4._expanded_pattern_list[1:]
patterns = pyfpgrowth.find_frequent_patterns(transactions, 3)
rules = pyfpgrowth.generate_association_rules(patterns, 0.5)
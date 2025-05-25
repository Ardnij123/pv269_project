#!/bin/python3.11
import argparse
import random
from math import log, log1p
from typing import Dict, Tuple
from heapq import nlargest
import itertools
from sys import stderr

HELP="""
This script provides a way of search for repeated sequence in fasta file.

The algorithm consists of building an automaton with k-meres from sequence
on nodes and probabilities of adding next character on edges. In the second
step the algorithm nondeterministically chooses a cycle in graph. Then it
searches the sequence for the subsequence defined by the cycle using algorithm
similiar to the one of `seqtk telo`.
"""

Kmere = str
AminoAcid = str
AA = AminoAcid


allowed_chars = {'A', 'C', 'G', 'T'}
unknowns = {'N'}
allowed_whitespace = {'\n'}

zeros = {char: 0 for char in allowed_chars}

rescale = {
        'log1p': log1p,
        'no-scale': (lambda x: x),
        }


parser = argparse.ArgumentParser(
    prog='repeats',
    description=HELP,
)

parser.add_argument('fasta_file', help='Bedpe file with sequences to search through')
parser.add_argument('-k', '--kmer_len', type=int, help='''Set length of k-meres''', default=20)
parser.add_argument('-t', '--abs-threshold', type=float, help='''Minimal number a k-mere has to be in the sequence to be considered to be searched for''', default=3)
parser.add_argument('-T', '--rel-threshold', type=float, help='''Only kmeres that are more frequent than REL_THRESHOLD * max frequency after rescaling''', default=0)
parser.add_argument('-S', '--scaling', type=str, help='''Add scaling of the frequences''', default='log1p', choices=rescale)

parser.add_argument('-m', '--max-drop', type=int, help='''Maximal drop in score to still count as one sequence''', default=200)
parser.add_argument('-i', '--insert-pen', type=float, help='''Penalty for insertion''', default=10)
parser.add_argument('-g', '--gap-pen', type=float, help='''Penalty for gap''', default=10)
parser.add_argument('-b', '--base-pen', type=float, help='''Penalty that is added in each step''', default=2)
parser.add_argument('-s', '--skip', type=int, help='''Skip first n bases of file''', default=0)
parser.add_argument('-e', '--exact-match', type=int, help='''Speed up the search by looking for exact matches of defined length at start of repetition''', default=7)

args = parser.parse_args()

fasta = args.fasta_file
k_len = args.kmer_len if args.kmer_len else 10


def stderrprint(string):
    print(string, file=stderr)


class CharNotAllowed(Exception):
    pass


"""
Generator of sequence from fasta file.

When reading character that is not allowed, raises:
    CharNotAllowed((sequence_id, char, kmere, position_in_seq, line))
Else yields tuple corresponding to:
    (sequence_id, added_char, kmere, position_in_seq, line)
"""
def fasta_reader(file, k=10, skip=0):
    with open(file, 'r') as fa:
        seq_id = ""
        kmere = ""
        pos = 0
        for line_no, line in enumerate(fa):
            if line[0] == '>':
                seq_id = line[1:-1]
                pos = 0
                continue

            for char in line:
                if char in unknowns:
                    kmere = ""
                    continue
                elif char in allowed_whitespace:
                    continue
                elif char not in allowed_chars:
                    raise CharNotAllowed((seq_id, char, kmere, pos, line_no))

                kmere = (kmere + char)[-k:]
                if skip <= 0:
                    yield (seq_id, char, kmere, pos, line_no)
                else:
                    skip -= 1
                pos += 1

class base_reader:
    def __init__(self, file):
        self.fa = open(file, 'r')
        self.seq_id = ""
        self.buffer = []
        self.offset = 0
        self.start = 0
        self.seq_offset = 0
        self.reader = self.file_read()
        self.seq_offsets = {}

    def __iter__(self):
        return self

    def __next__(self):
        if self.offset < len(self.buffer):
            char = self.buffer[self.offset]
            self.offset += 1
            return (self.seq_id, char)

        self.seq_id, char = next(self.reader)
        self.buffer.append(char)
        self.offset += 1
        return (self.seq_id, char)

    def file_read(self):
        for line in self.fa:
            if line[0] == '>':
                self.seq_id = line[1:-1]
                self.seq_offsets[self.seq_id] = self.start + self.offset
                continue

            for char in line:
                if char in unknowns:
                    kmere = ""
                    continue
                elif char in allowed_whitespace:
                    continue
                elif char not in allowed_chars:
                    raise CharNotAllowed((seq_id, char))

                yield (self.seq_id, char)

    def reset(self, start):
        assert start >= self.start
        if self.start + len(self.buffer) > start:
            self.buffer = self.buffer[(start-self.start):]
        else:
            for _ in range(start - self.start - len(self.buffer)):
                next(self.reader)
            self.buffer = []
        self.start = start
        self.offset = 0

    def get_offset(self, sequence):
        return self.seq_offsets[sequence]


aamapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def aamapper(aminoacid):
    return aamapping[aminoacid]


class Node:
    def __init__(self, idx, kmere):
        self.neighbors = [None, None, None, None]
        self._neighbors = []
        self.idx = idx
        self.kmere = kmere

    """
    def __str__(self):
        return self.kmere
    """


class Graph:
    def __init__(self, graph, suffix=1):
        self.kmere_mapping = dict()
        self.nodes = []
        for idx, kmere in enumerate(graph):
            node = Node(idx, kmere)
            self.kmere_mapping[kmere] = node
            self.nodes.append(node)

        for kmere in graph:
            for aminoacid, value in graph[kmere].items():
                self.assign(kmere, aminoacid, value)

        self.by_suffix = []
        self.by_suffix.append({"": self})
        for suf_len in range(1, suffix+1):
            _suffix = dict()
            for suf, nodes in self.by_suffix[-1].items():
                for base in allowed_chars:
                    _suffix[base+suf] = []
                for node in nodes:
                    _suffix[node.kmere[-suf_len:]].append(node)
            self.by_suffix.append(_suffix)

        self.bases = [[], [], [], []]
        for node in self.nodes:
            for base, value, other in node._neighbors:
                self.bases[base].append(node)

        self.__length = len(graph)

    def assign(self, kmere, aminoacid, value):
        aa = aamapper(aminoacid)
        other = (kmere + aminoacid)[-k_len:]
        # TODO: For some reason kmere mapping does not contain kmere sometimes
        # This might be due to N's and ends of sequences
        if other not in self.kmere_mapping:
            return
        self.kmere_mapping[kmere].neighbors[aa] = (value, self.kmere_mapping[other])
        self.kmere_mapping[kmere]._neighbors.append((aa, value, self.kmere_mapping[other]))

    def __len__(self):
        return self.__length

    def __iter__(self):
        yield from self.nodes


def generate_graph(fasta, k_len):
    graph: Dict[Kmere, Dict[AA, int]] = {"": zeros.copy()}

    old_kmere = ""
    for _, base, kmere, _, _ in fasta_reader(fasta, k_len):
        graph[old_kmere][base] += 1

        if kmere not in graph:
            graph[kmere] = zeros.copy()

        old_kmere = kmere
    return graph


def scale_graph(graph, scaling="log1p"):
    scaling = rescale[scaling]
    for kmere in graph:
        graph[kmere] = {base: scaling(count) for base, count in graph[kmere].items() if count > 0}
    return graph


def prune_graph(graph, abs_threshold=-1, rel_threshold=0):
    max_value = -1
    for kmere, edges in graph.items():
        for base, value in edges.items():
            max_value = max(max_value, value)

    threshold = max(max_value * rel_threshold / 100, abs_threshold)
    for kmere in graph:
        graph[kmere] = {base: value for base, value in graph[kmere].items() if value >= threshold}
    graph = {kmere: edges for kmere, edges in graph.items() if edges}
    return graph


# Legacy
def pick_cycle(graph):
    gone = {kmer: -1 for kmer in graph}

    kmer = list(graph)[random.randrange(len(graph))]
    sequence = ''

    while gone[kmer] == -1:
        gone[kmer] = len(sequence)
        out = graph[kmer]
        if not any(out.values()):
            raise Exception("Cycle not found, ended terminal node")

        total = sum(out.values())
        chance = random.randrange(total)
        for char, value in out.items():
            chance -= value
            if chance < 0:
                break

        sequence += char
        kmer += char
        if len(kmer) > k_len:
            kmer = kmer[-k_len:]
    return sequence


# Legacy
def sequence_loop(sequence):
    while True:
        yield from sequence


"""
Similiarity search algorithm using nondeterministic automaton with valuated
currently active states.

This algorithm works similiarly to the Needleman-Wunch algorithm in a sense
that the run of the algorithm can be rewritten as a matrix in which only
moves along directions composing main diagonal are allowed.
This algorithm however allows looping along the axis with queried sequence.

Pseudocode:

Init:
Graph G, V(G) - states of automaton, E(G) allowed transitions by reading sequence
    weights of edges correspond to increment of value by moving along the edge
Sequence corresponding to the graph

Parameters:
MaxDrop, InsertionPenalty, GapPenalty, BasePenalty

1) Values[0] := {(v: 0) for v in V(G)}; Starts[0] := {(v: 0) for v in V(G)};
   MaxValue := 0; MaxPosition = 0; MinPosition = 0; Position = 0

2) Do {
    Position := Position + 1
    Base := following base on input

    Set Value[Position][v] to maximum of:
        ( Value[Position-1][v] - InsertionPenalty )
        ( Value[Position][v-1] - GapPenalty ) v-1 is node with edge comming to v
        ( Value[Position-1][v-1] + weight of the transition edge ) if v-1 has correct transition

    Values[Position] := Values[Position] - BasePenalty
    Set Starts[Position] according to the states

    For each State v s.t. Values[Position][v] <= 0; Do
        Values[Position][v] = 0
        Starts[v] = Position
    Done

    Score := max(Values)
    If MaxValue < Score; Then
        MaxValue := Score
        MaxPosition := Position
        MinPosition := Starts[argmax(Values)]
    Fi

    Optionally: Remove states v from Values where Values[v] < MaxValue - MaxDrop

} While Score >= MaxValue - MaxDrop  

Return MinPosition, MaxPosition
"""
def single_search(sequence, graph, MaxDrop=200, InsertionPenalty=3, GapPenalty=3, BasePenalty=1,
                  offset=0, exact_start=3):
    max_value = -1
    max_position = offset+exact_start
    min_position = offset
    position = offset
    chrom = ""
    flood = 0
    batches = int(MaxDrop / GapPenalty) + 1

    try:
        start_seq = ""
        for _ in range(exact_start):
            chrom, base = next(sequence)
            start_seq += base
            position += 1
        next_states = [list(map(lambda x: (True, x), graph.by_suffix[exact_start][start_seq]))]

        if not next_states[0]:
            return (chrom, offset, offset+exact_start, 0, offset+exact_start)

        values = [(-1, 0) for _ in range(len(graph))]
        next_values = [(-1, 0) for _ in range(len(graph))]
        starts = [offset for _ in range(len(graph))]
        next_starts = [offset for _ in range(len(graph))]
        max_value = 0

        while any(next_states):
            current_states, next_states = next_states, [[] for _ in range(batches)]
            values, next_values = next_values, values
            starts, next_starts = next_starts, starts

            _chrom, base = next(sequence)
            if chrom != _chrom:
                break

            base = aamapper(base)
            position = position + 1

            _flood = flood + max(0, max_value - MaxDrop)


            for batch in current_states:
                for allow_insert, state in batch:
                    state_idx = state.idx
                    start, value = values[state_idx]
                    if value < _flood:
                        continue

                    if value <= 0:
                        state_start = position - 1
                    else:
                        state_start = starts[state_idx]

                    # Gaps + Correct base
                    gap_value = value - GapPenalty
                    gapped = (position, gap_value)
                    for nbase, increment, neigh in state._neighbors:
                        n_idx = neigh.idx
                        if nbase == base:
                            # Correct
                            new_value = value + increment
                            correct = (position, new_value)
                            if next_values[n_idx] < correct:
                                if new_value - flood > max_value:
                                    max_value = new_value - flood
                                    min_position = state_start
                                    max_position = position
                                    _flood = flood + max(0, max_value - MaxDrop)
                                if next_values[n_idx] < (position, -1):
                                    next_states[int((max_value+flood - new_value)/GapPenalty)].append((True, neigh))
                                next_values[n_idx] = correct
                                next_starts[n_idx] = state_start
                        else:
                            # Gap
                            if gapped > values[n_idx] and gap_value >= _flood:
                                values[n_idx] = gapped
                                current_states[int((max_value+flood - gap_value)/GapPenalty)].append((False, neigh))

                    # Insertion
                    if allow_insert and (new_value := value - InsertionPenalty) > _flood:
                        insertion = (position, new_value)
                        if next_values[state_idx] < insertion:
                            next_states[int((max_value+flood - new_value)/GapPenalty)].append((True, state))
                            next_values[state_idx] = insertion
                            next_starts[state_idx] = state_start

            flood += BasePenalty

        return chrom, min_position, max_position, max_value, position


    except StopIteration:
        return chrom, min_position, max_position, max_value, position


def repeats_search(fasta, graph, MinValue=200, MaxDrop=200, InsertionPenalty=3, GapPenalty=3, BasePenalty=1, fast_skip=True, skip=0, exact_start=3):
    stderrprint("# Starting search procedure")
    position = skip
    value = 0
    last_report = 0
    reader = base_reader(fasta)
    reader.reset(skip)

    while value >= 0:
        try:
            reader.reset(position)
        except StopIteration:
            return
        chrom, start, end, value, p_end = \
                single_search(reader, graph,
                              MaxDrop=MaxDrop,
                              InsertionPenalty=InsertionPenalty,
                              GapPenalty=GapPenalty,
                              BasePenalty=BasePenalty,
                              offset=position,
                              exact_start=exact_start)

        if value > MinValue:
            offset = reader.get_offset(chrom)
            _start = start - offset
            _end = end - offset

            last_report = end
            yield chrom, _start, _end, value

        if fast_skip:
            position = p_end
        else:
            position = end

        if last_report + 50000 <= position:
            last_report = position
            stderrprint(f"# Now at base: {position}")


stderrprint("# Generating graph of k-meres")
graph = generate_graph(fasta, k_len)
stderrprint("# Scaling graph of k-meres")
graph = scale_graph(graph, args.scaling)
stderrprint("# Pruning graph of k-meres")
graph = prune_graph(graph, args.abs_threshold, args.rel_threshold)
stderrprint("# Transforming graph of k-meres")
graph = Graph(graph, suffix=args.exact_match)


for chrom, start, end, value in repeats_search(
        fasta, graph,
        MaxDrop=args.max_drop,
        InsertionPenalty=args.insert_pen,
        GapPenalty=args.gap_pen,
        BasePenalty=args.base_pen,
        skip=args.skip,
        exact_start=args.exact_match):
    print(chrom, start, end, value)

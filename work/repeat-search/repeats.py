#!/bin/python3.11
import argparse
import random
from math import log, log1p

from typing import Dict, Tuple

Kmere = str
AminoAcid = str
AA = AminoAcid


HELP="""
This script provides a way of search for repeated sequence in fasta file.

The algorithm consists of building an automaton with k-meres from sequence
on nodes and probabilities of adding next character on edges. In the second
step the algorithm nondeterministically chooses a cycle in graph. Then it
searches the sequence for the subsequence defined by the cycle using algorithm
similiar to the one of `seqtk telo`.
"""

parser = argparse.ArgumentParser(
	prog='repeats',
	description=HELP,
)

rescale = {
        'log1p': log1p,
        'no-scale': (lambda x: x),
        }

parser.add_argument('fasta_file', help='Bedpe file with sequences to search through')
parser.add_argument('-k', '--kmer_len', type=int, help='''Set length of k-meres''')
parser.add_argument('-n', '--n_tries', type=int, help='''Try searching for repeats n-times''')
parser.add_argument('-B', '--out_bed', type=str, help='''Generate annotation of hits as a bed file''')
parser.add_argument('-t', '--abs-threshold', type=int, help='''Minimal number a k-mere has to be in the sequence to be considered to be searched for''', default=1)
parser.add_argument('-T', '--rel-threshold', type=float, help='''Only kmeres that are more frequent than REL_THRESHOLD * max frequency after rescaling''', default=0)
parser.add_argument('-S', '--scaling', type=str, help='''Add scaling of the frequences''', default='log1p', choices=rescale)

parser.add_argument('-m', '--max-drop', type=int, help='''Maximal drop in score to still count as one sequence''', default=200)
parser.add_argument('-i', '--insert-pen', type=float, help='''Penalty for insertion''', default=3)
parser.add_argument('-g', '--gap-pen', type=float, help='''Penalty for gap''', default=3)
parser.add_argument('-b', '--base-pen', type=float, help='''Penalty that is added in each step''', default=1)

args = parser.parse_args()

fasta = args.fasta_file
k_len = args.kmer_len if args.kmer_len else 10
allowed_chars = {'A', 'C', 'G', 'T'}
unknowns = {'N'}
allowed_whitespace = {'\n'}

zeros = {char: 0 for char in allowed_chars}


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
                seq_id = line[1:]
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

def base_reader(file, skip=0):
    with open(file, 'r') as fa:
        seq_id = ""

        for line in fa:
            if line[0] == '>':
                seq_id = line[1:]
                continue

            for char in line:
                if char in unknowns:
                    kmere = ""
                    continue
                elif char in allowed_whitespace:
                    continue
                elif char not in allowed_chars:
                    raise CharNotAllowed((seq_id, char))

                if skip <= 0:
                    yield (seq_id, char)
                skip -= 1



print("Generating graph of k-meres")
graph: Dict[Kmere, Dict[AA, int]] = {"": zeros.copy()}

old_kmere = ""
for _, base, kmere, _, _ in fasta_reader(fasta, k_len):
    graph[old_kmere][base] += 1

    if kmere not in graph:
        graph[kmere] = zeros.copy()

    old_kmere = kmere

scaling = rescale[args.scaling]
for kmere in graph:
    graph[kmere] = {base: scaling(count) for base, count in graph[kmere].items()}

max_value = -1
for kmere, edges in graph.items():
    for base, value in edges.items():
        max_value = max(max_value, value)

threshold = max(max_value * args.rel_threshold / 100, args.abs_threshold)
for kmere in graph:
    graph[kmere] = {base: value for base, value in graph[kmere].items() if value >= threshold}
graph = {kmere: edges for kmere, edges in graph.items() if edges}

print("Graph generated")


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
        ( Value[Position-1][v-1] - GapPenalty ) v-1 is node with edge comming to v
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
                  offset=0):
    max_value = -1
    max_position = 0
    min_position = 0
    position = 0
    chrom=""

    values = {kmere: 0 for kmere in graph}
    starts = {kmere: 0 for kmere in graph}

    try:
        for chrom, base in sequence:
            new_values = dict()
            new_starts = dict()
            position = position + 1

            for state in values:
                insertion = values[state] - InsertionPenalty
                if new_values.get(state, -1) < insertion:
                    new_values[state] = insertion
                    new_starts[state] = starts[state]

                for next_base, value in graph.get(state, dict()).items():
                    next_state = (state + next_base)[1:]
                    if next_base == base:
                        correct = values[state] + value
                        if new_values.get(next_state, -1) < correct:
                            new_values[next_state] = correct
                            new_starts[next_state] = starts[state]
                    else:
                        gap = values[state] - GapPenalty
                        if new_values.get(next_state, -1) < gap:
                            new_values[next_state] = gap
                            new_starts[next_state] = starts[state]
            
            values = {state: value - BasePenalty for state, value in new_values.items()}
            starts = new_starts

            for state in values:
                if values[state] <= 0:
                    values[state] = 0
                    starts[state] = position

            score, state = max([(values[state], state) for state in values], default=(-1, ""))
            if max_value < score:
                max_value = score
                max_position = position
                min_position = starts[state]

            if score < max_value - MaxDrop or score == -1:
                break

            values = {state: values[state] for state in values if values[state] >= max_value - MaxDrop}

        return chrom, min_position+offset, max_position+offset, max_value


    except StopIteration:
        return chrom, min_position+offset, max_position+offset, max_value


def repeats_search(fasta, graph, MinValue=200, MaxDrop=200, InsertionPenalty=3, GapPenalty=3, BasePenalty=1):
    position = 0
    value = 0
    last_report = 0

    while value >= 0:
        sequence = base_reader(fasta, skip=position)
        chrom, start, end, value = \
                single_search(sequence, graph, MaxDrop, InsertionPenalty, GapPenalty, BasePenalty, offset=position)

        if value > MinValue:
            last_report = end
            yield chrom, start, end, value

        position = end

        if last_report + 10000 < position:
            last_report = position
            print(f"Now at base: {position}")

        if position > 10000:
            return


for chrom, start, end, value in repeats_search(
        fasta, graph,
        MaxDrop=args.max_drop,
        InsertionPenalty=args.insert_pen,
        GapPenalty=args.gap_pen, BasePenalty=args.base_pen):
    print(chrom, start, end, value)

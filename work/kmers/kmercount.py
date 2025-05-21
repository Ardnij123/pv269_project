#!/bin/python3

from sys import argv
import math

file = argv[1]
window = 2

width = 180
all_lines = True


def result_table(results, width=width, show=100):
    results = sorted(results.items(), key=lambda x: x[1], reverse=True)
    results = results[:show]
    show = len(results)
    max_num = math.ceil(math.log10(results[0][1]))

    cols = width // (window + max_num + 3)
    rows = math.ceil(show / cols)

    for i in range(rows):
        for j in range(cols):
            try:
                slide, num = results[j*rows + i]
                print(slide, ("{:" + str(max_num) + "d}").format(num), end='   ')
            except IndexError:
                break
        print()


with open(file, 'r') as f:
    ch = ' '
    while len(ch) > 0:
        slides = dict()

        chrom = f.readline()[:-1]
        if len(chrom) > 0:
            print(f"Stats for {chrom}")  # header
        
        word = f.read(window)
        if len(word) == 0:
            break

        slides[word] = 1
        while (ch := f.read(1)) != '\n' and len(ch) > 0:
            word = word[1:] + ch
            slides[word] = slides.get(word, 0) + 1

        result_table(slides)

        if not all_lines:
            break

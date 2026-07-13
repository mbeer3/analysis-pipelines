#!/usr/bin/env python

import sys
import os.path
import argparse

VERSION = "3.2.1"


def byN(array, N):
    i = 0
    out = []
    while i < len(array):
        out.append(array[i])
        i += 1
        if len(out) == N:
            yield out
            out = []


def parseLibraryFile(fname):
    libraryDb = dict()
    dat = open(fname).read().strip().split()

    if len(dat) % 3 != 0:
        print >>sys.stderr, "ERROR: Library file malformed"
        print >>sys.stderr, "Number of fields =", len(dat)
        sys.exit(1)

    for (name, seq, libName) in byN(dat, 3):
        if seq in libraryDb:
            names = libraryDb[seq][3]
            names.append(name)
            libraryDb[seq] = [name, seq, libName, names]
        else:
            libraryDb[seq] = [name, seq, libName, [name]]

    return libraryDb


def parseCountFile(fname):
    countDb = dict()
    dat = open(fname).read().strip().split()

    if len(dat) % 2 != 0:
        print >>sys.stderr, "ERROR: Count file", fname, "malformed"
        print >>sys.stderr, "Number of fields =", len(dat)
        sys.exit(1)

    for (count, seq) in byN(dat, 2):
        countDb[seq] = int(count)

    return countDb


parser = argparse.ArgumentParser(prog="joinRNAi.py")
parser.add_argument("-o", "--output", help="Name of output joined file", required=False, default="join.txt")
parser.add_argument("libraryFile", nargs=1, help="Library file")
parser.add_argument("countFiles", nargs="+", help="One or more count files")
parser.add_argument("--version", action="version", version="%(prog)s " + VERSION)

args = parser.parse_args()

libraryFile = args.libraryFile[0]
library = parseLibraryFile(libraryFile)

counts = dict()

for countFile in args.countFiles:
    baseName = os.path.splitext(os.path.basename(countFile))[0]
    counts[baseName] = parseCountFile(countFile)

with open(args.output, "w") as fp:
    samples = sorted([x for x in counts])

    print >>fp, "\t".join(["NAME", "ALIAS", "SEQ", "LIB"] + samples)

    for seq in library:
        (name, _, libName, alias) = library[seq]

        if len(alias) == 1:
            aliasNames = ""
        else:
            aliasNames = ";".join(alias)

        out = [name, aliasNames, seq, libName]

        for si in samples:
            if seq in counts[si]:
                out.append(counts[si][seq])
            else:
                out.append(0)

        print >>fp, "\t".join(map(str, out))

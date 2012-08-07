#!/usr/bin/python

# Time-stamp: <2005-04-23 10:22:30 kasper>

import re
import string
import sys

class FastaIndex:
    "Class for indexing Fasta files"
    def __init__(self, filename=''):
        self.filename = filename
        self.file = open(filename, 'r')
        # See is there is pickeled version of this. If so return that instead

        regex = re.compile(">(\S+)")
        prevStart = 0
        start = 0
        sequenceStart = 0
        self.dir = {}
        prevId = ''
        while 1:
            line = self.file.readline()
            if not line:
                break
            match = regex.search(line)
            if match:
                id = match.group(1)
                start = self.file.tell() - len(line)
                sequenceStart = self.file.tell()
                if prevId:
                    self.dir[prevId] = [prevStart, sequenceStart, start - 1 - prevStart]
                prevId = id
                prevStart = start
                prevSequenceStart = sequenceStart
        # Dump the modules unless it was read in from disk

    def getSeq(self, id, start=-1, end=-1):
        offset, sequenceOffset, length = self.dir[id]
        self.file.seek(sequenceOffset)
        self.file.seek(start, 1)
        length = end - start
        seq = self.file.read(length)
        return seq
    
    def getLFasta(self, id):
        offset, sequenceOffset, length = self.dir[id]
        self.file.seek(offset)
        entry = self.file.read(length)








index = FastaIndex('/home/kasper/test/mRNA.fa')

lfa = index.getLFasta('NM_068129')
print lfa

entry = index.get('NM_068134', 0, 4)
print entry

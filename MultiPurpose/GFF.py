#!/usr/bin/python

# Time-stamp: <2005-08-09 11:31:43 kasper>

import re
import string

class GFF:
    def __init__(self, list):
        self.sanityCheck(list)
        self.seq = list[0]
        self.source = list[1]
        self.feature = list[2]
        self.start = int(list[3])
        self.end = int(list[4])
        self.score = list[5]
        self.strand = list[6]
        self.frame = list[7]
        if len(list) > 8 and list[8] != '.':
            self.attrib = {}
            for s in re.split("\s*;\s*", list[8]):
                pair = s.split("=")
#                 if len(pair[1:]) > 2:
#                     self.attrib[pair[0]] = pair[1:]
#                 else:
#                     self.attrib[pair[0]] = pair[1]
                if len(pair) == 2:
                    self.attrib[pair[0]] = pair[1]
        else:
            self.attrib = ''
    def __cmp__(self, other): 
        return cmp(self.seq, other.seq) \
               or cmp(self.start, other.start) \
               or cmp(self.end, other.end)
    def getline(self):
        list = [self.seq, self.source,
                self.feature, str(self.start),
                str(self.end), self.score,
                self.strand, self.frame]
        s = "\t".join(list)
        if self.attrib:
            list = ["=".join(x) for x in self.attrib.items()]
            s += "\t" + ";".join(list)
        s += "\n"
        return s
    def sanityCheck(self, list):
        assert len(list) == 8 or len(list) == 9, "Not appropriate nr of entries"
        assert str(list[3]).isdigit(), "Start not digits: " + str(list[3])
        assert str(list[4]).isdigit(), "End not digits: " + str(list[4])
    def overlaps(self, other):
        if self.seq == other.seq \
               and (self.start >= other.start and self.start <= other.end) \
               or self.end >= other.start and self.end <= other.end:
            return 1
        else:
            return 0

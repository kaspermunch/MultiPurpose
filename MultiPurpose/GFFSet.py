#!/usr/bin/python

# Time-stamp: <2005-04-14 21:41:29 kasper>

class GFFSet:
    def __init__(self, list=[]):
        self.list = list
        self.seq = list[0].seq
        self.start = list[0].start
        self.end = list[-1].end
        self.strand = list[0].strand
    def __getitem__(self, item):
        return self.list[item]
    def __delitem__(self, item):
        del self.list[item]
    def __cmp__(self, other): 
        return cmp(self.seq, other.seq) \
               or cmp(self.start, other.start) \
               or cmp(self.end, other.end)
    def __contains__(self, other):
        self.list.sort()
        start = self.list[0]
        end = self.list[-1]
        if selv.list.sequence == other.sequence \
           and (other.start >= start and other.start <= end)\
           or other.end >= start and other.end <= end:
            return 1
        else:
            return 0
    def __iter__(self):
        self.index = 0
        return self
    def next(self):
        if self.index < len(self.list):
            self.index += 1
            return self.list[self.index-1]
        else:
            raise StopIteration
    def add(self, other):
        self.list.append(other)
        self.list.sort()
    def len(self):
        return len(self.list)
    def overlaps(self, other):
        if self.seq == other.seq \
               and (self.start >= other.start and self.start <= other.end) \
               or self.end >= other.start and self.end <= other.end:
            return 1
        else:
            return 0


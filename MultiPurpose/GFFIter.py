#!/usr/bin/python

# Time-stamp: <2005-04-22 22:51:49 kasper>

import re
import string
import sys
from GFF import GFF

class GFFIter:
    "Iterator for iterating over GFF entries"
    def __init__(self, file=''):
        if file:
            self.file = file
            self.fh = open(file, 'r')
        else:
            self.fh = sys.stdin    
    def __iter__(self):
        return self
    def next(self):
        line = self.fh.readline()
        while re.match("^#|^\s$", line):
            line = self.fh.readline()
        line = line.strip()
        if line:
            list = line.split("\t")
            gff = GFF(list);        
            if gff:
                return gff
            else:
                raise StopIteration
        else:
            raise StopIteration

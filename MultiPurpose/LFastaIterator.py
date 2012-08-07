#!/usr/bin/python

# Time-stamp: <2005-11-02 13:10:20 kasper>

import re
import sys
from LFasta import LFasta

class LFastaIterator:
    "Iterator for iterating over Labels Fasta entries"
    def __init__(self, file=''):
        if file:
            self.file = file
            self.fh = open(file, 'r')
        else:
            self.fh = sys.stdin    
        self.saved_header = "";
    def __iter__(self):
        return self
    def next(self):
        entry = self.read_entry(self.fh)
        if entry:
            return entry
        else:
            raise StopIteration
    def read_entry(self, fh):
        header = ''
        seq = ''
        labels = ''
        pred = {}
        plpdata = ''
        while (1):
            if self.saved_header:
                line = self.saved_header
                self.saved_header = ''
            else:
                line = fh.readline()
            if not line:
                if header:
                    break
                else:
                    return 0
            line = line.strip()
            if not line:
                continue
            headerline = re.match('>((\S+).*)', line)
            if headerline and not header:
                header = headerline.group(1)
                id = headerline.group(2)
            elif headerline:
                self.saved_header = headerline.group(1)
                break
            elif header:
                comment = re.match('%', line)
                firstplp = re.match('# ' + id, line)
                secondplp = re.match('# {9}\S\s', line)
                plpline = re.match("\S\s\S\s+(?:\s+(?:(?i)(?:[+-]?)(?:(?=[0123456789]|[.])(?:[0123456789]*)(?:(?:[.])(?:[0123456789]{0,}))?)(?:(?:[Ee])(?:(?:[+-]?)(?:[0123456789]+))|))){2,}", line)
                if comment:
                    pass
                elif id and firstplp:
                    plpdata += line
                    continue
                elif secondplp:
                    plpdata += line
                    continue
                elif plpline:
                    plpdata += line
                    continue
                else:
                    labelline = re.match('#\s+(\S.*)', line)
                    predline = re.match('\?(\S)\s+(\S.*)', line)
                    whitespace = re.compile( '\s+')
                    if labelline:
                        str = labelline.group(1)
                        str = whitespace.sub('', str)
                        labels += str
                    elif predline:
                        prefix = predline.group(1)
                        str = predline.group(2)
                        str = whitespace.sub('', str)
                        if pred.has_key(prefix):
                            pred[prefix] += str
                        else:
                            pred[prefix] = str
                    else:
                        str = whitespace.sub('', line)
                        seq += str
        lfasta = LFasta(header=header, seq=seq, labels=labels, pred=pred)
#        lfasta = LFasta(header=header, seq=seq, labels=labels)
        return lfasta


# input = sys.argv[1]
# 
# output = sys.stdout
# 
# for lfa in LFastaIterator(input):
#     lfa.write(output)


# lfa = LFasta("ID", "DESC", "AAACCCTTT", "LLLLLLLLL")
# lfa = LFasta(seq="AAACCCTTT", labels="LLLLLLLLL", id="ID", desc="DESC")
# #lfa = LFasta("", "", "AAACCCTTT", "LLLLLLLLL", header="id desc")

# print lfa.id
# print lfa.desc
# print lfa.seq
# print lfa.labels



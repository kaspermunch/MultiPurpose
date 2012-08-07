#!/usr/bin/python

# Time-stamp: <2005-04-23 16:31:09 kasper>

import re
import sys
import string

            
class LFasta:
    "Holding a Labeled Fasta entry"
    def __init__(self, id="unknow id", desc="", header="", seq="", labels="", pred={}):
        assert header or id, 'No header or id'
        if header:
            self.header = header
            r = re.compile('(\S+)\s*(.*)')
            m = r.match(header)
            self.id = m.group(1)
            self.desc = m.group(2)
        else:
            self.id = id
            self.header = id
            if desc:
                self.desc = desc
                self.header = id + " " + desc
        assert len(seq) == len(labels), 'sequence and labels not equal length'
        self.seq = seq
        self.labels = labels
        self.pred = pred
    def length(self):
        return len(self.seq)
    def write(self, fh=""):
        linelen = 70
        entrystr = ''
        if not re.match(self.header, "^>"):
            entrystr += ">"
        entrystr += self.header + "\n"
        indentation = ""
        if self.labels:
            indentation = '   '
            labelprefix = '#  '
            length = self.length()
        for i in range(0, length, linelen):
            entrystr += indentation + self.seq[i:i+linelen] + "\n"
            if self.labels:
                entrystr += labelprefix + self.labels[i:i+linelen] + "\n"
            if self.pred:
                prefixes = self.pred.keys()
                prefixes.sort()
                for prefix in prefixes:
                    entrystr += "?" + prefix + " "
                    entrystr += self.pred[prefix][i:i+linelen] + "\n"
                entrystr += "\n"
#         if self.plpprint:
#             plpdata = self.plpdata
#             entrystr += "\n"
#             for line if plpdata:
#                 entrystr += line + "\n"
        fh.write(entrystr)

def LFastaGenerator(list):
    "Generates a LFasta entry form a string"
    header = ''
    seq = ''
    labels = ''
    pred = {}
    plpdata = ''
    assert list
    line = list.pop(0)
    headerline = re.match('>((\S+).*)', line)
    assert headerline
    header = headerline.group(1)
    id = headerline.group(2)
    while (1):
        if len(list):
            line = list.pop(0)
        else:
            break
        line = line.strip()
        if not line:
            continue
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
    # lfasta = LFasta(header=header, seq=seq, labels=labels)
    return lfasta

class LFastaIterator:
    "Iterator for iterating over Labels Fasta entries"
    Generator = LFasta
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
        entry = self.readEntry()
        if entry:
            return entry
        else:
            raise StopIteration
    def readEntry(self):
        list = []
        gotHeader = 0
        while (1):
            if self.saved_header:
                line = self.saved_header
                self.saved_header = ''
                gotHeader = 1
            else:
                line = self.fh.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue
            headerline = re.match('^>((\S+).*)', line)
            if headerline and not gotHeader:
                self.gotHeader = 1
                list.append(line)
            elif headerline:
                self.saved_header = headerline.group(1)
                break
            else:
                list.append(line)
        if len(list):
            lfasta = LFastaGenerator(list)
            return lfasta
        else:
            return 0

class LFastaIndex:
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


# index = FastaIndex('/home/kasper/test/mRNA.fa')
# 
# lfa = index.getLFasta('NM_068129')
# print lfa
# 
# entry = index.get('NM_068134', 0, 4)
# print entry




input = sys.argv[1]

output = sys.stdout

for lfa in LFastaIterator(input):
    lfa.write(output)


# lfa = LFasta("ID", "DESC", "AAACCCTTT", "LLLLLLLLL")
# lfa = LFasta(seq="AAACCCTTT", labels="LLLLLLLLL", id="ID", desc="DESC")
# #lfa = LFasta("", "", "AAACCCTTT", "LLLLLLLLL", header="id desc")

# print lfa.id
# print lfa.desc
# print lfa.seq
# print lfa.labels

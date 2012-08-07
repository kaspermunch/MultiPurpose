#!/usr/bin/python

# Time-stamp: <2006-05-02 14:09:31 kasper>

import re
import sys


class LFAstaGenerator:
    "Generates a LFasta entry form a string"
    def __init__(self, str):
        list = "\n".split(str)
        lfasta = read_entry(self, list)
        return lfasta
    def read_entry(self, list):
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
                line = list.pop([0])
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
    #     def revcomp(self):
    #         self.seq = revcomp(self.seq)
    #     def grep(self, regexp):
    #         pass
    def write(self, fh=False):
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

        if fh:
            sys.stdout.write(entrystr)
        else:
            fh.write(entrystr)

#     def grep(self, regex, flags=False):
# 
#         if flags:
#             regex = re.compile(regex, flags)
#         else:
#             regex = re.compile(regex)
# 
#         for m in regex.finditer(seq):
#             m.group(1)




# 
# 
# input = sys.argv[1]
# 
# output = sys.stdout
# 
# for lfa in LFastaIterator(input):
#     lfa.write(output)
# 

# lfa = LFasta("ID", "DESC", "AAACCCTTT", "LLLLLLLLL")
# lfa = LFasta(seq="AAACCCTTT", labels="LLLLLLLLL", id="ID", desc="DESC")
# #lfa = LFasta("", "", "AAACCCTTT", "LLLLLLLLL", header="id desc")

# print lfa.id
# print lfa.desc
# print lfa.seq
# print lfa.labels



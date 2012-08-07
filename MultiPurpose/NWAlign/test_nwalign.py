#!/usr/bin/env python
# standard modules
import sys
import string

# Biopython
from Bio import Fasta
from Bio.ParserSupport import AbstractConsumer

# alignment stuff
from NWAlign import NWAlign
from SubstitutionMatrix import SubstitutionMatrix
    
class FASTAHandler(AbstractConsumer):
    def __init__(self, seq_list):
        self.seq_list = seq_list
        
    def start_sequence(self):
        self.this_seq = ""
            
    def sequence(self, line):
        self.this_seq = self.this_seq + string.strip(line)

    def end_sequence(self):
        self.seq_list.append(self.this_seq)
        
def main():
    # create a substitution matrix
    sub_matrix = SubstitutionMatrix('blosum50')
    
    # set up for alignment
    aligner = NWAlign(sub_matrix)
    print "Testing a simple alignment..."
    seq1 = "HEAGAWGHEE"
    seq2 = "PAWHEAE"
    
    aligner.align(seq1, seq2)
    
    align1, align2 = aligner.get_optimal_alignment()
    score = aligner.get_optimal_score()
    
    print "Alignment Score:", score
    print align1.data
    print align2.data
    
    print "Testing a more complex alignment..."
    test_file = "PEPCarboxylase.fasta"
    
    print "Getting sequences from the file PEPCarboxylase.fasta..."
    seq_list = []
    scanner = Fasta._Scanner()
    handler = FASTAHandler(seq_list)
    file = open(test_file, 'r')
    scanner.feed(file, handler)
    scanner.feed(file, handler)
    #print seq_list
    
    print "Aligning sequences..."
    aligner = NWAlign(sub_matrix)
    aligner.align(seq_list[0][0:150], seq_list[1][0:150])
    
    align1, align2 = aligner.get_optimal_alignment()
    score = aligner.get_optimal_score()
    
    print "Alignment Score:", score
    line_width = 25
    current_position = 0
    current_position = current_position + line_width
    # pretty print the alignment
    while current_position < len(align1):
        print ""
        print align1.data[current_position - line_width:current_position]
        print align2.data[current_position - line_width:current_position]
        current_position = current_position + line_width
        
    # print whatever is left
    print ""
    print align1.data[current_position - line_width:len(align1) - 1]
    print align2.data[current_position - line_width:len(align2) - 1]
    
if __name__ == '__main__':
    sys.exit(main())

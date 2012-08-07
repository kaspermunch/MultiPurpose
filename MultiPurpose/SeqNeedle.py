#!/usr/bin/python

# Time-stamp: <2006-09-05 14:36:26 kasper>

# standard modules
import array

# Biopython classes
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq



import sys
import string

# Biopython
from Bio import Fasta
from Bio.ParserSupport import AbstractConsumer



DEFAULT_GAP_PENALTY = -8
GAP_CHAR = '-'


# Represent substitution matrices for alignments.
# 
# Substitution matrices are not actually represented as a matrix, but rather 
# as a dictionary, with the key being a tuple, and the value being the 
# substition value.
# So the entry for an A to R substitution in a BLOSUM50 would be represented 
# as {('A', 'R'): -2}. Since the matrix is symmetric, this lets us have less 
# objects in memory then if it was defined in a full matrix.


# import available matrices
from matrix_info import available_matrices

def get_available_matrices(self):
    """Return a list of the substitution matrices available."""
    return available_matrices


class SubstitutionMatrix:
    """Represent a substitution matrix of a given type."""
    
    def __init__(self, type):
        """Initialize with the type of substitution matrix to hold.
        
        Arguments:
        * type - A string with the type of the matrix. Use 
        'get_available_matrices' to find out the available types.
        """
        if type in available_matrices:
            exec 'from matrix_info import ' + type
            self.sub_matrix = eval(type)
            self.matrix_type = type
        else:
            raise IndexError \
              ("A %s substitution matrix is not available." % type)
        

        
    def get_sub_value(self, protein1, protein2):
        """Find the substitution value for two proteins.
        
        Since the substitution matrices are symmetric, the order of the
        two protiens passed is not important.
        
        Raises:
        * KeyError - If a substitution value for the two proteins is
        not found in the current matrix
        """
        try:
            sub_value = self.sub_matrix[(protein1, protein2)]
        # if we get one key error, try the proteins arranged in the 
        # opposite way
        except KeyError:
            sub_value = self.sub_matrix[(protein2, protein1)]
        
        return sub_value
 

class DPMatrix:
    """A generic Dynmaic Programming matrix.

    This is an abstract base class and should be inherited to provide a
    specific implementation of an algorithm.
    """
        
    def __init__(self, sub_matrix):
        """Initialize with a substitution matrix to use for alignments.
        
        Arguments:
        * sub_matrix - An initialized substitution matrix from the 
        'Substitution Matrix' class.
        """
        self.sub_matrix = sub_matrix
        self.dpmatrix = {}
        
    def fill_cell(self, column, row):
        pass
        
    def fill_matrix(self, sequence1, sequence2):
        """Fill the dpmatrix via a 'pull' recursion method."""
        self.seq1 = sequence1
        self.seq2 = sequence2
          
        last_cell = self.fill_cell(len(sequence1), len(sequence2))
        
        self.dpmatrix[(len(sequence1), len(sequence2))] = last_cell

            
class DPMatrixCell:
    """An individual cell in a DPMatrix.
    """
    
    def __init__(self, col_pos, row_pos, seq1item, seq2item):
        self.col_pos = col_pos
        self.row_pos = row_pos
        self.seq1item = seq1item
        self.seq2item = seq2item
        
        self.value = None
        self.parent_cell = None
        
    def set_value(self, value):
        self.value = value
        
    def set_parent(self, parent):
        self.parent_cell = parent
        
    def get_value(self):
        if self.value == None:
            raise IndexError('Value not set for this matrix cell.')
        
        return self.value
        
    def get_parent(self):
        
        return self.parent_cell


class Needle(DPMatrix):
    """Class representing a Needleman-Wunsch alignment algorithm.
    """
    
    def __init__(self, sub_matrix = None):
        DPMatrix.__init__(self, sub_matrix)
        
        self.gap_penalty = DEFAULT_GAP_PENALTY

    def set_gap_penalty(self, penalty):
        self.gap_penalty = penalty
        
    def fill_cell(self, column, row):
        """Fill a cell using the Needleman-Wunsch algorithm.
        
        In this algorithm we calculate the value of each cell i,j as the 
        maximum of:
            1. F(i-1, j-1) + s(xi, yj) : The value of cell i-1, j-1 plus
            the substitution matrix value of the two items (Amino Acids or
            Nucleotides) held in the cell.
            2. F(i, j-1) - d : The value of cell i, j-1 minus the gap penalty.
            3. F(i-1, j) - d : The value of cell i-1, j minus the gap penalty.
            
        Returns:
            A DPMatrixCell object that has it's information filled in.
        """
        # test to be sure we aren't falling off the matrix.
        assert (column >= 0)
        assert (row >= 0)
        
        # most special case, if we are at (0,0) in the matrix, this
        # is a zero and has no parent
        if (column == 0) and (row == 0):
            corner_cell = DPMatrixCell(column, row, "", "")
            corner_cell.set_value(0)
            return corner_cell
            
        # flags so we know which scores to calculate
        # default is to calculate them all.
        calculate1 = 1
        calculate2 = 1
        calculate3 = 1  
        # two other special cases, when we are the top row or in the
        # right-most column, then we only need to calculate a
        # single score.
        # Top Row -> only calculate score 3, and this will be our score
        if (row == 0) and (column != 0):
            calculate1 = 0
            calculate2 = 0
        # right-most column -> only calculate score 2 and that is
        # our score
        elif (column == 0) and (row != 0):
            calculate1 = 0
            calculate3 = 0
                
        # Now calculate the three scores
        # Score 1
        score1 = None
        while score1 is None and calculate1 == 1:
            try:
                score1 = self.dpmatrix[(column - 1, row - 1)].get_value() + \
                  self.sub_matrix.get_sub_value(self.seq1[column - 1], 
                    self.seq2[row - 1])
            # if we get either a key or index error the previous
            # cell item is not filled and needs to be filled
            except KeyError, IndexError:
                diag_cell = self.fill_cell(column - 1, row - 1)
                
                self.dpmatrix[(column - 1, row - 1)] = diag_cell
                
        # Score 2
        score2 = None
        while score2 is None and calculate2 == 1:
            try:
                score2 = self.dpmatrix[(column, row - 1)].get_value() + \
                  self.gap_penalty
            # if we get either a key or index error the previous cell
            # item is not filled and needs to be filled
            except KeyError, IndexError:
                upper_cell = self.fill_cell(column, row - 1)
                
                self.dpmatrix[(column, row - 1)] = upper_cell
                
        # Score 3
        score3 = None
        while score3 is None and calculate3 == 1:
            try:
                score3 = self.dpmatrix[(column - 1, row)].get_value() + \
                  self.gap_penalty
            # if we get either a key or index error the previous cell
            # item is not filled and needs to be filled
            except KeyError, IndexError:
                right_cell = self.fill_cell(column - 1, row)
                
                self.dpmatrix[(column - 1, row)] = right_cell
                
        # determine which score is highest of the three
        # first calculate the max of all the scores which are not none.
        max_list = []
        for n in [score1, score2, score3]:
            if n:
                max_list.append(n)
                
        max_score = max(max_list)
        
        # now set the appropriate value and previous item depending on 
        # which score is largest, and return the created cell
        if max_score == score1:
            new_cell = DPMatrixCell(column, row, self.seq1[column - 1], 
              self.seq2[row - 1])
            new_cell.set_value(max_score)
            new_cell.set_parent(self.dpmatrix[(column - 1, row - 1)])
            return new_cell
            
        elif max_score == score2:
            # two cases to create the new cell. 
            # 1. special case (right column)
            if score1 is None and score3 is None:
                new_cell = DPMatrixCell(column, row, "", self.seq2[row - 1])
            # 2. normal case (elsewhere in the matrix)
            else:
                new_cell = DPMatrixCell(column, row, self.seq1[column - 1], 
                  self.seq2[row - 1])
                
            new_cell.set_value(max_score)
            new_cell.set_parent(self.dpmatrix[(column, row - 1)])
            
            return new_cell
            
        elif max_score == score3:
            # two cases to create the new cell.
            # 1. speical case (top row)
            if score1 is None and score3 is None:
                new_cell = DPMatrixCell(column, row, self.seq1[column - 1], "")
            # 2. normal case (elsewhere in the matrix)
            else:
                new_cell = DPMatrixCell(column, row, self.seq1[column - 1], 
                  self.seq2[row - 1])
                
            new_cell.set_value(max_score)
            new_cell.set_parent(self.dpmatrix[(column - 1, row)])
            
            return new_cell
            
    def align(self, sequence1, sequence2):
        self.fill_matrix(sequence1, sequence2)
        
    def get_optimal_alignment(self):
        """Follow the traceback to get the optimal alignment."""
        # intialize the two sequences which will return the alignment
        align_seq1 = MutableSeq(array.array("c"), 
          Alphabet.Gapped(IUPAC.protein, GAP_CHAR))
        align_seq2 = MutableSeq(array.array("c"), 
          Alphabet.Gapped(IUPAC.protein, GAP_CHAR))
          
        # take care of the initial case with the bottom corner matrix
        # item
        current_cell = self.dpmatrix[(len(self.seq1), len(self.seq2))]
        align_seq1.append(current_cell.seq1item)
        align_seq2.append(current_cell.seq2item)
        
        next_cell = current_cell.get_parent()
        current_cell = next_cell
        next_cell = current_cell.get_parent()
        
        # keeping adding sequence until we reach (0, 0)
        while next_cell:
            # add the new sequence--three cases:
            # 1. Move up diaganolly, add a new seq1 and seq2 to the 
            # aligned sequences
            if ((next_cell.col_pos == current_cell.col_pos - 1) and
              (next_cell.row_pos == current_cell.row_pos - 1)):
                # print "case 1 -> seq1 %s, seq2 %s" % (
                # current_cell.seq1item, current_cell.seq2item)
                align_seq1.append(current_cell.seq1item)
                align_seq2.append(current_cell.seq2item)
            # 2. Move upwards, add a new seq2 and a gap in seq1
            elif ((next_cell.col_pos  == current_cell.col_pos) and
              (next_cell.row_pos == current_cell.row_pos - 1)):
                #print "case 2 -> seq2 %s" % current_cell.seq2item
                align_seq1.append(GAP_CHAR)
                align_seq2.append(current_cell.seq2item)
            # 3. Move to the right, add a new seq1 and a gap in seq2
            elif ((next_cell.col_pos == current_cell.col_pos - 1) and
              (next_cell.row_pos == current_cell.row_pos)):
                #print "case 3 -> seq1 % s" % current_cell.seq1item
                align_seq1.append(current_cell.seq1item)
                align_seq2.append(GAP_CHAR)
            
            # now move on to the next sequence
            current_cell = next_cell
            next_cell = current_cell.get_parent()
        
        # reverse the returned alignments since we are reading them in
        # backwards
        align_seq1.reverse()
        align_seq2.reverse()
        return align_seq1.toseq(), align_seq2.toseq()
            
        
    def get_optimal_score(self):
        # by definition, the bottom corner matrix item is the
        # optimal score for the alignment
        return self.dpmatrix[(len(self.seq1), len(self.seq2))].get_value()


        
if __name__ == '__main__':

    # create a substitution matrix
    sub_matrix = SubstitutionMatrix('pam120')
    
    # set up for alignment
    aligner = Needle(sub_matrix)
    print "Testing a simple alignment..."
    seq1 = "ATAGTAGGTACATCTTTAAGTATATTAATTCGTGCAGAATTAGGTCATCCTGGAACATTTATTGGAGATGATCAAATTTATAATGTAATTGTTACAGCACATGCTTTTATTATAATTTTTTTTATAGTAATACCTATTTTAATTGGTGGATTTGGAAATTGACTATTACCTTTAATATTAGGTGCCCCTGATATAGCTTTTCCTCGAATAAATAATATAAGTTTTTGATTATTACCCCCTTCTCTTACTCTTCTTCTTTCAAGTAGAATAGTAGAAAACGGAGCAGGTACAGGATGAACAGTTTATCCTCCACTTTCTTCAAGAATTGCTCATAGAGGTGCTTCTGTTGATTTAGCAATTTTTTCTCTTCATTTAGCAGGAATTTCTTCAATTTTAGGTTCAGTAAATTTTATTACTACAGCTATTAATATACGATCTAATGGTATTACGCTAGACCGAATACCTTTATTTGTATGATCAGTAGTAATTACAACTATCCTTCTTTTATTATCTTTACCTGTATTAGCTGGAGCAATTACAATATTATTAACTGATCGAAATTTAAATACATCTTTTTTTGACCCTGCAGGAGGAGGAGATCCAATTTTATACCAACATTTATTT"
    seq2 = "ATAGTAGGTACATCTTTAAGTATATTAATTCGTGCAGAATTAGGTCATCCTGGAACATTTATTGGAGATGATCAAATTTATAATGTAATTGTTACAGCACATGCTTTTATTATAATTTTTTTTATAGTAATACCTATTTTAATTGGTGGATTTGGAAATTGACTATTACCTTTAATATTAGGTGCCCCTGATATAGCTTTTCCTCGAATAAATAATATAAGTTTTTGATTATTACCCCCTTCTCTTACTCTTCTTCTTTCAAGTAGAATAGTAGAAAACGGAGCAGGTACAGGATGAACAGTTTATCCTCCACTTTCTTCAAGAATTGCTCATAGAGGTGCTTCTGTTGATTTAGCAATTTTTTCTCTTCATTTAGCAGGAATTTCTTCAATTTTAGGTTCAGTAAATTTTATTACTACAGCTATTAATATACGATCTAATGGTATTACGCTAGACCGAATACCTTTATTTGTATGATCAGTAGTAATTACAACTATCCTTCTTTTATTATCTTTACCTGTATTAGCTGGAGCAATTACAATATTATTAACTGATCGAAATTTAAATACATCTTTTTTTGACCCTGCAGGAGGAGGAGATCCAATTTTATACCAACATTTATTT"



    aligner.align(seq1, seq2)

    align1, align2 = aligner.get_optimal_alignment()
    score = aligner.get_optimal_score()
    
    print "Alignment Score:", score
    print align1.data
    print align2.data

#                     sub_matrix = SubstitutionMatrix('pam120')
#                     aligner = Needle(sub_matrix)                    
#                     aligner.align(seq1, seq2)
#                     aligned1, aligned2 = aligner.get_optimal_alignment()
#                     score = aligner.get_optimal_score()


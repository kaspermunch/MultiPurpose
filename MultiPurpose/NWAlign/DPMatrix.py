"""DPMatrix.py

This class models a dynamic progamming matrix for use in sequence 
alignment models. The design of this class is based on the description of 
dynamic programming matrices in:

    Durbin et al. 1998. Biological Sequence Analysis. Cambridge University 
    Press.
    
Classes:
* DPMatrix
* DPMatrixCell
"""

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

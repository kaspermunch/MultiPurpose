"""SubstitutionMatrix.py

Represent substitution matrices for alignments.

Substitution matrices are not actually represented as a matrix, but rather 
as a dictionary, with the key being a tuple, and the value being the 
substition value.
So the entry for an A to R substitution in a BLOSUM50 would be represented 
as {('A', 'R'): -2}. Since the matrix is symmetric, this lets us have less 
objects in memory then if it was defined in a full matrix.
"""

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

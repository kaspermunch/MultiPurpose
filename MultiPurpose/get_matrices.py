"""get_matrices.py

Fetch substitution matrices and write them as python dictionary objects 
(sparse matrices).
We fetch these remotely from a *very* useful page with info about a ton of 
matrices -> http://www.embl-heidelberg.de/~vogt/matrices/mlist1.html.
"""
# standard modules
import sys
import urllib
import string

##### constants
BASE_URL = "http://www.embl-heidelberg.de/~vogt/matrices/"
MATRIX_URLS = ['benner6.cmp', 'benner22.cmp', 'benner74.cmp',
'blosum100.cmp','blosum30.cmp', 'blosum35.cmp', 'blosum40.cmp',
'blosum45.cmp', 'blosum50.cmp', 'blosum55.cmp', 'blosum60.cmp',
'blosum62.cmp','blosum65.cmp','blosum70.cmp', 'blosum75.cmp',
'blosum80.cmp','blosum85.cmp','blosum90.cmp','blosum95.cmp',
'feng.cmp','fitch.cmp','genetic.cmp','gonnet.cmp','grant.cmp',
'ident.cmp', 'johnson.cmp','levin.cmp', 'mclach.cmp','miyata.cmp',
'nwsgappep.cmp', 'pam120.cmp','pam180.cmp','pam250.cmp', 'pam30.cmp',
'pam300.cmp', 'pam60.cmp','pam90.cmp', 'rao.cmp', 'risler.cmp',
'str.cmp']

def main():
    for url in MATRIX_URLS:
        sub_matrix = {}
        in_table = 0
        protein_top = []
        matrix_info = urllib.urlopen(BASE_URL + url)
        for line in matrix_info.readlines():
            strip_line = string.strip(line)
            # find the first line of the table
            if ((len(strip_line) != 0) and (strip_line[0] in string.uppercase) 
              and (strip_line[1] in string.whitespace)):
                protein_top = string.split(strip_line)
                # strip off the last item, whieh is "..."
                protein_top = protein_top[:-1]
                # set us as being in the table so we can start reading
                # the table values
                in_table = 1

            elif in_table == 1:
                line_info = string.split(line)
                # determine the start position in the top list of
                # proteins to  match up with the current line.
                # For instance, if we have
                # protein_top : A  B  C  D
                # line_info:       1  2  3 B
                # We want to start with the second item in the
                # protein_top list
                start_pos = (len(protein_top)) - (len(line_info) - 1)
                # now cycle through each item in the current line and
                # add the matrix info to our substitution matrix
                for n in range(len(line_info) - 1):
                    # find the key (the two proteins we are reading the
                    # substitution value for)
                    item_key = (protein_top[n + start_pos], 
                      line_info[len(line_info) - 1])
                    # convert the item into a string or float,
                    # depending on its type
                    try:
                        sub_value = string.atoi(line_info[n])
                    # if we get a ValueError then we've got a float
                    except ValueError:
                        sub_value = string.atof(line_info[n])
        
                    sub_matrix[item_key] = sub_value    

        matrix_info.close()
        # now that we've got the matrix info, we want to write it out
        # so we can save it in a file.
        print "#", BASE_URL + url
        print "%s = {" % url[:-4]
        counter = 0
        line_to_print = ""
        num_items = len(sub_matrix.keys())
        for n in sub_matrix.keys():
            line_to_print = line_to_print + str(n) + " : " + \
              str(sub_matrix[n]) + ", "
            counter = counter + 1
            
            # if we've reached the end, print everything out
            if counter == num_items:
                # strip off the trailing comma
                line_to_print = line_to_print[:-2]
                print line_to_print
            # every fourth item, we want to print the info
            elif (counter % 4) == 0:
                print line_to_print
                line_to_print = ""

            
        print "}"
                
        
    
if __name__ == '__main__':
    sys.exit(main())

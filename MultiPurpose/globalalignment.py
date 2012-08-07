#!/usr/bin/python

# Time-stamp: <2006-09-01 15:58:55 kasper>

# Determine the score of the optimal global alignment of two sequences
def globalAlignmentScore(s1, s2):

	# Initialize variables
	MATCH = 5
	MISMATCH = -4
	GAP = -6
	NUM_ROWS = len(s2)+1
	NUM_COLS = len(s1)+1

	# Create table and fill it with zeros
	table = createTable(NUM_ROWS, NUM_COLS)
		
	col = 0
	while (col < NUM_COLS):  # Fill in first row of table
		table[0][col] = col * GAP
		col = col + 1

	row = 0
	while (row < NUM_ROWS):  # Fill in first column of table
		table[row][0] = row * GAP
		row = row + 1
		
	row = 1
	while (row < NUM_ROWS):  # Fill in each row of table
		col = 1
		while (col < NUM_COLS):  # Fill in column for each row
			left = table[row][col-1] + GAP
			above = table[row-1][col] + GAP
			charAlign = 0
			if (s1[col-1] == s2[row-1]):
				charAlign = MATCH
			else:
				charAlign = MISMATCH
			diagonal = table[row-1][col-1] + charAlign
			table[row][col] = max3(left, above, diagonal)
			col = col + 1
		row = row + 1

	# Print out table (only useful for small tables)
	printTable(table)

	return table[NUM_ROWS-1][NUM_COLS-1]

# Create a table with the given number of rows and columns
def createTable(NUM_ROWS, NUM_COLS):
	table = []
	row = 0
	while (row < NUM_ROWS):
		table.append([])
		col = 0
		while (col < NUM_COLS):
			table[row].append(0)
			col = col + 1
		row = row + 1
	return table

# Print out a table (only useful for small tables)
def printTable(table):
	print ("Below is the filled in table")
	s = ""
	for i in range(0, len(table)):
		for j in range(0, len(table[0])):
			s = s + str(table[i][j]) + "\t"
		s = s + "\n"
	print s	

# Determine the maximum of 3 numbers
def max3(a, b, c):
	if ((a >= b) and (a >= c)): return a
	if ((b >= a) and (b >= c)): return b
	return c



########################################################
### End of functions ###################################
########################################################



# Use function to calculate global alignment score of two sequences
s1 = "AGCGTTA"
s2 = "AGCTTA"
#s2 = "ACGTGA"
optimalScore = globalAlignmentScore(s1, s2)
print s1
print s2
print "Global alignment score: ", optimalScore



import sys
import re
import string
import copy
from optparse import OptionParser
import bisect

import time



class Bisect(object):
    """
    Works on sorted non-overlapping index based intervals
    """
    def __init__(self, intervals, binIdx=[0,1], inverse=False):
        
        startIdx, endIdx = binIdx        
        self.intervals = intervals
        self.inverse = inverse
        self.starts = map(float, [x[startIdx] for x in intervals])
        self.ends = map(float, [x[endIdx] for x in intervals])

        assert self.starts == sorted(self.starts)
        assert self.ends == sorted(self.ends)

    def index_of_leftmost_overlapping(self, start):
        'Find index of leftmost end greater than start'   
        i = bisect.bisect_right(self.ends, start)
        if i != len(self.ends):
            assert self.ends[i] > start
            return i

    def index_of_rightmost_overlapping(self, end):
        'Find index of rightmost start less than end'
        i = bisect.bisect_left(self.starts, end)
        if i:
            assert self.starts[i-1] < end
            return i-1

    def overlappingIntervals(self, start, end):

        assert start <= end
        
        startIdx = self.index_of_leftmost_overlapping(start)
        endIdx = self.index_of_rightmost_overlapping(end)
        if None in [startIdx, endIdx] or endIdx < startIdx:
            return []
        else:
            if self.inverse:
                return [self.intervals[i][:] for i in range(len(self.intervals)) if i not in range(startIdx, endIdx+1)]
            else:
                return [self.intervals[i][:] for i in range(startIdx, endIdx+1)]



# FIXME: We don't need to find overlapping regions in collapse if we know that the set is
# already an ovoerlapping cluster. Then we can just return set.seq, set.start and set.end.

# FIXME: should take another look at cluster(), suspect redundant code there.

class Data(object):

    def __init__(self, lst, coordinateIdxs, linenr=None):
        self.linenr = linenr
        self.lst = lst
        self.nameIdx, self.startIdx, self.endIdx = coordinateIdxs
        self.name = self.lst[self.nameIdx]
        self.start = float(self.lst[self.startIdx])
        self.end = float(self.lst[self.endIdx])

        self.feature = '.'
        self.strand = '.'
        self.source = '.'
        self.score = '.'
        self.frame = '.'
        self.attrib = ''
        self.attribstr = ''

    def __iter__(self):
        self.i = 0
        return self

    def next(self):
        if self.i == len(self.lst):
            raise StopIteration
        else:
            self.i += 1
            return self.lst[self.i-1]

    def __cmp__(self, other):
        return cmp(self.name, other.name) \
               or cmp(self.start, other.start) \
               or cmp(self.end, other.end)

    def __str__(self):
        self.lst[self.nameIdx] = self.name
        self.lst[self.startIdx] = str(self.start)
        self.lst[self.endIdx] = str(self.end)
        return "\t".join(self.lst)

    def __repr__(self): 
        return str(self)

    def fields(self):
        self.lst[self.nameIdx] = self.name
        self.lst[self.startIdx] = str(self.start)
        self.lst[self.endIdx] = str(self.end)
        return self.lst
        
    def onlyCoordinates(self):
        self.feature = '.'
        self.strand = '.'
        self.source = '.'
        self.score = '.'
        self.frame = '.'

    def overlaps(self, other):
        if self.name == other.name \
               and (self.start >= other.start and self.start < other.end \
                    or self.end > other.start and self.end <= other.end \
                    or self.start < other.start and self.end > other.start):
            return 1
        else:
            return 0

    def nestedIn(self, other):
        if self.name == other.name and self.start >= other.start and self.end <= other.end:
            return True
        else:
            return False

    def endToEnd(self, other):
        if self.name == other.name and (self.start == other.end or self.end == other.start):
            return 1
        else:
            return 0

    def before(self,other):
        if self.name < other.name or self.name == other.name and self.end <= other.start:
            return 1
        else:
            return 0        

    def after(self,other):
        if self.name > other.name or self.name == other.name and self.start >= other.end:
            return 1
        else:
            return 0



class DataSet(object):
    "Class for holding a set of Data lines"

    def __init__(self, lst=[]):
        if len(lst):
            self.lst = lst
            self.name = lst[0].name
            self.start = lst[0].start # because the list is sorted
            ends = [item.end for item in self.lst]
            self.end = max(ends) # this could be lst[len(lst)].end
            self.strand = lst[0].strand
        else:
            self.lst = lst
            self.name = False
            self.start = False
            self.end = False
            self.strand = False
#         for i, j in zip(self.lst, sorted(self.lst)):
#             assert i.start == j.start, str(i) + " " +  str(j)

    def __str__(self):
        s = "%s %s %s\n" % (self.name, self.start, self.end)
        s += '\n'.join(map(str, self.lst))
        s += '\n'
        return s
       

    def __len__(self):
        return len(self.lst)

    def __getitem__(self, item):
        return self.lst[item]

    def __delitem__(self, item):
        del self.lst[item]

    def __iter__(self):
        self.index = 0
        return self

    def next(self):
        if self.index < len(self.lst):
            self.index += 1
            return self.lst[self.index-1]
        else:
            raise StopIteration

    def clear(self):
        self.lst = []
        self.name = False
        self.start = False
        self.end = False
        self.strand = False

    def add(self, other):
        self.lst.append(other)
        self.lst.sort()
        if not self.name:
            self.name = other.name
        else:
            assert self.name == other.name, "Tried to add a Data from a diffrent sequence " + other.name + " ne " + self.name
        if self.start:
            self.start = min(self.start, other.start)
        else:
            self.start = other.start
        if self.end:
            self.end = max(self.end, other.end)
        else:
            self.end = other.end
    def len(self):
        return len(self.lst)


class DataIterator(object):
    "Iterator for iterating over Data entries"

    def __init__(self, iterator, coordinateIdxs):

        self.iterator = iterator
        #self.fp = fp
        
        self.linenr = 0
        self.coordinateIdxs = coordinateIdxs

    def __iter__(self):
        return self

    def next(self):
        #line = self.fp.readline()
        line = self.iterator.next()

        # hack to allow lists as input:
        if isinstance(line, list) or isinstance(line, tuple):
            line = "\t".join(map(str, line))

        line = line.strip()
        if line:
            lst = line.split("\t")
            #coord = Data(lst);
            self.linenr += 1
            data = Data(lst, self.coordinateIdxs, self.linenr)        
            if data:
                return data
            else:
                raise StopIteration
        else:
            raise StopIteration        


class BufferedList(object):
    "Buffered Iterator for iterating over Data data list"

    def __init__(self, lst):
        self.lst = lst
        for i, j in zip(self.lst, sorted(self.lst, reverse=True)):
            assert i.start == j.start

    def __iter__(self):
        return self

    def next(self):
        data = False
        if len(self.lst):
            data = self.lst.pop()
            return data
        else:
            raise StopIteration

    def putback(self, data):
        self.lst.append(data)

#         for i, j in zip(self.lst, sorted(self.lst, reverse=True)):
#             if i.start != j.start:
#                 print "### error"
#                 print i
#                 print j
#                 assert 0


class BufferedIterator(object):
    "Buffered List for iterating over Data data input stream"
    def __init__(self, iterator):
        self.buffer = []
        self.iterator = iterator
    def __iter__(self):
        return self
    def next(self):
        data = False
        if len(self.buffer):
            data = self.buffer.pop()
        else:
            data = self.iterator.next()
        if data:
            return data
        else:
            raise StopIteration
    def putback(self, data):
        self.buffer.append(data)


class Mapping(object):

    def __init__(self, file1, file2, idx1=[0,1,2], idx2=[0,1,2], presorted=False, nonoverlapping=False):

        self.presorted = presorted
        self.nonoverlapping = nonoverlapping
        self.coordinateIdxs1 = idx1
        self.coordinateIdxs2 = idx2

        self.file1data = dict()
        self.file2data = dict()

        try:
            instream1 = open(file1, 'r')
            self.iter1 = iter(instream1)
        except:
            file1.seek(0)
            self.iter1 = iter(file1)
        try:
            instream2 = open(file2, 'r')
            self.iter2 = iter(instream2)
        except:
            file2.seek(0)
            self.iter2 = iter(file2)

    def __iter__(self):

        # Iterators to read in data lines:
        dataIterator1 = DataIterator(self.iter1, self.coordinateIdxs1)
        dataIterator2 = DataIterator(self.iter2, self.coordinateIdxs2)

        # Sort and make iterators:
        if not self.presorted:
            lst1 = []
            for data in dataIterator1:
                lst1.append(data)
                self.file1data[data.linenr] = data
            lst1.sort()
            lst1.reverse()
            lst2 = []
            # Make buffered iterator:
            self.file1 = BufferedList(lst1)
            for data in dataIterator2:
                lst2.append(data)
                self.file2data[data.linenr] = data
            lst2.sort()
            lst2.reverse()

            # Make buffered iterator:
            self.file2 = BufferedList(lst2)
        else:
            # Make buffered iterators:
            self.file1 = BufferedIterator(dataIterator1)
            self.file2 = BufferedIterator(dataIterator2)

        self.data = dict()

        # Initialise Sets to hold overlapping data lines within each file:
        cluster1 = self.nextCluster(self.file1)
        cluster2 = []
        while cluster1:

            # Get set of overlapping datas from file2:
            cluster2 = self.nextOverlapCluster(self.file2, cluster1)
            
            orig2 = copy.deepcopy(cluster2)

            if len(cluster2):

                overlap, nonoverlap = self.getMap(cluster1, cluster2)
                if self.nonoverlapping:
                    if nonoverlap:                        
                        for n in nonoverlap:
                            yield n
                else:
                    if overlap:
                        for o in overlap:
                            self.data[o[0]] = o[1:]
                            yield o
                        
                #self.printSet(overlap, sys.stdout)
                # We might need the data from cluster2 that extends
                # past the end of cluster1. So the collapsed regions
                # from file2 that we might need to look at again are
                # put back:
                for data in orig2[::-1]:
                    self.file2.putback(data)

            elif self.nonoverlapping:
                for l in cluster1:
                    yield l

            # Reinitialise:
            cluster1 = self.nextCluster(self.file1)
            cluster2 = []
            

    def printSet(self, set, out):
        "Prints a set of coord opjects."
        for data in set:
            out.write(str(data)+ "\n")


    def getMap(self, set1, set2):
        "Overlap of overlapping datas in two sorted sets of objects."

        overlap = []
        nonoverlap = []
        
        # where we start looking in set2 (it is the index where all
        # previous datas have ends before the start of the next set1 data):
        set2Start = 0

        # for data in set1:
        for i in range(len(set1)):
            # Indicates whether we have found the set2Start value:
            set2StartFound = False
            matches = list()
            matchIndexes = list()
            starts = list()
            ends = list()
            for j in range(set2Start, len(set2)):

                # Keep track of where to start in set2 next time:
                if not set2StartFound:
                    if i < len(set1) - 1 and set2[j].before(set1[i+1]):
                        set2Start = j
                    else:
                        set2StartFound = True                        

                # Check for overlap:
                if set1[i].overlaps(set2[j]):
                    # make an overlap data:
                    matches.append(set2[j].linenr)
                    matchIndexes.append(j)
                    starts.append(max(set1[i].start, set2[j].start))
                    ends.append(min(set1[i].end, set2[j].end))
# {{{                    

#                     # Keep track of where to start in set2 next time:
#                     if not set2StartFound:
#                         if i < len(set1) - 1 and set2[j].end < set1[i+1].start:
# 
# #set 1: ################
# 
# #set 2: ###################################  first one overlaps the one in set 1
#                            ######            the next one does not and j is incremented without set2start being incremented...
#                                       ################# ... so then it craches here....
# 
# 
# 
#                             #assert j == 0 or j == set2Start + 1, "Incremented by more than one."
#                             if not (j == 0 or j == set2Start + 1):
#                                 print set1
#                                 print set2
#                                 print j, set2Start + 1
#                                 sys.exit()
# 
#     #                        assert j == set2Start + 1, "Incremented by more than one."
#                             set2Start = j
#                         else:
#                             set2StartFound = 1                        

# }}}

                elif set1[i].before(set2[j]):
                    # Because the sets are sorted we are not going to find an overlap later on:
                    break

            if matches:
                #overlap.append([set1[i].linenr, matches, starts, ends])
                overlap.append([set1[i], [set2[x] for x in matchIndexes], starts, ends])

        return overlap, [set1[x] for x in [x for x in range(len(set1)) if x not in matchIndexes]]


    def nextCluster(self, buffIter):
        "Gets the next cluster of overlapping datas from the file1 stream."

        cluster = []

        for data in buffIter:

            # Get a cluster (often a singleton) of overlapping file1 datas:
            if not len(cluster):
                # If cluster1 is empty we add the line to it:
                cluster.append(data)
                continue
            elif data.overlaps(cluster[-1]):
                # If cluster1 is not empty and data overlaps the previous data:
                cluster.append(data)
                continue
            else:
                # If we don't have overlap we buffer the data for the next loop and go on:
                buffIter.putback(data)
                break

        cluster = DataSet(cluster)
        return cluster
        

    def nextOverlapCluster(self, buffIter, otherCluster):
        """Gets the next cluster of overlapping datas from the file2 stream
        that overlaps a cluster from the file1 stream."""

        cluster = []

        # Get a line from buffIter:
        for data in buffIter:

            if data.before(otherCluster):
                # Skip the data line if it is before the first line in set1:
                continue
            elif data.after(otherCluster):
                # Buffer and terminate loop if the data line is after the last line in set1:
                buffIter.putback(data)
                break
            else:
                # data overlaps set1 so we put it in set2:
                cluster.append(data)
                continue

        cluster = DataSet(cluster)

        return cluster


    def collapse(self, set):
        "Collapses overlapping datas in a sorted set of objects."
        i = 0
        j = 1
        collapsedSet = []

        new = copy.deepcopy(set)

        while j < len(new):
            while j < len(new) and new[i].overlaps(new[j]):
                new[i].end = max(new[i].end, new[j].end)
                j += 1
            new[i].onlyCoordinates()
            collapsedSet.append(new[i])
            i = j
            # if we are now at the last data we add it:
            if i == len(new) - 1:
                new[i].onlyCoordinates()
                collapsedSet.append(new[i])
                break
            j += 1
        if len(new) == 1:
            new[0].onlyCoordinates()
            collapsedSet.append(new[0])
        return DataSet(collapsedSet)



class Overlap(object):

    def __init__(self, file1, file2, coordinateIdxs1=[0,1,2], coordinateIdxs2=[0,1,2], presorted=False, nonoverlapping=False):

        self.presorted = presorted
        self.nonoverlapping = nonoverlapping
        self.coordinateIdxs1 = coordinateIdxs1
        self.coordinateIdxs2 = coordinateIdxs2

#         self.file1 = file1
#         self.file2 = file2
        self.file1data = dict()
        self.file2data = dict()

        try:
            instream1 = open(file1, 'r')
            iter1 = iter(instream1)
        except:
            iter1 = iter(file1)
        try:
            instream2 = open(file2, 'r')
            iter2 = iter(instream2)
        except:
            iter2 = iter(file2)

        # Iterators to read in data lines:
        dataIterator1 = DataIterator(iter1, coordinateIdxs1)
        dataIterator2 = DataIterator(iter2, coordinateIdxs2)

        # Sort and make iterators:
        if not presorted:
            lst1 = []
            for data in dataIterator1:
                lst1.append(data)
                self.file1data[data.linenr] = data
            lst1.sort()
            lst1.reverse()
            lst2 = []
            # Make buffered iterator:
            file1 = BufferedList(lst1)
            for data in dataIterator2:
                lst2.append(data)
                self.file2data[data.linenr] = data
            lst2.sort()
            lst2.reverse()

            # Make buffered iterator:
            file2 = BufferedList(lst2)
        else:
            # Make buffered iterators:
            file1 = BufferedIterator(dataIterator1)
            file2 = BufferedIterator(dataIterator2)

        self.data = dict()

        # Initialise Sets to hold overlapping data lines within each file:
        cluster1 = self.nextCluster(file1)
        cluster2 = []
        while cluster1:

            # Get set of overlapping datas from file2:
            cluster2 = self.nextOverlapCluster(file2, cluster1)
            
            orig2 = copy.deepcopy(cluster2)

            if len(cluster2):
                
                overlap = self.getMap(cluster1, cluster2)
                if overlap:
                    for o in overlap:
                        self.data[o[0]] = o[1:]
                    
                #self.printSet(overlap, sys.stdout)
                # We might need the data from cluster2 that extends
                # past the end of cluster1. So the collapsed regions
                # from file2 that we might need to look at again are
                # put back:
                for data in orig2[::-1]:
                    file2.putback(data)

            # Reinitialise:
            cluster1 = self.nextCluster(file1)
            cluster2 = []

        self.nonoverlappinglines = sorted(set(range(1,dataIterator1.linenr+1)).difference(set(self.data.keys())))


    def nonOverlappingLines(self):
        return [self.file1data[int(file1Line)] for file1Line in self.nonoverlappinglines]


    def __iter__(self):
#         return self
# 
#     def next(self):

        if self.presorted:
            print "You can't use the object as an iterator if input data is spedified as presorted"

#         if self.nonoverlapping:
#             for file1Line in self.nonoverlappinglines:
#                 yield self.file1data[int(file1Line)]
#         else:
        for file1Line in sorted(self.data.keys()):
            file2Lines = tuple([self.file2data[x] for x in self.data[file1Line][0]])
            yield OverlapInstance(file1Line,
                                  self.file1data[int(file1Line)],
                                  zip(*self.data[file1Line]+[file2Lines]), self.coordinateIdxs2)

#                 for file2Line, start, end in zip(*self.data[file1Line]):
#                     yield OverlapInstance(file1Line, start, end, self.file1data[int(file1Line)]), OverlapInstance(file2Line, start, end, self.file2data[int(file2Line)])
            

    def printSet(self, set, out):
        "Prints a set of coord opjects."
        for data in set:
            out.write(str(data)+ "\n")


    def getMap(self, set1, set2):
        "Overlap of overlapping datas in two sorted sets of objects."

        overlap = []

        # where we start looking in set2 (it is the index where all
        # previous datas have ends before the start of the next set1 data):
        set2Start = 0

        # for data in set1:
        for i in range(len(set1)):
            # Indicates whether we have found the set2Start value:
            set2StartFound = False
            matches = list()
            starts = list()
            ends = list()
            for j in range(set2Start, len(set2)):

                # Keep track of where to start in set2 next time:
                if not set2StartFound:
                    if i < len(set1) - 1 and set2[j].before(set1[i+1]):
                        set2Start = j
                    else:
                        set2StartFound = True                        

                # Check for overlap:
                if set1[i].overlaps(set2[j]):
                    # make an overlap data:
                    matches.append(set2[j].linenr)
                    starts.append(max(set1[i].start, set2[j].start))
                    ends.append(min(set1[i].end, set2[j].end))
# {{{                    
#                     # Keep track of where to start in set2 next time:
#                     if not set2StartFound:
#                         if i < len(set1) - 1 and set2[j].end < set1[i+1].start:
# 
# #set 1: ################
# 
# #set 2: ###################################  first one overlaps the one in set 1
#                            ######            the next one does not and j is incremented without set2start being incremented...
#                                       ################# ... so then it craches here....
# 
# 
# 
#                             #assert j == 0 or j == set2Start + 1, "Incremented by more than one."
#                             if not (j == 0 or j == set2Start + 1):
#                                 print set1
#                                 print set2
#                                 print j, set2Start + 1
#                                 sys.exit()
# 
#     #                        assert j == set2Start + 1, "Incremented by more than one."
#                             set2Start = j
#                         else:
#                             set2StartFound = 1                        
# }}}
                elif set1[i].before(set2[j]):
                    # Because the sets are sorted we are not going to find an overlap later on:
                    break

            if matches:
                overlap.append([set1[i].linenr, matches, starts, ends])

        return overlap


    def nextCluster(self, buffIter):
        "Gets the next cluster of overlapping datas from the file1 stream."

        cluster = []

        for data in buffIter:

            # Get a cluster (often a singleton) of overlapping file1 datas:
            if not len(cluster):
                # If cluster1 is empty we add the line to it:
                cluster.append(data)
                continue
            elif data.overlaps(cluster[-1]):
                # If cluster1 is not empty and data overlaps the previous data:
                cluster.append(data)
                continue
            else:
                # If we don't have overlap we buffer the data for the next loop and go on:
                buffIter.putback(data)
                break

        cluster = DataSet(cluster)
        return cluster
        

    def nextOverlapCluster(self, buffIter, otherCluster):
        """Gets the next cluster of overlapping datas from the file2 stream
        that overlaps a cluster from the file1 stream."""

        cluster = []

        # Get a line from buffIter:
        for data in buffIter:

            if data.before(otherCluster):
                # Skip the data line if it is before the first line in set1:
                continue
            elif data.after(otherCluster):
                # Buffer and terminate loop if the data line is after the last line in set1:
                buffIter.putback(data)
                break
            else:
                # data overlaps set1 so we put it in set2:
                cluster.append(data)
                continue

        cluster = DataSet(cluster)

        return cluster


    def collapse(self, set):
        "Collapses overlapping datas in a sorted set of objects."
        i = 0
        j = 1
        collapsedSet = []

        new = copy.deepcopy(set)

        while j < len(new):
            while j < len(new) and new[i].overlaps(new[j]):
                new[i].end = max(new[i].end, new[j].end)
                j += 1
            new[i].onlyCoordinates()
            collapsedSet.append(new[i])
            i = j
            # if we are now at the last data we add it:
            if i == len(new) - 1:
                new[i].onlyCoordinates()
                collapsedSet.append(new[i])
                break
            j += 1
        if len(new) == 1:
            new[0].onlyCoordinates()
            collapsedSet.append(new[0])
        return DataSet(collapsedSet)


class OverlapInstance(object):

    def __init__(self, linenr, line, lst, coordinateIdxs2):
        self.linenr = linenr
        self.line = line
        self.overlapCoordinates = lst
        self.nameIdx, self.startIdx, self.endIdx = coordinateIdxs2

    def __iter__(self):
        self.iterIdx = 0
        return self

    def next(self):
        if self.iterIdx < len(self.overlapCoordinates):
            i = self.iterIdx
            self.iterIdx += 1
            return self.overlapCoordinates[i]
#             linenr, start, end, line = self.overlapCoordinates[i]
#             return linenr, start, end, line.split()
        else:
            raise StopIteration

    def nrOverlaps(self):
        return len(self.overlapCoordinates)

#     def coordinates(self):
#         l = self.line.split()
#         return (l[self.nameIdx], l[self.startIdx], l[self.endIdx])

    
if __name__ == "__main__":
    file1, file2 = sys.argv[1], sys.argv[2]

    for overlapInstance in Overlap(file1, file2):
        print "line %d in file1: \"%s\"" % (overlapInstance.linenr, overlapInstance.line)
        for linenr, start, end, line in overlapInstance:
            print "overlaps line %d in file2: \"%s\" from %d to %d" % (linenr, line, start, end)
        print

    print "lines in file1 not overlapping anything in file2"
    for line in Overlap(file1, file2).nonOverlappingLines():
        print line
#     for line in Overlap(file1, file2, nonoverlapping=True):
#         print line



    for file1line, file2lines, starts, ends in Mapping(file1, file2, coordinateIdxs1=[0,3,4], coordinateIdxs2=[1,2,2], presorted=True):
        print file1line
        for l in file2lines:
            print '\t', l
        print

    for file1line in Mapping(file1, file2, coordinateIdxs1=[0,3,4], coordinateIdxs2=[0,3,4], presorted=True, nonoverlapping=True):
        print file1line

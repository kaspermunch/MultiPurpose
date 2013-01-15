
import os, sys, math, re, csv
import Overlap, StringIO

def midpointCoordinates(inp, outp, colIdx, header=False):
    """
    Replace a column of start coordinates with midpoint coordinates
    """
    inputFile = open(inp)
    outputFile = open(outp, 'w')
    if header:
        outputFile.write(inputFile.readline())
    prevPos = None
    for l in inputFile.xreadlines():
        lst = l.split()
        pos = int(lst[colIdx])
        if prevPos is not None:
            lst[colIdx] = str( int(prevPos + (pos - prevPos)/2.0) )
            print >>outputFile, "\t".join(lst)    
        prevPos = pos


def addEndCoordinates(inp, outp, colIdx, header=False):
    """
    Add end column to HapMap recombination map files.
    """
    inputFile = open(inp)
    outputFile = open(outp, 'w')

    if header:
        lst = inputFile.readline().split()
        l = lst[0:colIdx] + ["start", "end"]
        if len(lst) > colIdx+1:
            l += lst[colIdx+1:]
        print >>outputFile, "\t".join(l)
        
    prevPos = None
    for l in inputFile.xreadlines():
        lst = l.split()
        pos = lst[colIdx]
        if prevPos is not None:
            l = lst[0:colIdx] + [prevPos, pos]
            if len(lst) > colIdx+1:
                l += lst[colIdx+1:]
            print >>outputFile, "\t".join(l)    
            #print >>outputFile, "\t".join([prevPos, pos, lst[2]])    
        prevPos = pos

def categorizeBins(inputFileName, outputFileName, colIdx, percentiles = [0, 25, 50, 75], name='percentile', hasHeader=False):
    """
    Categorize bins from track assigning a level to each.
    """
    quantiles = [x/100.0 for x in percentiles]

    inputFile = open(inputFileName)
    if hasHeader:
        inputFile.readline()

    # get all values:
    values = list()
    for l in inputFile: 
        values.append(float(l.strip().split()[colIdx]))

    # find break values:
    breaks = list()
    for q in quantiles:
        breaks.append(sorted(values)[min(len(values)-1, int(len(values)*q))])

    # open again and add percentiles before printing:
    inputFile.close()
    inputFile = open(inputFileName) 

    with open(outputFileName, 'w') as f:
        if hasHeader:
            header = inputFile.readline().strip()
            print >>f, header + "\t" + name
        for l in inputFile:
            l = l.strip()
            x = float(l.split()[colIdx])
            print >>f, l + "\t%d" % int(quantiles[sum(x > b for b in breaks)-1] * 100)

    inputFile.close()

###################################################################################################
# summary stats
###################################################################################################

def means(buf, **kwargs):
    colIdx, binIdx, rmNaN = kwargs['colIdx'], kwargs['binIdx'], kwargs['rmNaN']
    assert type(colIdx) in [list, tuple]
    lst = list()
    for i in colIdx:
        m, v, s = threeMoments(buf, binIdx=binIdx, colIdx=i, rmNaN=rmNaN)
        result.append(m)
    return lst

#def sumOfSegmentLengths(buf, binSize, colIdx, binIdx, rmNaN=False):
#def sumOfSegmentLengths(buf, binIdx, **kwArgs): # kwargs to capture extra ignored args - to make it usable in smoothedSummaryStats()
def sumOfSegmentLengths(buf, **kwargs): # kwargs to capture extra ignored args - to make it usable in smoothedSummaryStats()
    binIdx = kwargs['binIdx']
    total = 0.0
    assert isinstance(binIdx, list) and len(binIdx) == 2
    for s in buf:
        l = s.split()
        total += float(l[binIdx[1]]) - float(l[binIdx[0]])
    return total


def poissonRateWithCI(buf, **kwargs): # kwargs to capture extra ignored args - to make it usable in smoothedSummaryStats()
    if 'scaling' in kwargs:
        scaling = kwargs['scaling']
    else:
        scaling = lambda x: 1
    binSize = kwargs['binSize']
    N = len(buf)
    L = float(binSize) * scaling(buf)
    rate = N / L
    if N < 2:
        lowerCI, upperCI = [float('nan')] * 2
    else:
        lowerCI = max(0, (1 - 1.96/math.sqrt(N-1)) * N / L)
        upperCI = (1 + 1.96/math.sqrt(N-1)) * N / L
    return [rate, lowerCI, upperCI]


def meanAndStdDev(buf, **kwargs):
    colIdx, rmNaN = kwargs['colIdx'], kwargs['rmNaN']

    if not buf:
        return float('nan'), float('nan')
#    buf = [float(x.split()[colIdx]) == float('nan') for x in buf]
    s2 = 0.0
    s = 0.0
    allIdentical = True
    first = float(buf[0].split()[colIdx])
    for l in buf:
        e = float(l.split()[colIdx])

        if rmNaN and math.isnan(e):
            continue 

        if e != first:
            allIdentical = False
        s += e
        s2 += e * e 

    if not buf:
        return float('nan'), float('nan')

    n = len(buf)
    u = s/n
    v = s2/n - u**2

    if allIdentical:
        return first, 0.0
    else:
        return u, math.sqrt(v)


#def threeMoments(buf, binSize, colIdx, binIdx, rmNaN=False):
def threeMoments(buf, **kwargs):
    colIdx, binIdx, rmNaN = kwargs['colIdx'], kwargs['binIdx'], kwargs['rmNaN']

    transformation = None
    if 'transformation' in kwargs:
        transformation = kwargs['transformation']
    
    if not buf:
        return float('nan'), float('nan'), float('nan')
    s2 = 0.0
    s = 0.0
    allIdentical = True

    if type(buf[0]) in (list, tuple):
        first = float(buf[0][colIdx])
    else:
        first = float(buf[0].split()[colIdx])
    values = list()
    weights = list()
    for l in buf:
        if type(l) in (list, tuple):
            lst = l
        else:
            lst = l.split()
        e = float(lst[colIdx])

        if transformation:
            e = transformation(e)

        if not (rmNaN and math.isnan(e)):
            values.append(e)
            if len(binIdx) == 2:
                weights.append(float(lst[binIdx[1]])-float(lst[binIdx[0]]))
                assert float(lst[binIdx[1]]) > float(lst[binIdx[0]])
        if e != first:
            allIdentical = False
    if allIdentical:
        if math.isnan(first):
            return float('nan'), float('nan'), float('nan')
        else:
            return first, 0.0, 0.0
    if len(binIdx) == 1:
        weights.append(1/float(len(values)))

    if not values:
        return float('nan'), float('nan'), float('nan')

    u = sum(k*x for k, x in zip(weights, values)) / sum(weights)

    v = sum(k*x**2 - u**2 for k, x in zip(weights, values)) / sum(weights)

    A = (sum(k*x**3 for k, x in zip(weights, values)) / sum(weights)) - (sum(k*x**2 + 2*u**3 for k, x in zip(weights, values)) * 3*u / sum(weights))

    s = A / math.sqrt(v)**3

    return u, v, s


def otherStat(buf, binSize, colIdx=[]):
    pass


def summaryStats(inputFileName, binIdx, stats, binSize, outputFileName=None, jumpSize=None, midPoint=False, hasHeader=False):
    """
    Assumes sorted input. Bins lines in a table on binIdxs using in windows of
    windowSize. If supplying more element in onebinIdxs these must be specified in order
    used in the sorting of the input file. Setting jumpSize larger than windowSize creates
    smoothening. The stats argument is a list of functions used to calculate stats for
    each set of binned rows. These functions are passed a list of the binned lines, the
    binSize and an optional colIdx list argument referring to relevant indexes in each
    line. Thay may return a float or a list/tuple of floats.
    """
    inputFile = open(inputFileName, 'r')
    if outputFileName:
        outputFile = open(outputFileName, 'w')
    else:
        outputFile = StringIO.StringIO()
    # maybe instead write a temp file and then read it back in and return as list...

    if not (isinstance(stats, list) or isinstance(stats, tuple)):
        stats = [stats]

    if not (isinstance(binIdx, list) or isinstance(binIdx, tuple)):
        binIdx = [binIdx]

    if hasHeader:
        header = inputFile.readline().strip().split()
        if isinstance(binIdx[0], str):
            binIdx[0] = header.index(binIdx[0])
        if len(binIdx) == 2 and isinstance(binIdx[1], str):
            binIdx[1] = header.index(binIdx[1])
    else:
        assert all(isinstance(binIdx[i], int) for i in range(len(binIdx)))

    if not jumpSize:
        jumpSize = binSize

    assert jumpSize <= binSize 

    def getInput():
        # read lines untill pos is not NA
        lst = None
        while not lst or not all(lst[b] != 'NA' for b in binIdx):
            if inputBuffer:
                l = inputBuffer.pop().strip()
            else:
                l = inputFile.readline().strip()
            if not l:
                return None, None
            else:
                lst = l.split()

        if len(binIdx) == 1:
            assert not midPoint, "don't use midPoint with one binby column"                
            return [float(lst[binIdx[0]])], l
        else:
            if midPoint:
                return [float(lst[binIdx[0]]) + (float(lst[binIdx[1]]) - float(lst[binIdx[0]]))/2.0], l              
            else:
                return [float(lst[x]) for x in binIdx], l

    def overlap(p):
        return p[0] < binStart + binSize
#         if len(p) > 1:
#             assert len(p) == 2
#             assert p[0] >= binStart, str(p) + " " +  str(binStart)
#             return not (p[1] <= binStart or p[0] >= binStart+binSize)
# #            return p[1] >= binStart and p[1] < binStart+binSize  or p[0] < binStart and p[1] >= binStart
#         else:
#             assert p[0] >= binStart, p[0]
#             return p[0] >= binStart and p[0] < binStart + binSize

    def checkForOverhang(p, l):
        if len(p) > 1 and p[1] > binStart+binSize:
            lst = l.split()

            assert p[0] == float(lst[binIdx[0]]) and p[1] == float(lst[binIdx[1]])

            lst[binIdx[0]] = p[0]
            lst[binIdx[1]] = binStart+binSize
            l1 = "\t".join(map(str, lst)) # subst positions

            lst[binIdx[0]] = binStart+binSize
            lst[binIdx[1]] = p[1]
            l2 = "\t".join(map(str, lst)) # subst positions
            
            assert l1 != l2
            assert l2 != l

            inputBuffer.append(l2)
            return (p[0], binStart+jumpSize), l1
        return p, l

    inputBuffer = list()
    binStart = 0
    buf = []
    pos, l = getInput()
    while pos:
        while overlap(pos):
            pos, l = checkForOverhang(pos, l)
            buf.append((pos, l))
            pos, l = getInput()
            if not pos: 
                break

        statsList = list()
        for s in stats:
            if len(buf):
                b = zip(*buf)[1] # get data lines in buf
            else:
                b = buf
            r = s(b)
            if isinstance(r, list) or isinstance(r, tuple):
                statsList.extend(r)
            else:
                statsList.append(r)
        assert statsList

        print >>outputFile, "\t".join(map(str, [binStart, binStart + binSize] + statsList))

        binStart += jumpSize

        while len(buf) and buf[0][0][-1] <= binStart:
            buf.pop(0)

    if not outputFileName:
        return [map(float, t.split('\t')) for t in outputFile.getvalue().strip().split('\n')]


def smoothedSummaryStats(inp, outp, colIdx, binIdx, binSize, windowSizes, stat, statResultMap, hasHeader=False):

    jumpSize = binSize

    midPointSet = None

    tableDict = dict()
    for w in windowSizes:

        def statFunc(buf):
            return stat(buf, binSize=w, colIdx=colIdx, binIdx=binIdx, rmNaN=True) # note: we remove nans generated for diff bins where sim rate is zero

        lst = summaryStats(inp, binIdx=binIdx, stats=statFunc, jumpSize=binSize, binSize=w, hasHeader=hasHeader)

        midPoints = [x[0]+(x[1]-x[0])/2.0 for x in lst]
        results = zip(*lst)
        if midPointSet is None:
            midPointSet = set(midPoints)
        else:
            assert midPointSet.issuperset(set(midPoints)), str(midPointSet.difference(midPoints))

        for k, v in statResultMap.items():
            tableDict["win.%d.%s" % (w, k)] = dict(zip(midPoints, results[v]))

    keys = sorted(tableDict.keys())
    midPoints = sorted(list(midPointSet))    
    
    with open(outp, 'w') as f:
        print >>f, "midPoint\t" + "\t".join(map(str, keys))
        for m in midPoints:
            l = []
            for w in keys:
                if m in tableDict[w]:
                    l.append(str(tableDict[w][m]))
                else:
                    l.append('NA')
            print >>f, str(m) + "\t" + "\t".join(l)


# def segmentOverlapToBins(inputFileName, outputFileName, binSize, nameIdx, startIdx, endIdx, tagPrefix=""):
#     """
#     Compute overlap of chromosomal segments to bins of defined width accross the chromosome.
#     """
#     inputFile = open(inputFileName)
#     outputFile = open(outputFileName, 'w')
# 
#     segmentList = list()
#     for l in inputFile:
#         lst = l.split()
#         segmentList.append([lst[nameIdx], int(lst[startIdx]), int(lst[endIdx])])
#         if segmentList:
#             assert lst[nameIdx] == segmentList[0][0] # check that all is same chromosome
# 
#     segmentList.sort()
#     segmentList = collapseCoordinates(segmentList, 0, 1, 2)
# 
#     name, start, end = 0, 1, 2
#     binStart = 0
#     overlapSum = 0
# 
#     for segment in segmentList:
#         if segment[end] is None:
#             assert int(segment[end]) > segment[end] # ?????
# 
#         while segment[start] >= binStart + binSize:
#             print >>outputFile, "%d\t%d" % (binStart+(binSize/2.0), overlapSum)
#             # start a new bin - possibly print remainder of current segment:
#             overlapSum = 0
#             binStart += binSize
# 
#         if segment[end] < binStart+binSize:
#             overlapSum += segment[end] - max(segment[start], binStart)
#             if not overlapSum >= 0:
#                 print segment[end], segment[start], binStart
#                 sys.exit()
#             assert overlapSum >= 0
#         else:
#             overlapSum += binStart+binSize - max(segment[start], binStart)
#             assert overlapSum >= 0
#             
#         while segment[end] > binStart + 2 * binSize:
#             print >>outputFile, "%d\t%d" % (binStart+(binSize/2.0), overlapSum)
#             # start a new bin:
#             overlapSum = binSize
#             binStart += binSize
# 
#         if segment[end] > binStart + binSize:
#             print >>outputFile, "%d\t%d" % (binStart+(binSize/2.0), overlapSum)
# 
#             # start a new bin:
#             overlapSum = segment[end] - (binStart + binSize)
#             assert overlapSum >= 0
#             binStart += binSize
# 
#     if overlapSum:
#         print >>outputFile, "%d\t%d" % (binStart+(binSize/2.0), overlapSum)


def computeBinDifferences(inputFileNames, outputFileName, colIdx1, colIdx2, normalize=False, relative=False, logRatio=False):
    """
    Compute log10 of ration between stats in bins across the chromosome.
    """
    assert len(inputFileNames) == 2

    track1FileName, track2FileName = inputFileNames
    track1File = open(track1FileName)
    track2File = open(track2FileName)

    outputFile = open(outputFileName, 'w')

    track1list = list()
    track2list = list()
    for l1 in track1File:
        l2 = track2File.readline()
#        if l2 is None:
        if not l2:
            break
        lst = l1.split()
        start1, end1, stat1 = float(lst[0]), float(lst[1]), float(lst[colIdx1])
        lst = l2.split()
        start2, end2, stat2 = float(lst[0]), float(lst[1]), float(lst[colIdx2])
        assert start1 == start2 and end1 == end2
        track1list.append((start1, end1, stat1))
        track2list.append((start2, end2, stat2))

    assert len(track1list) == len(track2list)

    if normalize:
        sum1 = sum(x for x in zip(*track1list)[2] if not math.isnan(x))
        track1list = [(x, y, z/float(sum1)) for x, y, z in track1list]
        sum2 = sum(x for x in zip(*track2list)[2] if not math.isnan(x))
        track2list = [(x, y, z/float(sum2)) for x, y, z in track2list]

#         print sum1, sum(x for x in zip(*track1list)[2] if not math.isnan(x))
#         print sum2, sum(x for x in zip(*track2list)[2] if not math.isnan(x))


    for i in range(len(track1list)):

        start1, end1, stat1 = track1list[i]
        start2, end2, stat2 = track2list[i]

        ratio = float(stat1) - float(stat2)

        if relative:
            if float(stat2):
                ratio /= float(stat2)
            else:
                ratio = float('nan')

        if logRatio:
            ratio = math.log10(ratio)

        ## # http://www.graphpad.com/FAQ/images/Ci%20of%20quotient.pdf
        ## # assuming that the count is large anough for a good normal approximation
        ## sem1 = sd1 / math.sqrt(n1)
        ## sem2 = sd2 / math.sqrt(n2)
        ## t = 1.96
        ## g = (t * sem2/stat2)**2
        ## if g >= 1:
        ##     ratio_SE = (ratio / (1-g)) * math.sqrt((1-g)*(sem1**2/stat1)+(sem2**2/stat2))
        ##     upperCI = ratio/(1-g) - t * ratio_SE
        ##     lowerCI = ratio/(1-g) + t * ratio_SE

        print >>outputFile, "\t".join([str(start1), str(end2), str(ratio)])


def computeInverseTrack(inputFileName, outputFileName, binIdx=[1,2]):

    startIdx, endIdx = binIdx
    prevEnd = None
    with open(outputFileName, 'w') as outputFile:
        with open(inputFileName, 'r') as inputFile:
            reader = csv.reader(inputFile, delimiter='\t')

            t = reader.next()
            start, end = float(t[startIdx]), float(t[endIdx])
            if start > 0:
                t[startIdx], t[endIdx] = str(0), str(start)
                print >>outputFile, "\t".join(t)
            prevEnd = end

            for t in reader:
                start, end = float(t[startIdx]), float(t[endIdx])
                t[startIdx], t[endIdx] = str(prevEnd), str(start)
                print >>outputFile, "\t".join(t)
                prevEnd = end

#         t[startIdx], t[endIdx] = str(prevEnd), str(start)
#         print >>outputFile, "\t".join(t)
# 
#         print prevEnd, float('inf')
                
def addEndCoordinates(inp, outp, colIdx, header=False):
    """
    Add end column to HapMap recombination map files.
    """
    inputFile = open(inp)
    outputFile = open(outp, 'w')

    if header:
        lst = inputFile.readline().split()
        l = lst[0:colIdx] + ["start", "end"]
        if len(lst) > colIdx+1:
            l += lst[colIdx+1:]
        print >>outputFile, "\t".join(l)
        
    prevPos = None
    for l in inputFile.xreadlines():
        lst = l.split()
        pos = lst[colIdx]
        if prevPos is not None:
            l = lst[0:colIdx] + [prevPos, pos]
            if len(lst) > colIdx+1:
                l += lst[colIdx+1:]
            print >>outputFile, "\t".join(l)    
            #print >>outputFile, "\t".join([prevPos, pos, lst[2]])    
        prevPos = pos

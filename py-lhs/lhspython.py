# lhspython.py
# cec24@phy.duke.edu
# ccs: july-8-2010
#
# generates npoints samples in ndim dimensions
# the samples should be generated as a random latin hypercube
# there are more complicated algorithms to generate an "optimal" sample
# but this is probably not worth the effort

from __future__ import division
import random

# just test in 1d to start, 
# LHS algorithm:
# divide the 1d region into npoints contiguous segments, 
# generate a random number in the segment
# this is the 1d lhs sample
# 
# in m'd we would generate m*1d lhs samples 
# we randomly combine m elements from the set of all samples in such
# a way that each sample covers a rooks move on the sample space, i.e in a 2d space if
# we pick one sample at i,j then 

def generate1dSample(npoints):
    """ returns a vector of n random samples 
    the samples are generated within equal segments of the unit interval
    """
    segmentSize = 1/npoints
    points = []
    for i in range(npoints):
        # this is the rolling lower bound of each segment
        segmentMin = i*segmentSize
    # now take a random location within this region
        points.append(segmentMin + (random.random()*segmentSize))
    #    print points[i]
    return points

def generateNdSample(ndim, npoints):
    """ generate a vector of ndimensional samples, of length npoints
    """
    # samples
    samples = []
    # rows holds ndim lists of nsamples each
    rowList = []
    for i in range(ndim):
    #    print(i)
        rowList.append(generate1dSample(npoints))

    # this is a shuffled list of indicies
    indexList = []
    for i in range(ndim):
        indexList.append(generateShuffledIndexList(npoints))

    for i in range(npoints):
        samples.append(makeSample(ndim, rowList, indexList))
        
    return samples

    
def makeSample(ndim, points, indexList):
    """ takes an ndim wide list of lists of nsamples (points) and 
    creates a single lhs sample 
    """
    # this is a list of the indices of points that we want access to
    sampleIndex = []    
    for i in range(ndim):
        sampleIndex.append(indexList[i].pop())
    #print sampleIndex

    sample = []

    for i in range(ndim):
        index = sampleIndex[i]
        # this is a sneaky way to deref a list to get to the sublist
        sample.append((points[i])[index])

    return sample

def generateShuffledIndexList(npoints):
    """ generate a list [1..npoints] and shuffle """
    indexList = range(npoints)
    random.shuffle(indexList)
    return(indexList)


if __name__ == "__main__":
    npoints = 100
    ndim = 10
    test = generateNdSample(ndim, npoints)
    
    for i in range(npoints):
        print
        for j in range(ndim):
            print (test[i])[j],"\t", 
        


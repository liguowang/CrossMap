#!/usr/bin/env python
'''
This python3 module was converted from python2.7 code using 2to3
'''


def tillingBed(chrName, chrSize, stepSize=10000):
    '''tilling whome genome into small sizes'''
    # tilling genome
    for start in range(0, chrSize, stepSize):
        end = start + stepSize
        if end < chrSize:
            yield (chrName, start, end)
        else:
            yield (chrName, start, chrSize)

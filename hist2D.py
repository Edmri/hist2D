#!python3
#! coding: utf-8

# import python package
from typing import Any
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib import cm as CM
import gzip
from operator import itemgetter
import random


##################    Extract data from files    ##################
def extractData(readFile, verbose):
    """
    Return a list of intervals (BED)

    Parameters:
    readFile (str) : Path to peaks (BED) file
    verbose  (bool): Enable verbose mode

    Return:
    (list[list[int | str]]): The peaks position
    """
    if verbose:
        print(f"extractData on {readFile}")
    peak_datas = []

    # Handle compresion
    if readFile.endswith(".gz"):
        with gzip.open(readFile, 'r') as bedF:
            for i in bedF:
                # Handle bytes as string
                line = i.decode("utf-8")
                line = line.strip()
                infos = line.split("\t")
                # genomic coord
                infos[1] = int(infos[1])
                infos[2] = int(infos[2])
                peak_datas.append(infos)
    else:
        with open(readFile, 'r') as bedF:
            for i in bedF:
                i = i.strip()
                infos = i.split("\t")
                # genomic coord
                infos[1] = int(infos[1])
                infos[2] = int(infos[2])
                peak_datas.append(infos)
    return(peak_datas)


def sortGenomeCoord(intervals, verbose):
    """
    Return sorted genomics coordinates

    Parameters:
    intervals (list[list[str | int]])  : Genomics intervals
    verbose   (bool)                   : Enbable verbose mode

    Return:
    (list[list[str | int]]): Sorted positions
    """
    if verbose:
        print("sortGenomeCoord")
    return(sorted(intervals, key=itemgetter(0, 1, 2)))


# get a dico with chromosome key and a list of peak/frag
def fileToDictList(readFile, verbose):
    """
    Return a per chromosome indexed dictionary

    Parameters:
    readFile    (str)   : Path to peak file
    verbose     (bool)  : Enbable verbose mode

    Return:
    (dict[str, list[str | int]]): The dna fragment positions
    """
    if verbose:
        print(f"fileToDictList on {readFile}")
    frag_datas = {}

    if readFile.endswith(".gz"):
        with gzip.open(readFile, 'r') as bedF:
            for i in bedF:
                line = i.decode("utf-8")
                infos = line.split("\t")
                chrNB = infos[0]
                infos[1] = int(infos[1])
                infos[2] = int(infos[2])
                if chrNB not in frag_datas.keys():
                    frag_datas[chrNB] = []
                    frag_datas[chrNB].append(infos)
                else:
                    frag_datas[chrNB].append(infos)

    else:
        with open(readFile, 'r') as bedF:
            for i in bedF:
                infos = i.split("\t")
                chrNB = infos[0]
                infos[1] = int(infos[1])
                infos[2] = int(infos[2])
                if chrNB not in datas.keys():
                    frag_datas[chrNB] = []
                    frag_datas[chrNB].append(infos)
                else:
                    frag_datas[chrNB].append(infos)

    for chrNB in frag_datas.keys():
        getInfo = False
        sortCoord = sortGenomeCoord(frag_datas[chrNB], getInfo)
        frag_datas[chrNB] = sortCoord

    return(frag_datas)


##################    2D histogram    ##################
def getCenterFragSize(sorted_intervals):
    """
    Return multiple information about fragments (peak, size, center, strand)

    Parameters:
    sorted_intervals (list[list[str | int]): A sorted list of intervals

    Return:
    (list[list[str | int]]): List of centers of peaks (and strand)
    """
    centerSize = []

    for i in sorted_intervals:
        # get the base position at the center fragment
        center = (int(i[1]) + int(i[2]) + 1) / 2

        # if the center is not an integer (i.e even size fragment) radomly
        # choose the base before or after the center position to avoid artefacts

        if (int(i[1]) + int(i[2]) + 1) % 2 != 0:
            center += random.choice([0.5, -0.5])

        # convert float in int format
        center = int("{0:.0f}".format(center))

        if(len(i)>= 4):
            strand = i[-1].strip()
            centerSize.append([i[0], center, strand])
        else:
            centerSize.append([i[0], center])
    return(centerSize)


# Return [[dist int, fagsize int],...]
def getDist(peakList, fragList, params, verbose):
    """
    Return distance between center's peaks and fragments and fragment's sizes

    Parameters:
    peakList (list[list[int | str]])       : The peaks position
    fragList (dict[str, list[str | int]])  : The dna fragment positions
    params   (dict[str, Any])                  : User paramters
    verbose  (bool)                  : Enbable verbose mode

    Return:
    (list[list[int]]): List of distance between center'peaks and fragments, and fragment's sizes
	"""
    if verbose:
        print("getDist")
    window = params["window"]
    minWin = -1 * window
    distGroup = []
    strand = ""
    # sort referenced peak by genomic coordinates
    sortpeakList = sortGenomeCoord(peakList, verbose)
    # get referenced peak center
    centerSize = getCenterFragSize(sortpeakList)
    # calcul of distance between peak list and dna fragments
    chrPeak = ""
    startIndx = 0
    for refPeak in centerSize:
        if chrPeak != refPeak[0]:
            startIndx = 0
        chrPeak = refPeak[0]
        midPeak = refPeak[1]
        if len(refPeak) == 3:
            strand = refPeak[2]
        if chrPeak in fragList.keys():
            maxEnd = len(fragList[chrPeak])
            distIdx = []
            for fragment in range(startIndx, maxEnd, 1):
                fragInfo = fragList[chrPeak][fragment]
                midFrag = (fragInfo[1] + fragInfo[2]) / 2
                frgSize = fragInfo[2] - fragInfo[1]
                midFrag = int("{0:.0f}".format(midFrag))
                distPeakFrag = midFrag - midPeak
                if distPeakFrag >= minWin and distPeakFrag <= window:
                    distIdx.append(fragment)
                    # change orientation when the peak is on downstream strand
                    if strand == "-":
                        distPeakFrag = distPeakFrag * -1
                    distGroup.append([int(distPeakFrag), int(frgSize)])
                if distPeakFrag > (window+250):
                    if distIdx:
                        startIndx = min(distIdx)
                    break
    return(distGroup)


def histTwoD(vplotList, params, verbose):
    """
    Return 2D histogram values from distance and size

    Parameters:
    vplotList (list[list[int]]) : List of distance between center'peaks and fragments, and fragment's sizes
    params    (dict[str, Any])            : User paramters
    verbose   (bool)            : Enbable verbose mode

    Return:
    xedges (ndarray, shape(nx+1,))  : The bin edges along the first dimension (distance)
    yedges (ndarray, shape(ny+1,))  : The bin edges along the second dimension (size)
    H      (ndarray, shape(nx, ny)) : The bi-dimensional histogram of samples x and y
	"""
    if verbose:
        print("conversion with histTwoD")
    window = params["window"]
    # draw selected fragments
    ## x : distance between peak's centrer and fragment's center
    x = [int(item[0]) for item in vplotList]
    ## y : fragment's size
    y = [int(item[1]) for item in vplotList]
    minY = min(y)
    maxY = max(y)
    delta = (maxY-minY)
    # H : The bi-dimensional histogram of samples x and y.
    #     Values in x are histogrammed along the first dimension and
    #     values in y are histogrammed along the second dimension.
    H, xedges, yedges = np.histogram2d(x, y, bins=(window, delta))
    # aplicated a ratio value on the intensity of each H elements to modify
    # signal intensity when drawing the 2D histogram
    if float(params["ratioFac"]) != 1.0:
        ratioFac = params["ratioFac"]
        H = [l*ratioFac for l in H]
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
    return(xedges, yedges, H)


def vplot(vplotList, params, verbose):
    """
    Save 2D histogram plot in png

    Parameters:
    vplotList (list[list[int]]) : List of distance between center'peaks and fragments, and fragment's sizes
    params    (dict[str, Any])            : User paramters
    verbose   (bool)            : Enbable verbose mode
	"""
    if verbose:
        print("vplot")
    #### init
    vName = params["vName"]
    xlab = params["xlabel"]
    ylab = params["ylabel"]
    ## Window size
    # get the min and the max values of the x-axis of the vplot
    window = params["window"]
    xmin = window * -1
    xmax = window
    #### 2D histogram
    # 2D histogram png params
    fig = plt.figure(figsize=(9.27, 8.78), dpi=100)
    grid = plt.GridSpec(4, 7)
    mainPlot = plt.subplot(grid[0:, :6])
    # get 2D histogram values
    xedges, yedges, H = histTwoD(vplotList, params, verbose)
    # set the min and the max values of the y-axis of the 2D histogram
    ymin = 0
    ymax = max(yedges)
    # Maximum value for the signal intensity
    if not params["colorLim"]:
        maxLim = np.max(H)
    else:
        maxLim = params["colorLim"]
    # pseudocolor plot of the 2D histogram in red colours
    histGraph = mainPlot.pcolormesh(xedges, yedges, H, vmin=0, vmax=maxLim, cmap=CM.Reds)
    # ADD colorbar and label on the 2D histogram
    cax = fig.add_axes([0.82, 0.110, 0.02, 0.770])
    cbar = fig.colorbar(histGraph, cax=cax)
    cbar.ax.set_ylabel('Counts')
    # RESIZE graph
    mainPlot.axis([xmin, xmax, ymin, ymax])
    mainPlot.set_xlabel(xlab)
    mainPlot.set_ylabel(ylab)
    #### save
    if vName:
        if verbose:
            print('Saving 2D histogram in png file')
        plt.savefig(vName, dpi=200)
    else:
        plt.show()
    plt.cla()
    plt.clf()
    plt.close()


##################  Control of the command line parameters  ##################
def positiveInt(v, tag, verbose):
    """
    Check if the value of the parameter is a positive integer (return a default value of 1)

    Parameters:
    v        (int)  : Value of the tag
    tag      (str)  : Tag name
    verbose  (bool) : Enbable verbose mode

    Return:
	"""
    if(verbose):
        print(f"positiveInt on {tag}: {v}")
    defaultVal = 1
    try:
        v = int(v)
    except ValueError:
        warnings.warn(f"Not integer type for {tag} value")
        warnings.warn("Value of 1 is used now")
    else:
        if v > 0:
            defaultVal = v
        else:
            warnings.warn("Negative integer not allowed, set value as 1")
    return(defaultVal)


def checkFloat(v, tag, verbose):
    """
    Check if the value of the parameter is a float

    Parameters:
    v        (int)  : Value of the tag
    tag      (str)  : Tag name
    verbose  (bool) : Enbable verbose mode

    Return:
    Raise an error if the value is not an float
	"""
    if(verbose):
        print(f"checkFloat on {tag} {v}")
    try:
        v = float(v)
    except ValueError:
        warnings.warn(f"Not float type for {tag} value")
        exit(1)
    return(v)


def checkOptions(options, params, verbose):
    """
    Returns a dico of parameters

    Parameters:
    options  (dict[str, Any]) : User options in a dico
    params   (dict[str, Any]) : Default options in a dico
    verbose  (bool) : Enbable verbose mode

    Return:
    (dico): all the user parameters with values
	"""
    if(verbose):
        print("checkOptions")
    for option in options.keys():
        k = option
        v = options[k]
        if k == "ratioFac" and v != "":
            params["ratioFac"] = float(v)
        if k == "colorLim" and v != "":
            v = checkFloat(v, k, verbose)
            params[k] = v
        if k == "windowSize" and v != "":
            v = positiveInt(v, k, verbose)
            params["window"] = int(v)
    return(params)


##################    CONFIG    ##################
# initiat config dico:
def defaultConf(verbose):
    """
    Returns a dico of paramters with default values

    Parameters:
    verbose  (bool) : Enbable verbose mode

    Return:
    (dict[str, str | int | float | bool]): default parameters with values
	"""
    if(verbose):
        print("defaultConf")
    params = {}
    params["vName"] = ""
    params["dict"] = False
    params["xlabel"] = "Distance from peak midpoints (bp)"
    params["ylabel"] = "Fragment size (bp)"
    return(params)


##################    command line    ##################
parser = argparse.ArgumentParser(description="Create a vplot png from fragment file and genomic region")
parser.add_argument("-n", help = "genomic region to plot (bed file)", required=True)
parser.add_argument("-p", help = "fragments sorted by genomic coordinate fragments (bed3)", required=True)
parser.add_argument("-o", help = "vplot file name (png)", required=True)
parser.add_argument("-w", help = "windows size valu (positive int)", required=True)
parser.add_argument("-l", help = "limit colour scale (default False, positive int)", default=False)
parser.add_argument("-r", help = "ratio to decrease the signal (default 1.0, float)", default=1.0, type=float)
parser.add_argument("-v", help = "verbose mode, true if used (default False)", action='store_true', default=False)
args = parser.parse_args()


### if parameters precise, get the value
nfrFile = args.n
peaks = args.p
vplotName = args.o
windowSize = args.w
verbose = args.v
colorLim = args.l
ratioFac = args.r

##################    run    ##################
confDico = {}
confDico["colorLim"] = colorLim
confDico["ratioFac"] = ratioFac
confDico["windowSize"] = windowSize

if verbose:
    print("Verbose mode activate")

## Congif
params = defaultConf(verbose)
params = checkOptions(confDico, params, verbose)
params["vName"] = vplotName

## Create vplot
if verbose:
    print("----> peak extraction")
nfrList = extractData(nfrFile, verbose)
if verbose:
    print("----> fragment extraction")
frgtList = fileToDictList(peaks, verbose)
if verbose:
    print("----> distance between fragments and peaks")
vplotList = getDist(nfrList, frgtList, params, verbose)
if verbose:
    print("----> 2D histogram")
vplot(vplotList, params, verbose)

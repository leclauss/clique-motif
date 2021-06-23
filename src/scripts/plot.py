import time
import matplotlib.pyplot as plt

from cliqueMotif import loadTS
from scripts.runBenchmark import getMetaInformation, parseOutput, run, algorithms
from scripts.getStats import getMotifPoints, getPointBasedScores


def plotMotif(algorithmName, tsPath, tsMetaPath, timeout=None, maxMemory=None):
    algorithm = algorithms[algorithmName]
    ts = loadTS(tsPath)
    info = getMetaInformation(tsMetaPath)
    length = int(info["length"])
    window = int(info["window"])
    radius = float(info["range"])
    motif = [int(index) for index in info["matchings"].split(",")]

    startTime = time.time()
    returnCode, output = run(algorithm(tsPath, length, window, radius), timeout=timeout, maxMemory=maxMemory)
    runtime = time.time() - startTime
    success = (returnCode == 0)
    print("Success:", success)
    if success:
        motifsFound, stats = parseOutput(output)
    else:
        motifsFound = []

    f1Score, precision, recall, _, index = getPointBasedScores(motif, motifsFound, window)
    bestMotifFound = motifsFound[index] if index != -1 else []
    print("Runtime:", "{:.3f}".format(runtime), "s", "\nF1-Score:", "{:.4f}".format(f1Score), "\nPrecision:",
          "{:.2%}".format(precision), "\nRecall:", "{:.2%}".format(recall))

    motifPoints = getMotifPoints(motif, window)
    foundPoints = getMotifPoints(bestMotifFound, window)

    truePositiveRanges = getRanges(motifPoints.intersection(foundPoints))
    falsePositiveRanges = getRanges(foundPoints.difference(motifPoints))
    falseNegativeRanges = getRanges(motifPoints.difference(foundPoints))

    #import matplotlib.dates as mdates
    #import datetime as dt
    #start = dt.datetime(2014,7,1,0,0,0)
    #end = dt.datetime(2015,2,1,0,0,0)
    #times = mdates.drange(start, end, dt.timedelta(minutes=30))

    fig, ax = plt.subplots()

    #plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%a, %Y-%m-%d'))  # %H:%M
    #plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=1))
    #plt.plot(times,ts)
    ax.plot(ts)
    #plt.gcf().autofmt_xdate()
    for start, end in truePositiveRanges:
        ax.axvspan(start - 0.5, end + 0.5, alpha=0.2, facecolor='green')
        #plt.axvspan(times[start], times[end], alpha=0.2, facecolor='green')
    for start, end in falsePositiveRanges:
        ax.axvspan(start - 0.5, end + 0.5, alpha=0.2, facecolor='green')
        #plt.axvspan(times[start], times[end], alpha=0.2, facecolor='green')
    for start, end in falseNegativeRanges:
        #ax.axvspan(start - 0.5, end + 0.5, alpha=0.2, facecolor='yellow')
        pass
        #plt.axvspan(times[start], times[end], alpha=0.2, facecolor='yellow')

    return bestMotifFound, ts, window


def getRanges(points):
    sortedList = sorted(points)
    if len(sortedList) == 0:
        return []
    startPoint = sortedList[0]
    lastPoint = sortedList[0]
    ranges = []
    for point in sortedList[1:]:
        if point > lastPoint + 1:
            ranges.append((startPoint, lastPoint))
            startPoint = point
        lastPoint = point
    ranges.append((startPoint, lastPoint))
    return ranges

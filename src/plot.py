import time
import matplotlib.pyplot as plt

from cliqueMotif import loadTS
from runBenchmark import getMetaInformation, getMotifPoints, getPointBasedScores, parseOutput, algorithms


def plotMotif(algorithm, tsPath, tsMetaPath):
    ts = loadTS(tsPath)
    info = getMetaInformation(tsMetaPath)
    length = int(info["length"])
    window = int(info["window"])
    radius = float(info["range"])
    motif = [int(index) for index in info["matchings"].split(",")]

    startTime = time.time()
    output = algorithms[algorithm](tsPath, length, window, radius)
    runtime = time.time() - startTime

    motifsFound, _ = parseOutput(output)
    f1Score, precision, recall, _, index = getPointBasedScores(motif, motifsFound, window)
    bestMotifFound = motifsFound[index] if index != -1 else []
    print("Runtime:", "{:.3f}".format(runtime), "s", "\nF1-Score:", "{:.4f}".format(f1Score), "\nPrecision:",
          "{:.2%}".format(precision), "\nRecall:", "{:.2%}".format(recall))

    motifPoints = getMotifPoints(motif, window)
    foundPoints = getMotifPoints(bestMotifFound, window)

    truePositiveRanges = getRanges(motifPoints.intersection(foundPoints))
    falsePositiveRanges = getRanges(foundPoints.difference(motifPoints))
    falseNegativeRanges = getRanges(motifPoints.difference(foundPoints))

    plt.clf()
    plt.cla()
    plt.close()
    plt.plot(ts)
    for start, end in truePositiveRanges:
        plt.axvspan(start - 0.5, end + 0.5, alpha=0.2, facecolor='green')
    for start, end in falsePositiveRanges:
        plt.axvspan(start - 0.5, end + 0.5, alpha=0.2, facecolor='red')
    for start, end in falseNegativeRanges:
        plt.axvspan(start - 0.5, end + 0.5, alpha=0.2, facecolor='yellow')
    plt.show()


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

import os
import re
import subprocess
import sys
import time
from pathlib import Path

from src.cliqueMotif import getTopMotif

path = Path(__file__).parent.absolute()

# settings
benchmarkPath = path / "benchmark"
algDir = path / "../tsgenerator/benchmark/algorithms"
radiusMultiplier = 1.0
algorithms = {
    "EMMA": lambda tsPath, _, window, radius: emma(tsPath, window, radius),
    "GrammarViz3.0": lambda tsPath, _, window, __: grammarViz(tsPath, window),
    "ScanMK": lambda tsPath, _, window, radius: scanMk(tsPath, window, radius),
    "SetFinder": lambda tsPath, _, window, radius: setFinder(tsPath, window, radius),
    "LearnMotifs": lambda tsPath, length, window, radius: learnMotifs(tsPath, length, window, radius),
    "ClusterMK": lambda tsPath, _, window, radius: clusterMk(tsPath, window, radius),
    "CliqueMotif": lambda tsPath, _, window, radius: cliqueMotif(tsPath, window, radius)
}


def emma(tsPath, window, radius):
    return getOutput(["java", "-jar", algDir / "EMMA.jar", tsPath,
                      str(window), str(2 * radius), "6", "4", "1.0"])


def grammarViz(tsPath, window):
    return getOutput(["java", "-jar", algDir / "GrammarViz3.jar", "-d", tsPath,
                      "-w", str(window), "-p", "6", "-a", "4"])


def scanMk(tsPath, window, radius):
    return getOutput(["java", "-jar", algDir / "ScanMK.jar", "d=" + str(tsPath),
                      "w=" + str(window), "k=1", "r=" + str(radius)])


def setFinder(tsPath, window, radius):
    return getOutput(["java", "-jar", algDir / "SetFinder.jar", "d=" + str(tsPath),
                      "w=" + str(window), "k=1", "r=" + str(radius)])


def learnMotifs(tsPath, length, window, radius):
    return getOutput(["java", "-jar", algDir / "LearnMotifs.jar", "dataSet=" + str(tsPath),
                      "ETA=0.1", "maxItr=1000", "numRandomRestarts=200", "alpha=2", "K=1", "pct=1",
                      "tsLength=" + str(length), "w=" + str(window), "t=" + str(radius)])


def clusterMk(tsPath, window, radius):
    return getOutput(["java", "-jar", algDir / "ClusterMK.jar", "d=" + str(tsPath),
                      "w=" + str(window), "k=1", "r=" + str(radius)])


def cliqueMotif(tsPath, window, radius):
    motif, stats = getTopMotif(window, radius, str(tsPath))
    output = ""
    for index in motif:
        output += str(index) + "\n"
    output += "end\n"
    for stat in stats:
        output += str(stat) + "\n"
    return output


def getOutput(args):
    return subprocess.run(args, capture_output=True, text=True).stdout


def main(args):
    runBenchmark(benchmarkPath)
    return 0


def runBenchmark(benchmarkDir):
    with open(benchmarkDir / "stats.csv", "w") as statsFile, open(benchmarkDir / "info.csv", "w") as infoFile:
        statsFile.write("Time Series Length,Method,Motif Shape,Noise,Motif Size,Window,Range,Algorithm,"
                        "Runtime,Precision,Recall,F1,Found Size\n")

        # walk benchmark directory recursively
        for directory, tsName, tsMetaName in getBenchmarkFiles(benchmarkDir):
            tsDir = Path(directory)
            tsPath = tsDir / tsName
            tsMetaPath = tsDir / tsMetaName

            # get meta information
            info = getMetaInformation(tsMetaPath)
            length = int(info["length"])
            window = int(info["window"])
            motifRange = float(info["range"])
            motif = [int(index) for index in info["matchings"].split(",")]

            for name, algorithm in algorithms.items():
                # start algorithm
                startTime = time.time()
                output = algorithm(tsPath, length, window, motifRange * radiusMultiplier)
                runtime = time.time() - startTime

                # save output and runtime
                with open(tsDir / ("output_" + name + ".txt"), "w") as outFile:
                    outFile.write(output)
                with open(tsDir / ("runtime_" + name + ".txt"), "w") as runtimeFile:
                    runtimeFile.write(str(runtime))

                # parse output, calculate measures and save stats
                motifsFound, infoLines = parseOutput(output)
                f1Score, precision, recall, foundSize = getPointBasedScores(motif, motifsFound, window)
                runInfo = info["length"] + "," + info["method"] + "," + info["type"] + "," + info["noise"] + "," + \
                          info["size"] + "," + info["window"] + "," + info["range"] + "," + name + ","
                statsFile.write(runInfo + str(runtime) + "," + str(precision) + "," + str(recall) + "," + str(f1Score)
                                + "," + str(foundSize) + "\n")
                if len(infoLines) > 0:
                    infoFile.write(runInfo + ','.join(infoLines) + "\n")


def getBenchmarkFiles(benchmarkDir):
    tsPattern = re.compile('time_series_\d+\.csv')
    tsMetaPattern = re.compile('time_series_meta_\d+\.csv')
    for root, dirNames, fileNames in os.walk(benchmarkDir):
        tsNames = list(filter(tsPattern.match, fileNames))
        tsMetaNames = list(filter(tsMetaPattern.match, fileNames))
        if len(tsNames) == 1 and len(tsMetaNames) == 1:
            tsName = tsNames[0]
            tsMetaName = tsMetaNames[0]
            yield root, tsName, tsMetaName


def getMetaInformation(tsMetaPath):
    info = {}
    with open(tsMetaPath) as metaFile:
        for line in metaFile:
            splits = line.split(",", 1)
            if len(splits) > 1:
                info[splits[0].strip()] = splits[1].strip()
    return info


def parseOutput(output):
    motifs = []
    currentMotif = []
    infoLines = []
    endFound = False
    for line in output.splitlines():
        if endFound:
            infoLines.append(line.strip())
        else:
            isNext = line.startswith("next")
            isEnd = line.startswith("end")
            if isNext or isEnd:
                if len(currentMotif) > 0:
                    motifs.append(currentMotif)
                    currentMotif = []
                if isEnd:
                    endFound = True
            else:
                currentMotif.append(int(line))
    if len(currentMotif) > 0:
        motifs.append(currentMotif)
    return motifs, infoLines


def getPointBasedScores(motif, motifsFound, window):
    # get relevant motif points
    motifPoints = getMotifPoints(motif, window)

    bestStats = (0.0, 0.0, 0.0, 0)  # f1Score, precision, recall, foundSize
    # calculate measures for each motif found
    for motifFound in motifsFound:
        foundPoints = getMotifPoints(motifFound, window)
        truePositives = len(motifPoints.intersection(foundPoints))
        falsePositives = len(foundPoints.difference(motifPoints))
        falseNegatives = len(motifPoints.difference(foundPoints))
        # calculate precision
        returnedPositives = truePositives + falsePositives
        precision = (truePositives / returnedPositives if returnedPositives > 0 else 0.0)
        # calculate recall
        actualPositives = truePositives + falseNegatives
        recall = (truePositives / actualPositives if actualPositives > 0 else 0.0)
        # calculate f1 score
        f1Score = (2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0)
        # save best f1 score
        if f1Score > bestStats[0]:
            bestStats = (f1Score, precision, recall, len(motifFound))
    return bestStats


def getMotifPoints(motif, window):
    motifPoints = set()
    for index in motif:
        motifPoints = motifPoints.union(set(range(index, index + window)))
    return motifPoints


if __name__ == '__main__':
    sys.exit(main(sys.argv))

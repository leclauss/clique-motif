import argparse
import sys
import ast
from pathlib import Path

from runBenchmark import getMetaInformation


def main():
    parser = argparse.ArgumentParser(description="Generate a stats file from the benchmark output.")
    parser.add_argument("inputPath", type=Path, help="relative path to the input benchmark file")
    parser.add_argument("outputPath", type=Path, help="relative path to the output file")
    parser.add_argument("timeout", type=float, help="timeout used for the benchmark")
    parser.add_argument("--sep", default=";", help="separator in the output csv file (default: ';')")
    parser.add_argument("--ground-truth", default=None, help="the algorithm whose output motif is ground truth "
                                                             "(default: use generator motif)")
    args = parser.parse_args()
    generateStatsFile(args.inputPath, args.outputPath, args.ground_truth, args.sep, args.timeout)
    return 0


def generateStatsFile(inputPath, outputPath, groundTruth, separator, timeout):
    headerStart = "TS Path"
    with open(inputPath, "r") as runFile:
        tsRuns = {}
        algorithms = []
        for line in runFile:
            if not line.startswith(headerStart):
                data = line.strip().split(";")
                algorithm = data[3]
                appendResult(tsRuns, tuple(data[0:3]), algorithm, data[4:])
                if algorithm not in algorithms:
                    algorithms.append(algorithm)

    with open(outputPath, "w") as outputFile:
        outputFile.write(separator.join(["Length", "Method", "Type", "Noise", "Size", "Window", "Range", ""]))
        for algorithm in algorithms:
            columns = ["Runtime", "Error", "Precision", "Recall", "F1", "Found_Size"]
            if algorithm == "CliqueMotif":
                columns += ["Nodes", "Edges", "File_Size", "Graph_Time", "Clique_Time"]
            algCols = []
            for column in columns:
                algCols.append(algorithm + "_" + column)
            outputFile.write(separator.join(algCols))
        outputFile.write("\n")

        for tsRun in tsRuns.keys():
            tsPath, tsMetaPath, radius = tsRun
            info = getMetaInformation(tsMetaPath)
            window = int(info["window"])
            line = [info["length"], info["method"], info["type"], info["noise"], info["size"], info["window"], radius]
            results = tsRuns[tsRun]
            if groundTruth is None:
                motif = [int(index) for index in info["matchings"].split(",")]
            else:
                motif = None
                if groundTruth in results.keys():
                    motifsFound = ast.literal_eval(results[groundTruth][1])
                    if motifsFound is not None and len(motifsFound) > 0:
                        motif = motifsFound[0]
            if motif is None:
                print("No ground truth motif found for", tsRun)

            for algorithm in algorithms:
                if algorithm not in results.keys():
                    print("No", algorithm, "result for", tsRun)
                    line += [""] * 6
                    if algorithm == "CliqueMotif":
                        line += [""] * 5
                else:
                    result = results[algorithm]
                    runtime = result[0]
                    motifsFound = ast.literal_eval(result[1])
                    stats = ast.literal_eval(result[2])
                    if motifsFound is None:
                        error = 1 if float(runtime) > timeout else 2
                        f1Score, precision, recall, lengthFound = 0.0, 0.0, 0.0, 0
                    else:
                        error = 0
                        if motif is None:
                            f1Score, precision, recall, lengthFound = "", "", "", ""
                        else:
                            f1Score, precision, recall, lengthFound, _ = getPointBasedScores(motif, motifsFound, window)
                    line += [runtime, str(error), str(precision), str(recall), str(f1Score), str(lengthFound)]
                    if algorithm == "CliqueMotif":
                        line += stats
            outputFile.write(separator.join(line) + "\n")


def appendResult(resultDict, tsRun, algorithm, result):
    if tsRun not in resultDict.keys():
        resultDict[tsRun] = {algorithm: result}
    else:
        resultDict[tsRun][algorithm] = result


def getPointBasedScores(motif, motifsFound, window):
    # get relevant motif points
    motifPoints = getMotifPoints(motif, window)

    bestStats = (0.0, 0.0, 0.0, 0, -1)  # f1Score, precision, recall, foundSize, index
    # calculate measures for each motif found
    for index, motifFound in enumerate(motifsFound):
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
            bestStats = (f1Score, precision, recall, len(motifFound), index)
    return bestStats


def getMotifPoints(motif, window):
    motifPoints = set()
    for index in motif:
        motifPoints = motifPoints.union(set(range(index, index + window)))
    return motifPoints


if __name__ == '__main__':
    sys.exit(main())

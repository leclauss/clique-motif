import os
import re
import subprocess
import sys
import time
import argparse
from pathlib import Path

from cliqueMotif import getTopMotif
from util import ParallelRun

path = Path(__file__).parent.absolute()

# algorithm settings
algDir = path / "../tsgenerator/benchmark/algorithms"
algorithms = {
    "CliqueMotif": lambda tsPath, _, window, radius, timeout: cliqueMotif(tsPath, window, radius, timeout),
    "EMMA": lambda tsPath, _, window, radius, timeout: emma(tsPath, window, radius, timeout),
    "GrammarViz3.0": lambda tsPath, _, window, __, timeout: grammarViz(tsPath, window, timeout),
    "ScanMK": lambda tsPath, _, window, radius, timeout: scanMk(tsPath, window, radius, timeout),
    "SetFinder": lambda tsPath, _, window, radius, timeout: setFinder(tsPath, window, radius, timeout),
    "ClusterMK": lambda tsPath, _, window, radius, timeout: clusterMk(tsPath, window, radius, timeout),
    "LearnMotifs": lambda tsPath, length, window, radius, timeout: learnMotifs(tsPath, length, window, radius, timeout)
}


def cliqueMotif(tsPath, window, radius, timeout):
    motif, stats = getTopMotif(window, radius, str(tsPath), timeout=timeout)
    return [motif] if motif is not None else None, stats


def emma(tsPath, window, radius, timeout):
    return parseOutput(run(["java", "-jar", algDir / "EMMA.jar", tsPath,
                            str(window), str(2 * radius), "6", "4", "1.0"], timeout=timeout)), None


def grammarViz(tsPath, window, timeout):
    return parseOutput(run(["java", "-jar", algDir / "GrammarViz3.jar", "-d", tsPath,
                            "-w", str(window), "-p", "6", "-a", "4"], timeout=timeout)), None


def scanMk(tsPath, window, radius, timeout):
    return parseOutput(run(["java", "-jar", algDir / "ScanMK.jar", "d=" + str(tsPath),
                            "w=" + str(window), "k=1", "r=" + str(radius)], timeout=timeout)), None


def setFinder(tsPath, window, radius, timeout):
    return parseOutput(run(["java", "-jar", algDir / "SetFinder.jar", "d=" + str(tsPath),
                            "w=" + str(window), "k=1", "r=" + str(radius)], timeout=timeout)), None


def clusterMk(tsPath, window, radius, timeout):
    return parseOutput(run(["java", "-jar", algDir / "ClusterMK.jar", "d=" + str(tsPath),
                            "w=" + str(window), "k=1", "r=" + str(radius)], timeout=timeout)), None


def learnMotifs(tsPath, length, window, radius, timeout):
    return parseOutput(run(["java", "-jar", algDir / "LearnMotifs.jar", "dataSet=" + str(tsPath),
                            "eta=0.1", "maxIter=1000", "numRandomRestarts=200", "alpha=2", "K=3", "pct=1",
                            "tsLength=" + str(length), "w=" + str(window), "t=" + str(radius ** 2)],
                           timeout=timeout)), None


def run(args, timeout=None):
    try:
        output = subprocess.run(args, stdout=subprocess.PIPE, timeout=timeout).stdout.decode("utf-8")
    except subprocess.TimeoutExpired:
        return None
    return output


def parseOutput(output):
    if output is None:
        return None
    motifs = []
    currentMotif = []
    for line in output.splitlines():
        isNext = line.startswith("next")
        isEnd = line.startswith("end")
        if isNext or isEnd:
            if len(currentMotif) > 0:
                motifs.append(currentMotif)
                currentMotif = []
            if isEnd:
                return motifs
        else:
            currentMotif.append(int(line))
    if len(currentMotif) > 0:
        motifs.append(currentMotif)
    return motifs


def main():
    parser = argparse.ArgumentParser(description="Run the latent motif discovery benchmark.")
    parser.add_argument("benchmarkPath", type=Path, help="relative path to the benchmark")
    parser.add_argument("outputFile", type=Path, help="relative path to the output file")
    parser.add_argument("--sep", default=";", help="separator in the output csv file (default: ';')")
    parser.add_argument("--algs", default=",".join(algorithms.keys()),
                        help="names of the algorithms to be tested, separated by ',' (default: "
                             "all available algorithms: " + ", ".join(algorithms.keys()) + ")")
    parser.add_argument("--radius", type=float, default=1.0, help="input radius multiplier (default: 1.0)")
    parser.add_argument("--threads", type=int, default=1, help="number of threads (default: 1)")
    parser.add_argument("--timeout", type=float, default=None, help="maximum time for a single time series in seconds"
                                                                    " (default: None)")
    args = parser.parse_args()
    return runBenchmark(args.benchmarkPath.absolute(), args.outputFile.absolute(), args.sep, args.algs.split(","),
                        args.radius, args.threads, args.timeout)


def runBenchmark(benchmarkPath, outFilePath, separator, algorithmNames, radiusMultiplier, threads, timeout):
    if outFilePath.exists():
        print("Error: output file", outFilePath, "already exists")
        return 1

    def output(*results):
        with open(outFilePath, "a") as outFile:
            outFile.write(separator.join(results) + "\n")

    # write header
    output(*["TS Path", "Meta Path", "Input Range", "Algorithm", "Runtime", "Found Motifs", "Stats"])

    # run benchmark parallelized
    pool = ParallelRun(__runAlgorithm, threads, log=print, outputFunction=output)

    # walk benchmark directory recursively
    for directory, tsName, tsMetaName in getBenchmarkFiles(benchmarkPath):
        tsDir = Path(directory)
        tsPath = tsDir / tsName
        tsMetaPath = tsDir / tsMetaName

        # get meta information
        info = getMetaInformation(tsMetaPath)
        length = int(info["length"])
        window = int(info["window"])
        motifRange = float(info["range"])
        inputRange = motifRange * radiusMultiplier

        for name in algorithmNames:
            pool.run((name, tsPath, tsMetaPath, length, window, inputRange, timeout))

    pool.join()
    return 0


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


def __runAlgorithm(name, tsPath, tsMetaPath, length, window, inputRange, timeout):
    algorithm = algorithms[name]

    # run algorithm
    startTime = time.time()
    motifsFound, stats = algorithm(tsPath, length, window, inputRange, timeout)
    runtime = time.time() - startTime

    # return results
    return [str(tsPath), str(tsMetaPath), str(inputRange), name, str(runtime), str(motifsFound), str(stats)]


if __name__ == '__main__':
    sys.exit(main())

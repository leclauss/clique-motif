import os
import re
import resource
import subprocess
import sys
import time
import argparse
from pathlib import Path

from util import ParallelRun

path = Path(__file__).parent.absolute()

# algorithm settings
algDir = path / "../tsgenerator/benchmark/algorithms"
algorithms = {
    "CliqueMotif": lambda tsPath, _, window, radius: cliqueMotif(tsPath, window, radius),
    "EMMA": lambda tsPath, _, window, radius: emma(tsPath, window, radius),
    "GrammarViz3.0": lambda tsPath, _, window, __: grammarViz(tsPath, window),
    "ScanMK": lambda tsPath, _, window, radius: scanMk(tsPath, window, radius),
    "SetFinder": lambda tsPath, _, window, radius: setFinder(tsPath, window, radius),
    "ClusterMK": lambda tsPath, _, window, radius: clusterMk(tsPath, window, radius),
    "LearnMotifs": lambda tsPath, length, window, radius: learnMotifs(tsPath, length, window, radius)
}


def cliqueMotif(tsPath, window, radius):
    return ["python3", "cliqueMotif.py", "-b", str(window), str(radius), str(tsPath)]


def emma(tsPath, window, radius):
    return ["java", "-jar", algDir / "EMMA.jar", str(tsPath), str(window), str(2 * radius), "6", "4", "1.0"]


def grammarViz(tsPath, window):
    return ["java", "-jar", algDir / "GrammarViz3.jar", "-d", str(tsPath), "-w", str(window), "-p", "6", "-a", "4"]


def scanMk(tsPath, window, radius):
    return ["java", "-jar", algDir / "ScanMK.jar", "d=" + str(tsPath), "w=" + str(window), "k=1", "r=" + str(radius)]


def setFinder(tsPath, window, radius):
    return ["java", "-jar", algDir / "SetFinder.jar", "d=" + str(tsPath), "w=" + str(window), "k=1", "r=" + str(radius)]


def clusterMk(tsPath, window, radius):
    return ["java", "-jar", algDir / "ClusterMK.jar", "d=" + str(tsPath), "w=" + str(window), "k=1", "r=" + str(radius)]


def learnMotifs(tsPath, length, window, radius):
    return ["java", "-jar", algDir / "LearnMotifs.jar", "dataSet=" + str(tsPath), "eta=0.1",
            "maxIter=1000", "numRandomRestarts=200", "alpha=2", "K=3", "pct=1",
            "tsLength=" + str(length), "w=" + str(window), "t=" + str(radius ** 2)]


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
    parser.add_argument("--memory", type=int, default=None, help="maximum memory per thread (default: None)")
    args = parser.parse_args()
    return runBenchmark(args.benchmarkPath.absolute(), args.outputFile.absolute(), args.sep, args.algs.split(","),
                        args.radius, args.threads, args.timeout, args.memory)


def runBenchmark(benchmarkPath, outFilePath, separator, algorithmNames, radiusMultiplier, threads, timeout, maxMemory):
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
            pool.run((name, tsPath, tsMetaPath, length, window, inputRange, timeout, maxMemory))

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


def __runAlgorithm(name, tsPath, tsMetaPath, length, window, inputRange, timeout, maxMemory):
    algorithm = algorithms[name]

    # run algorithm
    startTime = time.time()
    returnCode, output = run(algorithm(tsPath, length, window, inputRange), timeout=timeout, maxMemory=maxMemory)
    runtime = time.time() - startTime
    motifsFound, stats = parseOutput(output, infoOnly=(returnCode != 0))

    # return results
    return [str(tsPath), str(tsMetaPath), str(inputRange), name, str(runtime), str(motifsFound), str(stats)]


def run(args, timeout=None, maxMemory=None):
    def setLimits():
        resource.setrlimit(resource.RLIMIT_AS, (maxMemory, maxMemory))

    preFunc = None
    if maxMemory is not None:
        if args[0] == "java":
            args.insert(1, "-Xmx" + str(maxMemory))
        else:
            preFunc = setLimits

    if timeout is not None:
        args = ["timeout", str(timeout)] + args

    proc = subprocess.run(args, preexec_fn=preFunc, stdout=subprocess.PIPE)
    return proc.returncode, proc.stdout.decode("utf-8")


def parseOutput(output, infoOnly=False):
    motifs = []
    currentMotif = []
    infoLines = []
    endFound = infoOnly
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
    return (None if infoOnly else motifs), infoLines


if __name__ == '__main__':
    sys.exit(main())

import argparse
import itertools
import sys
import os
import subprocess
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="A latent motif discovery algorithm.")
    parser.add_argument("window", type=int, help="motif window")
    parser.add_argument("radius", type=float, help="motif radius")
    parser.add_argument("tsPath", type=Path, help="relative path to the time series")
    parser.add_argument("-v", "--verbose", action='store_true', help="relative path to the time series")
    parser.add_argument("-b", "--benchmark", action='store_true', help="relative path to the time series")
    args = parser.parse_args()
    log = print if args.verbose else doNothing

    motifIndices, stats = getTopMotif(args.window, args.radius, args.tsPath, log)
    if motifIndices is None:
        return 1

    log("top latent range motif size:", len(motifIndices))
    for index in motifIndices:
        print(index)
    if args.benchmark:
        print("end")
        for stat in stats:
            print(stat)
    return 0


def doNothing(*_, **__):
    pass


def loadTS(tsPath):
    with open(tsPath) as tsFile:
        ts = [float(line) for line in tsFile]
    return ts


def getTopMotif(windowSize, radius, tsPath, log=doNothing):
    # create distance graph
    log("creating graph... ", end="", flush=True)
    graphStartTime = time.time()
    mtx, nodeCount, edgeCount = createGraphSCAMP(tsPath, windowSize, 2 * radius)
    graphTime = time.time() - graphStartTime
    if mtx is None:
        log("failed (SCAMP error)")
        return None, (nodeCount, edgeCount, 0, graphTime, 0.0)
    log("done (" + str(graphTime) + " s)")
    mtxSize = len(mtx)
    log("nodes:", nodeCount, ", edges:", edgeCount, "(file size:", str(mtxSize), "B)")

    # find maximum clique
    log("running LMC... ", end="", flush=True)
    cliqueStartTime = time.time()
    motifIndices = findCliqueLMC(mtx, tsPath)
    cliqueTime = time.time() - cliqueStartTime
    if motifIndices is None:
        log("failed (LMC error)")
        return None, (nodeCount, edgeCount, mtxSize, graphTime, cliqueTime)
    log("done (" + str(cliqueTime) + " s)")

    return motifIndices, (nodeCount, edgeCount, mtxSize, graphTime, cliqueTime)


def createGraphSCAMP(tsPath, windowSize, radius, cpu_workers=1):
    # TODO multiple threads don't work (scamp output is not correct)
    ts = loadTS(tsPath)
    length = len(ts)
    nodes = length - windowSize + 1
    correlation_threshold = (1 - radius ** 2 / 2 / windowSize)
    # TODO replace subprocess call with pyscamp api (threshold parameter currently does not work):
    # pyscamp.selfjoin_knn(ts, windowSize, threshold=correlation_threshold)
    proc = subprocess.run(
        ["../algorithms/graph/scamp/build/SCAMP", "--window=" + str(windowSize), "--input_a_file_name=" + str(tsPath),
         "--no_gpu", "--num_cpu_workers=" + str(cpu_workers), "--threshold=" + str(correlation_threshold),
         "--output_a_file_name=/dev/null", "--output_a_index_file_name=/dev/null"], stdout=subprocess.PIPE)
    if proc.returncode != 0:
        return None, 0, 0

    edgeList = []
    for line in proc.stdout.decode("utf-8").splitlines():
        nums = line.split(" ")
        a = int(nums[0])
        b = int(nums[1])
        if b - a >= windowSize:
            edgeList.append(str(a + 1) + " " + str(b + 1))

    edgeCount = len(edgeList)
    mtxHeader = ["%%MatrixMarket matrix coordinate pattern symmetric", "%",
                 str(nodes) + " " + str(nodes) + " " + str(edgeCount)]
    mtx = "\n".join(itertools.chain(mtxHeader, edgeList))
    return mtx, nodes, edgeCount


def findCliqueLMC(mtx, tsPath):
    graphPath = Path(tsPath).parent.absolute() / ("distanceGraph-" + str(hash(tsPath))[:8] + ".mtx")
    os.mkfifo(graphPath)  # create named pipe (TODO does not work on Windows)
    try:
        # start LMC
        maxCliqueProcess = subprocess.Popen(["../algorithms/clique/lmc/LMC", graphPath], stdout=subprocess.PIPE)
        # write MTX twice (LMC needs to read twice)
        with open(graphPath, "w") as graphFile:
            graphFile.write(mtx)
            graphFile.write("\nR\n")
            graphFile.write(mtx)
        # wait for LMC
        stdout, _ = maxCliqueProcess.communicate()
        if maxCliqueProcess.returncode != 1:
            return None
    finally:
        os.remove(graphPath)  # remove pipe
    output = stdout.decode("utf-8")
    motifIndices = []
    for line in output.splitlines():
        if line.startswith("M "):
            indices = line.split(" ")
            for index in indices:
                try:
                    motifIndex = int(index) - 1
                    motifIndices.append(motifIndex)
                except ValueError:
                    pass
            break
    return motifIndices


if __name__ == "__main__":
    sys.exit(main())

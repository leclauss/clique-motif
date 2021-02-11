import sys
import os
import subprocess
import time
from pathlib import Path


def main():
    if len(sys.argv) < 4:
        print("Usage: python3 cliqueMotif.py windowSize radius path/to/ts.csv [-v]")
        return 1
    windowSize = int(sys.argv[1])
    radius = float(sys.argv[2])
    tsPath = sys.argv[3]
    log = print if "-v" in sys.argv else doNothing
    # TODO add timeout parameter, catch motifIndices=None below

    motifIndices = getTopMotif(windowSize, radius, tsPath, log)[0]

    log("top latent range motif size:", len(motifIndices))
    for index in motifIndices:
        print(index)
    return 0


def doNothing(*_, **__):
    pass


def loadTS(tsPath):
    with open(tsPath) as tsFile:
        ts = [float(line) for line in tsFile]
    return ts


def getTopMotif(windowSize, radius, tsPath, log=doNothing, timeout=None):
    # create named pipe
    graphPath = Path(tsPath).parent.absolute() / ("distanceGraph-" + str(hash(tsPath))[:8] + ".mtx")
    os.mkfifo(graphPath)

    try:
        # start maximum clique algorithm (pipe has to be opened for reading first)
        maxCliqueProcess = subprocess.Popen(["../algorithms/clique/lmc/LMC", graphPath], stdout=subprocess.PIPE)

        # create distance graph
        log("creating graph... ", end="", flush=True)
        startTime = time.time()
        try:
            nodeCount, edgeCount = createGraphSCAMP(tsPath, windowSize, 2 * radius, graphPath, timeout=timeout)
        except subprocess.TimeoutExpired:
            log("failed (timeout)")
            return None, (0, 0, max(0.0, timeout), 0.0)
        graphTime = time.time() - startTime
        log("done (" + str(graphTime) + " s)")
        log("nodes:", nodeCount, ", edges:", edgeCount)

        # find maximum clique
        log("running LMC... ", end="", flush=True)
        startTime = time.time()
        try:
            stdout, _ = maxCliqueProcess.communicate(timeout=(timeout - graphTime) if timeout is not None else None)
        except subprocess.TimeoutExpired:
            log("failed (timeout)")
            return None, (nodeCount, edgeCount, graphTime, max(0.0, timeout - graphTime))
    finally:
        os.remove(graphPath)
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

    cliqueTime = time.time() - startTime
    log("done (" + str(cliqueTime) + " s)")

    return motifIndices, (nodeCount, edgeCount, graphTime, cliqueTime)


def createGraphSCAMP(tsPath, windowSize, radius, graphPath, cpu_workers=1, timeout=None):
    # TODO multiple threads don't work (scamp output is not correct)
    ts = loadTS(tsPath)
    length = len(ts)
    nodes = length - windowSize + 1
    correlation_threshold = (1 - radius ** 2 / 2 / windowSize)
    # TODO replace subprocess call with pyscamp api (threshold parameter currently does not work):
    # pyscamp.selfjoin_knn(ts, windowSize, threshold=correlation_threshold)
    output = subprocess.run(
        ["../algorithms/graph/scamp/build/SCAMP", "--window=" + str(windowSize), "--input_a_file_name=" + tsPath,
         "--no_gpu", "--num_cpu_workers=" + str(cpu_workers), "--threshold=" + str(correlation_threshold),
         "--output_a_file_name=/dev/null", "--output_a_index_file_name=/dev/null"], stdout=subprocess.PIPE,
        timeout=timeout).stdout.decode("utf-8")

    mtx = ""
    edgeCount = 0
    for line in output.splitlines():
        nums = line.split(" ")
        a = int(nums[0])
        b = int(nums[1])
        if abs(a - b) >= windowSize:
            edgeCount += 1
            mtx += "\n" + str(a + 1) + " " + str(b + 1)
    mtx = "%%MatrixMarket matrix coordinate pattern symmetric\n%\n" + str(nodes) + " " + str(nodes) \
          + " " + str(edgeCount) + mtx
    with open(graphPath, "w") as graphFile:
        graphFile.write(mtx)
    return nodes, edgeCount


if __name__ == "__main__":
    sys.exit(main())

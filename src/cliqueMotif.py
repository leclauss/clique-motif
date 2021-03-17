import argparse
import itertools
import signal
import sys
import os
import subprocess
import time
from pathlib import Path
from threading import Lock


def main():
    parser = argparse.ArgumentParser(description="A latent motif discovery algorithm.")
    parser.add_argument("window", type=int, help="motif window")
    parser.add_argument("radius", type=float, help="motif radius")
    parser.add_argument("tsPath", type=Path, help="relative path to the time series")
    parser.add_argument("-v", "--verbose", action='store_true', help="relative path to the time series")
    parser.add_argument("-b", "--benchmark", action='store_true', help="relative path to the time series")
    args = parser.parse_args()
    log = print if args.verbose else doNothing

    motifIndices, stats = getTopMotif(args.window, args.radius, args.tsPath, log, True)

    if motifIndices is not None:
        log("top latent range motif size:", len(motifIndices))
        for index in motifIndices:
            print(index)
        if args.benchmark:
            print("end")
    if args.benchmark:
        for stat in stats:
            print(stat)

    return 0 if motifIndices is not None else 1


def doNothing(*_, **__):
    pass


def loadTS(tsPath):
    with open(tsPath) as tsFile:
        ts = [float(line) for line in tsFile]
    return ts


def getTopMotif(windowSize, radius, tsPath, log=doNothing, handleSignals=False):
    signalHandler = SignalHandler() if handleSignals else None

    # create distance graph
    log("creating graph... ", end="", flush=True)
    graphStartTime = time.time()
    mtx, nodeCount, edgeCount = createGraphSCAMP(tsPath, windowSize, 2 * radius, _signalHandler=signalHandler)
    graphTime = time.time() - graphStartTime
    if mtx is None:
        log("failed (SCAMP error)")
        return None, (nodeCount, edgeCount, 0, graphTime, 0.0)
    log("done (" + str(graphTime) + " s)")
    mtxSize = len(mtx)
    log("nodes:", nodeCount, ", edges:", edgeCount, "(file size:", str(mtxSize), "B)")

    if edgeCount == 0:
        return [], (nodeCount, edgeCount, mtxSize, graphTime, 0.0)

    # find maximum clique
    log("running LMC... ", end="", flush=True)
    cliqueStartTime = time.time()
    motifIndices = findCliqueLMC(mtx, tsPath, _signalHandler=signalHandler)
    cliqueTime = time.time() - cliqueStartTime
    if motifIndices is None:
        log("failed (LMC error)")
    else:
        log("done (" + str(cliqueTime) + " s)")

    return motifIndices, (nodeCount, edgeCount, mtxSize, graphTime, cliqueTime)


def createGraphSCAMP(tsPath, windowSize, radius, cpu_workers=1, _signalHandler=None):
    # TODO multiple threads don't work (scamp output is not correct)
    ts = loadTS(tsPath)
    length = len(ts)
    nodes = length - windowSize + 1
    correlation_threshold = (1 - radius ** 2 / 2 / windowSize)

    # TODO replace subprocess call with pyscamp api (threshold parameter currently does not work):
    # pyscamp.selfjoin_knn(ts, windowSize, threshold=correlation_threshold)
    proc = startSubprocess(["../algorithms/graph/scamp/build/SCAMP", "--window=" + str(windowSize),
                            "--input_a_file_name=" + str(tsPath), "--no_gpu", "--num_cpu_workers=" + str(cpu_workers),
                            "--threshold=" + str(correlation_threshold), "--output_a_file_name=/dev/null",
                            "--output_a_index_file_name=/dev/null"], _signalHandler)
    if proc is None:
        return None, 0, 0

    returnCode, output = waitSubprocess(proc, _signalHandler)
    if returnCode != 0:
        return None, 0, 0

    edgeList = []
    for line in output.splitlines():
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


def findCliqueLMC(mtx, tsPath, _signalHandler=None):
    graphPath = Path(tsPath).parent.absolute() / ("distanceGraph-" + str(hash(tsPath))[:8] + ".mtx")
    os.mkfifo(graphPath)  # create named pipe (TODO does not work on Windows)
    try:
        # start LMC
        maxCliqueProcess = startSubprocess(["../algorithms/clique/lmc/LMC", graphPath], _signalHandler)
        if maxCliqueProcess is None:
            return None
        # write MTX twice (LMC needs to read twice)
        try:
            with open(graphPath, "w") as graphFile:
                graphFile.write(mtx)
                graphFile.write("\nR\n")
                graphFile.write(mtx)
        except BrokenPipeError:
            return None
        finally:
            # wait for LMC
            returnCode, output = waitSubprocess(maxCliqueProcess, _signalHandler)
        if returnCode != 1:
            return None
    finally:
        os.remove(graphPath)  # remove pipe
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


class SignalHandler:
    def __init__(self):
        self.lock = Lock()
        self.stopped = False
        self.processes = []
        signal.signal(signal.SIGINT, self.exitGracefully)
        signal.signal(signal.SIGTERM, self.exitGracefully)

    def getLock(self):
        return self.lock

    def exitGracefully(self, _signum, _frame):
        with self.lock:
            self.stopped = True
            for proc in self.processes:
                proc.kill()

    def isStopped(self):
        return self.stopped

    def addSubprocess(self, proc):
        self.processes.append(proc)

    def removeSubprocess(self, proc):
        self.processes.remove(proc)


def startSubprocess(args, signalHandler=None):
    if signalHandler is None:
        proc = subprocess.Popen(args, stdout=subprocess.PIPE)
    else:
        with signalHandler.getLock():
            if signalHandler.isStopped():
                return None
            proc = subprocess.Popen(args, stdout=subprocess.PIPE)
            signalHandler.addSubprocess(proc)
    return proc


def waitSubprocess(proc, signalHandler=None):
    stdout, _ = proc.communicate()
    if signalHandler is not None:
        with signalHandler.getLock():
            signalHandler.removeSubprocess(proc)
    return proc.returncode, stdout.decode("utf-8")


if __name__ == "__main__":
    sys.exit(main())

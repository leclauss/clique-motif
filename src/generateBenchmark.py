import subprocess
import sys
from pathlib import Path
import shutil

path = Path(__file__).parent.absolute()

# settings
tsgeneratorPath = path / "../tsgenerator/tsgenerator/build/tsgenerator"
benchmarkPath = path / "benchmark"

repeat = 250
generator = "latent motif"
defaultLength = 4000
defaultWindow = 30
defaultSize = 9
defaultNoise = 2.0
defaultRangeFactor = 1.0
defaultMethod = "boundedNormalRandomWalk"
defaultType = "box"

lengths = [2000, 4000, 8000, 16000, 32000, 64000]
windows = []  # [30, 35, 40, 45]
sizes = []  # [3, 5, 7, 9]
noises = []  # [2.0, 2.25, 2.5, 2.75]
rangeFactors = []  # [1.0, 1.5, 2.0, 2.5]
methods = []  # ["boundedNormalRandomWalk", "linearRandomWalk", "piecewiseLinearRandom", "splineRepeated"]
motifTypes = []  # ["box", "semicircle", "positiveflank", "sine"]

defaultDelta = 1.0
defaultStep = 1.0
defaultHeight = 10.0
defaultTimes = 3
defaultMaximum = 20.0


def main(args):
    if Path.is_dir(benchmarkPath):
        if "-c" in args or "--clean" in args:
            print("Cleaning Benchmark Directories")
            shutil.rmtree(benchmarkPath)
        else:
            print("Error: benchmark directory already exists.\nRun with --clean to remove it.")
            return 1

    print("Building Benchmark")
    Path.mkdir(benchmarkPath)
    generateBenchmark()
    return 0


def generateBenchmark():
    for method in methods:
        createTimeSeries("variableMethod", method=method)

    for motifType in motifTypes:
        createTimeSeries("variableType", motifType=motifType)

    for noise in noises:
        createTimeSeries("variableNoise", noise=noise)

    for length in lengths:
        createTimeSeries("variableLength", length=length)

    for size in sizes:
        createTimeSeries("variableSize", size=size)

    for window in windows:
        createTimeSeries("variableWindow", window=window)

    for rangeFactor in rangeFactors:
        createTimeSeries("variableRangeFactor", factor=rangeFactor)


def createTimeSeries(folderName, rep=repeat, gen=generator, length=defaultLength, window=defaultWindow,
                     size=defaultSize, noise=defaultNoise, factor=defaultRangeFactor, method=defaultMethod,
                     motifType=defaultType):
    outDir = folderName + "/method_" + method + "_type_" + motifType + "_noise_" + str(noise) \
             + "_length_" + str(length) + "_size_" + str(size) + "_window_" + str(window) + "_factor_" + str(factor)
    outPath = benchmarkPath / outDir
    Path.mkdir(outPath, parents=True)

    delta = defaultDelta
    step = defaultStep
    height = defaultHeight
    times = defaultTimes
    maximum = defaultMaximum

    if method in ["splineRepeated", "piecewiseLinearRandom"]:
        times = 10
        delta = 20.0
        step = 20.0
    elif method in ["linearRandomWalk", "boundedLinearRandomWalk"]:
        delta = 5.0
        step = 10.0

    for i in range(rep):
        args = [tsgeneratorPath,
                "--generator", gen,
                "--length", str(length),
                "--window", str(window),
                "--delta", str(delta),
                "--noise", str(noise),
                "--type", motifType,
                "--size", str(size),
                "--height", str(height),
                "--step", str(step),
                "--times", str(times),
                "--Method", method,
                "--maximum", str(maximum)]
        print("Generating time series " + str(i) + " with", *args[1:])
        while subprocess.run(args, cwd=outPath).returncode != 0:
            print("retry")


if __name__ == '__main__':
    sys.exit(main(sys.argv))

import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path
import shutil

path = Path(__file__).parent.absolute()

# settings
tsgeneratorPath = path / "../tsgenerator/tsgenerator/build/tsgenerator"
benchmarkPath = path / "benchmark"
threads = 8

repeat = 1
generator = "latent motif"
defaultLength = 4000
defaultWindow = 30
defaultSize = 9
defaultNoise = 2.0
defaultMethod = "boundedNormalRandomWalk"
defaultType = "box"

lengths = [2000, 4000, 8000, 16000, 32000, 64000]  # [4000, 5000, 6000, 7000]
sizes = []  # [3, 5, 7, 9]
sizesPercent = [0.0016, 0.0024, 0.0032, 0.004]  # size proportional to length
windows = [20, 30, 40, 50]  # [30, 35, 40, 45]
noises = [1.5, 2.0, 2.5, 3.0]  # [2.0, 2.25, 2.5, 2.75]
methods = ["boundedNormalRandomWalk", "linearRandomWalk", "piecewiseLinearRandom", "splineRepeated"]
motifTypes = ["box", "semicircle", "positiveflank", "sine"]

defaultDelta = 1.0
defaultStep = 1.0
defaultHeight = 10.0
defaultTimes = 3
defaultMaximum = 20.0


def main(args):
    if benchmarkPath.is_dir():
        if "-c" in args or "--clean" in args:
            print("Cleaning Benchmark Directories")
            shutil.rmtree(benchmarkPath)
        else:
            print("Error: benchmark directory already exists.\nRun with --clean to remove it.")
            return 1
    exhaustive = True if ("-e" in args or "--exhaustive" in args) else False

    print("Building Benchmark")
    benchmarkPath.mkdir(parents=True)
    generateBenchmark(exhaustive)
    return 0


def generateBenchmark(exhaustive):
    poolArgs = []

    if exhaustive:
        for length in lengths:
            for size in getAllSizes(length):
                for window in windows:
                    for noise in noises:
                        for method in methods:
                            for motifType in motifTypes:
                                poolArgs.append(("", repeat, generator, length, window, size, noise, method, motifType))
    else:
        for length in lengths:
            poolArgs.append(("variableLength", repeat, generator, length, defaultWindow, defaultSize,
                             defaultNoise, defaultMethod, defaultType))

        for size in getAllSizes(defaultLength):
            poolArgs.append(("variableSize", repeat, generator, defaultLength, defaultWindow, size,
                             defaultNoise, defaultMethod, defaultType))

        for window in windows:
            poolArgs.append(("variableWindow", repeat, generator, defaultLength, window, defaultSize,
                             defaultNoise, defaultMethod, defaultType))

        for noise in noises:
            poolArgs.append(("variableNoise", repeat, generator, defaultLength, defaultWindow, defaultSize,
                             noise, defaultMethod, defaultType))

        for method in methods:
            poolArgs.append(("variableMethod", repeat, generator, defaultLength, defaultWindow, defaultSize,
                             defaultNoise, method, defaultType))

        for motifType in motifTypes:
            poolArgs.append(("variableType", repeat, generator, defaultLength, defaultWindow, defaultSize,
                             defaultNoise, defaultMethod, motifType))

    with Pool(threads) as threadPool:
        threadPool.starmap(createTimeSeries, poolArgs)


def getAllSizes(length):
    return sorted(set([round(sizePercent * length) for sizePercent in sizesPercent] + sizes))


def createTimeSeries(folderName="", rep=repeat, gen=generator, length=defaultLength, window=defaultWindow,
                     size=defaultSize, noise=defaultNoise, method=defaultMethod, motifType=defaultType):
    if folderName != "":
        folderName += "/"
    outDir = folderName + "method_" + method + "_type_" + motifType + "_noise_" + str(noise) \
             + "_length_" + str(length) + "_size_" + str(size) + "_window_" + str(window)
    outPath = benchmarkPath / outDir
    outPath.mkdir(parents=True)

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

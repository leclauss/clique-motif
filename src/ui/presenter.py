import math
from pathlib import Path
from PyQt5.QtCore import QObject, pyqtSignal, QThread
from scipy import stats

from scripts.getStats import getMotifPoints
from scripts.plot import getRanges
from ui.model import Model
from ui.view import View
from cliqueMotif import loadTS, getTopMotif


class Presenter:
    def __init__(self, view: View, model: Model):
        self.view = view
        self.model = model

        self.runThread = None
        self.worker = None

    def openFile(self):
        fileName = self.view.openFileNameDialog()
        self.loadData(fileName)

    def loadData(self, filePath):
        if not self.model.running and filePath is not None:  # TODO feedback if running; check if it is a file
            self.model.tsPath = filePath
            self.model.ts = loadTS(filePath)  # TODO catch exceptions

            self.view.plotTs(self.model.ts, Path(filePath).name)
            self.view.clearMotifPlot()
            self.view.clearMotifTable()

            self.view.setProgress(0, "ready")
            self.view.setRunnable(True)

    def setRunning(self, running):
        self.view.setRunning(running)
        self.model.running = running

    def run(self):
        if not self.model.running:
            # run
            self.setRunning(True)
            self.view.setProgress(0, "started")
            # create thread
            self.runThread = QThread()
            self.worker = Worker(self.model.tsPath, range(20, 50), 0.97)
            self.worker.moveToThread(self.runThread)
            # setup
            self.runThread.started.connect(self.worker.run)
            self.worker.finished.connect(self.runThread.quit)
            self.worker.finished.connect(self.worker.deleteLater)
            self.runThread.finished.connect(self.runThread.deleteLater)
            # connect signals and slots
            self.worker.progress[int].connect(lambda value: self.view.setProgress(value=value))
            self.worker.progress[str].connect(lambda text: self.view.setProgress(text=text))
            self.worker.progress[int, str].connect(lambda value, text: self.view.setProgress(value, text))
            self.worker.motif.connect(lambda motif, window: self.addMotif(motif, window))
            self.runThread.finished.connect(self.finishRun)
            # start thread
            self.runThread.start()
        else:
            # cancel
            self.setRunning(False)
            self.view.setProgress(text="cancelled")
            self.runThread.join()  # TODO stop

    def addMotif(self, motif, window):
        motifId = self.view.addMotifRow(window, 0.97, len(motif))
        self.model.motifs[motifId] = (motif, window)

    def finishRun(self):
        self.view.setProgress(100, "finished")
        self.setRunning(False)

    def showMotif(self):
        motifId = self.view.getSelectedMotif()
        if motifId is not None:
            motif, window = self.model.motifs[motifId]

            ranges = getRanges(getMotifPoints(motif, window))
            self.view.showRanges(ranges)

            subsequences = []
            mean = [0] * window
            for i in motif:
                subsequence = stats.zscore(self.model.ts[i:i + window])
                subsequences.append(subsequence)
                mean += subsequence
            mean /= len(motif)
            self.view.plotMotif(mean, subsequences)

    def save(self):
        pass  # TODO


class Worker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal([int], [str], [int, str])
    motif = pyqtSignal(list, int)

    def __init__(self, tsPath, windows, correlation):
        super().__init__()
        self.tsPath = tsPath
        self.windows = windows
        self.correlation = correlation

    def run(self):
        runs = len(self.windows)
        run = 0
        for window in self.windows:
            radius = math.sqrt((1 - self.correlation) * window / 2)
            motif, _ = getTopMotif(window, radius, self.tsPath)
            if len(motif) > 0:
                self.motif.emit(motif, window)
            run += 1
            self.progress[int].emit(int(100.0 * run / runs))
        self.finished.emit()

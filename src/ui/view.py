from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt

from ui.window import Ui_MainWindow


# main window
class View(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # create ts plot canvas
        self.figureTs = plt.figure()
        self.axesTs = self.figureTs.add_subplot(111)
        self.canvasTs = FigureCanvas(self.figureTs)
        self.toolbarTs = NavigationToolbar(self.canvasTs, self)
        self.ui.frameTsPlot.layout().addWidget(self.canvasTs)
        self.ui.frameTsPlot.layout().addWidget(self.toolbarTs)
        self.drawnRanges = []
        self.clearTsPlot()

        # create motif plot canvas
        self.figureMotif = plt.figure()
        self.axesMotif = self.figureMotif.add_subplot(111)
        self.canvasMotif = FigureCanvas(self.figureMotif)
        self.toolbarMotif = NavigationToolbar(self.canvasMotif, self)
        self.ui.frameMotifPlot.layout().addWidget(self.canvasMotif)
        self.ui.frameMotifPlot.layout().addWidget(self.toolbarMotif)
        self.clearMotifPlot()

        # setup motif table
        self.ui.tableMotifs.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        self.ui.tableMotifs.setColumnHidden(0, True)

        # create popup menu
        widgetAction = QtWidgets.QWidgetAction(self.ui.toolButton)
        widgetAction.setDefaultWidget(self.ui.widgetAdvanced)
        widgetMenu = Menu(self.ui.toolButton)
        widgetMenu.addAction(widgetAction)
        self.ui.toolButton.setMenu(widgetMenu)

    # methods for file drag and drop
    @staticmethod
    def checkDragDropEvent(event):
        urls = event.mimeData().urls()
        if urls and len(urls) == 1 and urls[0].scheme() == 'file':
            return True
        else:
            return False

    def dragEnterEvent(self, event):
        if self.checkDragDropEvent(event):
            event.acceptProposedAction()

    def dragMoveEvent(self, event):
        if self.checkDragDropEvent(event):
            event.acceptProposedAction()

    @staticmethod
    def dropEventFunction(event, function):
        if View.checkDragDropEvent(event):
            filepath = str(event.mimeData().urls()[0].path())
            function(filepath)

    # connect to the presenter
    def connect(self, presenter):
        self.ui.buttonLoad.clicked.connect(presenter.openFile)
        self.dropEvent = lambda e: self.dropEventFunction(e, presenter.loadData)
        self.ui.buttonRun.clicked.connect(presenter.run)
        self.ui.buttonSave.clicked.connect(presenter.save)
        self.ui.tableMotifs.itemSelectionChanged.connect(presenter.showMotif)

    def openFileNameDialog(self):
        title = ""
        fileTypes = "CSV (*.csv *.txt);;All Files (*)"
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self, title, "", fileTypes)
        if fileName:
            return fileName
        else:
            return None

    def addMotifRow(self, length, correlation, occurrences):
        row = self.ui.tableMotifs.rowCount()
        rowContents = [row, length, correlation, occurrences]
        self.ui.tableMotifs.insertRow(row)
        for i, content in enumerate(rowContents):
            item = QtWidgets.QTableWidgetItem()
            item.setData(Qt.DisplayRole, content)
            self.ui.tableMotifs.setItem(row, i, item)
        return row

    def getSelectedMotif(self):
        selectedRow = self.ui.tableMotifs.currentRow()
        if selectedRow is None or selectedRow < 0:
            return None
        return self.ui.tableMotifs.item(selectedRow, 0).data(Qt.DisplayRole)

    def clearMotifTable(self):
        self.ui.tableMotifs.setRowCount(0)

    def setProgress(self, value=None, text=None):
        if text is not None:
            self.ui.progressBar.setFormat("%p% (Status: " + text + ")")
        if value is not None:
            self.ui.progressBar.setValue(value)

    def setRunnable(self, runnable):
        self.ui.buttonRun.setEnabled(runnable)

    def setRunning(self, running):
        text = "Cancel" if running else "Run"
        self.ui.buttonRun.setText(text)
        self.ui.buttonLoad.setEnabled(not running)

    def drawTs(self):
        self.figureTs.tight_layout()
        self.canvasTs.draw()

    def plotTs(self, ts, path=None):
        self.axesTs.cla()
        self.axesTs.plot(ts, '-')
        self.drawTs()
        text = "Time Series"
        if path is not None:
            text += " (" + str(path) + ")"
        self.ui.labelTs.setText(text)

    def showRanges(self, ranges):
        for r in self.drawnRanges:
            r.remove()
        self.drawnRanges = []
        for start, end in ranges:
            self.drawnRanges.append(self.axesTs.axvspan(start - 0.5, end + 0.5, alpha=0.2, facecolor='green'))
        self.drawTs()

    def clearTsPlot(self):
        self.drawnRanges = []
        self.axesTs.cla()
        self.drawTs()

    def drawMotif(self):
        self.figureMotif.tight_layout()
        self.canvasMotif.draw()

    def plotMotif(self, motif, subsequences=None):
        if subsequences is None:
            subsequences = []
        self.axesMotif.cla()
        for subsequence in subsequences:
            self.axesMotif.plot(subsequence, '-', color="lightgray")
        self.axesMotif.plot(motif, '-', linewidth=3)
        self.drawMotif()

    def clearMotifPlot(self):
        self.axesMotif.cla()
        self.drawMotif()


# drop-down menu class
class Menu(QtWidgets.QMenu):
    def showEvent(self, event):
        if self.isVisible():
            button = self.parentWidget()
            if button is not None:
                pos = button.mapToGlobal(button.rect().bottomRight())
                self.move(pos - self.rect().topRight())
        super().showEvent(event)

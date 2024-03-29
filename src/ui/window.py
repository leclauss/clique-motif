# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/home/sp/work/clique-motif/src/ui/window.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1500, 1000)
        MainWindow.setAcceptDrops(True)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("/home/sp/work/clique-motif/src/ui/icon128.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.setWindowIcon(icon)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.widgetMenu = QtWidgets.QWidget(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.widgetMenu.sizePolicy().hasHeightForWidth())
        self.widgetMenu.setSizePolicy(sizePolicy)
        self.widgetMenu.setObjectName("widgetMenu")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.widgetMenu)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.buttonLoad = QtWidgets.QPushButton(self.widgetMenu)
        self.buttonLoad.setMinimumSize(QtCore.QSize(100, 30))
        self.buttonLoad.setObjectName("buttonLoad")
        self.horizontalLayout.addWidget(self.buttonLoad)
        self.buttonRun = QtWidgets.QPushButton(self.widgetMenu)
        self.buttonRun.setEnabled(False)
        self.buttonRun.setMinimumSize(QtCore.QSize(120, 30))
        self.buttonRun.setObjectName("buttonRun")
        self.horizontalLayout.addWidget(self.buttonRun)
        self.progressBar = QtWidgets.QProgressBar(self.widgetMenu)
        self.progressBar.setMinimumSize(QtCore.QSize(0, 30))
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")
        self.horizontalLayout.addWidget(self.progressBar)
        self.toolButton = QtWidgets.QToolButton(self.widgetMenu)
        self.toolButton.setMinimumSize(QtCore.QSize(30, 30))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.toolButton.setFont(font)
        self.toolButton.setStyleSheet("QToolButton::menu-indicator{width:0px;}")
        self.toolButton.setPopupMode(QtWidgets.QToolButton.InstantPopup)
        self.toolButton.setObjectName("toolButton")
        self.horizontalLayout.addWidget(self.toolButton)
        self.verticalLayout.addWidget(self.widgetMenu)
        self.widgetAdvanced = QtWidgets.QWidget(self.centralwidget)
        self.widgetAdvanced.setObjectName("widgetAdvanced")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.widgetAdvanced)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.widget = QtWidgets.QWidget(self.widgetAdvanced)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.widget.sizePolicy().hasHeightForWidth())
        self.widget.setSizePolicy(sizePolicy)
        self.widget.setObjectName("widget")
        self.gridLayout = QtWidgets.QGridLayout(self.widget)
        self.gridLayout.setObjectName("gridLayout")
        self.comboBox = QtWidgets.QComboBox(self.widget)
        self.comboBox.setObjectName("comboBox")
        self.gridLayout.addWidget(self.comboBox, 0, 1, 1, 1)
        self.checkBox = QtWidgets.QCheckBox(self.widget)
        self.checkBox.setObjectName("checkBox")
        self.gridLayout.addWidget(self.checkBox, 1, 1, 1, 1)
        self.label = QtWidgets.QLabel(self.widget)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.horizontalSlider = QtWidgets.QSlider(self.widget)
        self.horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider.setObjectName("horizontalSlider")
        self.gridLayout.addWidget(self.horizontalSlider, 2, 1, 1, 1)
        self.horizontalSlider_2 = QtWidgets.QSlider(self.widget)
        self.horizontalSlider_2.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider_2.setObjectName("horizontalSlider_2")
        self.gridLayout.addWidget(self.horizontalSlider_2, 3, 1, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.widget)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 2, 0, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.widget)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 3, 0, 1, 1)
        self.horizontalLayout_2.addWidget(self.widget)
        self.verticalLayout.addWidget(self.widgetAdvanced)
        self.widgetResults = QtWidgets.QWidget(self.centralwidget)
        self.widgetResults.setObjectName("widgetResults")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.widgetResults)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout()
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.widgetTsPlot = QtWidgets.QWidget(self.widgetResults)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(2)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.widgetTsPlot.sizePolicy().hasHeightForWidth())
        self.widgetTsPlot.setSizePolicy(sizePolicy)
        self.widgetTsPlot.setObjectName("widgetTsPlot")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.widgetTsPlot)
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_4.setSpacing(0)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.labelTs = QtWidgets.QLabel(self.widgetTsPlot)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.labelTs.sizePolicy().hasHeightForWidth())
        self.labelTs.setSizePolicy(sizePolicy)
        self.labelTs.setStyleSheet("padding: 8px 10px 7px 10px;\n"
"background-color: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #ffffff, stop:1 #e0e0e0);\n"
"border: 1px solid #b8b8b8;\n"
"border-bottom: 0px;")
        self.labelTs.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.labelTs.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.labelTs.setObjectName("labelTs")
        self.verticalLayout_4.addWidget(self.labelTs)
        self.frameTsPlot = QtWidgets.QFrame(self.widgetTsPlot)
        self.frameTsPlot.setStyleSheet("background-color: #ffffff;")
        self.frameTsPlot.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frameTsPlot.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frameTsPlot.setObjectName("frameTsPlot")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.frameTsPlot)
        self.verticalLayout_5.setContentsMargins(5, 5, 5, 5)
        self.verticalLayout_5.setSpacing(0)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.verticalLayout_4.addWidget(self.frameTsPlot)
        self.verticalLayout_7.addWidget(self.widgetTsPlot)
        self.widgetMotifPlot = QtWidgets.QWidget(self.widgetResults)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.widgetMotifPlot.sizePolicy().hasHeightForWidth())
        self.widgetMotifPlot.setSizePolicy(sizePolicy)
        self.widgetMotifPlot.setObjectName("widgetMotifPlot")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.widgetMotifPlot)
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.setSpacing(0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_5 = QtWidgets.QLabel(self.widgetMotifPlot)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_5.sizePolicy().hasHeightForWidth())
        self.label_5.setSizePolicy(sizePolicy)
        self.label_5.setStyleSheet("padding: 8px 10px 7px 10px;\n"
"background-color: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #ffffff, stop:1 #e0e0e0);\n"
"border: 1px solid #b8b8b8;\n"
"border-bottom: 0px;")
        self.label_5.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.label_5.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.label_5.setObjectName("label_5")
        self.verticalLayout_2.addWidget(self.label_5)
        self.frameMotifPlot = QtWidgets.QFrame(self.widgetMotifPlot)
        self.frameMotifPlot.setStyleSheet("background-color: #ffffff;")
        self.frameMotifPlot.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frameMotifPlot.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frameMotifPlot.setObjectName("frameMotifPlot")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout(self.frameMotifPlot)
        self.verticalLayout_6.setContentsMargins(5, 5, 5, 5)
        self.verticalLayout_6.setSpacing(0)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.verticalLayout_2.addWidget(self.frameMotifPlot)
        self.verticalLayout_7.addWidget(self.widgetMotifPlot)
        self.horizontalLayout_3.addLayout(self.verticalLayout_7)
        self.widget1 = QtWidgets.QWidget(self.widgetResults)
        self.widget1.setObjectName("widget1")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.widget1)
        self.verticalLayout_3.setContentsMargins(-1, 0, 0, 0)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.tableMotifs = QtWidgets.QTableWidget(self.widget1)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tableMotifs.sizePolicy().hasHeightForWidth())
        self.tableMotifs.setSizePolicy(sizePolicy)
        self.tableMotifs.setStyleSheet("QHeaderView::section {padding: 6px 10px;}")
        self.tableMotifs.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.tableMotifs.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.tableMotifs.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.tableMotifs.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.tableMotifs.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.tableMotifs.setVerticalScrollMode(QtWidgets.QAbstractItemView.ScrollPerPixel)
        self.tableMotifs.setHorizontalScrollMode(QtWidgets.QAbstractItemView.ScrollPerPixel)
        self.tableMotifs.setObjectName("tableMotifs")
        self.tableMotifs.setColumnCount(4)
        self.tableMotifs.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.tableMotifs.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableMotifs.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableMotifs.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableMotifs.setHorizontalHeaderItem(3, item)
        self.tableMotifs.horizontalHeader().setHighlightSections(False)
        self.tableMotifs.horizontalHeader().setSortIndicatorShown(True)
        self.tableMotifs.verticalHeader().setVisible(False)
        self.verticalLayout_3.addWidget(self.tableMotifs)
        self.buttonSave = QtWidgets.QPushButton(self.widget1)
        self.buttonSave.setObjectName("buttonSave")
        self.verticalLayout_3.addWidget(self.buttonSave)
        self.horizontalLayout_3.addWidget(self.widget1)
        self.verticalLayout.addWidget(self.widgetResults)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1500, 22))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "CliqueMotif"))
        self.buttonLoad.setText(_translate("MainWindow", "Load Data"))
        self.buttonRun.setText(_translate("MainWindow", "Run"))
        self.progressBar.setFormat(_translate("MainWindow", "%p% (Status: no data)"))
        self.toolButton.setText(_translate("MainWindow", "☰"))
        self.checkBox.setText(_translate("MainWindow", "Auto set parameters"))
        self.label.setText(_translate("MainWindow", "Algorithm"))
        self.label_2.setText(_translate("MainWindow", "Length"))
        self.label_3.setText(_translate("MainWindow", "Correlation"))
        self.labelTs.setText(_translate("MainWindow", "Time Series"))
        self.label_5.setText(_translate("MainWindow", "Motif"))
        self.tableMotifs.setSortingEnabled(True)
        item = self.tableMotifs.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "ID"))
        item = self.tableMotifs.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Length"))
        item = self.tableMotifs.horizontalHeaderItem(2)
        item.setText(_translate("MainWindow", "Correlation"))
        item = self.tableMotifs.horizontalHeaderItem(3)
        item.setText(_translate("MainWindow", "Occurrences"))
        self.buttonSave.setText(_translate("MainWindow", "Save Motifs"))

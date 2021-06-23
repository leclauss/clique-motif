import sys
from PyQt5 import QtWidgets

from ui.model import Model
from ui.presenter import Presenter
from ui.view import View


def main():
    app = QtWidgets.QApplication(sys.argv)
    view = View()
    model = Model()
    presenter = Presenter(view, model)
    view.connect(presenter)
    view.show()
    sys.exit(app.exec_())

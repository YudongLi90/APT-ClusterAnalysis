import os
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

import ui.test as main

import sys
import requests

class MainDialog(QMainWindow, main.Ui_MainWindow):

    def __init__(self, parent=None):
        super(MainDialog, self).__init__(parent)
        self.setupUi(self)
        self.pushButton_loadAPT.clicked.connect(self.button_loadAPT)

    
    @pyqtSlot(name="pushButton_loadAPT")
    def button_loadAPT(self):
        fname = QFileDialog.getOpenFileName(self,'Open APT/POS file', os.getcwd(),"Files (*.APT or *.POS)")
        self.lineEdit_filename.setText(fname[0])


if __name__ == "__main__":

    app = QApplication(sys.argv)
    form = MainDialog()
    form.show()
    app.exec_()
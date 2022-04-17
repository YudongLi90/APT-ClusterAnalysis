from ast import Pass
import os
from PyQt5.QtWidgets import QApplication, QWidget, QMainWindow, QFileDialog, QGridLayout
from PyQt5 import uic
from PyQt5.QtCore import pyqtSignal, pyqtSlot, Qt, QThread, QObject, QPointF
import pyqtgraph as pg
import pyqtgraph.opengl as gl

import ui.MainWindow as MainWindow

import sys
import requests

from aptca.core.APTData import APTData

class MainWindow(QMainWindow):
    load_data_ready = pyqtSignal(object)
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        print(os.getcwd())
        ui_file = './aptca/GUI/ui/MainWindow.ui'
        self.ui = uic.loadUi(ui_file)

        self.ui.pushButton_openAPT.clicked.connect(self.button_openAPT)
        self.ui.pushButton_openRRNG.clicked.connect(self.button_openRRNG)
        self.ui.pushButton_loadData.clicked.connect(self.button_loadData)
        self.ui.listWidget_ionSelection.itemClicked.connect(self.show_selection)
        self.ui.pushButton_filterIon.clicked.connect(self.button_filterIon)
        self.ui.pushButton_showAtoms.clicked.connect(self.button_showAtoms)
        
    @pyqtSlot(name="pushButton_openAPT")
    def button_openAPT(self):
        fname = QFileDialog.getOpenFileName(self,'Open APT/POS file', os.getcwd(),"Files (*.APT or *.POS)")
        self.ui.lineEdit_pos_filename.setText(fname[0])

    @pyqtSlot(name="pushButton_openRRNG")
    def button_openRRNG(self):
        fname = QFileDialog.getOpenFileName(self,'Open RRNG file', os.getcwd(),"Files (*.RRNG)")
        self.ui.lineEdit_RRNG_filename.setText(fname[0])

    @pyqtSlot(name="pushButton_loadData")
    def button_loadData(self):
        apt_fname = self.ui.lineEdit_pos_filename.text()
        RRNG_fname = self.ui.lineEdit_RRNG_filename.text()
        self.aptData = APTData()
        self.aptData.load_data_fromfile(apt_file=apt_fname,rrng_file=RRNG_fname)
        for ion in self.aptData.ions:
            self.ui.listWidget_ionSelection.addItem(ion)

    def show_selection(self):
        items = self.ui.listWidget_ionSelection.selectedItems()
        x = []
        for i in range(len(items)):
            x.append(str(self.ui.listWidget_ionSelection.selectedItems()[i].text()))
        for xi in x:
            self.ui.listView_ionSelected.addItem(xi)
        self.ions_selected = x
    
    @pyqtSlot(name="pushButton_filterIon")
    def button_filterIon(self):
        self.selectedAtomCloud = APTData.filter_ions(self.ions_selected)

    @pyqtSlot(name="pushButton_showAtoms")
    def button_showAtoms(self):
        Pass


    





if __name__ == "__main__":

    app = QApplication(sys.argv)
    form = MainWindow()
    form.show()
    app.exec_()
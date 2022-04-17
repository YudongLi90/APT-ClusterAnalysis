from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot
from aptca.core.APTData import APTData
from aptca.GUI.QtLogger import QtLogger
from aptca.utils.logging import log

class DataLoader(QObject):
    data_ready = pyqtSignal(object)
    finished = pyqtSignal()
    def __init__(self,apt_fname, RRNG_fname):
        super().__init__()
        self.apt_fname = apt_fname
        self.RRNG_fname = RRNG_fname
        self.aptData = APTData()
    
    @pyqtSlot(name = 'load_data')
    def load_data(self):
        log.info("Start reading files:")
        log.info(f"Position data: {self.apt_fname}")
        log.info(f"RRNG data: {self.RRNG_fname}")
        try:
            self.aptData.load_data_fromfile(apt_file=self.apt_fname,rrng_file=self.RRNG_fname)
            self.data_ready.emit(self.aptData)
        except:
            log.info("Cannot load file.")
            self.data_ready.emit(None)
        self.finished.emit()



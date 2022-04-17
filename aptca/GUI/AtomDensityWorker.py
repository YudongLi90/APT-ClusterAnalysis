from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot
from aptca.core.APTData import APTData
from aptca.GUI.QtLogger import QtLogger
from aptca.utils.logging import log

class AtomDensityWorker(QObject):
    data_ready = pyqtSignal(object)
    finished = pyqtSignal()
    def __init__(self,aptdata: APTData):
        super().__init__()
        if isinstance(aptdata, APTData):
            self.aptdata = aptdata
        else:
            raise ValueError("Please provide correct APTData object.")
    
    @pyqtSlot(name = 'calc_atomdensity')
    def calc_atomdensity(self):
        log.info("Start computing density:")
        try:
            self.aptdata.atomcloud.gen_atom_density_distribution(plot=False)
            self.data_ready.emit(self.aptdata)
        except:
            log.info("Error occured computing atom density.")
            self.data_ready.emit(None)
        self.finished.emit()
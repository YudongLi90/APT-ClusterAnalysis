from os.path import abspath, join, dirname
import sys
from PyQt5 import QtWidgets, QtCore
import logging

class QtLogger(logging.Handler,QtCore.QObject):
    appendPlainText = QtCore.pyqtSignal(str)
    def __init__(self, QPlainTextEdit):
        super().__init__()
        QtCore.QObject.__init__(self)
        self.widget = QPlainTextEdit
        self.widget.setReadOnly(True)
        self.appendPlainText.connect(self.widget.appendPlainText)
        
    
    def get_logger(self):
        _logger_form = logging.Formatter('%(asctime)s %(levelname)s %(funcName)s(%(lineno)d) %(message)s',
                                  datefmt='%d/%m/%Y %H:%M:%S')
        _log_path = join(abspath(dirname(__file__)), "..", "aptca.log")
        _file_handler = logging.FileHandler(_log_path,mode='w')
        _file_handler.setFormatter(_logger_form)
        _file_handler.setLevel(logging.DEBUG)
        log = logging.getLogger("aptca")
        log.addHandler(_file_handler)
        log.addHandler(self)
        log.setLevel(logging.INFO)
        return log

    def emit(self, record):
        msg = self.format(record)
        self.appendPlainText.emit(msg)
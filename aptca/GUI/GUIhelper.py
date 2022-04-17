from typing import List
from PyQt5 import QtCore,QtWidgets

def update_slider_and_label_range(slider:QtWidgets.QSlider,
                                     labelmin: QtWidgets.QLabel,
                                     labelmax: QtWidgets.QLabel,
                                     range: List):
    slider.setMinimum(range[0])
    slider.setMaximum(range[1])
    labelmin.setText(f"{range[0]:.2f}")
    labelmax.setText(f"{range[1]:.2f}")
    #TODO update location range min max

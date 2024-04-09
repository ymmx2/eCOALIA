__author__ = 'Maxime'
# -*- coding: utf-8 -*-

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtGui import *
from PyQt6.QtCore import *
from PyQt6.QtWidgets import *
import os
import sys
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib as mpl
mpl.rcParams['path.simplify'] = True
mpl.rcParams['path.simplify_threshold'] = 1.0
import matplotlib.style as mplstyle
mplstyle.use('fast')
from matplotlib.ticker import MultipleLocator
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg  import FigureCanvasQTAgg as FigureCanvas
import numpy as np




class LineEdit(QLineEdit):
    KEY = Qt.Key.Key_Return

    def __init__(self, *args, **kwargs):
        QLineEdit.__init__(self, *args, **kwargs)
        QREV = QRegularExpressionValidator(QRegularExpression("[+-]?\\d*[\\.]?\\d+"))
        QREV.setLocale(QLocale(QLocale.Language.English))
        self.setValidator(QREV)


class EEG_Viewer(QMainWindow):
    def __init__(self, parent=None ):
        super(EEG_Viewer, self).__init__()
        self.parent = parent

        #######################################
        self.centralWidget = QWidget()
        self.setCentralWidget(self.centralWidget)
        self.mainVBOX_param_scene = QVBoxLayout()
        # self.setMinimumHeight(700)
        # self.setMinimumWidth(800)
        self.mascene = EEG_plot(self)


        self.paramPlotV = QVBoxLayout()
        self.horizontalSliders  = QScrollBar(Qt.Orientation.Horizontal)
        self.horizontalSliders.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.horizontalSliders.valueChanged.connect(self.update_plot)
        self.horizontalSliders.setMinimum(0)
        self.horizontalSliders.setMaximum(1)



        self.paramPlot = QHBoxLayout()
        l_gain = QLabel('Gain')
        self.e_gain = LineEdit('1')
        l_win = QLabel('Window')
        self.e_win = LineEdit('10')
        l_spacing = QLabel('vertical spacing')
        self.e_spacing = LineEdit('1')
        l_linewidth = QLabel('linewidth')
        self.e_linewidth = LineEdit('1')
        self.e_gain.returnPressed.connect(self.update_plot)
        self.e_win.returnPressed.connect(self.udpate_plot_plus_slider)
        self.e_spacing.returnPressed.connect(self.update_plot)
        self.e_linewidth.returnPressed.connect(self.update_plot)

        self.paramPlot.addWidget(l_gain)
        self.paramPlot.addWidget(self.e_gain)
        self.paramPlot.addWidget(l_win)
        self.paramPlot.addWidget(self.e_win)
        self.paramPlot.addWidget(l_spacing)
        self.paramPlot.addWidget(self.e_spacing)
        self.paramPlot.addWidget(l_linewidth)
        self.paramPlot.addWidget(self.e_linewidth)


        self.paramPlotV.addWidget(self.horizontalSliders)
        self.paramPlotV.addLayout(self.paramPlot)




        self.mainVBOX_param_scene.addWidget(self.mascene)
        self.mainVBOX_param_scene.addLayout(self.paramPlotV)

        self.centralWidget.setLayout(self.mainVBOX_param_scene)


    def updateslider(self):
        self.horizontalSliders.setMinimum(0)
        self.horizontalSliders.setMaximum(int(np.ceil(self.t[-1]/int(self.e_win.text()) )-1))
        self.horizontalSliders.setPageStep(1)
        self.horizontalSliders.update()
        # self.horizontalSliders.setPageStep(self.horizontalSliders.width()/ (int(self.Sigs_dict['t'][-1]/10)+1))

    def udpate_plot_plus_slider(self):
        self.updateslider()
        self.mascene.update()

    def update_plot(self):
        self.mascene.update()

    def update(self, Sigs_dict,Sigs_Color,LFP_Names,t):
        self.Sigs_dict = Sigs_dict
        self.LFP_Names = LFP_Names
        self.Sigs_Color = Sigs_Color
        self.t = t
        self.updateslider()
        self.mascene.update()


class EEG_plot(QGraphicsView):
    def __init__(self, parent=None):
        super(EEG_plot, self).__init__(parent)
        self.parent = parent
        self.setStyleSheet("border: 0px")
        self.scene = QGraphicsScene(self)
        self.setScene(self.scene)
        self.figure = Figure(facecolor='white')
        self.canvas = FigureCanvas(self.figure)
        # self.toolbar = NavigationToolbar(self.canvas, self)

        self.widget = QWidget()
        self.widget.setLayout(QVBoxLayout())
        self.widget.layout().setContentsMargins(0, 0, 0, 0)
        self.widget.layout().setSpacing(0)
        self.scroll = QScrollArea(self.widget)
        self.scroll.setWidget(self.canvas)


        self.axes = self.figure.add_subplot(111)
        self.axes.set_xlabel("Time (s)")

        # self.canvas.setGeometry(0, 0, 1500, 500)
        layout = QVBoxLayout()
        # layout.addWidget(self.toolbar)
        layout.addWidget(self.scroll)
        self.setLayout(layout)


    def update(self):
        Sigs_dict = self.parent.Sigs_dict
        Sigs_Color = self.parent.Sigs_Color
        LFP_Names = self.parent.LFP_Names
        t = self.parent.t
        Fs= int(1/(t[1]-t[0]))

        win_num = self.parent.horizontalSliders.value()


        self.figure.clear()
        self.figure.subplots_adjust(left=0.1, bottom=0.01, right=1, top=1, wspace=0.0, hspace=0.0)
        self.axes = self.figure.add_subplot(1, 1, 1)
        gain = float(self.parent.e_gain.text())
        win= float(self.parent.e_win.text())
        spacing = float(self.parent.e_spacing.text())
        linewidth = float(self.parent.e_linewidth.text())
        ts = int(win*(win_num) * Fs)
        te = ts + int(win * Fs)
        if te > len(t):
            te=len(t)

        for i, key in enumerate(LFP_Names):
            self.axes.plot(t[ts:te], gain*(Sigs_dict[i,ts:te]-np.mean(Sigs_dict[i,ts:te]))+i*spacing,color=Sigs_Color[i], linewidth=linewidth  )



        self.axes.set_yticks(np.arange(len(LFP_Names[:]))*spacing)
        self.axes.set_yticklabels(LFP_Names)
        minorLocator = MultipleLocator(1)
        minorLocator.MAXTICKS = 10000
        self.axes.xaxis.set_minor_locator(minorLocator)
        majorLocator = MultipleLocator(999999999)
        self.axes.xaxis.set_major_locator(majorLocator)
        self.axes.xaxis.grid(which='both', color='#B0B0B0', linestyle='-', linewidth=0.5)

        self.axes.autoscale(enable=True, axis='both', tight=True)

        # self.figure.set_size_inches(len(LFP_Names[:])*10, self.parent.width()/8)
        self.canvas.setGeometry(0, 0, self.parent.width()-100, int((self.parent.height()-100)*spacing))
        self.canvas.draw_idle()
        # self. canvas.show()

    def resizeEvent(self, event):
        self.update()
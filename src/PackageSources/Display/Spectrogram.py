__author__ = "Yochum Maxime"
__copyright__ = "No Copyright"
__credits__ = ["Yochum Maxime"]
__email__ = "maxime.yochum@univ-rennes1.fr"
__status__ = "Prototype"

from PyQt6.QtGui import *
from PyQt6.QtCore import *
from PyQt6.QtWidgets import *

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
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
import numpy as np
from scipy import signal


class Spectrogram_Viewer(QGraphicsView):
    def __init__(self, parent=None):
        super(Spectrogram_Viewer, self).__init__(parent)

        self.parent = parent
        self.setStyleSheet("border: 0px")
        self.scene = QGraphicsScene(self)
        self.setScene(self.scene)
        self.figure = Figure(facecolor='white')  # Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)


        self.axes = self.figure.add_subplot(211)
        self.axes.set_xlabel("Time (s)")

        self.canvas.setGeometry(0, 0, 1500, 500)
        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        self.canvas.show()
        self.newpopup=[]


    def update(self, LFPs=None, Names=None ,Fs=None , plot1D2D=None, cut=None , Fmax=None, Fseg=None , Colors=None):
        N = LFPs.shape[0]
        tp = np.arange(LFPs.shape[1]) / Fs

        cutint = int(cut * Fs)
        if cutint >= len(tp):
            cutint = len(tp) - Fs
        tp = tp[cutint:]
        LFPs = LFPs[:, cutint:]

        if plot1D2D:
            nb_line = N
            nb_column = 1
        else:
            sqrt_N = int(np.sqrt(N))
            nb_line = sqrt_N
            nb_column = int(np.ceil(N / sqrt_N))

        self.figure.clear()
        self.figure.subplots_adjust(left=0.03, bottom=0.02, right=0.99, top=0.95, wspace=0.0, hspace=0.0)

        maxval = np.max(LFPs)
        minval = np.min(LFPs)
        if minval <= maxval:
            minval = maxval - 0.01


        for l in np.arange(nb_line):
            for c in np.arange(nb_column):
                idx = (l) * nb_column + c + 1
                if idx <= N:
                    if c == 0 and l == 0:
                        ax1 = self.figure.add_subplot(nb_line, nb_column, idx)

                        f, t, Sxx = signal.spectrogram(LFPs[idx - 1, :], Fs, nperseg=int(Fs * Fseg),
                                                       noverlap=int(Fs / 8), nfft=int(Fs * Fseg))
                        maxvect = np.max(Sxx, axis=0)
                        Sxx = Sxx / maxvect
                        colormesh = ax1.pcolormesh(t, f[:np.where(f > Fmax)[0][0]],
                                                   Sxx[:np.where(f > Fmax)[0][0], :])

                        ax1.text(0.1, 0.9,   Names[idx - 1], color=Colors[idx - 1], weight='bold',
                                 ha='center', va='center', transform=ax1.transAxes,
                                 bbox=dict(facecolor=(0.2, 0.2, 0.2), alpha=0.3, edgecolor=Colors[idx - 1],
                                           boxstyle="round"))

                        ax1.ignore_existing_data_limits = True
                        ax1.update_datalim(colormesh.get_datalim(ax1.transData))
                        ax1.autoscale_view()
                    else:
                        ax = self.figure.add_subplot(nb_line, nb_column, idx, sharex=ax1)

                        f, t, Sxx = signal.spectrogram(LFPs[idx - 1, :], Fs, nperseg=int(Fs * Fseg),
                                                       noverlap=int(Fs / 8), nfft=int(Fs * Fseg))
                        maxvect = np.max(Sxx, axis=0)
                        Sxx = Sxx / maxvect
                        colormesh = ax.pcolormesh(t, f[:np.where(f > Fmax)[0][0]],
                                                   Sxx[:np.where(f > Fmax)[0][0], :])

                        ax.text(0.1, 0.9,   Names[idx - 1], color=Colors[idx - 1], weight='bold',
                                ha='center', va='center', transform=ax.transAxes,
                                bbox=dict(facecolor=(0.2, 0.2, 0.2), alpha=0.3, edgecolor=Colors[idx - 1],
                                          boxstyle="round"))

                        ax.ignore_existing_data_limits = True
                        ax.update_datalim(colormesh.get_datalim(ax.transData))
                        ax.autoscale_view()
        self.canvas.draw_idle()
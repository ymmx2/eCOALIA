__author__ = "Yochum Maxime"
__copyright__ = "No Copyright"
__credits__ = ["Yochum Maxime"]
__email__ = "maxime.yochum@univ-rennes1.fr"
__status__ = "Prototype"



from scipy import signal
import numpy as np


def iir_band_filter_EEG(ite_data, fs, ftype, order=None, lowcut=None, highcut=None,
                        zerophase=None ):
    fe = fs / 2.0
    if not lowcut == '' and not highcut == '':
        wn = [lowcut / fe, highcut / fe]
    elif not lowcut == '':
        wn = lowcut / fe
    elif not highcut == '':
        wn = highcut / fe


    z, p, k = signal.iirfilter(order, wn, btype=ftype, ftype="butter", output="zpk" )
    sos = signal.zpk2sos(z, p, k)
    ite_data = signal.sosfilt(sos, ite_data)
    return ite_data



def signalfilter_EEG(Sigs , fs, ftype='bandpass', order=None, lowcut=None, highcut=None  ):

    if not lowcut == '':
        if lowcut <= 0:
            lowcut = 1 / fs
    if not highcut == '':
        if highcut >= fs / 2:
            highcut = fs / 2 - 1

    for idx_lfp  in range(Sigs.shape[0]):
        Sigs[idx_lfp,:] = iir_band_filter_EEG(Sigs[idx_lfp,:], fs,  ftype, order=order, lowcut=lowcut, highcut=highcut )
    return Sigs
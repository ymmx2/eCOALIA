__author__ = "Yochum Maxime"
__copyright__ = "No Copyright"
__credits__ = ["Yochum Maxime"]
__email__ = "maxime.yochum@univ-rennes1.fr"
__status__ = "Prototype"

# Practical Example 1: running from a previously saved simulation file

#make import
from PyQt6.QtWidgets import QApplication
from PackageSources.Model import Cortex_Model_NeoNMM
from PackageSources.Computation.Loading import LoadSimul, Save_Simulation
from PackageSources.Computation.Generate_Signal import Plot_Generate_ParamEvol, Plot_Generate_Stim_signal
from PackageSources.Computation.Filter import signalfilter_EEG
from PackageSources.Display.EEG_Viewer import EEG_Viewer
from PackageSources.Display.Spectrogram import Spectrogram_Viewer
from PackageSources.Computation.Classes import stim_sig, ParamEvolClass
import sys
import numpy as np


def main():
    
    #  create the model
    Model = Cortex_Model_NeoNMM.Cortex(Nb_NMM=1)
    
    #load a file
    SaveFile_Name = r'SaveFiles/1NMM_alpha.txt'
    Model, List_Stim, List_ParamEvol = LoadSimul(FilePath=SaveFile_Name, Model=Model)

    # Parameter of numerical integration
    Fs = 1024
    T = 10

    # Numerical integration 
    LFPs, tp, Pulses, PSPs, ESs = Model.Compute_Time(T, Fs, Stim = List_Stim,
                                                   List_ParamEvol = List_ParamEvol,
                                                   Pre_Post = False)

    #Fitering
    LFPs = signalfilter_EEG(LFPs, Fs, ftype='bandpass', order=3, lowcut=1, highcut=80 )
    
    app = QApplication(sys.argv)
    
    # Display LFP
    ex0 = EEG_Viewer()
    ex0.setWindowTitle('LFPs')
    ex0.update(LFPs, Model.popColor, Model.popName, tp)
    ex0.showMaximized()
    
    # Display Firing Rates
    import matplotlib.colors as colors
    colors_list = list(colors._colors_full_map.values())
    ex1 = EEG_Viewer( )
    ex1.setWindowTitle('Pulses')
    ex1.update(Pulses[0], colors_list[0:Pulses.shape[1]], Model.get_Pulse_Names(), tp )
    ex1.showMaximized()
    
    # Display membrane perturbation
    ex2 = EEG_Viewer( )
    ex2.setWindowTitle('Vm(t)')
    ex2.update(ESs[0], Model.popColor, Model.popName, tp)
    ex2.showMaximized()

    # Display spectrogram
    ex3 = Spectrogram_Viewer( )
    ex3.setWindowTitle('Spectrogram_Viewer')
    ex3.update(LFPs=LFPs,Names=Model.popName, Fs=Fs, plot1D2D=False, cut=1 , Fmax=50, Fseg=0.5,  Colors=Model.popColor)
    ex3.showMaximized()

    sys.exit(app.exec())


if __name__ == '__main__':
    main()


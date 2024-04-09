__author__ = "Yochum Maxime"
__copyright__ = "No Copyright"
__credits__ = ["Yochum Maxime"]
__email__ = "maxime.yochum@univ-rennes1.fr"
__status__ = "Prototype"

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
    #create the model
    Model = Cortex_Model_NeoNMM.Cortex(Nb_NMM=3)
    
    List_Stim = []
    List_ParamEvol = []

    # #load a file
    # SaveFile_Name = r'SaveFiles/1NMM_alpha.txt'
    # Model, List_Stim, List_ParamEvol = LoadSimul(FilePath=SaveFile_Name, Model=Model)

    # Modify model parameters in a static manner
    Model.pop.A[0] = 10 # modify one parameter for NMM[0]
    
    # Modify model parameters in a dynamic manner 
    # Modify_Add_rem  the List_ParamEvol
    myParamEvol = ParamEvolClass()
    myParamEvol.NMM = [0] # NMM number where the evolution will be applyed
    myParamEvol.typeinterp = 'linear'   # The type of interpolation that will be done inbetween specified time points
                                        # choose among the scipy innterp1D list: ‘linear’, ‘nearest’, ‘nearest-up’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’, or ‘next’. ‘zero’, ‘slinear’, ‘quadratic’ and ‘cubic’
    myParamEvol.Name = 'A' # The name of the parameter that will be updated
    myParamEvol.time = [0,10,20, 50] # time points
    myParamEvol.val  = [0,5,5, 2] # value of the parameter for each time points

    List_ParamEvol.append(myParamEvol)
    
    # Modify_Add_rem  the List_Stim
    # 'Sinus','Square', 'Sawtooth','Triangle'  : need A, f0 and phi
    # 'Ramp' : need A, Aend
    # 'Sinus Ramp' : need A, Aend, f0 and phi
    # 'Square Pulse Biphasic' : need A, f0, phi, pulsewidth
    # 'Constant' : need A
    # 'Chirp' : need A, f0, f1 and phi
    # 'Chirp Pulse' : need A, f0, f1, phi, pulsewidth
    # 'Chirp Pulse Biphasic': need A, f0, f1, phi, pulsewidth
    # 'Square Pulse Rand': need A, f0 and pulsewidth
    # 'Square Pulse Biphasic Rand' :  need A, f0 and pulsewidth
    # myStim = stim_sig()
    # myStim.kind = 'Sinus' #kind of stimulation
    # myStim.times = 0. # Time the stimulation will start, every stimulation need this value
    # myStim.timee = 1. # Time the stimulation will end, every stimulation need this value
    # myStim.A = 1. # magnitude of the stimulation, every stimulation need this value
    # myStim.Aend = 2# magnitude of the stimulation at the end for ramp stimualtion
    # myStim.f0 = 10. # Frequency of the stimulation
    # myStim.f1 = 30 # Frequency of the stimulation at the end for ramp stimualtion
    # myStim.phi = 0. # Phase of the stimulation in rd
    # myStim.pulsewidth = 0.001 # length of the stimulation pulse
    # myStim.pop = [0] # Label of the NMM where the stimulation will be applyed
    # myStim.Fs = 1024 # Frequency for ploting the stimulation signals
    # List_Stim.append(myStim)

    myStim = stim_sig()
    myStim.kind = 'Square Pulse Biphasic' #kind of stimulation
    myStim.times = 0. # Time the stimulation will start, every stimulation need this value
    myStim.timee = 10. # Time the stimulation will end, every stimulation need this value
    myStim.A = 1. # magnitude of the stimulation, every stimulation need this value
    myStim.f0 = 2. # Frequency of the stimulation
    myStim.phi = 0.0 # Phase of the stimulation in rd
    myStim.pulsewidth = 0.001 # length of the stimulation pulse
    myStim.pop = [0] # Label of the NMM where the stimulation will be applyed
    myStim.Fs = 1024 # Frequency for ploting the stimulation signals
    List_Stim.append(myStim)
    
    
    myStim = stim_sig()
    myStim.kind = 'Chirp' #kind of stimulation
    myStim.times = 0. # Time the stimulation will start, every stimulation need this value
    myStim.timee = 20. # Time the stimulation will end, every stimulation need this value
    myStim.A = 1. # magnitude of the stimulation, every stimulation need this value
    myStim.f0 = 1. # Frequency of the stimulation
    myStim.f1 = 10. # Frequency of the stimulation
    myStim.phi = 0.1 # Phase of the stimulation in rd
    myStim.pop = [1, 2] # Label of the NMM where the stimulation will be applyed
    myStim.Fs = 1024 # Frequency for ploting the stimulation signals
    List_Stim.append(myStim)
    
    
    
    # # Save a simulation file
    # SaveFile_Name = r'SaveFiles/single_node_simulation.txt'
    # Save_Simulation(fileName=SaveFile_Name,stim=List_Stim,evol=List_ParamEvol ,model=Model)

    # Parameter of numerical integration
    Fs = 1024
    T = 50
    Pre_Post = False

    # Numerical integration 
    LFPs, tp, Pulses, PSPs, ESs = Model.Compute_Time(T, Fs, Stim = List_Stim,
                                                   List_ParamEvol = List_ParamEvol,
                                                   Pre_Post = Pre_Post)

    #Fitering
    LFPs = signalfilter_EEG(LFPs, Fs, ftype='bandpass', order=3, lowcut=1, highcut=80 )
    
    # Display LFP
    app = QApplication(sys.argv)
    ex = EEG_Viewer( )
    ex.setWindowTitle('LFPs')
    ex.update(LFPs, Model.popColor, Model.popName, tp)
    ex.showMaximized()

    ex4 = Spectrogram_Viewer( )
    ex4.setWindowTitle('Spectrogram_Viewer')
    ex4.update(LFPs=LFPs,Names=Model.popName, Fs=Fs, plot1D2D=False, cut=1 , Fmax=50, Fseg=0.5,  Colors=Model.popColor)
    ex4.showMaximized()
    
    Plot_Generate_ParamEvol(List_ParamEvol)

    Plot_Generate_Stim_signal(List_Stim, model=Model)

    sys.exit(app.exec())


if __name__ == '__main__':
    main()


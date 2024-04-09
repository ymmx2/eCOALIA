__author__ = "Yochum Maxime"
__copyright__ = "No Copyright"
__credits__ = ["Yochum Maxime"]
__email__ = "maxime.yochum@univ-rennes1.fr"
__status__ = "Prototype"

# Practical example 2: simulating a seizure-like activity 
#make import
from PyQt6.QtWidgets import QApplication
from PackageSources.Model import Cortex_Model_NeoNMM
from PackageSources.Computation.Loading import LoadSimul, get_electrode, Save_Simulation, LoadLeadfield
from PackageSources.Computation.Generate_Signal import Plot_Generate_ParamEvol, Plot_Generate_Stim_signal
from PackageSources.Computation.Filter import signalfilter_EEG
from PackageSources.Display.EEG_Viewer import EEG_Viewer
from PackageSources.Display.Mesh3DView import Mesh_SimpleView
from PackageSources.Display.Spectrogram import Spectrogram_Viewer
from PackageSources.Computation.Classes import stim_sig, ParamEvolClass
import sys
import numpy as np


def main():
    #create the model
    Model = Cortex_Model_NeoNMM.Cortex(Nb_NMM=67)

    #load a file
    SaveFile_Name = r'SaveFiles/67NMM_spikewave_2.txt'
    Model,List_Stim, List_ParamEvol = LoadSimul(FilePath=SaveFile_Name, Model=Model)

    #parameter of the simulation
    Fs = 1024
    T = 40
    Pre_Post = False

    ### Modify parameter statically
    Model.pop.CI1I1[24] = 0 # modify one parameter for NMM[24]
    
    ### Modify model parameters dynamically:
    myParamEvol = ParamEvolClass()
    myParamEvol.NMM = [24] # NMM number where the evolution will be applied
    myParamEvol.Name = 'CPI2' # The name of the parameter that will be updated
    myParamEvol.time = [0.0, 10.0, 11.0, 15.0, 30.0] # time points
    myParamEvol.val  = [55.0, 35.0, 10.0, 10.0, 55.0]	 # value of the parameter for each time points
    myParamEvol.typeinterp = 'linear'   # The type of interpolation that will be done in between specified time points
    List_ParamEvol.append(myParamEvol)
    
    myParamEvol = ParamEvolClass()
    myParamEvol.NMM = [24] # NMM number where the evolution will be applied
    myParamEvol.Name = 'Ba' # The name of the parameter that will be updated
    myParamEvol.time = [0.0, 10.0, 11.0, 25.0, 30.0]	 # time points
    myParamEvol.val  = [19.0, 19.0, 2.0, 5.0, 19.0]	 # value of the parameter for each time points
    myParamEvol.typeinterp = 'linear'   # The type of interpolation that will be done in between specified time points
    List_ParamEvol.append(myParamEvol)

    # Save a simulation file
    SaveFile_Name = r'SaveFiles/67NMM_interictal_to_ictal.txt'
    Save_Simulation(fileName=SaveFile_Name, stim=List_Stim, evol=List_ParamEvol, model=Model)

    # Integrate delayed system equations 
    LFPs, tp, Pulses, PPSs, ESs = Model.Compute_Time_with_delay(T, Fs, 
                                                  Stim = List_Stim,
                                                   List_ParamEvol = List_ParamEvol,
                                                   Pre_Post = False)


    # Forward problem
    
    # Choose the montage (number of electrodes) and corresponding 
    montage = 32
    if montage == 21:
        FileName_Leadfield = r"Ressources/LeadField_66x21_DESIKAN_RL.mat"
        FileName_Electrode = r"Ressources/EEG_ElectrodeNames_21_DESIKAN_66_RL.txt"
    elif montage ==32:
        FileName_Leadfield = r"Ressources/LeadField_66x32_DESIKAN_RL.mat"
        FileName_Electrode = r"Ressources/EEG_ElectrodeNames_32_DESIKAN_66_RL.txt"
    elif montage ==65:
        FileName_Leadfield = r"Ressources/LeadField_66x65_DESIKAN_RL.mat"
        FileName_Electrode = r"Ressources/EEG_ElectrodeNames_65_DESIKAN_66_RL.txt"
    elif montage ==110:
        FileName_Leadfield = r"Ressources/LeadField_66x110_DESIKAN_RL.mat"
        FileName_Electrode = r"Ressources/EEG_ElectrodeNames_110_DESIKAN_66_RL.txt"
    elif montage ==131:
        FileName_Leadfield = r"Ressources/LeadField_66x131_DESIKAN_RL.mat"
        FileName_Electrode = r"Ressources/EEG_ElectrodeNames_131_DESIKAN_66_RL.txt"
    elif montage ==200:
        FileName_Leadfield = r"Ressources/LeadField_66x200_DESIKAN_RL.mat"
        FileName_Electrode = r"Ressources/EEG_ElectrodeNames_200_DESIKAN_66_RL.txt"
    elif montage ==256:
        FileName_Leadfield = r"Ressources/LeadField_66x256_DESIKAN_RL.mat"
        FileName_Electrode = r"Ressources/EEG_ElectrodeNames_256_DESIKAN_66_RL.txt"
    elif montage ==257:
        FileName_Leadfield = r"Ressources/LeadField_66x257_DESIKAN_RL.mat"
        FileName_Electrode = r"Ressources/EEG_ElectrodeNames_257_DESIKAN_66_RL.txt"


    LeadField = LoadLeadfield(FileName = FileName_Leadfield)

    EEG_Names, EEG_Color = get_electrode(filename=FileName_Electrode)

    #Apply the lead field with the correct sources
    EEG = np.dot(LeadField,  LFPs[1:,:])

    # Create a QT app for the displays
    app = QApplication(sys.argv)
    
    # Display LFP signals
    ex1 = EEG_Viewer( )
    ex1.setWindowTitle('LFPs')
    ex1.update(LFPs, Model.popColor, Model.popName, tp)
    ex1.showMaximized()
    
    # Display LFP spectrogram 
    ex2 = Spectrogram_Viewer( )
    ex2.setWindowTitle('Spectrogram_Viewer')
    ex2.update(LFPs=LFPs, Names=Model.popName, Fs=Fs, 
               plot1D2D=False, cut=1 , Fmax=50, Fseg=0.5,  Colors=Model.popColor)
    ex2.showMaximized()
    
    # Display the 3D anatomical model with the source activity and electrodes
    ex3 = Mesh_SimpleView(parent=app, LFPs= LFPs[1:,:] , Fs=Fs, 
                          Names=Model.popName,Colors=Model.popColor, FileName = FileName_Electrode)
    ex3.showMaximized()

    # Display EEG signals
    ex4 = EEG_Viewer( )
    ex4.setWindowTitle('EEG')
    ex4.update(EEG, EEG_Color, EEG_Names, tp)
    ex4.showMaximized()

    # Display dynamic parameters
    Plot_Generate_ParamEvol(List_ParamEvol)
    
    # Display stimulation signal if any
    Plot_Generate_Stim_signal(List_Stim, model=Model)

    sys.exit(app.exec())

if __name__ == '__main__':
    main()


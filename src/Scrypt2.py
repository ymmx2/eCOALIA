__author__ = "Yochum Maxime"
__copyright__ = "No Copyright"
__credits__ = ["Yochum Maxime"]
__email__ = "maxime.yochum@univ-rennes1.fr"
__status__ = "Prototype"

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
    Model = Cortex_Model_NeoNMM.Cortex(Nb_NMM=66)

    #load a file
    SaveFile_Name = r'SaveFiles/67NMM_resting_alpha_4.txt'
    Model,List_Stim, List_ParamEvol = LoadSimul(FilePath=SaveFile_Name, Model=Model )

    #parameter of the simulation
    Fs = 1024
    T = 10
    Pre_Post = False


    # Modify_Add_rem  the List_ParamEvol
    myParamEvol = ParamEvolClass()
    myParamEvol.NMM = [0,1,2] # NMM number where the evolution will be applyed
    myParamEvol.typeinterp = 'linear'   # The type of interpolation tha twill be done inbetween specified time points
                                        # choose among the scipy innterp1D list: ‘linear’, ‘nearest’, ‘nearest-up’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’, or ‘next’. ‘zero’, ‘slinear’, ‘quadratic’ and ‘cubic’
    myParamEvol.Name = 'A' # The name of the parameter that will be updated
    myParamEvol.time = [0,10,20] # time points
    myParamEvol.val  = [0,5,2] # value of the parameter for each time points
    List_ParamEvol.append(myParamEvol)

    #Sans delai
    LFPs, tp, Pulses, PPSs, ESs = Model.Compute_Time(T, Fs, Stim = List_Stim,
                                                   List_ParamEvol = List_ParamEvol,
                                                   Pre_Post = Pre_Post)


    #Avec delais
    # LFPs, tp, Pulses, PPSs, ESs = Model.Compute_Time_with_delay(T, Fs, Stim = List_Stim,
    #                                                List_ParamEvol = List_ParamEvol,
    #                                                Pre_Post = Pre_Post)


    #Fitering
    LFPs = signalfilter_EEG(LFPs,Fs,ftype='bandpass', order=3, lowcut=1, highcut=80 )
    # ftype =  'bandpass' 'lowpass' 'highpass' 'bandstop'

    # Problème direct
    montage = 21
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

    EEG_Names,EEG_Color = get_electrode(filename=FileName_Electrode)


    #Apply the leadfield with the correct sources
    # EEG = np.dot(LeadField, np.delete(LFPs, (47), axis=0)[14:-1,:])
    # EEG = EEG / np.max(np.abs(EEG))
    EEG = np.dot(LeadField,  LFPs[1:,:])
    #Display LFP


    app = QApplication(sys.argv)
    ex = EEG_Viewer( )
    ex.setWindowTitle('LFPs')
    ex.update(LFPs, Model.popColor, Model.popName, tp)
    ex.showMaximized()

    ex4 = Spectrogram_Viewer( )
    ex4.setWindowTitle('Spectrogram_Viewer')
    ex4.update(LFPs=LFPs,Names=Model.popName, Fs=Fs, plot1D2D=False, cut=0 , Fmax=100, Fseg=2,  Colors=Model.popColor)
    ex4.showMaximized()

    Plot_Generate_ParamEvol(List_ParamEvol)

    Plot_Generate_Stim_signal(List_Stim, model=Model)


    # ex2 = Mesh_SimpleView(parent=app, LFPs=np.delete(LFPs, (47), axis=0)[14:-1,:], Fs=Fs, Names=Model.popName,Colors=Model.popColor, FileName = FileName_Electrode)
    ex2 = Mesh_SimpleView(parent=app, LFPs= LFPs[1:,:] , Fs=Fs, Names=Model.popName,Colors=Model.popColor, FileName = FileName_Electrode)
    ex2.showMaximized()

    ex3 = EEG_Viewer( )
    ex3.setWindowTitle('EEG')
    ex3.update(EEG, EEG_Color, EEG_Names, tp)
    ex3.showMaximized()

    sys.exit(app.exec())
    plt.show()

if __name__ == '__main__':
    main()


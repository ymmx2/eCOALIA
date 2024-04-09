__author__ = 'Maxime'
from PackageSources.Model import Model_NeoNMM
from PyQt6 import QtCore, QtWidgets
import numpy as np
from PackageSources.Computation.Generate_Signal import  Generate_ParamEvol_signals,Generate_Stim_signal
import copy
import time

class SenderObject(QtCore.QObject):
    something_happened = QtCore.pyqtSignal(float)

def CM_1D_unidimensional(n):
    I = np.eye(n)
    return np.roll(I, -1, axis=1)

def msg_cri(s):
    msg = QtWidgets.QMessageBox()
    msg.setIcon(QtWidgets.QMessageBox.Icon.Critical)
    msg.setText(s)
    msg.setWindowTitle(" ")
    msg.setStandardButtons(QtWidgets.QMessageBox.StandardButton.Ok)
    msg.exec()

class Cortex():#1 principal cell 3 interneurons
    updateTime = SenderObject()
    def __init__(self,ConnectivityMat=[],DelayMat=[],Nb_NMM=None):

        self.listlabel=[ "P -> P", "P -> I1", "P -> I2", "P -> I3", "P -> I4",
                        'Delay']
        self.listCMnames=['CM_P_P','CM_P_I1','CM_P_I2','CM_P_I3','CM_P_I4',
                          'DelayMat']
        self.modelName = 'Model_NeoNMM.pop_Coalia'

        if Nb_NMM==None:
            self.Nb_NMM = ConnectivityMat[0].shape[0]
        else:
            self.Nb_NMM = Nb_NMM

        self.Nb_NMM_m1 = self.Nb_NMM -1

        if not ConnectivityMat ==[]:
            self.CM_P_P     =ConnectivityMat[0]
            self.CM_P_I1    =ConnectivityMat[1]
            self.CM_P_I2    =ConnectivityMat[2]
            self.CM_P_I3    =ConnectivityMat[3]
            self.CM_P_I4    =ConnectivityMat[4]
            
        else:
            self.CM_P_P     =CM_1D_unidimensional(self.Nb_NMM)
            self.CM_P_I1    =np.zeros((self.Nb_NMM,self.Nb_NMM))
            self.CM_P_I2    =np.zeros((self.Nb_NMM,self.Nb_NMM))
            self.CM_P_I3    =np.zeros((self.Nb_NMM,self.Nb_NMM))
            self.CM_P_I4    =np.zeros((self.Nb_NMM,self.Nb_NMM))
            


        self.list_of_CM = [self.CM_P_P ,self.CM_P_I1 ,self.CM_P_I2 ,self.CM_P_I3, self.CM_P_I4  ]

        if not DelayMat ==[]:
            self.DelayMat = DelayMat
        else:
            self.DelayMat = np.zeros((self.Nb_NMM,self.Nb_NMM))
        # for Desikan
        Regions = ['TH', 'r_BSTS', 'r_CAC', 'r_CMF', 'r_CUN', 'r_ENT', 'r_FP', 'r_FUS',
         'r_IP', 'r_IT', 'r_ISTC', 'r_LOCC', 'r_LOF', 'r_LING', 'r_MOF', 'r_MT',
         'r_PARC', 'r_PARH', 'r_POPE', 'r_PORB', 'r_PTRI', 'r_PCAL', 'r_PSTC', 'r_PC',
         'r_PREC', 'r_PCUN', 'r_RAC', 'r_RMF', 'r_SF', 'r_SP', 'r_ST', 'r_SMAR',
         'r_TP', 'r_TT', 'l_BSTS', 'l_CAC', 'l_CMF', 'l_CUN', 'l_ENT', 'l_FP',
         'l_FUS', 'l_IP', 'l_IT', 'l_ISTC', 'l_LOCC', 'l_LOF', 'l_LING', 'l_MOF',
         'l_MT', 'l_PARC', 'l_PARH', 'l_POPE', 'l_PORB', 'l_PTRI', 'l_PCAL', 'l_PSTC',
         'l_PC', 'l_PREC', 'l_PCUN', 'l_RAC', 'l_RMF', 'l_SF', 'l_SP', 'l_ST',
         'l_SMAR', 'l_TP', 'l_TT']

        self.pop = Model_NeoNMM.pop_Coalia()
        self.ODE_solver = getattr(self.pop, Model_NeoNMM.get_ODE_solver()[0])
        self.pop.NbNMMs = self.Nb_NMM
        self.pop.init_vector()
        self.pop.init_vector_param()

        self.popName = []
        self.popColor = []
        for idx in range(self.Nb_NMM):
            if idx<len(Regions):
                self.popName.append(Regions[idx])
            else:
                self.popName.append('C')
            if idx ==0:
                # self.popName.append('T')
                self.popColor.append('#912f19')
            else:
                # self.popName.append('C')
                self.popColor.append('#193c91')

        self.Nb_pulses = len(self.get_Pulse_Names())
        self.Nb_ppss = len(self.get_PPS_Names())
        self.Nb_ESs = len(self.get_ExtraSigs_Names())

        self.dt = 1/1024
        self.Stop = False

    def set_connectivityMat(self,list_of_CM):
        self.CM_P_P     = list_of_CM[0]
        self.CM_P_I1    = list_of_CM[1]
        self.CM_P_I2    = list_of_CM[2]
        self.CM_P_I3    = list_of_CM[3]
        self.CM_P_I4    = list_of_CM[4]
        
        self.list_of_CM = [self.CM_P_P ,self.CM_P_I1 ,self.CM_P_I2 ,self.CM_P_I3,self.CM_P_I4  ]


    def set_DelayMat(self,DelayMat):
        self.DelayMat = DelayMat

    def set_pop_param(self,index = 0, param = 'A' ,val = 5):
        # if len(param)==1:
        setattr(self.pop[index],param ,val)
        # else:
        #     for

    def Compute_Time(self, T, Fs, Stim=[], List_ParamEvol=[] , Pre_Post=True):
        # Original CM and Delay
        # self.initvectors()
        list_of_CM_org = copy.deepcopy(self.list_of_CM)
        DelayMat_org = copy.deepcopy(self.DelayMat)

        self.pop.init_vector()

        self.pop.CM_P_P = self.CM_P_P
        self.pop.CM_P_I1 = self.CM_P_I1
        self.pop.CM_P_I2 = self.CM_P_I2
        self.pop.CM_P_I3 = self.CM_P_I3
        self.pop.CM_P_I4 = self.CM_P_I4
        
        self.pop.DelayMat = self.DelayMat

        self.pop.NonNullMat()

        tp = np.linspace(0, T, int(T * Fs) + 1)
        self.dt = tp[1] - tp[0]

        self.pop.dt = self.dt

        LFPs = np.zeros((self.Nb_NMM, len(tp)))
        Pulses = np.zeros((self.Nb_NMM, self.Nb_pulses, len(tp)))
        PPSs = np.zeros((self.Nb_NMM, self.Nb_ppss, len(tp)))
        ESs = np.zeros((self.Nb_NMM, self.Nb_ESs, len(tp)))

        self.pop.Nb_NMM_m1 = 1
        # diviser par le numbre ded NMM -1 comme le modèle de Siouar
        # self.Nb_NMM_m1 = self.Nb_NMM - 1
        # if self.Nb_NMM_m1 <= 0:
        #     self.Nb_NMM_m1 = 1


        ystim, t = Generate_Stim_signal(Stim, model=self.pop, Fs=Fs)
        if not ystim == []:
            taille = int(T * Fs) + 1
            if t[-1] >= T:
                ystim = ystim[:, 0:taille]
            elif t[-1] < T:
                x = np.zeros((ystim.shape[0], taille))
                x[:, :ystim.shape[1]] = ystim
                ystim = x

        yparamEvol = Generate_ParamEvol_signals(List_ParamEvol, T + 1 / Fs, Fs)
        try:
            for pevol in yparamEvol:
                val = getattr(self.pop, pevol[0])
                for nmm in pevol[1]:
                    val[nmm] = pevol[2][0]
                setattr(self.pop, pevol[0], val)
        except:
            self.Stop = True
            msg_cri('Unable to generate parameter evolution.\nPlease check that all evolution are correct.')
            return LFPs, tp, Pulses, PPSs, ESs

        self.pop.Pre_Post = Pre_Post

        length = len(tp)
        for i in range(length):
            if not yparamEvol == {}:
                for pevol in yparamEvol:
                    val = getattr(self.pop, pevol[0])
                    for nmm in pevol[1]:
                        val[nmm] = pevol[2][i]
                    setattr(self.pop, pevol[0], val)

            if not ystim == []:
                self.pop.Stim = ystim[:, i]



            # self.pop.rk4()  # self.do_steptime()
            self.ODE_solver()
            # self.get_vector_Pulse()
            self.pop.apply_connectivity_Mat()
            # self.apply_external_input()
            LFPs[:, i] = self.pop.LFPoutput
            Pulses[:, 0, i] = self.pop.OutSigmoidEXC
            Pulses[:, 1, i] = self.pop.OutSigmoidI1
            Pulses[:, 2, i] = self.pop.OutSigmoidI2
            Pulses[:, 3, i] = self.pop.OutSigmoidI3
            Pulses[:, 4, i] = self.pop.OutSigmoidI4
            Pulses[:, 5, i] = self.pop.OutSigmoidEXC1
                      
            PPSs[:, 0, i] = self.pop.PPSEXC
            PPSs[:, 1, i] = self.pop.PPSI1
            PPSs[:, 2, i] = self.pop.PPSI2basal
            PPSs[:, 3, i] = self.pop.PPSI2apical
            PPSs[:, 4, i] = self.pop.PPSI3
            PPSs[:, 5, i] = self.pop.PPSI4
            PPSs[:, 6, i] = self.pop.PPSEXC1
            PPSs[:, 7, i] = self.pop.PPSExtI_P_P
            ESs[:, 0, i] = self.pop.EEGoutput
                        
            if i % (Fs // 10) == 0:
                if i > 0:
                    self.updateTime.something_happened.emit(i / length)
                if self.Stop == True:
                    self.set_connectivityMat(list_of_CM_org)
                    self.DelayMat = DelayMat_org
                    return LFPs, tp, Pulses, PPSs, ESs
        self.set_connectivityMat(list_of_CM_org)
        self.DelayMat = DelayMat_org
        return LFPs, tp, Pulses, PPSs, ESs

    def Compute_Time_with_delay(self, T, Fs, Stim=[], List_ParamEvol=[] , Pre_Post=True):
        # self.initvectors()
        # Original CM and Delay
        list_of_CM_org = copy.deepcopy(self.list_of_CM)
        DelayMat_org = copy.deepcopy(self.DelayMat)

        self.pop.init_vector()

        self.pop.CM_P_P = self.CM_P_P
        self.pop.CM_P_I1 = self.CM_P_I1
        self.pop.CM_P_I2 = self.CM_P_I2
        self.pop.CM_P_I3 = self.CM_P_I3
        self.pop.CM_P_I4 = self.CM_P_I4

        self.pop.DelayMat = self.DelayMat

        self.pop.NonNullMat()

        tp = np.linspace(0, T, int(T * Fs) + 1)
        self.dt = tp[1] - tp[0]

        self.pop.dt = self.dt

        LFPs = np.zeros((self.Nb_NMM, len(tp)))
        Pulses = np.zeros((self.Nb_NMM, self.Nb_pulses, len(tp)))
        PPSs = np.zeros((self.Nb_NMM, self.Nb_ppss, len(tp)))
        ESs = np.zeros((self.Nb_NMM, self.Nb_ESs, len(tp)))

        self.pop.Nb_NMM_m1 = 1
        # diviser par le numbre ded NMM -1 comme le modèle de Siouar
        # self.pop.Nb_NMM_m1 = self.Nb_NMM - 1
        # if self.pop.Nb_NMM_m1 <= 0:
        #     self.pop.Nb_NMM_m1 = 1

        # self.init_vectors()
        # get delay matrix
        self.pop.convert_delay_in_index()


        ystim, t = Generate_Stim_signal(Stim, model=self.pop, Fs=Fs)
        if not ystim == []:
            taille = int(T * Fs) + 1
            if t[-1] >= T:
                ystim = ystim[:, 0:taille]
            elif t[-1] < T:
                x = np.zeros((ystim.shape[0], taille))
                x[:, :ystim.shape[1]] = ystim
                ystim = x

        yparamEvol = Generate_ParamEvol_signals(List_ParamEvol, T + 1 / Fs, Fs)
        try:
            for pevol in yparamEvol:
                val = getattr(self.pop, pevol[0])
                for nmm in pevol[1]:
                    val[nmm] = pevol[2][0]
                setattr(self.pop, pevol[0], val)
        except:
            self.Stop = True
            msg_cri('Unable to generate parameter evolution.\nPlease check that all evolution are correct.')
            return LFPs, tp, Pulses, PPSs, ESs

        # set position for stim (pre or post sigmoid
        self.pop.Pre_Post = Pre_Post

        length = len(tp)
        for i in range(length):
            if not yparamEvol == {}:
                for pevol in yparamEvol:
                    val = getattr(self.pop, pevol[0])
                    for nmm in pevol[1]:
                        val[nmm] = pevol[2][i]
                    setattr(self.pop, pevol[0], val)

            if not ystim == []:
                self.pop.Stim = ystim[:, i]


            # self.pop.rk4()  # self.do_steptime()
            self.ODE_solver()
            self.pop.compute_pulse_delayed()
            self.pop.apply_connectivity_Mat_delay()
            # self.apply_external_input()
            LFPs[:, i] = self.pop.LFPoutput
            Pulses[:, 0, i] = self.pop.OutSigmoidEXC
            Pulses[:, 1, i] = self.pop.OutSigmoidI1
            Pulses[:, 2, i] = self.pop.OutSigmoidI2
            Pulses[:, 3, i] = self.pop.OutSigmoidI3
            Pulses[:, 4, i] = self.pop.OutSigmoidI4
            Pulses[:, 5, i] = self.pop.OutSigmoidEXC1
                      
            PPSs[:, 0, i] = self.pop.PPSEXC
            PPSs[:, 1, i] = self.pop.PPSI1
            PPSs[:, 2, i] = self.pop.PPSI2basal
            PPSs[:, 3, i] = self.pop.PPSI2apical
            PPSs[:, 4, i] = self.pop.PPSI3
            PPSs[:, 5, i] = self.pop.PPSI4
            PPSs[:, 6, i] = self.pop.PPSEXC1
            PPSs[:, 7, i] = self.pop.PPSExtI_P_P
            ESs[:, 0, i] = self.pop.EEGoutput
            
            ESs[:, 0, i] = self.pop.EEGoutput
            if i % (Fs // 10) == 0:
                if i > 0:
                    self.updateTime.something_happened.emit(i / length)
                if self.Stop == True:
                    self.set_connectivityMat(list_of_CM_org)
                    self.DelayMat = DelayMat_org
                    return LFPs, tp, Pulses, PPSs, ESs
        self.set_connectivityMat(list_of_CM_org)
        self.DelayMat = DelayMat_org
        return LFPs, tp, Pulses, PPSs, ESs

    def Compute_EEEG(self, PPS, NNM_Pos, Normals, EEG_pos):

        from scipy.spatial import distance
        EEG = np.zeros((len(EEG_pos), PPS.shape[2]))

        for i in range(EEG_pos.shape[0]):
            Res = np.zeros(EEG.shape[1])

            for n in range(3):
                pos = EEG_pos[i, :]
                NNM_Pos_n = NNM_Pos + Normals * self.pop.sources[n][1]
                Distance_from_electrode = distance.cdist([pos, pos], NNM_Pos_n, 'euclidean')[0, :]
                # Distance_from_electrode_square = Distance_from_electrode * Distance_from_electrode
                # U = (Points - pos) / Distance_from_electrode[:, None]
                # U = np.dot(Points - pos,Normals) / np.linalg.norm(pos)
                U = np.zeros(NNM_Pos.shape[0])
                for k in range(NNM_Pos.shape[0]):
                    vect_projectOnto = NNM_Pos[k] - pos
                    projection = (vect_projectOnto * np.dot(Normals[k], vect_projectOnto) / np.dot(vect_projectOnto,
                                                                                                   vect_projectOnto))
                    norm_vect_projectOnto = np.linalg.norm(vect_projectOnto)
                    norm_projection = np.linalg.norm(projection + vect_projectOnto)
                    U[k] = norm_vect_projectOnto - norm_projection

                if n == 0:
                    for k in range(NNM_Pos.shape[0]):
                        Res = Res + ((U[k] * PPS[k,6,:]) / ((4 * np.pi * self.pop.sigma)  * Distance_from_electrode[k])) # PPSEXC1
                        Res = Res + ((U[k] * PPS[k,0,:]) / ((4 * np.pi * self.pop.sigma)  * Distance_from_electrode[k])) # PPSEXC
                elif n == 1 :
                    for k in range(NNM_Pos.shape[0]):
                        Res = Res - ((U[k] * PPS[k,2,:]) / ((4 * np.pi * self.pop.sigma)  * Distance_from_electrode[k])) # sstbasal
                        Res = Res + ((U[k] * PPS[k,3,:]) / ((4 * np.pi * self.pop.sigma)  * Distance_from_electrode[k])) # sstapical
                        Res = Res + ((U[k] * PPS[k,5,:]) / ((4 * np.pi * self.pop.sigma)  * Distance_from_electrode[k])) # NGFC/RELN
                        Res = Res - ((U[k] * PPS[k,7,:]) / ((4 * np.pi * self.pop.sigma)  * Distance_from_electrode[k])) # cortical excitatiory input

                elif n == 2 :
                    for k in range(NNM_Pos.shape[0]):
                        Res = Res - ((U[k] * PPS[k,1,:]) / ((4 * np.pi * self.pop.sigma)  * Distance_from_electrode[k])) # PPSI1

            EEG[i, :] = Res
        return EEG

    def get_NMM_Variables(self):
        return Model_NeoNMM.get_Variable_Names()

    def get_Pulse_Names(self):
        return Model_NeoNMM.get_Pulse_Names()

    def get_PPS_Names(self):
        return Model_NeoNMM.get_PPS_Names()

    def get_ExtraSigs_Names(self):
        return Model_NeoNMM.get_ExtraSigs_Names()

    def get_ODE_solver_Time(self):
        return Model_NeoNMM.get_ODE_solver_Time()

    def get_NMM_pop(self):
        return Model_NeoNMM.pop_Coalia()

    def get_z_Model_GUI(self):
        return Model_NeoNMM.get_z_Model_GUI()


if __name__ == '__main__':
    C = Cortex(Nb_NMM=10)
    C.Compute_Time(10, 1000, Stim=[], List_ParamEvol=[] , Pre_Post=True)



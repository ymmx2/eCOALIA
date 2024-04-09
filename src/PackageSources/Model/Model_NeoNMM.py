import numpy as np
import random
from PackageSources.Model.Model_NeoNMM_GUI import z_Model_GUI

def get_z_Model_GUI():
    return z_Model_GUI

def get_LFP_Name():
    return ["LFPoutput"]

def get_Pulse_Names():
    return [ 'OutSigmoidEXC', 'OutSigmoidI1',  'OutSigmoidI2', 'OutSigmoidI3' , 'OutSigmoidI4', 'OutSigmoidEXC1']

def get_PPS_Names():
    return ['PPSEXC' ,'PPSI1', 'PPSI2basal', 'PPSI2apical', 'PPSI3' , 'PPSI4', 'PPSEXC1', 'PPSExtI_P_P']

def get_ExtraSigs_Names():
    return ['Vm']

def get_ODE_solver():
    return ["Eul_Maruyama"]

def get_ODE_solver_Time():
    return ["Eul_Time_Maruyama"]

def get_Variable_Names():
    return ['A','Bb', 'Ba', 'G','D','R',
            'a1','a2','bb1','bb2', 'ba1', 'ba2', 'g1','g2','d1','d2','r1','r2',
            'CPP','CP1P','CI1P','CI2aP', 'CI2bP','CI4P',
            'CPP1',
            'CPI1','CI1I1','CI1bI1','CI2I1','CI4I1',
            'CPI1b','CI1I1b',
            'CPI2','CI3I2','CI4I2',
            'CPI3','CI2I3','CI4I3',
            'CI2I4','CI4I4','CI4bI4',
            'CI4I4b',
            'Pv0','I1v0','I2v0','I3v0','I4v0',
            'Pe0','I1e0','I2e0','I3e0','I4e0',
            'Pr0','I1r0','I2r0','I3r0','I4r0',
            'Pm','I1m','I2m','I3m','I4m',
            'Ps','I1s','I2s','I3s','I4s',
            'Pcoef','I1coef','I2coef','I3coef','I4coef',
            'k_P','k_Pp','k_I1','k_I2','k_I3','k_I4']




class pop_Coalia:
    def __init__(self,):
        self.dt = 1./2048. #initial value for dt
        self.NbODEs = 38 # number of ODE equations to solve
        self.NbNMMs = 1 # number of NMM in the model
        self.Nb_NMM_m1=1 # possibility to divide connectivity matrices by a factor, not used here.
        self.init_vector() # call init_vector
        self.init_vector_param() # call init_vector_param
        self.Pre_Post = True # if True apply stim before sigmoid function, if False apply stim after sigmoid function

    def init_vector_simple(self):
        self.dydx = np.zeros((self.NbODEs,self.NbNMMs)) # matrix for the first derivated of RK4
        self.dydx1 = np.zeros((self.NbODEs,self.NbNMMs)) # matrix for the second derivated of RK4
        self.dydx2 = np.zeros((self.NbODEs,self.NbNMMs)) # matrix for the third derivated of RK4
        self.dydx3 = np.zeros((self.NbODEs,self.NbNMMs)) # matrix for the fourth derivated of RK4
        self.y    =np.zeros((self.NbODEs,self.NbNMMs)) # matrix for the original model state of RK4
        self.yt    =np.zeros((self.NbODEs,self.NbNMMs)) # matrix for the temporary model state of RK4

    def init_vector(self):
        self.dydx = np.zeros((self.NbODEs,self.NbNMMs))# matrix for the first derivated of RK4
        self.dydx1 = np.zeros((self.NbODEs,self.NbNMMs))# matrix for the second derivated of RK4
        self.dydx2 = np.zeros((self.NbODEs,self.NbNMMs))# matrix for the third derivated of RK4
        self.dydx3 = np.zeros((self.NbODEs,self.NbNMMs))# matrix for the fourth derivated of RK4
        self.y    =np.zeros((self.NbODEs,self.NbNMMs))# matrix for the original model state of RK4
        self.yt    =np.zeros((self.NbODEs,self.NbNMMs))# matrix for the temporary model state of RK4
        self.ExtI_P_P      = np.zeros((self.NbNMMs)) #vector that receive external P to P inputs
        self.ExtI_P_I1     = np.zeros((self.NbNMMs)) #vector that receive external P to I1 inputs
        self.ExtI_P_I2     = np.zeros((self.NbNMMs)) #vector that receive external P to I2 inputs
        self.ExtI_P_I3     = np.zeros((self.NbNMMs)) #vector that receive external P to I3 inputs
        self.ExtI_P_I4     = np.zeros((self.NbNMMs)) #vector that receive external P to I4 inputs
        self.Stim = np.zeros((self.NbNMMs)) #vector that have stimulation to apply
        self.bruitP = np.zeros((self.NbNMMs)) #vector for the noise in P
        self.bruitI1 = np.zeros((self.NbNMMs))#vector for the noise in I1
        self.bruitI2 = np.zeros((self.NbNMMs))#vector for the noise in I2
        self.bruitI3 = np.zeros((self.NbNMMs))#vector for the noise in I3
        self.bruitI4 = np.zeros((self.NbNMMs))#vector for the noise in I4

        self.EEGoutput = np.zeros((self.NbNMMs)) #vector for the LFPs
        self.LFPoutput = np.zeros((self.NbNMMs)) #vector for the LFPs

        self.OutSigmoidEXC = np.zeros((self.NbNMMs))#vector for the Firing rate of corresponding pop
        self.OutSigmoidI1 = np.zeros((self.NbNMMs))#vector for the Firing rate of corresponding pop
        self.OutSigmoidI1b = np.zeros((self.NbNMMs))#vector for the Firing rate of corresponding pop
        self.OutSigmoidI2 = np.zeros((self.NbNMMs))#vector for the Firing rate of corresponding pop
        self.OutSigmoidEXC1 = np.zeros((self.NbNMMs))#vector for the Firing rate of corresponding pop
        self.OutSigmoidI3 = np.zeros((self.NbNMMs))#vector for the Firing rate of corresponding pop
        self.OutSigmoidI4 = np.zeros((self.NbNMMs))#vector for the Firing rate of corresponding pop
        self.OutSigmoidI4b = np.zeros((self.NbNMMs))#vector for the Firing rate of corresponding pop
        self.PPSEXC = np.zeros((self.NbNMMs))#vector for the PSP of corresponding pop
        self.PPSI1 = np.zeros((self.NbNMMs))#vector for the PSP of corresponding pop
        self.PPSI1b = np.zeros((self.NbNMMs))#vector for the PSP of corresponding pop
        self.PPSI2apical = np.zeros((self.NbNMMs))#vector for the PSP of corresponding pop
        self.PPSI2basal = np.zeros((self.NbNMMs))#vector for the PSP of corresponding pop
        self.PPSExtI_P_P = np.zeros((self.NbNMMs))#vector for the PSP of corresponding pop
        
        self.PPSEXC1 = np.zeros((self.NbNMMs))#vector for the PSP of corresponding pop
        self.PPSI3 = np.zeros((self.NbNMMs))#vector for the PSP of corresponding pop
        self.PPSI4 = np.zeros((self.NbNMMs))#vector for the PSP of corresponding pop
        self.PPSI4b = np.zeros((self.NbNMMs))#vector for the PSP of corresponding pop

        self.CM_P_P = np.zeros((self.NbNMMs, self.NbNMMs)) # connectivity matrix between P and P
        self.CM_P_I1 = np.zeros((self.NbNMMs, self.NbNMMs)) # connectivity matrix between P and I1
        self.CM_P_I2 = np.zeros((self.NbNMMs, self.NbNMMs)) # connectivity matrix between P and I2
        self.CM_P_I3 = np.zeros((self.NbNMMs, self.NbNMMs)) # connectivity matrix between P and I3
        self.CM_P_I4 = np.zeros((self.NbNMMs, self.NbNMMs)) # connectivity matrix between P and I4
        self.DelayMat = np.zeros((self.NbNMMs, self.NbNMMs))# delay matrix
        self.Delay_in_index_Mat = np.zeros((self.NbNMMs, self.NbNMMs), dtype=np.int32)# delay used to comput a buffer delay matrix
        self.EmptyMat = np.full((5), True)# vector to know if some conncectivity matric are all zeros to save computation time (not apply dot product is the matrix is all zero)
        self.history_pulseP = np.zeros((self.NbNMMs, 1))# delay in number of sample
        self.max_delay_index = 0# size for the buffer delay matrix
        self.pulseP_delay_Mat = np.zeros((self.NbNMMs, self.NbNMMs))# buffer delay matrix

    def init_vector_param(self):
        self.A = np.ones((self.NbNMMs)) * 5 # initialisation of all parameters values
        self.Bb = np.ones((self.NbNMMs)) * 30.
        self.Ba = np.ones((self.NbNMMs)) * 30.
        self.G = np.ones((self.NbNMMs)) * 30.
        self.D = np.ones((self.NbNMMs)) * 5.
        self.R = np.ones((self.NbNMMs)) * 5.
        self.a1 = np.ones((self.NbNMMs)) * 100.
        self.a2 = np.ones((self.NbNMMs)) * 100.
        self.bb1 = np.ones((self.NbNMMs)) * 30.
        self.bb2 = np.ones((self.NbNMMs)) * 30.
        self.ba1 = np.ones((self.NbNMMs)) * 30.
        self.ba2 = np.ones((self.NbNMMs)) * 30.
        self.g1 = np.ones((self.NbNMMs)) * 150.
        self.g2 = np.ones((self.NbNMMs)) * 150.
        self.d1 = np.ones((self.NbNMMs)) * 20.
        self.d2 = np.ones((self.NbNMMs)) * 20.
        self.r1 = np.ones((self.NbNMMs)) * 5.
        self.r2 = np.ones((self.NbNMMs)) * 5.

        self.CPP = np.ones((self.NbNMMs)) * 0.
        self.CP1P = np.ones((self.NbNMMs)) * 100.
        self.CI1P = np.ones((self.NbNMMs)) * 50.
        self.CI2aP = np.ones((self.NbNMMs)) * 20.
        self.CI2bP = np.ones((self.NbNMMs)) * 20.
        self.CI4P = np.ones((self.NbNMMs)) * 0.
        self.CPP1 = np.ones((self.NbNMMs)) * 135.
        self.CPI1 = np.ones((self.NbNMMs)) * 50.
        self.CI1I1 = np.ones((self.NbNMMs)) * 200.
        self.CI1bI1 = np.ones((self.NbNMMs)) * 0
        self.CI2I1 = np.ones((self.NbNMMs)) * 13.5
        self.CI4I1 = np.ones((self.NbNMMs)) * 0.
        self.CPI1b = np.ones((self.NbNMMs)) * 0.
        self.CI1I1b = np.ones((self.NbNMMs)) * 0.
        self.CPI2 = np.ones((self.NbNMMs)) * 20.
        self.CI3I2 = np.ones((self.NbNMMs)) * 20.
        self.CI4I2 = np.ones((self.NbNMMs)) * 0.
        self.CPI3 = np.ones((self.NbNMMs)) * 0.
        self.CI2I3 = np.ones((self.NbNMMs)) * 20.
        self.CI4I3 = np.ones((self.NbNMMs)) * 0.
        self.CI2I4 = np.ones((self.NbNMMs)) * 0.
        self.CI4I4 = np.ones((self.NbNMMs)) * 0.
        self.CI4bI4 = np.ones((self.NbNMMs)) * 0.
        self.CI4I4b = np.ones((self.NbNMMs)) * 0.

        self.Pv0 = np.ones((self.NbNMMs)) * 6.
        self.I1v0 = np.ones((self.NbNMMs)) * 6.
        self.I2v0 = np.ones((self.NbNMMs)) * 6.
        self.I3v0 = np.ones((self.NbNMMs)) * 6.
        self.I4v0 = np.ones((self.NbNMMs)) * 6.
        self.Pe0 = np.ones((self.NbNMMs)) * 5.
        self.I1e0 = np.ones((self.NbNMMs)) * 5.
        self.I2e0 = np.ones((self.NbNMMs)) * 5.
        self.I3e0 = np.ones((self.NbNMMs)) * 5.
        self.I4e0 = np.ones((self.NbNMMs)) * 5.
        self.Pr0 = np.ones((self.NbNMMs)) * 0.56
        self.I1r0 = np.ones((self.NbNMMs)) * 0.56
        self.I2r0 = np.ones((self.NbNMMs)) * 0.56
        self.I3r0 = np.ones((self.NbNMMs)) * 0.56
        self.I4r0 = np.ones((self.NbNMMs)) * 0.56
        self.Pm = np.ones((self.NbNMMs)) * 120
        self.I1m = np.ones((self.NbNMMs)) * 0.
        self.I2m = np.ones((self.NbNMMs)) * 0.
        self.I3m = np.ones((self.NbNMMs)) * 0.
        self.I4m = np.ones((self.NbNMMs)) * 0.
        self.Ps = np.ones((self.NbNMMs)) * 3.0
        self.I1s = np.ones((self.NbNMMs)) * 0.
        self.I2s = np.ones((self.NbNMMs)) * 0.
        self.I3s = np.ones((self.NbNMMs)) * 0.
        self.I4s = np.ones((self.NbNMMs)) * 0.
        self.Pcoef = np.ones((self.NbNMMs)) * 1.
        self.I1coef = np.ones((self.NbNMMs)) * 0.
        self.I2coef = np.ones((self.NbNMMs)) * 0.
        self.I3coef = np.ones((self.NbNMMs)) * 0.
        self.I4coef = np.ones((self.NbNMMs)) * 0.


        self.k_P = np.ones((self.NbNMMs)) * 1.
        self.k_Pp = np.ones((self.NbNMMs)) * 1.
        self.k_I1 = np.ones((self.NbNMMs)) * 1.
        self.k_I2 = np.ones((self.NbNMMs)) * 1.
        self.k_I3 = np.ones((self.NbNMMs)) * 1.
        self.k_I4 = np.ones((self.NbNMMs)) * 1.
        
        
        self.k_P = np.ones((self.NbNMMs)) * 1.
        self.k_Pp = np.ones((self.NbNMMs)) * 1.
        self.k_I1 = np.ones((self.NbNMMs)) * 1.
        self.k_I2 = np.ones((self.NbNMMs)) * 1.
        self.k_I3 = np.ones((self.NbNMMs)) * 1.

        # self.Thresh = np.ones((self.NbNMMs)) * 6.
        # self.Pv0_UD = np.ones((self.NbNMMs)) * 6.
        # self.Pe0_UD = np.ones((self.NbNMMs)) * 5.
        # self.Pr0_UD = np.ones((self.NbNMMs)) * 0.56
        # self.Pv0_used = np.ones((self.NbNMMs)) * 6.
        # self.Pe0_used = np.ones((self.NbNMMs)) * 5.
        # self.Pr0_used = np.ones((self.NbNMMs)) * 0.56

        # self.sources = np.vstack((np.array([0., 2]),np.array([0, 0.5]),np.array([0, 0.2])))#
        self.sources = np.array([[0., 2.],[0., 0.5],[0., 0.2]], dtype="float64")
        self.position = np.array([10., 2], dtype="float64")
        #self.position = np.array([1., 2], dtype=np.float64)
        self.sigma = 4e-4
        self.compute_distance()
        

    def compute_distance(self):
        self.distance = np.sqrt((self.sources[:, 0] - self.position[0]) ** 2 + (self.sources[:, 1]  - self.position[1]) ** 2)

    def random_seeded(self,seed):
        random.seed(int(seed))

    def sigmP(self,v):
        return  self.Pe0/(1+np.exp( self.Pr0*(self.Pv0-v)))

    def sigmI1(self,v):
        return  self.I1e0/(1+np.exp( self.I1r0*(self.I1v0-v)))

    def sigmI2(self,v):
        return  self.I2e0/(1+np.exp( self.I2r0*(self.I2v0-v)))

    def sigmI3(self,v):
        return  self.I3e0/(1+np.exp( self.I3r0*(self.I3v0-v)))

    def sigmI4(self,v):
        return  self.I4e0/(1+np.exp( self.I4r0*(self.I4v0-v)))

    def noiseP(self):
        for i in range(self.NbNMMs):
            self.bruitP[i] = self.Pcoef[i] * np.random.normal(self.Pm[i],self.Ps[i])

    def noiseI1(self):
        for i in range(self.NbNMMs):
            self.bruitI1[i] = self.I1coef[i] * np.random.normal(self.I1m[i],self.I1s[i])

    def noiseI2(self):
        for i in range(self.NbNMMs):
            self.bruitI2[i] = self.I2coef[i] * np.random.normal(self.I2m[i],self.I2s[i])

    def noiseI3(self):
        for i in range(self.NbNMMs):
            self.bruitI3[i] = self.I3coef[i] * np.random.normal(self.I3m[i],self.I3s[i])

    def noiseI4(self):
        for i in range(self.NbNMMs):
            self.bruitI4[i] = self.I4coef[i] * np.random.normal(self.I4m[i], self.I4s[i])

    def gainFactor_fun(self,v1, v2):
        tmax = v1*0.
        a = np.logical_and(v1 == v2, np.logical_not(v1 == 0))
        tmax[a] = 1/v1[a]
        a = np.logical_and(np.logical_not(v1 == v2), np.logical_not(v1 == 0))
        tmax[a] = np.log(v2[a] / v1[a]) / (v2[a] - v1[a])
        return np.exp(v1 * tmax - 1) * v2

    def PSP(self, y0, y1, y2, V, v1, v2):
        return (V * self.gainFactor_fun(v1, v2) * y0 - (v1 + v2) * y2 - v1 * v2 * y1)


    def Eul_Maruyama(self):
        
        self.bruitP  = self.Pcoef  * self.Pm
        self.bruitI1 = self.I1coef * self.I1m
        self.bruitI2 = self.I2coef * self.I2m
        self.bruitI3 = self.I3coef * self.I3m
        self.bruitI4 = self.I4coef * self.I4m

        self.dydx1=self.derivT()
        self.y += (self.dydx1 * self.dt)

        gainFactor = self.A * self.gainFactor_fun(self.a1, self.a2)

        self.y[17] += gainFactor * self.Pcoef * self.Ps * np.random.normal(loc=0.0, scale=np.sqrt(self.dt),size=self.NbNMMs)
        self.y[19] += gainFactor * self.I1coef * self.I1s * np.random.normal(loc=0.0, scale=np.sqrt(self.dt),size=self.NbNMMs)
        self.y[21] += gainFactor * self.I2coef * self.I2s * np.random.normal(loc=0.0, scale=np.sqrt(self.dt),size=self.NbNMMs)
        self.y[23] += gainFactor * self.I3coef * self.I3s * np.random.normal(loc=0.0, scale=np.sqrt(self.dt),size=self.NbNMMs)
        self.y[25] += gainFactor * self.I4coef * self.I4s * np.random.normal(loc=0.0, scale=np.sqrt(self.dt),size=self.NbNMMs)
        
                
        # # up-down
        # self.Pv0_used[self.OutSigmoidI2 > self.Thresh ] = self.Pv0_UD[self.OutSigmoidI2 > self.Thresh ]*1
        # self.Pv0_used[self.OutSigmoidI2 <= self.Thresh ] = self.Pv0[self.OutSigmoidI2 <= self.Thresh ]*1

        # self.Pe0_used[self.OutSigmoidI2 > self.Thresh ] = self.Pe0_UD[self.OutSigmoidI2 > self.Thresh ]*1
        # self.Pe0_used[self.OutSigmoidI2 <= self.Thresh ] = self.Pe0[self.OutSigmoidI2 <= self.Thresh ]*1

        # self.Pr0_used[self.OutSigmoidI2 > self.Thresh ] = self.Pr0_UD[self.OutSigmoidI2 > self.Thresh ]*1
        # self.Pr0_used[self.OutSigmoidI2 <= self.Thresh ] = self.Pr0[self.OutSigmoidI2 <= self.Thresh ]*1
        self.compute_LFP()


    def Eul_Time_Maruyama(self,N,stim):
        self.init_vector()
        lfp = np.zeros((N,self.NbNMMs))

        if np.sum(np.abs(stim)) == 0:
            for k in range(N):
                self.Eul_Maruyama()
                lfp[k,:]= self.EEGoutput
        else:
            for k in range(N):
                self.Stim = stim[:, k]
                self.Eul_Maruyama()
                lfp[k,:]= self.EEGoutput

        return lfp


    def derivT(self , ):

        self.EEGoutput = (  self.CP1P * self.y[8,:]
                          + self.CPP * self.y[0,:]
                          - self.CI1P * self.y[2,:]
                          - self.CI2bP * self.y[6,:]
                          - self.CI4P * self.y[12,:]
                          - self.CI2aP * self.y[36,:] 
                          + self.y[16,:]
                          + self.y[26,:]) # comput LFP signal of the P Population
        self.OutSigmoidEXC = self.sigmP(self.k_P * self.Stim * self.Pre_Post
                                        +self.EEGoutput ) # compute sigmoid of the LFP signal and apply stimulation if applyed before sigmoide


        self.dydx[0,:] = self.y[1,:]
        self.dydx[1,:] =  self.PSP(self.k_P * self.Stim * (not self.Pre_Post)
                                 + self.OutSigmoidEXC,self.y[0,:],self.y[1,:], self.A, self.a1, self.a2)# compute H transfert function on Firing rate P signal and apply stimulation if applyed after sigmoide


        self.OutSigmoidI1 = self.sigmI1( self.k_I1 * self.Stim * self.Pre_Post
                                         + self.CPI1 * self.y[0,:]
                                         - self.CI1bI1 * self.y[4,:]
                                         - self.CI1I1 * self.y[2,:]
                                         - self.CI2I1 * self.y[6,:]
                                         - self.CI4I1 * self.y[12,:]
                                         + self.y[18,:]
                                         +(self.y[28,:] ))# compute sigmoid on the firing rate of I1 population and apply stimulation if applyed before sigmoide
        self.dydx[2,:] = self.y[3,:]
        self.dydx[3,:] =  self.PSP(self.k_I1 * self.Stim * (not self.Pre_Post)
                                 + self.OutSigmoidI1 ,self.y[2,:],self.y[3,:], self.G, self.g1, self.g2)# compute H transfert function on Firing rate I1 signal and apply stimulation if applyed after sigmoide

        self.OutSigmoidI1b = self.sigmI1(   self.CPI1b*self.y[0,:]
                                          - self.CI1I1b*self.y[2,:])# compute sigmoid on the firing rate of I1b population and apply stimulation if applyed before sigmoide
        self.dydx[4,:] = self.y[5,:]
        self.dydx[5,:] =  self.PSP(self.OutSigmoidI1b,self.y[4,:],self.y[5,:], self.G, self.g1, self.g2)# compute H transfert function on Firing rate I1b signal and apply stimulation if applyed after sigmoide

        # Equations for inhibitory interneurons I2 (inputs) SST_BASAL
        self.OutSigmoidI2 = self.sigmI2(  self.k_I2 * self.Stim * self.Pre_Post
                                        + self.CPI2*self.y[0,:]
                                        - self.CI3I2*self.y[10,:]
                                        - self.CI4I2*self.y[12,:]
                                        + self.y[20,:]
                                        + self.y[30,:] )# compute sigmoid on the firing rate of I2 population and apply stimulation if applyed before sigmoide
        self.dydx[6,:] = self.y[7,:]
        self.dydx[7,:] =  self.PSP(self.k_I2 * self.Stim * (not self.Pre_Post)
                                 + self.OutSigmoidI2 ,self.y[6,:],self.y[7,:], self.Bb, self.bb1, self.bb2)# compute H transfert function on Firing rate I2 signal and apply stimulation if applyed after sigmoide


        # Equations for Pyramidal cells P1 (inputs)
        self.OutSigmoidEXC1 =  self.sigmP(self.k_Pp * self.Stim * self.Pre_Post
                                        +self.CPP1*self.y[0,:])# compute sigmoid on the firing rate of Pp population and apply stimulation if applyed before sigmoide
        self.dydx[8,:]  =  self.y[9,:]
        self.dydx[9,:]  =   self.PSP(self.k_Pp * self.Stim * (not self.Pre_Post)
                                 + self.OutSigmoidEXC1,self.y[8,:],self.y[9,:], self.A, self.a1, self.a2)# compute H transfert function on Firing rate Pp signal and apply stimulation if applyed after sigmoide


        # Equations for inhibitory interneurons I3 (inputs) VIP
        self.OutSigmoidI3 = self.sigmI3(self.k_I3 * self.Stim * self.Pre_Post
                                        + self.CPI3 * self.y[0,:]
                                        - self.CI2I3*self.y[6,:]
                                        - self.CI4I3*self.y[12,:]
                                        + self.y[22,:]
                                        + self.y[32,:] )# compute sigmoid on the firing rate of I3 population and apply stimulation if applyed before sigmoide
        self.dydx[10,:] = self.y[11,:]
        self.dydx[11,:] =  self.PSP(self.k_I3 * self.Stim * (not self.Pre_Post)
                                  + self.OutSigmoidI3 ,self.y[10,:],self.y[11,:], self.D, self.d1, self.d2)# compute H transfert function on Firing rate I3 signal and apply stimulation if applyed after sigmoide


        # Equations for inhibitory interneurons I4 (inputs) RLN
        self.OutSigmoidI4 = self.sigmI4(self.k_I4 * self.Stim * self.Pre_Post
                                        - self.CI2I4 * self.y[6, :]
                                        - self.CI4I4 * self.y[12, :]
                                        - self.CI4bI4 * self.y[14, :]
                                        + self.y[24, :]
                                        + self.y[34, :])# compute sigmoid on the firing rate of I4 population and apply stimulation if applyed before sigmoide
        self.dydx[12, :] = self.y[13, :]
        self.dydx[13, :] = self.PSP(self.k_I4 * self.Stim * (not self.Pre_Post)
                                    + self.OutSigmoidI4, self.y[12, :], self.y[13, :], self.R, self.r1, self.r2)# compute H transfert function on Firing rate I4 signal and apply stimulation if applyed after sigmoide


        # Equations for inhibitory interneurons I4b (inputs) RLNb

        self.OutSigmoidI4b = self.sigmI4(self.k_I4 * self.Stim * self.Pre_Post
                                        - self.CI4I4b * self.y[12, :] )# compute sigmoid on the firing rate of I4b population and apply stimulation if applyed before sigmoide
        self.dydx[14, :] = self.y[15, :]
        self.dydx[15, :] = self.PSP(self.k_I4 * self.Stim * (not self.Pre_Post)
                                    + self.OutSigmoidI4b, self.y[14, :], self.y[15, :], self.R, self.r1, self.r2)# compute H transfert function on Firing rate I4b signal and apply stimulation if applyed after sigmoide


   
        ###################################

        # Equations for input noise Pyramidal
        self.dydx[16,:] = self.y[17,:]
        self.dydx[17,:] = self.PSP(self.bruitP,self.y[16,:],self.y[17,:], self.A, self.a1, self.a2)
        # Equations for input noise I1
        self.dydx[18,:] = self.y[19,:]
        self.dydx[19,:] = self.PSP(self.bruitI1,self.y[18,:],self.y[19,:],self.A, self.a1, self.a2)
        # Equations for input noise I2
        self.dydx[20,:] = self.y[21,:]
        self.dydx[21,:] = self.PSP(self.bruitI2,self.y[20,:],self.y[21,:],self.A, self.a1, self.a2)
        # Equations for input noise I3
        self.dydx[22,:] = self.y[23,:]
        self.dydx[23,:] = self.PSP(self.bruitI3,self.y[22,:],self.y[23,:],self.A, self.a1, self.a2)
        # Equations for input noise I4
        self.dydx[24,:] = self.y[25,:]
        self.dydx[25,:] = self.PSP(self.bruitI4,self.y[24,:],self.y[25,:],self.A, self.a1, self.a2)

        # afférences

        self.dydx[26,:] = self.y[27,:]
        self.dydx[27,:] = self.PSP(self.ExtI_P_P,self.y[26,:],self.y[27,:],self.A, self.a1, self.a2)

        self.dydx[28,:] = self.y[29,:]
        self.dydx[29,:] = self.PSP(self.ExtI_P_I1,self.y[28,:],self.y[29,:],self.A, self.a1, self.a2)

        self.dydx[30,:] = self.y[31,:]
        self.dydx[31,:] = self.PSP(self.ExtI_P_I2,self.y[30,:],self.y[31,:],self.A, self.a1, self.a2)

        self.dydx[32,:] = self.y[33,:]
        self.dydx[33,:] = self.PSP(self.ExtI_P_I3,self.y[32,:],self.y[33,:],self.A, self.a1, self.a2)

        self.dydx[34,:] = self.y[35,:]
        self.dydx[35,:] = self.PSP(self.ExtI_P_I4,self.y[34,:],self.y[35,:],self.A, self.a1, self.a2)


        ##################################
        # Equations for interneurons I2 (inputs) SST_APICAL
        self.dydx[36,:]  =  self.y[37,:]
        self.dydx[37,:]  =   self.PSP(self.k_I2 * self.Stim * (not self.Pre_Post)
                                 + self.OutSigmoidI2, self.y[36,:],self.y[37,:], self.Ba, self.ba1, self.ba2)# compute H transfert function on Firing rate Pp signal and apply stimulation if applyed after sigmoide


        self.PPSEXC = self.CPP*self.y[0,:]  # get PSP of P population on P 
        self.PPSI1 = self.CI1P*self.y[2,:]     # get PSP of I1 population on P
        self.PPSI1b = self.y[4,:]    # get PSP of I1b population 
        self.PPSI2basal = self.CI2bP * self.y[6,:]     # get PSP of I2 basal projection on P
        self.PPSEXC1 = self.CP1P*self.y[8,:]   # get PSP of Pp population on P
        self.PPSI3 = self.y[10,:]    # get PSP of I3 population
        self.PPSI4 = self.CI4P*self.y[12,:]    # get PSP of I4 population on P
        self.PPSI4b = self.y[14,:]   # get PSP of I4b population
        self.PPSI2apical = self.CI2aP*self.y[36,:] # get PSP of I2 apical projection on P 
        self.PPSExtI_P_P = self.y[26,:] + self.y[16,:]# get PSP of External Input on PYR

        return self.dydx+0. #needed to creat a copy


    def NonNullMat(self): # compute is the matrix are all zeros or not
        if np.max(self.CM_P_P)==0.:
            self.EmptyMat[0]=False
        else:
            self.EmptyMat[0]=True
        if np.max(self.CM_P_I1)==0.:
            self.EmptyMat[1]=False
        else:
            self.EmptyMat[1]=True
        if np.max(self.CM_P_I2)==0.:
            self.EmptyMat[2]=False
        else:
            self.EmptyMat[2]=True
        if np.max(self.CM_P_I3)==0.:
            self.EmptyMat[3]=False
        else:
            self.EmptyMat[3]=True
        if np.max(self.CM_P_I4)==0.:
            self.EmptyMat[4]=False
        else:
            self.EmptyMat[4]=True

    def apply_connectivity_Mat(self):# Apply the connectivity matrix
        if self.EmptyMat[0]:
            self.ExtI_P_P	= np.dot(self.CM_P_P  / self.Nb_NMM_m1   ,self.OutSigmoidEXC)
        if self.EmptyMat[1]:
            self.ExtI_P_I1 	= np.dot(self.CM_P_I1 / self.Nb_NMM_m1   ,self.OutSigmoidEXC)
        if self.EmptyMat[2]:
            self.ExtI_P_I2 	= np.dot(self.CM_P_I2 / self.Nb_NMM_m1   ,self.OutSigmoidEXC)
        if self.EmptyMat[3]:
            self.ExtI_P_I3 	= np.dot(self.CM_P_I3 / self.Nb_NMM_m1   ,self.OutSigmoidEXC)
        if self.EmptyMat[4]:
            self.ExtI_P_I4 	= np.dot(self.CM_P_I4 / self.Nb_NMM_m1   ,self.OutSigmoidEXC)

    def apply_connectivity_Mat_delay(self):# Apply the connectivity matrix if their exist a delay (use the buffer)
        if self.EmptyMat[0]:
            self.ExtI_P_P  = np.sum(self.CM_P_P  / self.Nb_NMM_m1 * self.pulseP_delay_Mat, axis=1)
        if self.EmptyMat[1]:
            self.ExtI_P_I1 = np.sum(self.CM_P_I1 / self.Nb_NMM_m1 * self.pulseP_delay_Mat, axis=1)
        if self.EmptyMat[2]:
            self.ExtI_P_I2 = np.sum(self.CM_P_I2 / self.Nb_NMM_m1 * self.pulseP_delay_Mat, axis=1)
        if self.EmptyMat[3]:
            self.ExtI_P_I3 = np.sum(self.CM_P_I3 / self.Nb_NMM_m1 * self.pulseP_delay_Mat, axis=1)
        if self.EmptyMat[4]:
            self.ExtI_P_I4 = np.sum(self.CM_P_I4 / self.Nb_NMM_m1 * self.pulseP_delay_Mat, axis=1)


    def convert_delay_in_index(self):# compute the delay in sample from its value in second
        self.Delay_in_index_Mat = (self.DelayMat / self.dt).astype(np.int32)
        self.max_delay_index = np.max(self.Delay_in_index_Mat)+1
        self.history_pulseP = np.zeros((self.NbNMMs, self.max_delay_index))

    def compute_pulse_delayed(self):# GEt signal at previous time by using the buffer.
        if np.any(self.EmptyMat[:5]):
            self.history_pulseP[:, 1:] = self.history_pulseP[:, 0:-1]
            self.history_pulseP[:, 0]=self.OutSigmoidEXC
            for idx2 in range(self.NbNMMs):
                for idx in range(self.NbNMMs):
                    self.pulseP_delay_Mat[idx2, idx]  = self.history_pulseP[idx,self.Delay_in_index_Mat[idx2,idx]]

    ## commented to change the configuration to get BRE
    # def compute_LFP(self): 
    #     self.LFPoutput = (1 / (4 * np.pi * self.sigma)) * \
    #                      (self.PPSEXC1 * (1 / self.distance[0] - 1 / self.distance[1] ) + # par la pop P'
    #                       self.PPSEXC * (1 / self.distance[2] - 1 / self.distance[1] ) + # par la pop P
    #                       self.PPSI2basal * (1 / self.distance[1] - 1 / self.distance[0] ) +
    #                       self.PPSI2apical * ( 1 / self.distance[0]  - 1 / self.distance[1]) +
    #                       self.PPSI1 * (1 / self.distance[1] - 1 / self.distance[0] )+
    #                       self.PPSI4 * ( 1 / self.distance[0]  - 1 / self.distance[1])+
    #                       self.PPSExtI_P_P * ( 1 / self.distance[1]  - 1 / self.distance[0]))
                         
                         
    def compute_LFP(self): 
        self.LFPoutput =  (self.PPSEXC1 * np.abs(self.distance[0]- self.distance[1] ) + # par la pop P'
                          - self.PPSEXC * np.abs(self.distance[2] - self.distance[1] ) + # par la pop P
                          - self.PPSI2basal * np.abs(self.distance[1] - self.distance[0] ) +
                          self.PPSI2apical * np.abs(self.distance[0] - self.distance[1]) +
                          -self.PPSI1 * np.abs(self.distance[1] - self.distance[0] )+
                          self.PPSI4 * np.abs(self.distance[0]  - self.distance[1])+
                          - self.PPSExtI_P_P * np.abs(self.distance[1]  - self.distance[0]))


                         
if __name__ == '__main__':
    pop = pop_Coalia()
    pop.NbNMMs = 1
    pop.init_vector()
    pop.init_vector_param()
    lfp = pop.Eul_Time_Maruyama(100)
    
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(lfp[:,0])
    plt.ylabel('LFP')

    plt.show()
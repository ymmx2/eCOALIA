__author__ = 'Maxime'

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import copy
from scipy import interpolate


def Fs_max(List_Stim):
    Fs = List_Stim[0].Fs
    for Stim in List_Stim:
        if Stim.Fs>Fs:
            Fs = Stim.Fs
    return Fs

def Te_max(List_Stim):
    Te = List_Stim[0].timee
    for Stim in List_Stim:
        if Stim.timee>Te:
            Te = Stim.timee
    return Te


def Generate_Stim_signal(List_Stim, model = [], Fs=[]):
    if List_Stim == []:
        #msg_cri('No stimulation signal has been set yet.\nPlease add at least one.')
        return [],[]
    else:
        nbStim = len(List_Stim)
        #extract param



        if Fs == []:
            Fs = Fs_max(List_Stim)

        Te = Te_max(List_Stim)


        t = np.arange(Te * Fs) / Fs

        taille =len(t)
        x = np.zeros((nbStim,taille))
        for idx,Stim in enumerate(List_Stim):
            A = Stim.A
            Aend = Stim.Aend
            f0 = Stim.f0
            f1 = Stim.f1
            phi = Stim.phi
            ts = int(Stim.times * Fs)
            te = int(Stim.timee * Fs)
            type = Stim.kind
            pulseW=Stim.pulsewidth

            RiseTC = Stim.RiseTC
            FallTC = Stim.FallTC
            RiseTW = Stim.RiseTW
            FallTW = Stim.FallTW

            if type ==  'Sinus':
                x[idx,ts:te] = A * np.sin (2 * np.pi * f0 * t[ts:te] + phi)
            elif type ==  'Square':
                x[idx,ts:te] = A * signal.square(2 * np.pi * f0 * t[ts:te] + phi, duty= 1/2)
            elif type ==  'Sawtooth':
                x[idx,ts:te] = A *signal.sawtooth(2 * np.pi * f0 * t[ts:te] + phi, 0)
            elif type ==  'Triangle':
                x[idx,ts:te] = A *signal.sawtooth(2 * np.pi * f0 * t[ts:te]+phi, 0.5)
            elif type ==  'Ramp':
                x[idx,ts:te] = A +  (t[ts:te] - t[ts])/(Stim.timee-Stim.times) * (Aend-A)
            elif type ==  "Sinus Ramp":
                x[idx,ts:te] = (A +  (t[ts:te] - t[ts])/(Stim.timee-Stim.times) * (Aend-A)) * np.sin (2 * np.pi * f0 * t[ts:te] + phi)
            elif type ==  'Square Pulse':
                Fspulse = int(round(pulseW * Fs))
                s=signal.sawtooth(2 * np.pi * f0 * t[ts:te] + phi, 0)
                xp = np.zeros(len(s))
                ds =  np.diff(s)
                for i in range(len(ds)):
                    if ds[i]>0:
                        xp[i:i+Fspulse+1]=1*A
                x[idx,ts:te] = xp
            elif type ==  'Square Pulse Biphasic':
                Fspulse = int(round(pulseW * Fs/2))
                s=signal.sawtooth(2 * np.pi * f0 * t[ts:te] + phi, 0)
                xp = np.zeros(len(s))
                ds =  np.diff(s)
                for i in range(len(ds)):
                    if ds[i]>0:
                        xp[i:i+Fspulse]=1*A
                        xp[i+Fspulse:i+2*Fspulse]=-1*A
                x[idx,ts:te] = xp


            elif type ==  'Constant':
                x[idx,ts:te] = np.ones(x[idx,ts:te].shape)*A
            elif type == 'Chirp':
                x[idx,ts:te] = A * signal.chirp(t[ts:te]-t[ts], f0=f0, f1=f1, t1=t[te-1]-t[ts], method='logarithmic')

            elif type == 'Chirp Pulse':
                Fspulse = int(round(pulseW * Fs ))
                s= signal.chirp(t[ts:te]-t[ts], f0=f0, f1=f1, t1=t[te-1]-t[ts], method='logarithmic')
                s[s>=0] = 1
                s[s<0] = 0
                xp = np.zeros(len(s))
                ds =  np.diff(s)
                for i in range(len(ds)):
                    if ds[i]>0:
                        xp[i:i+Fspulse]=A
                if len(xp)>len(s):
                    xp = xp[:len(s)]
                x[idx,ts:te] = xp


            elif type == 'Chirp Pulse Biphasic':
                Fspulse = int(round(pulseW * Fs/2))
                s= signal.chirp(t[ts:te]-t[ts], f0=f0, f1=f1, t1=t[te-1]-t[ts], method='logarithmic')
                s[s>=0] = 1
                s[s<0] = 0
                xp = np.zeros(len(s))
                ds =  np.diff(s)
                for i in range(len(ds)):
                    if ds[i]>0:
                        xp[i:i+Fspulse]=A
                        xp[i+Fspulse:i+2*Fspulse]=-A
                if len(xp)>len(s):
                    xp = xp[:len(s)]
                x[idx,ts:te] = xp
            elif type == "Square Pulse Rand":
                Fspulse = int(round(pulseW * Fs/2))
                xp = np.zeros(len(t[ts:te]))
                temps = (te - ts) / Fs
                Nb = int(f0 * temps)
                pos = np.random.randint(0, te - ts, Nb)
                for i in pos:
                    if i + Fspulse + 1 < te:
                        xp[i:i + Fspulse + 1] = xp[i:i + Fspulse + 1] + 1 * A

                x[idx,ts:te] = xp

            elif type == "Square Pulse Biphasic Rand":
                Fspulse = int(round(pulseW * Fs/2))
                xp = np.zeros(len(t[ts:te]))
                temps = (te - ts) / Fs
                Nb = int(f0 * temps)
                pos = np.random.randint(0, te - ts, Nb)
                for i in pos:
                    if i+2*Fspulse < te:
                        xp[i:i + Fspulse + 1] = xp[i:i + Fspulse + 1] + 1 * A
                        xp[i + Fspulse:i + 2 * Fspulse+1] = xp[i:i + Fspulse + 1] - 1 * A

                x[idx,ts:te] = xp

        nbpop =  model.NbNMMs
        y = np.zeros((nbpop,taille))#1 signal par pop [CA1,CA3,GD,Sub, CED,CES]
        for idx,Stim in enumerate(List_Stim):
            for pop in Stim.pop:
                i_pop = int(pop)
                y[i_pop,:] += x[idx,:]

        return y,t

def Plot_Generate_Stim_signal(List_Stim, model):
    y, t = Generate_Stim_signal(List_Stim, model=model.pop)
    if not y == []:
        N = y.shape[0]

        sqrt_N = int(np.sqrt(N))
        nb_line = sqrt_N
        nb_column = int(np.ceil(N / sqrt_N))

        fig = plt.figure()

        maxval= np.max(y)
        minval= np.min(y)
        for l in np.arange(nb_line):
            for c in np.arange(nb_column):
                idx = (l)*nb_column + c +1
                if idx <= N:
                    if c==0 and l ==0:
                        ax1 =plt.subplot(nb_line,nb_column,idx)
                        plt.plot(t,y[idx-1,:])
                        ax1.set_ylim(minval, maxval)
                        ax1.text(0.1, 0.9,str(idx-1), color='white', weight='bold', ha='center', va='center', transform=ax1.transAxes,
                                 bbox=dict(facecolor=(0.2, 0.2, 0.2), edgecolor='black',boxstyle="round"))
                    else:
                        ax =plt.subplot(nb_line,nb_column,idx, sharex=ax1)
                        plt.plot(t,y[idx-1,:])
                        ax.set_ylim(minval, maxval)
                        ax.text(0.1, 0.9,str(idx-1), color='white', weight='bold', ha='center', va='center', transform=ax.transAxes,
                                 bbox=dict(facecolor=(0.2, 0.2, 0.2), edgecolor='black',boxstyle="round"))

        fig.show()




def Plot_Generate_ParamEvol(List_ParamEvol):
    if not List_ParamEvol == []:
        te = 0
        for pevol in List_ParamEvol:
            idx = sorted(range(len(pevol.time)), key=lambda k: pevol.time[k])
            if pevol.time[idx[-1]]> te:
                te = pevol.time[idx[-1]]


        fig = plt.figure()
        ax = plt.subplot(111)
        for pevol in List_ParamEvol:
            idx = sorted(range(len(pevol.time)), key=lambda k: pevol.time[k])
            t = [pevol.time[i] for i in idx]
            v = [pevol.val[i] for i in idx]

            if t[idx[-1]]< te:
                t.append(te)
                v.append(v[idx[-1]])
            if t[idx[0]]> 0.:
                t.insert(0, 0.)
                v.insert(0, v[0])
            # idx = sorted(range(len(pevol.time)), key=lambda k: pevol.time[k])
            plt.plot(t ,v,'ko')

            f = interpolate.interp1d(t,v ,kind=pevol.typeinterp)
            xnew = np.arange(np.min(t), np.max(t), 1./1024)
            plt.plot(xnew,f(xnew), label=pevol.Name +':'+str(pevol.NMM).replace('[','').replace(']',''))

        plt.ylabel('Parameter value')
        plt.xlabel('Time [s]')
        plt.legend()
        fig.show()

def Generate_ParamEvol_signals(List_ParamEvolORG, T=5, Fs=[]):
    List_ParamEvol =  copy.deepcopy(List_ParamEvolORG)
    if List_ParamEvol == []:
        #msg_cri('No stimulation signal has been set yet.\nPlease add at least one.')
        return {}
    else:
        dt = 1/Fs
        Sig=[]
        for pevol in List_ParamEvol:
            idx = sorted(range(len(pevol.time)), key=lambda k: pevol.time[k])
            t = [pevol.time[i] for i in idx]
            v = [pevol.val[i] for i in idx]

            if t[idx[-1]]< T:
                t.append(T)
                v.append(v[idx[-1]])
            if t[idx[0]]> 0.:
                t.insert(0, 0.)
                v.insert(0, v[0])

            f = interpolate.interp1d(t,v ,kind=pevol.typeinterp)
            xnew = np.arange(np.min(t), np.max(t), dt)
            Sig.append((pevol.Name,pevol.NMM,f(xnew)))
        return Sig










if __name__ == '__main__':
    import matplotlib.pyplot as plt


    T=1
    Fs = 1000
    samples = int(T*Fs)
    t = np.arange(samples) / Fs
    A_Up = 2
    A_Down = 1
    RiseTC = 0.001
    FallTC =  0.0015
    RiseTW = 0.01
    FallTW = 0.02

    # t_Rise = np.arange(int(T_Up*Fs)) / Fs
    # y_Rise = (1 - np.exp(-t_Rise/TC_Rise)) * A_Up
    # t_Decay = np.arange(int(T_Down*Fs)) / Fs
    # y_Decay = 1 - np.exp(-t_Decay/T_Down) * A_Down
    #
    #
    #
    # y_Rise2 = np.exp(-t_Rise / TC_Rise)
    # y_Rise2 = y_Rise2 - np.min(y_Rise2)
    # y_Rise2 = (y_Rise2 / np.max(y_Rise2)) * (np.max(y_Rise) + np.max(y_Decay)) - np.max(y_Decay)
    #
    # y_Decay = y_Decay - np.max(y_Decay)
    # y = np.hstack((y_Rise, y_Rise2, y_Decay))

    t_Rise = np.arange(int(RiseTW * Fs)) / Fs
    y_Rise = (1 - np.exp(-t_Rise / RiseTC)) * A_Up
    t_Decay = np.arange(int(FallTW * Fs)) / Fs
    y_Decay = 1 - np.exp(-t_Decay / FallTC) * A_Down
    y_Decay = y_Decay - np.max(y_Decay)

    y_Rise2 = np.exp(-t_Rise / RiseTC)
    y_Rise2 = y_Rise2 - np.min(y_Rise2)
    y_Rise2 = (y_Rise2 / np.max(y_Rise2)) * (np.max(y_Rise) - np.min(y_Decay)) + np.min(y_Decay)

    y = np.hstack((y_Rise[:-1], y_Rise2[:-1], y_Decay))


    plt.figure()
    plt.plot(t_Rise,y_Rise)
    plt.plot(t_Rise[-1]+t_Decay,y_Decay)
    plt.plot(y)

    plt.show()



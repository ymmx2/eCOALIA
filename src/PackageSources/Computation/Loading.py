__author__ = "Yochum Maxime"
__copyright__ = "No Copyright"
__credits__ = ["Yochum Maxime"]
__email__ = "maxime.yochum@univ-rennes1.fr"
__status__ = "Prototype"

import os
import numpy as np
import scipy.io
from PackageSources.Computation.Classes import ParamEvolClass, stim_sig
from PackageSources.Model import Cortex_Model_NeoNMM


def get_electrode(filename=None):
    with open(filename, "r") as f:
        rows = f.readlines()
        Regions_elec = []
        for row in rows:
            stripped_line = row.strip()
            line_list = stripped_line.split()
            Regions_elec.append(line_list[0])

    EEG_Color = ["#000000"] * len(Regions_elec)
    return Regions_elec,EEG_Color

def LoadSimul(FilePath = "", Model= None):
    stim, evol, model, modelname, cm, delay = Load_Simulation(fileName=FilePath, listcm=Model.listCMnames)

    if not stim == None:
        List_Stim = stim
    else:
        List_Stim = []
    if not evol == None:
        List_ParamEvol = evol
    else:
        List_ParamEvol = []
    if not model == None:
        if modelname == Model.modelName:
            names = list(model.keys())
            n_mod = len(model[names[0]])

            Model = Cortex_Model_NeoNMM.Cortex(Nb_NMM=n_mod)

            names = list(model.keys())
            for n in names:
                if n == "Name":
                    Model.popName = model[n]
                elif n == "Color":
                    Model.popColor = model[n]
                else:
                    try:
                        if np.issubdtype(type(getattr(Model.pop, n)[0]), np.integer):
                            setattr(Model.pop, n, np.array(model[n], dtype=np.int32))
                        else:
                            setattr(Model.pop, n, np.array(model[n]))
                    except:
                        pass


    if not cm == None:
        Model.set_connectivityMat(cm)
        Model.set_DelayMat(delay)

    return Model, List_Stim, List_ParamEvol




def Load_Simulation(fileName, listcm):
    stim=None
    evol=None
    model=None
    modelname=None
    cm=None
    delay=None
    with  open(fileName  , 'r') as f:
        num_characters = os.fstat(f.fileno()).st_size
        line =f.readline()
        i=0
        while (i<num_characters) or i==0:
            i=f.tell()
            if  "Stim_info::" in line:
                stim,line=read_stim(f)
            if  "Evol_info::" in line:
                evol,line=read_evol(f)
            if  "Model_info::" in line:
                model,modelname, line=read_model(f )
            if  "CM_info::" in line:
                cm,delay,line = read_cm(f,listcm)
            line =f.readline()
        f.close()
        return stim, evol , model, modelname, cm, delay


def read_stim(f):
    line =f.readline()
    nbstim = int(line.split('=')[-1])
    stim=[]
    for i in range(nbstim):
        stim.append(stim_sig())
    line =f.readline()
    while not("::" in line or line == ''):
        if not (line =='' or line == "\n"):
            lsplit = line.split("\t")
            name = lsplit[0]
            vals = lsplit[1:]

            if any(([x in name for x in ["pop"]])):
                for i in range(nbstim):
                    listval = [x.strip() for x in eval(vals[i])]
                    setattr(stim[i],name,listval)
            if any(([x in name for x in ["kind"]])):
                for i in range(nbstim):
                    setattr(stim[i],name,vals[i])
            elif any(([x in name for x in ["Fs"]])):
                for i in range(nbstim):
                    setattr(stim[i],name,int(vals[i]))
            elif any(([x in name for x in ["Add_Edit","dVm"]])):
                pass
            else:
                for i in range(nbstim):
                    try:
                        setattr(stim[i],name,float(vals[i]))
                    except:
                        pass
        line =f.readline()
    return stim,line

def read_evol(f):
    line =f.readline()
    nbevol = int(line.split('=')[-1])
    evol=[]
    for i in range(nbevol):
        evol.append(ParamEvolClass())
    line =f.readline()
    while not("::" in line or line == ''):
        if not (line =='' or line == "\n"):
            lsplit = line.split("\t")
            name = lsplit[0]
            vals = lsplit[1:]
            if any(([x in name for x in ["Name","typeinterp"]])):
                for i in range(nbevol):
                    setattr(evol[i],name,vals[i])
            elif any(([x in name for x in ["Add_Edit"]])):
                pass
            elif any(([x in name for x in ["NMM"]])):
                for i in range(nbevol):
                    val=vals[i].replace("[","")
                    val=val.replace("]","")
                    setattr(evol[i],name,[int(x) for x in val.split(",")])
            else:
                for i in range(nbevol):
                    val=vals[i].replace("[","")
                    val=val.replace("]","")
                    setattr(evol[i],name,[float(x) for x in val.split(",")])
        line =f.readline()
    return evol,line


def read_model(f):
    line =f.readline()
    modelname = line.split('=')[-1]
    modelname=modelname.replace(" ", "")
    modelname=modelname.replace("\n", "")
    line =f.readline()
    nbmodel = int(line.split('=')[-1])
    model={}
    line =f.readline()
    while not("::" in line or line == ''):
        if not (line =='' or line == "\n"):
            lsplit = line.split("\t")
            name = lsplit[0]
            vals = lsplit[1:nbmodel+1]
            try:
                vals=  [float(x) for x in vals]
            except:
                pass
            model[name] = vals
        line =f.readline()
    return model,modelname,line


def read_cm(f,listcm):

    count = -1

    line =f.readline()
    nbcm = int(line.split('=')[-1])

    list_of_CM=[]
    for i in range(len(listcm)-1):
        list_of_CM.append(np.zeros((nbcm,nbcm)))

    Delay= np.zeros((nbcm,nbcm))


    line =f.readline()
    l = 0
    while not(line == ''):
        if not (line =='' or line == "\n"):
            if any(([x in line for x in listcm])):
                count += 1
                line =f.readline()
                l = 0

            lsplit = line.split("\t")[:nbcm]
            lsplit=  [float(x) for x in lsplit]

            if count < (len(listcm)-1):
                list_of_CM[count][l,:]=np.array(lsplit)
                l +=1
            else:
                Delay[l,:]=np.array(lsplit)
                l +=1


        line =f.readline()
    return list_of_CM,Delay,line


def Save_Simulation(fileName=None,stim=None,evol=None ,model=None):
    cm=model.list_of_CM
    delay = model.DelayMat
    listCMnames=model.listCMnames
    f = open(fileName , 'w')
    if  stim:
        write_stim(f,stim)
    if  evol:
        write_evol(f,evol)
    if  model.pop:
        write_model(f,model,model.get_NMM_Variables(),model.modelName)
    if  cm:
        write_cm(f,listCMnames,cm,delay)
    f.close()

def write_stim(f,stim):
    f.write("Stim_info::\n")
    stim_dict = dict(stim[0].__dict__.items())
    name = list(stim_dict.keys())
    nbstim = len(stim)
    f.write("Nb_stim = " + str(nbstim) + "\n")
    for n in name:
        f.write(n+"\t")
        for i in range(nbstim):
            val = getattr(stim[i],n)
            f.write(str(val)+"\t")
        f.write("\n")

def write_evol(f,evol):
    f.write("Evol_info::\n")
    evol_dict = dict(evol[0].__dict__.items())
    name = list(evol_dict.keys())
    nbevol = len(evol)
    f.write("Nb_evol = " + str(nbevol) + "\n")
    for n in name:
        f.write(n+"\t")
        for i in range(nbevol):
            val = getattr(evol[i],n)
            f.write(str(val)+"\t")
        f.write("\n")

def write_model(f,model,listvar,modelName):
    f.write("Model_info::\n")
    f.write("Model_Name = " + modelName +"\n")
    nbmodel = model.pop.NbNMMs
    f.write("Nb_NMM = " + str(nbmodel) + "\n")
    for n in listvar:
        f.write(n+"\t")
        val = getattr(model.pop,n)
        for i in range(nbmodel):
            f.write(str(val[i])+"\t")
        f.write("\n")
    if model.popName:
        f.write('Name' + "\t")
        for n in model.popName:
            f.write(n + "\t")
        f.write("\n")
    if model.popColor:
        f.write('Color' + "\t")
        for n in model.popColor:
            f.write(n + "\t")
        f.write("\n")

def write_cm(f,listcm,cm,delay):
    f.write("CM_info::\n")
    nbcm = len(listcm)
    raw,col = cm[0].shape
    f.write("Nb_NMM = " + str(raw) + "\n")
    for n in range(nbcm):
        f.write(listcm[n]+":\n")
        if n<(nbcm-1):
            for l in range(raw):
                for c in range(col):
                    f.write(str('%E' % cm[n][l,c])+"\t")
                f.write("\n")
        else:
            for l in range(raw):
                for c in range(col):
                    f.write(str('%E' % delay[l,c])+"\t")
                f.write("\n")

def LoadLeadfield(FileName):
    mat = scipy.io.loadmat(FileName)
    return np.array(mat[[l for l in list(mat.keys()) if not '__' in l][0]])
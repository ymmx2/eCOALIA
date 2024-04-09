__author__ = "Yochum Maxime"
__copyright__ = "No Copyright"
__credits__ = ["Yochum Maxime"]
__email__ = "maxime.yochum@univ-rennes1.fr"
__status__ = "Prototype"


class stim_sig():
    def __init__(self):
        self.kind = 'Sinus'
        self.times = 0.
        self.timee = 1.
        self.A = 1.
        self.Aend = ''
        self.f0 = 10.
        self.f1 = ''
        self.phi =0.
        self.pulsewidth =0.001
        self.pop = []
        self.Add_Edit = False #Add = False, Edit = True
        self.dVm = []
        self.Fs = 1024
        self.RiseTC = 0.001
        self.FallTC = 0.0015
        self.RiseTW = 0.01
        self.FallTW = 0.02

class ParamEvolClass():
    def __init__(self):
        self.NMM =  [0]
        self.Name = 'A'
        self.time = []
        self.val = []
        self.typeinterp = 'linear'
        self.Add_Edit = False
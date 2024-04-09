from PyQt6.QtGui import *
from PyQt6.QtCore import *
from PyQt6.QtWidgets import *
import os
import sys

class z_Model_GUI(QMainWindow):
    Mod_OBJ = pyqtSignal(list)
    Close_OBJ = pyqtSignal()
    def __init__(self, parent= None, model=None, cortex = None ):
        super(z_Model_GUI, self).__init__()
        self.model = model
        self.cortex = cortex

        self.isclosed = False

        self.left = 0
        self.top = 0
        self.width = 1500
        self.height = 800
        self.height_per_line = 20
        self.Width_per_Edit = 50
        self.Width_param = 400
        self.resize(self.width,self.height)
        # self.setGeometry(self.left, self.top, self.width, self.height)
        # self.setFixedWidth(self.width)
        # self.setFixedHeight(self.height)

        #update button
        self.B_Apply = QPushButton('Apply')
        self.B_Apply.setFixedSize(self.Width_param,30)
        self.B_Apply.clicked.connect(self.B_apply_fun)
        #population
        GB_Population = QGroupBox('Population')
        GB_Population.setFixedWidth(self.Width_param)
        layout_range = QVBoxLayout()
        self.pop_No = QComboBox()
        self.curI = 0
        if self.model == None:
            for idx in range(len(self.cortex.popName)):
                self.pop_No.addItem(str(idx) + '_' + self.cortex.popName[idx])
                self.pop_No.currentIndexChanged.connect(self.pop_No_fun)
        else:
            self.pop_No.setEnabled(False)
        layout_range.addWidget(self.pop_No)
        GB_Population.setLayout(layout_range)

        if self.model == None:
            self.model = self.cortex.pop

            #Synaptic Gain (mV)
        labelGroup = "Synaptic Gain (mV)"
        label = ['A : excitatory (Glu)','Bb : inhibitory (SSTb, GABAa)','Ba : inhibitory (SSTa, GABAa)','G : inhibitory (PV, GABAa)', 'D : inhibitory (VIP, GABAa)','R : inhibitory (NGFC, GABAb)']
        edit  = ['0','0','0','0','0', '0']
        GB_Syn_Gain = QGroupBox(labelGroup)
        GB_Syn_Gain.setFixedWidth(self.Width_param)
        layout_range = QVBoxLayout()
        # layout.setFixedHeight(  height_per_line* len(label))
        grid = QGridLayout()
        layout_range.addLayout(grid)
        Edit_List = []
        for idx in range(len(label)):
            Label = QLabel(label[idx])
            Label.setFixedHeight(self.height_per_line)
            Edit = LineEdit(edit[idx])
            Edit.setFixedHeight(self.height_per_line)
            Edit.setFixedWidth(self.Width_per_Edit)
            grid.addWidget(Label, idx, 0)
            grid.addWidget(Edit, idx, 1, Qt.AlignmentFlag.AlignLeft)
            Edit_List.append(Edit)
        GB_Syn_Gain.setLayout(layout_range)
        self.Edit_Syn_Gain = Edit_List

        # Synaptic Gain (mV)
        labelGroup = "PSP Constant Rates (s-1)"
        label3 =['AMPA: ', 'SSTb, GABAa,s : ', 'SSTa, GABAa,s : ', 'PV, GABAa,f ', 'VIP, GABAa : ', 'NGFC, GABAb']
        edit = ['0', '0', '0', '0', '0', '0']
        label2 = ['a2', 'bb2', 'ba2',  'g2', 'd2', 'r2']
        edit2 = ['0', '0', '0', '0', '0,' '0', '0']
        GB_PSP_CR = QGroupBox(labelGroup)
        GB_PSP_CR.setFixedWidth(self.Width_param)
        layout_range = QVBoxLayout()
        # layout.setFixedHeight(  height_per_line* len(label))
        grid = QGridLayout()
        layout_range.addLayout(grid)
        Edit_List = []
        for idx in range(len(label)):
            Label3 = QLabel(label3[idx])
            Label3.setFixedHeight(self.height_per_line)
            Label = QLabel(label[idx])
            Label.setFixedHeight(self.height_per_line)
            Edit = LineEdit(edit[idx])
            Edit.setFixedHeight(self.height_per_line)
            Edit.setFixedWidth(self.Width_per_Edit)
            grid.addWidget(Label3, idx, 0)
            grid.addWidget(Label, idx, 1)
            grid.addWidget(Edit, idx, 2, Qt.AlignmentFlag.AlignLeft)
            Edit_List.append(Edit)

            Label = QLabel(label2[idx])
            Label.setFixedHeight(self.height_per_line)
            Edit = LineEdit(edit2[idx])
            Edit.setFixedHeight(self.height_per_line)
            Edit.setFixedWidth(self.Width_per_Edit)
            grid.addWidget(Label, idx, 3)
            grid.addWidget(Edit, idx, 4, Qt.AlignmentFlag.AlignLeft)
            Edit_List.append(Edit)
        GB_PSP_CR.setLayout(layout_range)
        self.Edit_PSP_CR = Edit_List

        # Sigmoid function
        labelGroup = "Sigmoid Functions"
        label = [' ','v0 (inflection point)', 'e0 (2x(firing at inflection))', 'r0 (slope at inflection)']
        edit = [['PC', 'I1','I2b/a', 'I3','I4'],
                ['0', '0', '0', '0', '0'],
                ['0', '0', '0', '0', '0'],
                ['0', '0', '0', '0', '0']]
        GB_Sig = QGroupBox(labelGroup)
        GB_Sig.setFixedWidth(self.Width_param)
        GB_Sig.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout_range = QVBoxLayout()
        # layout.setFixedHeight(  height_per_line* len(label))
        grid = QGridLayout()
        layout_range.addLayout(grid)
        Edit_List = []
        for idx in range(len(label)):
            Label = QLabel(label[idx])
            Label.setFixedHeight(self.height_per_line)
            Label.setFixedWidth(self.Width_per_Edit*3)
            grid.addWidget(Label, idx , 0)
            if idx ==0:
                for i,ed in enumerate(edit[idx]):
                    Edit = QLabel(ed)
                    Edit.setFixedHeight(self.height_per_line)
                    Edit.setFixedWidth(self.Width_per_Edit )
                    grid.addWidget(Edit, idx, i+1)
            else:
                for i, ed in enumerate(edit[idx]):
                    Edit = LineEdit(ed)
                    Edit.setFixedHeight(self.height_per_line)
                    Edit.setFixedWidth(self.Width_per_Edit )
                    grid.addWidget(Edit, idx, i+1)
                    Edit_List.append(Edit)
        GB_Sig.setLayout(layout_range)
        self.Edit_Sig = Edit_List
        GB_Sig.setAlignment(Qt.AlignmentFlag.AlignTop)

        # Noises
        labelGroup = "Noise Inputs p(t)"
        label = [' ', 'mean', 'std', 'coef']
        edit = [['PC', 'I1', 'I2b/a', 'I3', 'I4'],
                ['0', '0', '0', '0', '0'],
                ['0', '0', '0', '0', '0'],
                ['0', '0', '0', '0', '0']]
        GB_Noise = QGroupBox(labelGroup)
        GB_Noise.setFixedWidth(self.Width_param)
        GB_Noise.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout_range = QVBoxLayout()
        # layout.setFixedHeight(  height_per_line* len(label))
        grid = QGridLayout()
        layout_range.addLayout(grid)
        Edit_List = []
        for idx in range(len(label)):
            Label = QLabel(label[idx])
            Label.setFixedHeight(self.height_per_line)
            Label.setFixedWidth(self.Width_per_Edit * 3)
            grid.addWidget(Label, idx, 0)
            if idx == 0:
                for i, ed in enumerate(edit[idx]):
                    Edit = QLabel(ed)
                    Edit.setFixedHeight(self.height_per_line)
                    Edit.setFixedWidth(self.Width_per_Edit)
                    grid.addWidget(Edit, idx, i + 1)
            else:
                for i, ed in enumerate(edit[idx]):
                    Edit = LineEdit(ed)
                    Edit.setFixedHeight(self.height_per_line)
                    Edit.setFixedWidth(self.Width_per_Edit)
                    grid.addWidget(Edit, idx, i + 1)
                    Edit_List.append(Edit)
        GB_Noise.setLayout(layout_range)
        self.Edit_Noise= Edit_List
        GB_Noise.setAlignment(Qt.AlignmentFlag.AlignTop)
        #Graph
        script_dir = os.path.dirname(__file__)
        oImage = QImage(os.path.join(script_dir, "Model_NeoNMM_GUI.png"))
        img_w = oImage.width()
        img_h = oImage.height()
        r_w = img_w / (self.width - self.Width_param)
        r_h = img_h / self.height
        factor = 1
        if (r_w <= 1 and r_h <= 1) or (r_w > 1 and r_h > 1):
            if r_w >= r_h:
                oImage = oImage.scaledToWidth(self.width - self.Width_param,Qt.SmoothTransformation)
                factor = r_w
            else:
                oImage = oImage.scaledToHeight(self.height,Qt.SmoothTransformation)
                factor = r_h
        elif r_w <= 1 and r_h > 1:
            oImage = oImage.scaledToHeight(self.height,Qt.SmoothTransformation)
            factor = r_h
        elif r_w > 1 and r_h <= 1:
            oImage = oImage.scaledToWidth(self.width - self.Width_paramv)
            factor = r_w
        
        self.Graph = self.graph(self,oImage,factor)
        self.Graph.setFixedSize(oImage.width(),oImage.height())

        #stim + up down layout
        self.Edit_stim_k_updown_l = QHBoxLayout()
        # Stim
        labelGroup = "Stim k factor"
        label = ['kP','kPp', 'kI1 ', 'kI2a/b', 'kI3', 'kI4']
        edit = ['0', '0', '0', '0', '0', '0']
        GB_stim_k = QGroupBox(labelGroup)
        GB_stim_k.setFixedWidth(self.Width_param//2)
        layout_range = QVBoxLayout()
        # layout.setFixedHeight(  height_per_line* len(label))
        grid = QGridLayout()
        layout_range.addLayout(grid)
        Edit_List = []
        for idx in range(len(label)):
            Label = QLabel(label[idx])
            Label.setFixedHeight(self.height_per_line)
            Edit = LineEdit(edit[idx])
            Edit.setFixedHeight(self.height_per_line)
            Edit.setFixedWidth(self.Width_per_Edit)
            grid.addWidget(Label, idx, 0)
            grid.addWidget(Edit, idx, 1, Qt.AlignmentFlag.AlignLeft)
            Edit_List.append(Edit)
        GB_stim_k.setLayout(layout_range)
        self.Edit_stim_k = Edit_List

        # # Updown
        # labelGroup = "Up down"
        # label = ['Thresh','Pv0_UD','Pe0_UD','Pr0_UD']
        # edit = ['0', '0', '0', '0']
        # GB_Up_down = QGroupBox(labelGroup)
        # GB_Up_down.setFixedWidth(self.Width_param//2)
        # layout_range = QVBoxLayout()
        # # layout.setFixedHeight(  height_per_line* len(label))
        # grid = QGridLayout()
        # layout_range.addLayout(grid)
        # Edit_List = []
        # for idx in range(len(label)):
        #     Label = QLabel(label[idx])
        #     Label.setFixedHeight(self.height_per_line)
        #     Edit = LineEdit(edit[idx])
        #     Edit.setFixedHeight(self.height_per_line)
        #     Edit.setFixedWidth(self.Width_per_Edit)
        #     grid.addWidget(Label, idx, 0)
        #     grid.addWidget(Edit, idx, 1, Qt.AlignLeft)
        #     Edit_List.append(Edit)
        # GB_Up_down.setLayout(layout_range)
        # self.Edit_Up_down = Edit_List

        self.Edit_stim_k_updown_l.addWidget(GB_stim_k)
        # self.Edit_stim_k_updown_l.addWidget(GB_Up_down)


        self.paramlayout = QVBoxLayout()
        self.paramlayout.addWidget(self.B_Apply)
        self.paramlayout.addWidget(GB_Population)
        self.paramlayout.addWidget(GB_Syn_Gain)
        self.paramlayout.addWidget(GB_PSP_CR)
        self.paramlayout.addWidget(GB_Sig)
        self.paramlayout.addWidget(GB_Noise)
        self.paramlayout.addLayout(self.Edit_stim_k_updown_l)
        self.paramlayout.setAlignment(Qt.AlignmentFlag.AlignTop)


        self.CentralWidget = QWidget()
        self.globallayout = QHBoxLayout()
        self.globallayout.addLayout(self.paramlayout)
        self.globallayout.addWidget(self.Graph)

        self.CentralWidget.setLayout(self.globallayout)
        self.setCentralWidget(self.CentralWidget)

        self.update_values()

    def pop_No_fun(self):
        # self.model = self.cortex.pop[self.pop_No.currentIndex()]
        self.curI = self.pop_No.currentIndex()
        self.update_values()

    def B_apply_fun(self):
        if not self.model==None:
            edit = ['A', 'Bb', 'Ba','G', 'D', 'R']
            for i, ed in enumerate(self.Edit_Syn_Gain):
                getattr(self.model, edit[i])[self.curI] = float(ed.text())

            edit =['a1', 'a2', 'bb1', 'bb2', 'ba1', 'ba2', 'g1', 'g2', 'd1', 'd2', 'r1', 'r2']
            for i, ed in enumerate(self.Edit_PSP_CR):
                getattr(self.model, edit[i])[self.curI] = float(ed.text())

            edit = ['Pv0', 'I1v0', 'I2v0','I3v0','I4v0','Pe0', 'I1e0', 'I2e0','I3e0','I4e0','Pr0', 'I1r0', 'I2r0','I3r0','I4r0']

            for i, ed in enumerate(self.Edit_Sig):
                getattr(self.model, edit[i])[self.curI] = float(ed.text())

            edit = ['Pm', 'I1m', 'I2m','I3m','I4m','Ps', 'I1s', 'I2s','I3s','I4s','Pcoef', 'I1coef','I2coef','I3coef','I4coef']

            for i, ed in enumerate(self.Edit_Noise):
                getattr(self.model, edit[i])[self.curI] = float(ed.text())

            edit = ['k_P', 'k_Pp','k_I1', 'k_I2', 'k_I3', 'k_I4']
            for i, ed in enumerate(self.Edit_stim_k):
                getattr(self.model, edit[i])[self.curI] = float(ed.text())

            for i,ed in enumerate(self.Graph.edit_):
                getattr(self.model, self.Graph.C[i])[self.curI] = float(ed.text())

            # edit = ['Thresh','Pv0_UD','Pe0_UD','Pr0_UD']
            # for i, ed in enumerate(self.Edit_Up_down):
            #     getattr(self.model, edit[i])[self.curI] = float(ed.text())

            for i, ed in enumerate(self.Graph.edit_):
                getattr(self.model, self.Graph.C[i])[self.curI] = float(ed.text())

            if not self.cortex == None:
                self.cortex.pop = self.model
                self.Mod_OBJ.emit([self.cortex])
            else:
                self.Mod_OBJ.emit([self.model])
                self.Close_OBJ.emit()
                self.isclosed = True
                self.close()
            msgBox = QMessageBox()
            msgBox.setText("The model has been updated.")
            msgBox.exec()


    def closeEvent(self, event):
        self.Close_OBJ.emit()
        self.isclosed = True
        self.close()


    def update_values(self):

        edit = [str(getattr(self.model, 'A')[self.curI]), str(getattr(self.model, 'Bb')[self.curI]), str(getattr(self.model, 'Ba')[self.curI]), str(getattr(self.model, 'G')[self.curI]), str(getattr(self.model, 'D')[self.curI]), str(getattr(self.model, 'R')[self.curI])]
        for i, ed in enumerate(self.Edit_Syn_Gain):
            ed.setText(edit[i])

        edit = [str(getattr(self.model, 'a1')[self.curI]),str(getattr(self.model, 'a2')[self.curI]),
                str(getattr(self.model, 'bb1')[self.curI]),str(getattr(self.model, 'bb2')[self.curI]),
                str(getattr(self.model, 'ba1')[self.curI]),str(getattr(self.model, 'ba2')[self.curI]),
                str(getattr(self.model, 'g1')[self.curI]),str(getattr(self.model, 'g2')[self.curI]),
                str(getattr(self.model, 'd1')[self.curI]),str(getattr(self.model, 'd2')[self.curI]),
                str(getattr(self.model, 'r1')[self.curI]),str(getattr(self.model, 'r2')[self.curI])]
        for i, ed in enumerate(self.Edit_PSP_CR):
            ed.setText(edit[i])

        edit = [str(getattr(self.model, 'Pv0')[self.curI]), str(getattr(self.model, 'I1v0')[self.curI]), str(getattr(self.model, 'I2v0')[self.curI]), str(getattr(self.model, 'I3v0')[self.curI]), str(getattr(self.model, 'I4v0')[self.curI]),
                str(getattr(self.model, 'Pe0')[self.curI]), str(getattr(self.model, 'I1e0')[self.curI]), str(getattr(self.model, 'I2e0')[self.curI]), str(getattr(self.model, 'I3e0')[self.curI]), str(getattr(self.model, 'I4e0')[self.curI]),
                str(getattr(self.model, 'Pr0')[self.curI]), str(getattr(self.model, 'I1r0')[self.curI]), str(getattr(self.model, 'I2r0')[self.curI]), str(getattr(self.model, 'I3r0')[self.curI]), str(getattr(self.model, 'I4r0')[self.curI])]

        for i, ed in enumerate(self.Edit_Sig):
            ed.setText(edit[i])

        
        edit = [ str(getattr(self.model, 'Pm')[self.curI]), str(getattr(self.model, 'I1m')[self.curI]), str(getattr(self.model, 'I2m')[self.curI]),
                 str(getattr(self.model, 'I3m')[self.curI]),str(getattr(self.model, 'I4m')[self.curI]),
                 str(getattr(self.model, 'Ps')[self.curI]), str(getattr(self.model, 'I1s')[self.curI]), str(getattr(self.model, 'I2s')[self.curI]),
                 str(getattr(self.model, 'I3s')[self.curI]),str(getattr(self.model, 'I4s')[self.curI]),
                 str(getattr(self.model, 'Pcoef')[self.curI]), str(getattr(self.model, 'I1coef')[self.curI]), str(getattr(self.model, 'I2coef')[self.curI]),
                 str(getattr(self.model, 'I3coef')[self.curI]),str(getattr(self.model, 'I4coef')[self.curI]) ]


        for i, ed in enumerate(self.Edit_Noise):
            ed.setText(edit[i])

        edit = [str(getattr(self.model, 'k_P')[self.curI]), str(getattr(self.model, 'k_Pp')[self.curI]),str(getattr(self.model, 'k_I1')[self.curI]),
                str(getattr(self.model, 'k_I2')[self.curI]),str(getattr(self.model, 'k_I3')[self.curI]),str(getattr(self.model, 'k_I4')[self.curI])]
        for i, ed in enumerate(self.Edit_stim_k):
            ed.setText(edit[i])

        # edit = [str(getattr(self.model, 'Thresh')[self.curI]), str(getattr(self.model, 'Pv0_UD')[self.curI]),
        #         str(getattr(self.model, 'Pe0_UD')[self.curI]),
        #         str(getattr(self.model, 'Pr0_UD')[self.curI])]
        # for i, ed in enumerate(self.Edit_Up_down):
        #     ed.setText(edit[i])

        self.Graph.updatevalues()

    class graph(QWidget):

        def __init__(self, parent = None, oImage=None, factor=None):
            super().__init__()
            self.parent = parent
            # self.setGeometry(300,300,300,220)
            self.setAutoFillBackground(True)

            # sImage = oImage.scaled(QSize(self.width(), self.height()))  # resize Image to widgets size
            # sImage = oImage.scaled(QSize(width, height))  # resize Image to widgets size

            palette = QPalette()
            palette.setBrush(QPalette.ColorRole.Window,QBrush(oImage))  # 10 = Windowrole
            self.setPalette(palette)

            self.C = ["CPP", "CP1P", "CI1P", "CI2bP", "CI2aP", "CI4P", "CPI1", "CI1bI1", "CI2I1", "CI4I1", "CPI1b","CI1I1", "CI1I1b", "CPI2",
                      "CI3I2", "CI4I2", "CPP1", "CPI3", "CI2I3", "CI4I3","CI2I4","CI4I4","CI4bI4","CI4I4b"]
            v = []
            for c in self.C:
                v.append(getattr(self.parent.model,c)[self.parent.curI])
            pos = [[638.0, 571.0], #CPP 
                    [710.0, 395.0], #CP1P
                    [466.0, 561.0], #CI1P 
                    [600.0, 425.0], #CI2bP
                    [595.0, 110.0], #CI2aP
                    [483.0, 35.0], #CI4P
                    [480.0, 505.0], #CPI1
                    [276.0, 608.0], #CI1bI1
                    [370.0, 443.0], # CI2I1
                    [300.0, 350.0], #CI4I1
                    [451.0, 623.0], # CPI1b
                    [233.0, 507.0], #CI1I1
                    [179.0, 417.0], # CI1I1b
                    [485.0, 455.0], # CPI2
                    [516.0, 287.0], #CI3I2
                    [380.0, 300.0], #CI4I2
                    [783.0, 509.0], #CPP1
                    [870.0, 480.0], #CPI3
                    [800.0, 350.0], #CI2I3
                    [843.0, 220.0], #CI4I3
                    [380.0, 170.0], # CI2I4
                    [321.0, 35.0], #CI4I4
                    [225.0, 169.0], #CI4bI4 
                    [207.0, 35.0] ] #CI4I4b

            self.edit_ = []
            self.label_ = []
            for i, c in enumerate(self.C):
                self.label_.append(QLabel(c, self))
                self.label_[i].setAlignment(Qt.AlignmentFlag.AlignCenter)
                # self.label_[i].setFixedSize(50,20)
                self.label_[i].setFont(QFont("Times", 14, weight=QFont.Weight.Bold))
                self.edit_.append(LineEdit(str(v[i]), self))
                self.edit_[i].setFixedSize(50,20)
                l_w = self.label_[i].fontMetrics().boundingRect(self.label_[i].text()).width() / 2
                l_h = self.label_[i].height() / 2
                e_w = self.edit_[i].width() / 2
                e_h = self.edit_[i].height() / 2
                self.label_[i].move(int((pos[i][0] / factor) - l_w), int((pos[i][1] / factor) - l_h - e_h * 1.5))
                self.edit_[i].move(int((pos[i][0] / factor) - e_w), int((pos[i][1] / factor) - e_h))

            # self.CP1P  self.CI1P self.CI2P self.CI3aP
            # self.CPI1  self.CI1aI1  self.CI2I1  self.CPI1a  self.CI1I1a
            # self.CPI2  self.CI3I2
            # self.CPP1
            # self.CI2I3
        def updatevalues(self):
            for i, c in enumerate(self.C):
                self.edit_[i].setText(str(getattr(self.parent.model, c)[self.parent.curI]))



class LineEdit(QLineEdit):
    KEY = Qt.Key.Key_Return
    def __init__(self, *args, **kwargs):
        QLineEdit.__init__(self, *args, **kwargs)
        QREV = QRegularExpressionValidator(QRegularExpression("[+-]?\\d*[\\.]?\\d+"))
        QREV.setLocale(QLocale(QLocale.Language.English))
        self.setValidator(QREV)

def main():
    # from Cortex_Model_Siouar_TargetGains_Pyramidalconnection import Cortex
    from Model_NeoNMM import pop_Coalia
    pop=pop_Coalia()
    # C = Cortex(Nb_NMM=2)
    app = QApplication(sys.argv)
    app.setStyle('Windows')
    ex = z_Model_GUI(app,model = pop)
    # ex = Model_GUI(app,cortex = C)
    ex.setWindowTitle('window')
    ex.show()
    sys.exit(app.exec( ))


if __name__ == '__main__':
    main()
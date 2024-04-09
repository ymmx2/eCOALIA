__author__ = "Yochum Maxime"
__copyright__ = "No Copyright"
__credits__ = ["Yochum Maxime"]
__email__ = "maxime.yochum@univ-rennes1.fr"
__status__ = "Prototype"

from PyQt6.QtGui import *
from PyQt6.QtCore import *
from PyQt6.QtWidgets import *
import numpy as np
import sys
import copy
import scipy.io as sio
import os
import vtk
import scipy.io
from scipy import signal
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from matplotlib import cm


class LineEdit_Int(QLineEdit):
    KEY = Qt.Key.Key_Return
    def __init__(self, *args, **kwargs):
        QLineEdit.__init__(self, *args, **kwargs)
        QREV = QRegularExpressionValidator(QRegularExpression("[+-]?\\d+"))
        QREV.setLocale(QLocale(QLocale.Language.English))
        self.setValidator(QREV)

class LineEdit(QLineEdit):
    KEY = Qt.Key.Key_Return
    def __init__(self, *args, **kwargs):
        QLineEdit.__init__(self, *args, **kwargs)
        QREV = QRegularExpressionValidator(QRegularExpression("[+-]?\\d*[\\.]?\\d+"))
        QREV.setLocale(QLocale(QLocale.Language.English))
        self.setValidator(QREV)

class Mesh_SimpleView(QMainWindow):
    def __init__(self,parent=None, LFPs=None, Fs=None, Names=None,Colors=None,FileName=None):
        super(Mesh_SimpleView, self).__init__()
        if getattr(sys, 'frozen', False):
            self.application_path = os.path.dirname(sys.executable)
        elif __file__:
            self.application_path = os.path.dirname(__file__)
        self.setWindowTitle("Source View")
        self.parent = parent
        self.LFPs = LFPs
        self.Fs = Fs
        self.Names = Names
        self.Colors = Colors
        self.power = np.zeros((len(self.Names),2))

        self.centralWidget = QWidget()
        self.setCentralWidget(self.centralWidget)
        self.mainHBOX_param_scene = QVBoxLayout()


        self.get_mesh_path(filename_vtk=r'Ressources/Cortex_Collin_15002.vtk')
        self.get_parcellisation_path(filename_mat=r'Ressources/scout_Desikan_kiliany_66_RL.mat')
        self.get_electrode_path(filename_mat=FileName)


        layout_parambandeau = QHBoxLayout()






        band_group = QButtonGroup()
        self.display_names_CB = QCheckBox('display régions names')
        self.rWithout = QRadioButton("Without [0Hz-")
        self.rDelta = QRadioButton("Delta [0-4Hz[")
        self.rTheta = QRadioButton("Theta [4-8Hz[")
        self.rAlpha = QRadioButton("Alpha [8-12Hz[")
        self.rBeta = QRadioButton("Beta [12-30Hz[")
        self.rgamma = QRadioButton("gamma [30Hz-")
        self.rgamma.setChecked(True)
        self.rAmplitude_Mean = QRadioButton("Amplitude_Mean")
        self.rAmplitude_Max = QRadioButton("Amplitude_Max")
        self.rAmplitude_RMS = QRadioButton("Amplitude_RMS")
        band_group.addButton(self.rWithout)
        band_group.addButton(self.rDelta)
        band_group.addButton(self.rTheta)
        band_group.addButton(self.rAlpha)
        band_group.addButton(self.rBeta)
        band_group.addButton(self.rgamma)
        band_group.addButton(self.rAmplitude_Mean)
        band_group.addButton(self.rAmplitude_Max)
        band_group.addButton(self.rAmplitude_RMS)
        layout_type = QHBoxLayout()
        grid = QGridLayout()
        layout_type.addLayout(grid)
        l = 0
        c = 0
        grid.addWidget(self.display_names_CB,l,c)
        l += 1
        grid.addWidget(QLabel('Freq based'),l,c)
        l += 1
        grid.addWidget(self.rWithout,l,c)
        l += 1
        grid.addWidget(self.rDelta,l,c)
        l += 1
        grid.addWidget(self.rTheta,l,c)
        l += 1
        grid.addWidget(self.rAlpha,l,c)
        l += 1
        grid.addWidget(self.rBeta,l,c)
        l += 1
        grid.addWidget(self.rgamma,l,c)
        l = 0
        c = 1
        grid.addWidget(QLabel('Amp based'),l,c)
        l += 1
        grid.addWidget(self.rAmplitude_Mean,l,c)
        l += 1
        grid.addWidget(self.rAmplitude_Max,l,c)
        l += 1
        grid.addWidget(self.rAmplitude_RMS,l,c)

        layout_win =QHBoxLayout()
        window_l = QLabel('Window (s)')
        self.window_LE = LineEdit('1')
        layout_win.addWidget(window_l)
        layout_win.addWidget(self.window_LE)
        self.UpdateColor_PB = QPushButton('Update Color')
        self.UpdateColor_PB.clicked.connect(self.UpdateColor_PB_click)
        layout = QVBoxLayout()
        layout.addLayout(layout_type)
        layout.addLayout(layout_win)
        layout.addWidget(self.UpdateColor_PB)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)


        layout_parambandeau.addLayout(layout)


        self.MV_plot = VTK_MeshView(self)
        layout_opacity = QVBoxLayout()
        Opacity_L = QLabel()
        Opacity_L.setFixedSize(QSize(20,20))
        pixic = QPixmap(os.path.join(self.application_path,'icons', r'Ressources/icons/opacity.png'))
        Opacity_L.setPixmap(QPixmap.fromImage(QImage(pixic).scaled(20, 20, Qt.AspectRatioMode.KeepAspectRatio)))
        Opacity_L.repaint()
        self.Opacity_S = QSlider()
        self.Opacity_S.setRange(0,100)
        self.Opacity_S.setValue(100)
        self.Opacity_S.setSingleStep(1)
        layout_opacity.addWidget(Opacity_L)
        layout_opacity.addWidget(self.Opacity_S)
        layout_brightness = QVBoxLayout()
        brightnessL = QLabel()
        brightnessL.setFixedSize(QSize(20,20))
        pixic2 =QPixmap(os.path.join(self.application_path, 'icons', 'Ressources/icons/brightness.png'))
        brightnessL.setPixmap(QPixmap.fromImage(QImage(pixic2).scaled(20, 20, Qt.AspectRatioMode.KeepAspectRatio)))
        brightnessL.repaint()
        self.Brightness_S = QSlider()
        self.Brightness_S.setRange(0, 100)
        self.Brightness_S.setValue(0)
        self.Brightness_S.setSingleStep(1)
        layout_brightness.addWidget(brightnessL)
        layout_brightness.addWidget(self.Brightness_S)
        layout_gammacor = QVBoxLayout()
        gammacorL = QLabel()
        gammacorL.setFixedSize(QSize(20, 20))
        pixic3 = QPixmap(os.path.join(self.application_path, 'icons', 'Ressources/icons/gamma.png'))
        img3 = QImage(pixic3).scaled(20, 20, Qt.AspectRatioMode.KeepAspectRatio)
        gammacorL.setPixmap(QPixmap.fromImage(img3))
        gammacorL.repaint()
        self.Gammacor_S = QSlider()
        self.Gammacor_S.setRange(0, 100)
        self.Gammacor_S.setValue(10)
        self.Gammacor_S.setSingleStep(1)
        layout_gammacor.addWidget(gammacorL)
        layout_gammacor.addWidget(self.Gammacor_S)
        self.gamma=1.

        layout_threshold = QVBoxLayout()
        thresholdL = QLabel()
        thresholdL.setFixedSize(QSize(20, 20))
        pixic3 = QPixmap(os.path.join(self.application_path, 'icons', 'Ressources/icons/threshold.png'))
        img3 = QImage(pixic3).scaled(20, 20, Qt.AspectRatioMode.KeepAspectRatio)
        thresholdL.setPixmap(QPixmap.fromImage(img3))
        thresholdL.repaint()
        self.Threshold_S = RangeSlider()
        self.Threshold_S.setRangeLimit(0, 1000)
        self.Threshold_S.setRange(0,1000)
        layout_threshold.addWidget(thresholdL)
        layout_threshold.addWidget(self.Threshold_S)

        mesh_l = QHBoxLayout()
        mesh_l.setContentsMargins(0,0,0,0)
        mesh_l.addWidget(self.MV_plot)
        mesh_l.addLayout(layout_opacity)
        mesh_l.addLayout(layout_brightness)
        mesh_l.addLayout(layout_gammacor)
        # mesh_l.addLayout(layout_threshold_min)
        # mesh_l.addLayout(layout_threshold_max)
        mesh_l.addLayout(layout_threshold)
        self.threshold_max = 1.
        self.threshold_min = 0.


        lecture_l = QHBoxLayout()
        self.colormap = QComboBox()

        self.colormap.addItems(['viridis', 'plasma', 'inferno', 'magma', 'cividis',
                                'Greys_r', 'Purples_r', 'Blues_r', 'Greens_r', 'Oranges_r', 'Reds_r',
                                'YlOrBr_r', 'YlOrRd_r', 'OrRd_r', 'PuRd_r', 'RdPu_r', 'BuPu_r',
                                'GnBu_r', 'PuBu_r', 'YlGnBu_r', 'PuBuGn_r', 'BuGn_r', 'YlGn_r',
                                'binary_r', 'gist_yarg_r', 'gist_gray', 'gray', 'bone', 'pink',
                                'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
                                'hot', 'afmhot', 'gist_heat', 'copper',
                                'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
                                'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'
                                ])
        self.colormap.setFixedWidth(100)
        self.colormap.currentIndexChanged.connect(self.colorize_power)
        self.myslider = StyledTextScrollBar()
        self.Play_PB = QPushButton('Play')
        self.Play_PB.setFixedWidth(50)
        self.Rec_RB = QCheckBox('Rec')
        self.Rec_RB.setFixedWidth(50)
        self.Rec_RB.setVisible(False)
        self.timer = None
        lecture_l.addWidget(self.colormap)
        lecture_l.addWidget(self.myslider)
        lecture_l.addWidget(self.Play_PB)
        lecture_l.addWidget(self.Rec_RB)

        self.myslider.valueChanged.connect( self.timevalueChanged )
        self.Opacity_S.valueChanged.connect( self.OpacityvalueChanged )
        self.Brightness_S.valueChanged.connect( self.BrightnessvalueChanged )
        self.Gammacor_S.valueChanged.connect(self.Gammacor_S_click)
        self.Threshold_S.valueChanged.connect(self.Threshold_S_click)
        self.Play_PB.clicked.connect(self.Play_PB_click)

        wiget1 = QWidget()
        wiget1.setLayout(layout_parambandeau)
        wiget2 = QWidget()
        wiget2.setLayout(mesh_l)
        self.mainsplitter = QSplitter(Qt.Orientation.Vertical)
        self.mainsplitter.addWidget(wiget1)
        self.mainsplitter.addWidget(wiget2)
        self.mainsplitter.setStretchFactor(0, 1)
        self.mainsplitter.setStretchFactor(1, 2)
        self.mainsplitter.setSizes([1500, 1500 * 3 ])

        self.mainHBOX_param_scene.addWidget(self.mainsplitter)
        # self.mainHBOX_param_scene.addLayout(layout_parambandeau)
        # self.mainHBOX_param_scene.addLayout(mesh_l)
        self.mainHBOX_param_scene.addLayout(lecture_l)
        self.centralWidget.setLayout(self.mainHBOX_param_scene)
        self.setMinimumSize(400, 400)



    def OpacityvalueChanged(self):
        if self.MV_plot.IGotMesh:
            self.MV_plot.set_opacity(percent =self.Opacity_S.value()/100)

    def BrightnessvalueChanged(self):
        if self.MV_plot.IGotMesh:
            self.MV_plot.set_brightness(percent=self.Brightness_S.value() / 100)

    def Gammacor_S_click(self):
        if self.MV_plot.IGotMesh:
            self.gamma = self.Gammacor_S.value()/10
            self.colorize_power()
            self.timevalueChanged()

    def Threshold_S_click(self):
        if self.MV_plot.IGotMesh:
            self.threshold_min, self.threshold_max = self.Threshold_S.getRange()
            self.threshold_min /= 1000
            self.threshold_max /= 1000
            if self.threshold_max < self.threshold_min:
                self.threshold_max, self.threshold_min = self.threshold_min, self.threshold_max
            self.colorize_power()
            self.timevalueChanged()

    def Play_PB_click(self):
        if self.MV_plot.IGotMesh:
            if self.timer is None:
                window = float(self.window_LE.text())
                if window < 0.01:
                    window = 0.01
                delay = int(window * 1000)
                if self.Rec_RB.isChecked():
                    fileName = QFileDialog.getSaveFileName(caption='Movie file', filter="AVI (*.avi)")
                    if (fileName[0] == []):
                        return
                    self.imageFilter = vtk.vtkWindowToImageFilter()
                    self.imageFilter.SetInput(self.MV_plot.vtkWidget.GetRenderWindow())
                    self.imageFilter.SetInputBufferTypeToRGB()
                    self.imageFilter.ReadFrontBufferOff()
                    self.imageFilter.Update()
                    self.writer = vtk.vtkAVIWriter()
                    self.writer.SetInputConnection(self.imageFilter.GetOutputPort())
                    self.writer.SetQuality(1)
                    self.writer.SetRate(int(1 / window))
                    self.writer.SetCompressorFourCC("DIB")
                    self.writer.SetFileName(fileName[0])
                    self.writer.Start()
                else:
                    self.imageFilter = None
                    self.writer = None

                self.timer = QTimer()
                self.timer.setInterval(delay)
                self.timer.timeout.connect(self.fromPlay)
                self.timer.start()
                self.Play_PB.setText('Pause')
            else:
                self.timer.stop()
                self.timer = None
                if not self.writer is None:
                    self.writer.End()
                    self.imageFilter = None
                self.writer = None
                self.Play_PB.setText('Play')

    def fromPlay(self):
        self.myslider.setValue(self.myslider.value() + 1)
        if not self.imageFilter is None:
            self.imageFilter.Modified()
            self.writer.Write()
        if self.myslider.value() >= self.myslider.maximum():
            self.timer.stop()
            self.timer = None
            if not self.writer is None:
                self.writer.End()
            self.imageFilter = None
            self.writer = None
            self.Play_PB.setText('Play')


    def timevalueChanged(self, play = False):
        if self.Colors is not None:
            timeint = self.myslider.value()
            window = float(self.window_LE.text())
            if window < 0.01:
                window = 0.01
            self.myslider.setSliderText(text=str("%.2f" % (window*timeint)) + ' s')

            couleurs = self.Colors[:,timeint,:]
            self.MV_plot.update_colors(couleurs,self.Vertex_correspondance)
            self.MV_plot.Rendering()

    def get_mesh_path(self,filename_vtk=None ):
        reader = vtk.vtkGenericDataObjectReader()
        reader.SetFileName(filename_vtk)
        reader.Update()
        # get the brain mesh
        self.GetOutput = reader.GetOutput()
        self.points = np.array(self.GetOutput.GetPoints().GetData())


    def get_parcellisation_path(self, filename_mat=None):
        mat = scipy.io.loadmat(filename_mat)
        Scouts = mat['Scouts'][0]
        dtypes = Scouts[0].dtype
        names = list(dtypes.names)
        self.patch = []
        self.Regions = []
        self.RegionsCoord = []

        for i in range(Scouts.shape[0]):
            self.RegionsCoord.append(self.points[Scouts[i]['Seed'][0][0]-1])
            self.Regions.append(Scouts[i]['Label'][0])
            region = {}
            for n in names:
                if not n == 'Handles':
                    region[n] = Scouts[i][n][0]
            self.patch.append(region)
        self.RegionsCoord = np.array(self.RegionsCoord)
        self.Vertex_correspondance = []
        for i in range(self.points.shape[0]):
            c = 0
            for i_p, p in enumerate(self.patch):
                if i + 1 in p['Vertices']:
                    c = i_p
                    break
            self.Vertex_correspondance.append(c)

    def get_electrode_path(self,  filename_mat=None ):
        with open(filename_mat, "r") as f:
            rows = f.readlines()
            self.Regions_elec = []
            self.E_coord = []
            for row in rows:
                stripped_line = row.strip()
                line_list = stripped_line.split()
                self.Regions_elec.append(line_list[0])
                self.E_coord.append([float(item) for item in line_list[1:]])
        self.E_coord = np.array(self.E_coord)




    def UpdateColor_PB_click(self):

        lfp = np.array(self.LFPs)

        window = float(self.window_LE.text())
        if window < 1/self.Fs:
            window = 1/self.Fs
        taille = lfp.shape[1]
        mini = 0
        segment = int(window*self.Fs)
        maxi = int(taille / segment)
        setp = 1
        self.slider_update(mini,maxi-1,setp)
        nbregion = len(self.Regions)

        if self.rDelta.isChecked() or self.rTheta.isChecked() or self.rAlpha.isChecked() or self.rBeta.isChecked() or self.rgamma.isChecked() or self.rWithout.isChecked():
            f,t,y = signal.spectrogram(lfp, fs=self.Fs, nperseg=segment,noverlap=0)
            if self.rDelta.isChecked():
                fmin = 0
                fmax = 4
            elif self.rTheta.isChecked():
                fmin = 4
                fmax = 8
            elif self.rAlpha.isChecked():
                fmin = 8
                fmax = 12
            elif self.rBeta.isChecked():
                fmin = 12
                fmax = 30
            elif self.rgamma.isChecked():
                fmin = 30
                fmax = f[-1]
            elif self.rWithout.isChecked():
                fmin = 0
                fmax = f[-1]

            self.power = np.sum(y[:, np.bitwise_and(f >= fmin , f < fmax), :], axis=1)
            self.power = self.power / np.max(self.power)
        elif self.rAmplitude_Mean.isChecked() or self.rAmplitude_Max.isChecked() or self.rAmplitude_RMS.isChecked():
            self.power = np.zeros((nbregion,maxi))
            for i in range(maxi):
                if self.rAmplitude_Mean.isChecked():
                    self.power[:,i] = np.mean(lfp[:, i * segment:(i + 1) * segment],axis=1)
                elif self.rAmplitude_Max.isChecked():
                    self.power[:, i] = np.max(lfp[:, i * segment:(i + 1) * segment], axis=1)
                elif self.rAmplitude_RMS.isChecked():
                    self.power[:, i] = np.sqrt(np.mean(lfp[:, i * segment:(i + 1) * segment] ** 2, axis=1))
            self.power = self.power / np.max(self.power)

        self.colorize_power()
        self.MV_plot.draw_Mesh(self.GetOutput, self.points, self.RegionsCoord,  Regions_elec=self.Regions_elec,E_coord= self.E_coord)
        self.timevalueChanged()
        self.OpacityvalueChanged()
        self.BrightnessvalueChanged()

    def colorize_power(self):
        cmap = cm.get_cmap(self.colormap.currentText())
        power = copy.deepcopy(self.power)
        power[power>self.threshold_max] = 1#self.threshold_max
        power[power<self.threshold_min] = 0#self.threshold_min
        self.Colors = np.uint8(cmap(np.power(power,self.gamma)) * 255)[:, :, :3]


    def slider_update(self,min,max,interval):
        self.myslider.setMinimum(min)
        self.myslider.setMaximum(max)
        self.myslider.setPageStep(interval)

class VTK_MeshView(QMainWindow):
    def __init__(self, parent=None):
        super(VTK_MeshView, self).__init__(parent)
        self.parent=parent
        self.frame = QFrame()
        self.vl = QVBoxLayout()
        self.vtkWidget = QVTKRenderWindowInteractor(self.frame)
        self.vl.addWidget(self.vtkWidget)

        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(.1, .1, .1)
        # self.ren.SetUseDepthPeeling(1)
        # self.ren.SetMaximumNumberOfPeels(100)
        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
        self.style1 = vtk.vtkInteractorStyleTrackballCamera()
        self.style1.SetDefaultRenderer(self.ren)
        self.iren.SetInteractorStyle(self.style1)

        self.ren.ResetCamera()

        self.frame.setLayout(self.vl)
        self.setCentralWidget(self.frame)

        self.show()
        self.iren.Initialize()
        self.iren.Start()

        self.IGotMesh = False
        self.gamma = 1

    def set_opacity(self,percent=None):
        # try:
        self.actor.GetProperty().SetOpacity(percent)
        self.Rendering()
        # except:
        #     pass

    def set_brightness(self,percent=None):
        self.actor.GetProperty().SetAmbient(percent)
        self.Rendering()

    def set_gamma(self,gamma=1):
        self.gamma = gamma


    def update_colors(self,couleurs=None,Vertex_correspondance=None):
        Colors = vtk.vtkUnsignedCharArray()
        Colors.SetNumberOfComponents(3)
        Colors.SetName("Colors")
        for i in range(self.points.shape[0]):
            if couleurs is None:
                Colors.InsertNextTuple3(255,255,255)
            else:
                c = couleurs[Vertex_correspondance[i],:]
                Colors.InsertNextTuple3(c[0], c[1], c[2])
        self.GetOutput.GetPointData().SetScalars(Colors)
        self.GetOutput.Modified()

    def draw_Mesh(self,GetOutput,points,patch, Regions_elec ,E_coord  ):
        self.GetOutput = GetOutput
        self.points = points
        self.patch = patch
        self.Regions_elec = Regions_elec
        self.E_coord = E_coord
        self.ren.RemoveAllViewProps()
        self.update_colors((np.random.random((len(self.patch),3))*255).astype(np.uint8),self.parent.Vertex_correspondance)
        vtkpoints = vtk.vtkPoints()
        for p in points:
            vtkpoints.InsertNextPoint(p)

        mapper = vtk.vtkPolyDataMapper()
        self.actor = vtk.vtkActor()
        mapper.SetInputData(self.GetOutput)
        self.actor.SetMapper(mapper)
        self.actor.GetProperty().SetPointSize(20)

        self.ren.AddActor(self.actor)


        #electrode shpere
        self.shperes = []
        self.shperes_m = []
        self.shperes_a = []

        for i in range(len(self.Regions_elec)):
            source = vtk.vtkSphereSource()
            source.SetCenter(self.E_coord[i,0],self.E_coord[i,1],self.E_coord[i,2])
            source.SetRadius(4.0)

            self.shperes.append(source)
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(source.GetOutputPort())
            self.shperes_m.append(mapper)
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            self.shperes_a.append(actor)
            self.ren.AddActor(actor)

        if self.parent.display_names_CB.isChecked():
            for i in range(len(self.Regions_elec)):
                pos = self.E_coord[i]
                name = self.Regions_elec[i]
                textActor = vtk.vtkBillboardTextActor3D()
                textActor.SetInput(name)
                pos2 = pos + pos/np.linalg.norm(pos) * 10
                textActor.SetPosition(pos2[0], pos2[1], pos2[2])
                textActor.GetTextProperty().SetFontSize(20)
                textActor.GetTextProperty().SetJustificationToCentered()
                # textActor.SetDisplayOffset(-20, -20)
                textActor.PickableOff()
                self.ren.AddActor(textActor)


        self.set_center()
        self.Rendering()
        self.IGotMesh = True

    def Rendering(self):
        self.iren.ReInitialize()
        self.iren.GetRenderWindow().Render()


    def set_center(self):
        self.ren.ResetCamera()

    def namearegion(self,index):
        if self.IGotMesh:
            coord = vtk.vtkCoordinate()
            coord.SetCoordinateSystemToWorld()
            x, y, z = self.parent.RegionsCoord[index]
            coord.SetValue(x, y, z)
            x, y = coord.GetComputedDisplayValue(self.ren)
            print(x, y )
            self.style1.addregiontex(x, y, index)
            self.Rendering()



class RangeSlider(QSlider):

    valueChanged = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(Qt.Orientation.Vertical)
        self.first_position = 0
        self.second_position = 1000

        self.opt = QStyleOptionSlider()
        self.opt.minimum = 0
        self.opt.maximum = 1000
        self.opt.orientation = Qt.Orientation.Vertical

        self.setTickPosition(QSlider.TickPosition.NoTicks)
        self.setTickInterval(1)

        self.setSizePolicy(
            QSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Expanding, QSizePolicy.ControlType.Slider)
        )

    def setRangeLimit(self, minimum: int, maximum: int):
        self.opt.minimum = minimum
        self.opt.maximum = maximum

    def setRange(self, start: int, end: int):
        self.first_position = start
        self.second_position = end

    def getRange(self):
        return (self.first_position, self.second_position)

    def setTickPosition(self, position: QSlider.TickPosition):
        self.opt.tickPosition = position

    def setTickInterval(self, ti: int):
        self.opt.tickInterval = ti

    def paintEvent(self, event: QPaintEvent):

        painter = QPainter(self)

        # Draw rule
        self.opt.initFrom(self)
        self.opt.rect = self.rect()
        self.opt.sliderPosition = 0
        self.opt.subControls = QStyle.SubControl.SC_SliderGroove | QStyle.SubControl.SC_SliderTickmarks

        #   Draw GROOVE
        self.style().drawComplexControl(QStyle.ComplexControl.CC_Slider, self.opt, painter)

        #  Draw INTERVAL

        color = self.palette().color(self.palette().ColorGroup.All,QPalette.ColorRole.Highlight)
        color.setAlpha(160)
        painter.setBrush(QBrush(color))
        painter.setPen(Qt.PenStyle.NoPen)

        self.opt.sliderPosition = self.first_position
        x_left_handle = (
            self.style()
                .subControlRect(QStyle.ComplexControl.CC_Slider, self.opt, QStyle.SubControl.SC_SliderHandle)
                .bottom()
        )

        self.opt.sliderPosition = self.second_position
        x_right_handle = (
            self.style()
                .subControlRect(QStyle.ComplexControl.CC_Slider, self.opt, QStyle.SubControl.SC_SliderHandle)
                .top()
        )

        groove_rect = self.style().subControlRect(
            QStyle.ComplexControl.CC_Slider, self.opt, QStyle.SubControl.SC_SliderGroove
        )

        selection = QRect(
            groove_rect.x()-4,
            x_left_handle,
            groove_rect.height()+4,
            x_right_handle - x_left_handle
        ).adjusted(-1, 1, 1, -1)

        painter.drawRect(selection)

        # Draw first handle

        self.opt.subControls = QStyle.SubControl.SC_SliderHandle
        self.opt.sliderPosition = self.first_position
        self.style().drawComplexControl(QStyle.ComplexControl.CC_Slider, self.opt, painter)

        # Draw second handle
        self.opt.sliderPosition = self.second_position
        self.style().drawComplexControl(QStyle.ComplexControl.CC_Slider, self.opt, painter)

    def mousePressEvent(self, event: QMouseEvent):

        self.opt.sliderPosition = self.first_position
        self._first_sc = self.style().hitTestComplexControl(
            QStyle.ComplexControl.CC_Slider, self.opt, event.pos(), self
        )

        self.opt.sliderPosition = self.second_position
        self._second_sc = self.style().hitTestComplexControl(
            QStyle.ComplexControl.CC_Slider, self.opt, event.pos(), self
        )

    def mouseMoveEvent(self, event: QMouseEvent):

        distance = self.opt.maximum - self.opt.minimum

        pos = self.style().sliderValueFromPosition(
            0, distance, event.pos().y(), self.rect().height()
        )

        if self._first_sc == QStyle.SubControl.SC_SliderHandle:
            if pos <= self.second_position:
                self.first_position = pos
                self.update()
                self.valueChanged.emit()
                return

        if self._second_sc == QStyle.SubControl.SC_SliderHandle:
            if pos >= self.first_position:
                self.second_position = pos
                self.update()
                self.valueChanged.emit()


class TextScrollBarStyle(QProxyStyle):
    def drawComplexControl(self, control, option, painter, widget):
        # call the base implementation which will draw anything Qt will ask
        super().drawComplexControl(control, option, painter, widget)
        # check if the control type is that of a scroll bar and its orientation matches
        if control == QStyle.ComplexControl.CC_ScrollBar and option.orientation == Qt.Orientation.Horizontal:
            # the option is already provided by the widget's internal paintEvent;
            # from this point on, it's almost the same as explained above, but
            # setting the pen might be required for some styles
            painter.setPen(widget.palette().color(widget.palette().ColorGroup.Active,QPalette.ColorRole.WindowText))
            margin = self.frameMargin(widget)
            subPageRect = self.subControlRect(control, option, QStyle.SubControl.SC_ScrollBarSubPage, widget)
            painter.drawText(subPageRect.adjusted(margin, 0, 0, 0), Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignVCenter, widget.preText)

            addPageRect = self.subControlRect(control, option, QStyle.SubControl.SC_ScrollBarAddPage, widget)
            painter.drawText(addPageRect.adjusted(0, 0, -margin, 0), Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignVCenter, widget.postText)

            sliderRect = self.subControlRect(control, option, QStyle.SubControl.SC_ScrollBarSlider, widget)
            painter.drawText(sliderRect, Qt.AlignmentFlag.AlignCenter, widget.sliderText)

    def frameMargin(self, widget):
        # a helper function to get the default frame margin which is usually added
        # to widgets and sub widgets that might look like a frame, which usually
        # includes the slider of a scrollbar
        option = QStyleOptionFrame()
        option.initFrom(widget)
        return self.pixelMetric(QStyle.PixelMetric.PM_DefaultFrameWidth, option, widget)

    def subControlRect(self, control, option, subControl, widget):
        rect = super().subControlRect(control, option, subControl, widget)
        if (control == QStyle.ComplexControl.CC_ScrollBar and isinstance(widget, StyledTextScrollBar) and
            option.orientation == Qt.Orientation.Horizontal and subControl == QStyle.SubControl.SC_ScrollBarSlider):
                # get the *default* groove rectangle (the space in which the slider can move)
                grooveRect = super().subControlRect(control, option, QStyle.SubControl.SC_ScrollBarGroove, widget)
                # if the slider has text, ensure that the width is large enough to show it
                width = max(rect.width(), widget.sliderWidth + self.frameMargin(widget))
                # compute the position of the slider according to the current value
                pos = self.sliderPositionFromValue(widget.minimum(), widget.maximum(),
                    widget.sliderPosition(), grooveRect.width() - width)
                # return the new rectangle
                return QRect(grooveRect.x() + pos, rect.y(), width, rect.height())
        return rect

    def hitTestComplexControl(self, control, option, pos, widget):
        if control == QStyle.ComplexControl.CC_ScrollBar:
            sliderRect = self.subControlRect(control, option, QStyle.SubControl.SC_ScrollBarSlider, widget)
            if pos in sliderRect:
                return QStyle.SubControl.SC_ScrollBarSlider
        return super().hitTestComplexControl(control, option, pos, widget)

class StyledTextScrollBar(QScrollBar):
    def __init__(self, sliderText='', preText='', postText=''):
        super().__init__(Qt.Orientation.Horizontal)
        self.setStyle(TextScrollBarStyle())
        self.preText = preText
        self.postText = postText
        self.sliderText = sliderText
        self.sliderTextMargin = 2
        self.sliderWidth = self.fontMetrics().horizontalAdvance(sliderText) + self.sliderTextMargin + 2

    def setPreText(self, text):
        self.preText = text
        self.update()

    def setPostText(self, text):
        self.postText = text
        self.update()

    def setSliderText(self, text):
        self.sliderText = text
        self.sliderWidth = self.fontMetrics().horizontalAdvance(text) + self.sliderTextMargin + 2

    def setSliderTextMargin(self, margin):
        self.sliderTextMargin = margin
        self.sliderWidth = self.fontMetrics().horizontalAdvance(self.sliderText) + margin + 2
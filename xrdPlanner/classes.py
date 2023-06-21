import os
import sys
import json
import numpy as np
import pyqtgraph as pg
import Dans_Diffraction as dif
from PyQt6 import QtWidgets, QtCore, QtGui
from pyFAI import calibrant
import xrdPlanner.resources

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # set path home
        self.path = os.path.dirname(__file__)
        # add an icon
        self.setWindowIcon(QtGui.QIcon(':/icons/xrdPlanner.png'))

        # Drag-and-Drop cif-file
        #  dropEvent()
        #  - check dropped file is a cif
        #  get_cif_reference()
        #  - use gemmi to get cell, centring and crystal  system from cif
        #  - use pyFAI get_d_spacings() to create contours
        self.setAcceptDrops(True)
        
        # set path to settings file
        self.file_dump = os.path.join(self.path, 'settings.json')
        # save parameters to file
        self.init_par()

        # experimental darkmode?
        if self.plo.darkmode:
            self.init_darkmode()
        
        # set window color
        pg.setConfigOptions(background=self.plo.plot_bg_color, antialias=True)

        ###########
        # CONTOUR #
        ###########
        # generate contour levels
        self.cont_geom_num = np.linspace(self.plo.conic_tth_min, self.plo.conic_tth_max, self.plo.conic_tth_num)

        # get colormap
        self.cont_cmap = pg.colormap.get(self.plo.conic_colormap)
        
        # reverse the colormap useful to increase visibility in darkmode
        if self.plo.reverse_cmap:
            self.cont_cmap.reverse()
        
        # What standards should be available as reference
        # The d spacings will be imported from pyFAI
        self.geo.ref_library = calibrant.names()
        # dict to store custom reference data
        self.geo.ref_custom = {}
        self.geo.ref_custom_hkl = {}
        
        # get the detector specs
        # - update: overwrite existing file after load
        # - reset: overwrite existing file with defaults
        self.detectors = self.get_det_library(update=self.plo.update_det_bank, reset=self.plo.reset_det_bank)

        # pick current detector
        self.det = self.get_specs_det(self.detectors, self.geo.det_type, self.geo.det_size)

        # translate unit for plot title
        self.geo.unit_names = ['2\U0001D6F3 [\u00B0]', 'd [\u212B\u207B\u00B9]', 'q [\u212B]', 'sin(\U0001D6F3)/\U0001D706 [\u212B]']
        if self.geo.unit >= len(self.geo.unit_names):
            print(f'Error: Valid geo.unit range is from 0 to {len(self.geo.unit_names)-1}, geo.unit={self.geo.unit}')
            raise SystemExit

        # init the hkl tooltip
        font = QtGui.QFont()
        font.setPixelSize(self.plo.conic_hkl_label_size)
        font.setBold(True)
        QtWidgets.QToolTip.setFont(font)

        ###########
        # GENERAL #
        ###########

        # menubar is displayed within the main window on Windows
        # so we need to make space for it
        # no idea about other OS, if there are issues fix them here
        if sys.platform == 'win32':
            self.offset_win32 = self.menuBar().height() - int(round(self.plo.slider_margin/2, 0))
        else:
            self.offset_win32 = 0
        
        # remove window margins
        self.setContentsMargins(0, 0, 0, 0)
        
        # define grid layout
        self.layout = QtWidgets.QGridLayout()
        self.layout.setContentsMargins(0, self.plo.slider_margin, 0, 0)
        
        # make a widget, set the layout
        centralwidget = QtWidgets.QWidget()
        centralwidget.setLayout(self.layout)
        centralwidget.setContentsMargins(0, 0, 0, 0)
        self.setCentralWidget(centralwidget)

        # add the plot to the layout
        self.ax = pg.plot()
        # added to avoid the error:
        # qt.pointer.dispatch: skipping QEventPoint(id=0 ts=0 pos=0,0 scn=482.023,246.011
        # gbl=482.023,246.011 Released ellipse=(1x1 ∡ 0) vel=0,0 press=-482.023,-246.011
        # last=-482.023,-246.011 Δ 482.023,246.011) : no target window
        self.ax.viewport().setAttribute(QtCore.Qt.WidgetAttribute.WA_AcceptTouchEvents, False)
        self.layout.addWidget(self.ax)

        # populate the menus with detectors, references and units
        self.init_menus()
        
        # initialize the screen
        self.init_screen()

        # add the slider frame
        self.sliderWidget = SliderWidget(self, self.geo, self.plo, self.lmt)

    def init_darkmode(self):
        # set window color
        self.plo.conic_label_fill = '#000000'    # [str]    Contour label fill color
        # - reference contour section - 
        self.plo.conic_ref_color = '#202020'     # [color]  Reference contour color
        self.plo.det_module_color = '#EEEEEE'    # [color]  Detector module border color
        self.plo.det_module_fill = '#EEEEEE'     # [color]  Detector module background color
        # - general section - 
        self.plo.plot_bg_color = '#000000'       # [str]    Plot background color
        self.plo.unit_label_color = '#808080'    # [str]    Label color
        self.plo.unit_label_fill = '#000000'     # [str]    Label fill color
        # - slider section - 
        self.plo.slider_border_color = '#202020' # [str]    Slider window border color
        self.plo.slider_bg_color = '#AA303030'   # [str]    Slider window background color
        self.plo.slider_bg_hover = '#303030'     # [str]    Slider window hover color
        self.plo.slider_label_color = '#808080'  # [str]    Slider window label color
        
    def init_screen(self):
        # init the plot for contours and beam center
        self.ax.setAspectLocked()
        #remove axes
        self.ax.getPlotItem().hideAxis('bottom')
        self.ax.getPlotItem().hideAxis('left')
        # disable pan/zoom
        self.ax.setMouseEnabled(x=False, y=False)
        # disable right click  context menu
        self.ax.setMenuEnabled(False)
        # hide the autoscale button
        self.ax.hideButtons()

        # enable debug mode functions
        if self.plo.set_debug:
            # disable pan/zoom
            self.ax.setMouseEnabled(x=True, y=True)
            # disable right click  context menu
            self.ax.setMenuEnabled(True)
        
        # container for contour lines
        self.patches = {'beamcenter':None, 'conic':[], 'reference':[], 'labels':[]}

        # add beam center scatter plot
        self.patches['beamcenter'] = pg.ScatterPlotItem()
        self.ax.addItem(self.patches['beamcenter'])

        # add empty plot per reference contour line
        for i in range(self.plo.conic_ref_num):
            ref = pg.PlotCurveItem(useCache=True)
            self.ax.addItem(ref)
            self.patches['reference'].append(ref)
            self.patches['reference'][i].setClickable(True, width=self.plo.conic_ref_linewidth)
            self.patches['reference'][i].sigClicked.connect(self.show_tooltip)
            self.patches['reference'][i].name = None
        
        # define font to label conics
        font = QtGui.QFont()
        font.setPixelSize(self.plo.conic_label_size)
        font.setBold(True)

        # add empty plot per contour line
        for i in range(self.plo.conic_tth_num):
            curve = pg.PlotCurveItem(useCache=True)
            self.ax.addItem(curve)
            self.patches['conic'].append(curve)
            temp_label = pg.TextItem(anchor=(0.5,0.5), fill=pg.mkBrush(self.plo.conic_label_fill))
            temp_label.setFont(font)
            self.patches['labels'].append(temp_label)
            self.ax.addItem(temp_label)
        
        # build detector modules
        self.build_detector()

        # figure out proper plot dimensions
        # limit the axis x and y
        # get proper dimensions
        # fix the window size
        # resize the window
        # resize the plot
        self.resize_window()

        # add unit label
        self.add_unit_label()

        # create cones and draw contour lines
        self.update_screen()
        self.set_window_title()

    def init_menus(self):
        menuBar = QtWidgets.QMenuBar()
        self.setMenuBar(menuBar)

        menu_det = menuBar.addMenu('Detector')
        group_det = QtGui.QActionGroup(self)
        group_det.setExclusive(True)

        # menu Detectors
        for d in self.detectors:
            d = d.upper()
            d_menu = QtWidgets.QMenu(d, self)
            d_menu.setStatusTip('')
            menu_det.addMenu(d_menu)
            for s in self.detectors[d]['size']:
                s = s.upper()
                det_action = QtGui.QAction(s, self, checkable=True)
                self.set_menu_action(det_action, self.change_detector, d, s)
                d_menu.addAction(det_action)
                group_det.addAction(det_action)
                if d == self.geo.det_type.upper() and s == self.geo.det_size.upper():
                    det_action.setChecked(True)
        
        self.menu_ref = menuBar.addMenu('Reference')
        self.group_ref = QtGui.QActionGroup(self)
        self.group_ref.setExclusive(True)
        
        # menu Reference: add None
        ref_action = QtGui.QAction('None', self, checkable=True)
        self.set_menu_action(ref_action, self.change_reference, 'None')
        self.menu_ref.addAction(ref_action)
        self.group_ref.addAction(ref_action)
        if 'None' == self.geo.reference:
            ref_action.setChecked(True)
        
        # menu Reference: add pyFAI library
        self.sub_menu_pyFAI = QtWidgets.QMenu('pyFAI', self)
        self.sub_menu_pyFAI.setStatusTip('')
        self.menu_ref.addMenu(self.sub_menu_pyFAI)
        for ref_name in self.geo.ref_library:
            ref_action = QtGui.QAction(ref_name, self, checkable=True)
            self.set_menu_action(ref_action, self.change_reference, ref_name)
            self.sub_menu_pyFAI.addAction(ref_action)
            self.group_ref.addAction(ref_action)
            if ref_name == self.geo.reference:
                ref_action.setChecked(True)

        # menu Reference: add Custom
        self.sub_menu_custom = QtWidgets.QMenu('Custom', self)
        self.sub_menu_custom.setStatusTip('Drag and Drop *.cif files.')
        self.menu_ref.addMenu(self.sub_menu_custom)
        
        # menu Units
        menu_unit = menuBar.addMenu('Units')
        group_unit = QtGui.QActionGroup(self)
        group_unit.setExclusive(True)
        for unit_index, unit_name in enumerate(self.geo.unit_names):
            unit_action = QtGui.QAction(unit_name, self, checkable=True)
            self.set_menu_action(unit_action, self.change_units, unit_index)
            menu_unit.addAction(unit_action)
            group_unit.addAction(unit_action)
            if unit_index == self.geo.unit:
                unit_action.setChecked(True)

    def add_unit_label(self):
        font = QtGui.QFont()
        font.setPixelSize(self.plo.unit_label_size)
        self.unit_label = pg.TextItem(anchor=(0.0,0.0), color=self.plo.unit_label_color, fill=self.plo.unit_label_fill)
        self.unit_label.setText(self.geo.unit_names[self.geo.unit])
        self.unit_label.setFont(font)
        self.ax.addItem(self.unit_label)
        self.unit_label.setPos(-self.plo.xdim, self.plo.ydim)

    def resize_window(self):
        # figure out proper plot dimensions
        self.plo.xdim = (self.det.hms * self.det.hmn + self.det.pxs * self.det.hgp * self.det.hmn + self.det.cbh)/2
        self.plo.ydim = (self.det.vms * self.det.vmn + self.det.pxs * self.det.vgp * self.det.vmn + self.det.cbh)/2
        
        # limit the axis x and y
        self.ax.setXRange(-self.plo.xdim, self.plo.xdim, padding=0, update=True)
        self.ax.setYRange(-self.plo.ydim, self.plo.ydim, padding=0, update=True)

        # get proper dimensions
        width = int(np.ceil(self.plo.plot_size * self.plo.xdim / self.plo.ydim))
        height = self.plo.plot_size + self.plo.slider_margin//2 + self.offset_win32

        # fix the window size
        if self.plo.plot_size_fixed:
            self.setMaximumHeight(height)
            self.setMinimumHeight(height)
            self.setMaximumWidth(width)
            self.setMinimumWidth(width)

        # resize the window
        self.resize(width, height)

    def set_menu_action(self, action, target, *args):
        action.triggered.connect(lambda: target(*args))

    def set_window_title(self):
        if self.geo.reference == 'None':
            self.setWindowTitle(self.det.name)
        else:
            self.setWindowTitle(f'{self.det.name} - {self.geo.reference}')

    def change_detector(self, det_name, det_size):
        self.det = self.get_specs_det(self.detectors, det_name, det_size)
        self.ax.clear()
        self.init_screen()
        self.sliderWidget.center_frame()

    def change_units(self, unit_index):
        self.geo.unit = unit_index
        self.unit_label.setText(self.geo.unit_names[unit_index])
        self.draw_contours()

    def change_reference(self, ref_name):
        self.geo.reference = ref_name
        self.get_reference()
        self.draw_reference()

    def get_cif_reference(self, fpath):
        xtl = dif.Crystal(fpath)
        # :return xval: arrray : x-axis of powder scan (units)
        # :return inten: array : intensity values at each point in x-axis
        # :return reflections: (h, k, l, xval, intensity) array of reflection positions, grouped by min_overlap
        xval, inten, reflections = xtl.Scatter.powder(scattering_type='xray', units='dspace', powder_average=True, min_overlap=0.02)
        # reject low intensities: based on median or mean?
        # median is always around unity -> useless
        # mean rejects many, add adjustable multiplicator?
        used = reflections[reflections[:,4] > reflections[:,4].max() * self.plo.conic_ref_min_int]
        # sort by intensity -> ascending -> flip
        ordered = used[used[:, 4].argsort()][::-1]
        # pick the strongest
        ordered = ordered[:self.plo.conic_ref_num]
        # assign dspacing
        self.plo.cont_ref_dsp = ordered[:,3]
        # hkl -> integer
        # cast hkl array to list of tuples (for easy display)
        hkl = ordered[:,:3].astype(int)
        if self.plo.conic_hkl_show_int:
            self.plo.cont_ref_hkl = list(zip(hkl[:,0], hkl[:,1], hkl[:,2], ordered[:,4].astype(int)))
        else:
            self.plo.cont_ref_hkl = list(zip(hkl[:,0], hkl[:,1], hkl[:,2]))

        self.geo.reference = os.path.basename(fpath)
        self.geo.ref_custom[self.geo.reference] = self.plo.cont_ref_dsp
        self.geo.ref_custom_hkl[self.geo.reference] = self.plo.cont_ref_hkl

        ref_action = QtGui.QAction(self.geo.reference, self, checkable=True)
        self.set_menu_action(ref_action, self.change_reference, self.geo.reference)
        self.sub_menu_custom.addAction(ref_action)
        self.group_ref.addAction(ref_action)
        ref_action.setChecked(True)
        
        # update window title
        self.set_window_title()

        self.draw_reference()

    def get_reference(self):
        if self.geo.reference in self.geo.ref_library:
            # get the d spacings for the calibrtant from pyFAI
            self.plo.cont_ref_dsp = np.array(calibrant.get_calibrant(self.geo.reference).get_dSpacing()[:self.plo.conic_ref_num])
            self.plo.cont_ref_hkl = None
        elif self.geo.reference in self.geo.ref_custom:
            # get custom d spacings
            self.plo.cont_ref_dsp = self.geo.ref_custom[self.geo.reference]
            self.plo.cont_ref_hkl = self.geo.ref_custom_hkl[self.geo.reference]
        else:
            # set all d-spacings to -1
            self.plo.cont_ref_dsp = np.zeros(self.plo.conic_ref_num)
            self.plo.cont_ref_hkl = None
        # update window title
        self.set_window_title()

    def get_specs_geo(self):
        ######################
        # Setup the geometry #
        ######################
        geo = container()
        geo.det_type = 'EIGER2' # [str]  Pilatus3 / Eiger2
        geo.det_size = '4M'     # [str]  300K 1M 2M 6M / 1M 4M 9M 16M
        geo.ener = 21.0         # [keV]  Beam energy
        geo.dist = 75.0         # [mm]   Detector distance
        geo.yoff = 0.0          # [mm]   Detector offset (vertical)
        geo.xoff = 0.0          # [mm]   Detector offset (horizontal)
        geo.rota = 25.0         # [deg]  Detector rotation
        geo.tilt = 0.0          # [deg]  Detector tilt
        geo.unit = 1            # [0-3]  Contour legend
                                #          0: 2-Theta
                                #          1: d-spacing
                                #          2: q-space
                                #          3: sin(theta)/lambda
        geo.reference = 'None'  # [str]  Plot reference contours
                                #          pick from pyFAI
        return geo

    def get_specs_plo(self):
        ################
        # Plot Details #
        ################
        plo = container()
        # - geometry contour section - 
        plo.conic_tth_min = 5               # [int]    Minimum 2-theta contour line
        plo.conic_tth_max = 150             # [int]    Maximum 2-theta contour line
        plo.conic_tth_num = 30              # [int]    Number of contour lines
        plo.beamcenter_marker = 'o'         # [marker] Beam center marker
        plo.beamcenter_size = 6             # [int]    Beam center size
        plo.conic_linewidth = 3.0           # [float]  Contour linewidth
        plo.conic_label_size = 14           # [int]    Contour label size
        plo.conic_label_fill = '#FFFFFF'    # [str]    Contour label fill color
        plo.conic_colormap = 'viridis'      # [cmap]   Contour colormap
        # - reference contour section - 
        plo.conic_ref_color = '#DCDCDC'     # [color]  Reference contour color
        plo.conic_ref_linewidth = 12.0      # [float]  Reference contour linewidth
        plo.conic_ref_num = 100             # [int]    Number of reference contours
        plo.conic_ref_min_int = 0.01        # [int]    Minimum display intensity (cif)
        plo.conic_hkl_label_size = 14       # [int]    Font size of hkl tooltip
        plo.conic_hkl_show_int = False      # [bool]   Include intensity in hkl tooltip
        # - module section - 
        plo.det_module_alpha = 0.20         # [float]  Detector module alpha
        plo.det_module_width = 1            # [float]  Detector module border width
        plo.det_module_color = '#404040'    # [color]  Detector module border color
        plo.det_module_fill = '#404040'     # [color]  Detector module background color
        # - general section - 
        plo.conic_steps = 100               # [int]    Conic resolution
        plo.plot_size = 768                 # [int]    Plot size, px
        plo.plot_size_fixed = True          # [int]    Fix window size
        plo.plot_bg_color = '#FFFFFF'       # [str]    Plot background color
        plo.unit_label_size = 16            # [int]    Label size, px
        plo.unit_label_color = '#808080'    # [str]    Label color
        plo.unit_label_fill = '#FFFFFF'     # [str]    Label fill color
        # - slider section - 
        plo.slider_margin = 12              # [int]    Slider frame top margin
        plo.slider_border_width = 1         # [int]    Slider frame border width
        plo.slider_border_radius = 1        # [int]    Slider frame border radius (px)
        plo.slider_border_color = '#808080' # [str]    Slider frame border color
        plo.slider_bg_color = '#AAC0C0C0'   # [str]    Slider frame background color
        plo.slider_bg_hover = '#C0C0C0'     # [str]    Slider frame hover color
        plo.slider_label_size = 14          # [int]    Slider frame label size
        plo.slider_column_width = 75        # [int]    Slider label column width
        plo.slider_label_color = '#000000'  # [str]    Slider frame label color
        plo.enable_slider_ener = True       # [bool]   Show energy slider
        plo.enable_slider_dist = True       # [bool]   Show distance slider
        plo.enable_slider_rota = True       # [bool]   Show rotation slider
        plo.enable_slider_yoff = True       # [bool]   Show vertical offset slider
        plo.enable_slider_xoff = True       # [bool]   Show horizontal offset slider
        plo.enable_slider_tilt = True       # [bool]   Show tilt slider
        # - update/reset - 
        plo.update_settings = True          # [bool]   Update settings file after load
        plo.update_det_bank = True          # [bool]   Update detector bank after load
        plo.reset_settings = False          # [bool]   Reset settings file
        plo.reset_det_bank = False          # [bool]   Reset detector bank
        # - debug/testing -
        plo.set_debug = False               # [bool]   Debug mode
        plo.darkmode = False                # [bool]   Darkmode
        plo.reverse_cmap = False            # [bool]   Reverse colormap

        return plo

    def get_specs_lmt(self):
        ##########
        # Limits #
        ##########
        lmt = container()
        lmt.ener_min =  5.0    # [float] Energy minimum [keV]
        lmt.ener_max =  100.0  # [float] Energy maximum [keV]
        lmt.ener_stp =  1.0    # [float] Energy step size [keV]
        lmt.dist_min =  40.0   # [float] Distance minimum [mm]
        lmt.dist_max =  1000.0 # [float] Distance maximum [mm]
        lmt.dist_stp =  1.0    # [float] Distance step size [mm]
        lmt.xoff_min = -150.0  # [float] Horizontal offset minimum [mm]
        lmt.xoff_max =  150.0  # [float] Horizontal offset maximum [mm]
        lmt.xoff_stp =  1.0    # [float] Horizontal offset step size [mm]
        lmt.yoff_min = -250.0  # [float] Vertical offset minimum [mm]
        lmt.yoff_max =  250.0  # [float] Vertical offset maximum [mm]
        lmt.yoff_stp =  1.0    # [float] Vertical offset step size [mm]
        lmt.rota_min = -60.0   # [float] Rotation minimum [deg]
        lmt.rota_max =  60.0   # [float] Rotation maximum [deg]
        lmt.rota_stp =  1.0    # [float] Rotation step size [deg]
        lmt.tilt_min = -25.0   # [float] Tilt minimum [deg]
        lmt.tilt_max =  25.0   # [float] Tilt maximum [deg]
        lmt.tilt_stp =  1.0    # [float] Tilt step size [deg]
        
        return lmt

    def get_defaults(self):
        # load the defaults
        # geo: geometry and detector specs
        self.geo = self.get_specs_geo()
        # plo: plot details
        self.plo = self.get_specs_plo()
        # lmt: geometry limits
        self.lmt = self.get_specs_lmt()

    def get_specs_det(self, detectors, det_type, det_size):
        det_type = det_type.upper()
        det_size = det_size.upper()

        if det_type not in detectors.keys():
            print(f'Unknown detector type: {det_type}')
            print(f'Current databank entries: {", ".join(detectors.keys())}.')
            raise SystemExit
        
        if det_size not in detectors[det_type]['size'].keys():
            print(f'Unknown detector type/size combination: {det_type}/{det_size}')
            print(f'Current {det_type} databank sizes: {", ".join(detectors[det_type]["size"].keys())}.')
            raise SystemExit
        
        det = container()
        det.hms = detectors[det_type]['hms']
        det.vms = detectors[det_type]['vms']
        det.pxs = detectors[det_type]['pxs']
        det.hgp = detectors[det_type]['hgp']
        det.vgp = detectors[det_type]['vgp']
        det.cbh = detectors[det_type]['cbh']
        det.hmn, det.vmn = detectors[det_type]['size'][det_size]
        det.name = f'{det_type} {det_size}'

        return det

    def get_det_library(self, update=True, reset=False):
        ###########################
        # Detector Specifications #
        ###########################
        detectors = dict()
        ###############################
        # Specifications for Pilatus3 #
        ###############################
        detectors['PILATUS3'] = {
            'hms' : 83.8,    # [mm]  Module size (horizontal)
            'vms' : 33.5,    # [mm]  Module size (vertical)
            'pxs' : 172e-3,  # [mm]  Pixel size
            'hgp' : 7,       # [pix] Gap between modules (horizontal)
            'vgp' : 17,      # [pix] Gap between modules (vertical)
            'cbh' : 0,       # [mm]  Central beam hole
            'size' : {'300K':(1,3),'1M':(2,5),'2M':(3,8),'6M':(5,12)},
            }
        ###############################
        # Specifications for Pilatus4 #
        ###############################
        # Note: These are probably not
        # the correct PILATUS4 specs
        # and are only meant to play
        # around!
        detectors['PILATUS4'] = {
            'hms' : 75.0,    # [mm]  Module size (horizontal)
            'vms' : 39.0,    # [mm]  Module size (vertical)
            'pxs' : 150e-3,  # [mm]  Pixel size
            'hgp' : 19,      # [pix] Gap between modules (horizontal)
            'vgp' : 6,       # [pix] Gap between modules (vertical)
            'cbh' : 0,       # [mm]  Central beam hole
            'size' : {'1M':(2,4),'2M':(3,6),'4M':(4,8)}
            }
        
        #############################
        # Specifications for Eiger2 #
        #############################
        detectors['EIGER2'] = {
            'hms' : 77.1,    # [mm]  Module size (horizontal)
            'vms' : 38.4,    # [mm]  Module size (vertical)
            'pxs' : 75e-3,   # [mm]  Pixel size
            'hgp' : 38,      # [pix] Gap between modules (horizontal)
            'vgp' : 12,      # [pix] Gap between modules (vertical)
            'cbh' : 0,       # [mm]  Central beam hole
            'size' : {'1M':(1,2),'4M':(2,4),'9M':(3,6),'16M':(4,8)},
            }
        
        #############################
        # Specifications for MPCCD #
        #############################
        detectors['MPCCD'] = {
            'hms' : 51.2,    # [mm]  Module size (horizontal)
            'vms' : 25.6,    # [mm]  Module size (vertical)
            'pxs' : 50e-3,   # [mm]  Pixel size
            'hgp' : 18,      # [pix] Gap between modules (horizontal)
            'vgp' : 27,      # [pix] Gap between modules (vertical)
            'cbh' : 3,       # [mm]  Central beam hole
            'size' : {'4M':(2,4)},
            }

        ##############################
        # Specifications for RAYONIX #
        ##############################
        detectors['RAYONIX'] = {
            'hms' : 75.0,   # [mm]  Module size (horizontal)
            'vms' : 75.0,   # [mm]  Module size (vertical)
            'pxs' : 39e-3,  # [mm]  Pixel size
            'hgp' : 0,      # [pix] Gap between modules (horizontal)
            'vgp' : 0,      # [pix] Gap between modules (vertical)
            'cbh' : 0,      # [mm]  Central beam hole
            'size' : {'MX225-HS':(3,3),'MX300-HS':(4,4)},
            }
        
        #############################
        # Specifications for PHOTON #
        #############################
        detectors['PHOTON-III'] = {
            'hms' : 100.0,  # [mm]  Module size (horizontal)
            'vms' : 70.0,   # [mm]  Module size (vertical)
            'pxs' : 50e-3,  # [mm]  Pixel size
            'hgp' : 0,      # [pix] Gap between modules (horizontal)
            'vgp' : 0,      # [pix] Gap between modules (vertical)
            'cbh' : 0,      # [mm]  Central beam hole
            'size' : {'7':(1,1),'14':(1,2),'28':(2,2)},
            }
        
        #############################
        # Specifications for PHOTON #
        #############################
        detectors['PHOTON-II'] = {
            'hms' : 100.0,  # [mm]  Module size (horizontal)
            'vms' : 70.0,   # [mm]  Module size (vertical)
            'pxs' : 50e-3,  # [mm]  Pixel size
            'hgp' : 0,      # [pix] Gap between modules (horizontal)
            'vgp' : 0,      # [pix] Gap between modules (vertical)
            'cbh' : 0,      # [mm]  Central beam hole
            'size' : {'7':(1,1),'14':(1,2)},
            }
        
        # make file dump
        det_db = os.path.join(self.path, 'detector_db.json')
        if not os.path.exists(det_db) or reset:
            with open(det_db, 'w') as wf:
                json.dump(detectors, wf, indent=4)
        else:
            with open(det_db, 'r') as of:
                for key, vals in json.load(of).items():
                    detectors[key] = vals
        
        if update:
            with open(det_db, 'w') as wf:
                json.dump(detectors, wf, indent=4)
        
        return detectors

    def build_detector(self):
        # build detector modules
        # beam position is between the modules (even) or at the center module (odd)
        # determined by the "+det.hmn%2" part
        for i in range(-self.det.hmn//2+self.det.hmn%2, self.det.hmn-self.det.hmn//2):
            for j in range(-self.det.vmn//2+self.det.vmn%2, self.det.vmn-self.det.vmn//2):
                # - place modules along x (i) and y (j) keeping the gaps in mind ( + (det.hgp*det.pxs)/2)
                # - the " - ((det.hms+det.hgp*det.pxs)/2)" positions the origin (the beam) at the center of a module
                #   and "det.hmn%2" makes sure this is only active for detectors with an odd number of modules
                # - define sets of panels that collectively move to realize a central hole offset for MPCCD detectors
                #   that are used at SACLA/SPring-8:
                #   x = (...) + (det.cbh/2)*(2*(j&det.vmn)//det.vmn-1)
                #   y = (...) + (det.cbh/2)*(1-2*(i&det.hmn)//det.hmn)
                # - negative values of det.cbh for 'clockwise' offset order
                origin_x = i * (self.det.hms + self.det.hgp * self.det.pxs) \
                             - ((self.det.hms + self.det.hgp * self.det.pxs)/2) * (self.det.hmn % 2) \
                             + (self.det.hgp * self.det.pxs)/2 \
                             + (self.det.cbh/2) * (2*(j & self.det.vmn) // self.det.vmn-1)
                origin_y = j * (self.det.vms + self.det.vgp * self.det.pxs) \
                             - ((self.det.vms + self.det.vgp * self.det.pxs)/2) * (self.det.vmn%2) \
                             + (self.det.vgp * self.det.pxs)/2 \
                             + (self.det.cbh/2) * (1-2*(i & self.det.hmn) // self.det.hmn)
                # add the module
                rect_item = QtWidgets.QGraphicsRectItem(origin_x, origin_y,  self.det.hms, self.det.vms)
                rect_item.setPen(pg.mkPen(color = self.plo.det_module_color, width = self.plo.det_module_width))
                rect_item.setBrush(pg.mkBrush(color = self.plo.det_module_fill))
                rect_item.setOpacity(self.plo.det_module_alpha)
                self.ax.addItem(rect_item)

    def draw_contours(self):
        # calculate the offset of the contours resulting from yoff and rotation
        # shift the grid to draw the cones, to make sure the contours are drawn
        # within the visible area
        _comp_shift = -(self.geo.yoff + np.tan(np.deg2rad(self.geo.rota))*self.geo.dist)
        # update beam center
        self.patches['beamcenter'].setData([self.geo.xoff],[_comp_shift],
                                            symbol = self.plo.beamcenter_marker,
                                            size = self.plo.beamcenter_size,
                                            brush = pg.mkBrush(self.cont_cmap.map(0, mode='qcolor')),
                                            pen = pg.mkPen(None))
        for _n, _ttd in enumerate(self.cont_geom_num):
            self.patches['conic'][_n].setVisible(False)
            self.patches['labels'][_n].setVisible(False)
            # current fraction for colormap
            _f = _n/len(self.cont_geom_num)

            # convert theta in degrees to radians
            theta = np.deg2rad(_ttd)

            # convert theta in degrees to radians
            # for some reason I defined it negative some time ago
            # now there's no turning back!
            omega = -np.deg2rad(self.geo.tilt + self.geo.rota)

            # calculate the conic section corresponding to the theta angle
            # :returns False is conic is outside of visiblee area
            x, y, label_pos = self.calc_conic(omega, theta, steps=self.plo.conic_steps)
            if x is False or len(x) == 0:
                continue

            # plot the conic section
            self.patches['conic'][_n].setData(x, y, pen=pg.mkPen(self.cont_cmap.map(_f, mode='qcolor'), width=self.plo.conic_linewidth))
            self.patches['conic'][_n].setVisible(True)

            # Conversion factor keV to Angstrom: 12.398
            # sin(t)/l: np.sin(Theta) / lambda -> (12.398/geo_energy)
            stl = np.sin(theta/2)/(12.398/self.geo.ener)

            # d-spacing: l = 2 d sin(t) -> 1/2(sin(t)/l)
            dsp = 1/(2*stl)

            # prepare the values in the different units / labels
            #_units = {0:np.rad2deg(theta), 1:_dsp, 2:_stl*4*np.pi, 3:_stl}
            _units = {0:_ttd, 1:dsp, 2:stl*4*np.pi, 3:stl}
            self.patches['labels'][_n].setPos(self.geo.xoff, label_pos)
            self.patches['labels'][_n].setText(f'{_units[self.geo.unit]:.2f}', color=self.cont_cmap.map(_f, mode='qcolor'))
            self.patches['labels'][_n].setVisible(True)

    def draw_reference(self):
        # plot reference contour lines
        # standard contour lines are to be drawn
        for _n in range(self.plo.conic_ref_num):
            self.patches['reference'][_n].setVisible(False)
            # number of d-spacings might be lower than the maximum number of allowed contours
            if _n < len(self.plo.cont_ref_dsp):
                _d = self.plo.cont_ref_dsp[_n]
                # None adds a list of zeros
                # catch those here
                if _d <= 0:
                    continue
                # lambda = 2 * d * sin(theta)
                # 2-theta = 2 * (lambda / 2*d)
                # lambda -> (12.398/geo_energy)
                lambda_d = (12.398/self.geo.ener) / (2*_d)
                if lambda_d > 1.0:
                    continue
                
                # get theta
                theta = 2 * np.arcsin(lambda_d)
                
                # convert theta in degrees to radians
                # for some reason I defined it negative some time ago
                # now there's no turning back!
                omega = -np.deg2rad(self.geo.tilt + self.geo.rota)
                
                # calculate the conic section corresponding to the theta angle
                # :returns False is conic is outside of visiblee area
                x, y, _ = self.calc_conic(omega, theta, steps=self.plo.conic_steps)
                if x is False:
                    continue

                # plot the conic section
                self.patches['reference'][_n].setData(x, y, pen=pg.mkPen(self.plo.conic_ref_color, width=self.plo.conic_ref_linewidth))
                self.patches['reference'][_n].setVisible(True)
                
                # if hkl are available
                # put them in the proper container for the contour
                # so indexing gets it right
                if self.plo.cont_ref_hkl:
                    self.patches['reference'][_n].name = self.plo.cont_ref_hkl[_n]
                else:
                    self.patches['reference'][_n].name = None
    
    def calc_conic(self, omega, theta, steps=100):
        # skip drawing smaller/larger +-90 deg contours
        # reject overlap of the 'backscattering'
        # -> limitation of the current implementation
        if theta > np.pi/2 + abs(omega):
            return False, False, False

        # y axis offset of the cone center
        dy_cone = self.geo.dist * np.tan(omega)
        # change in 'r', the length of the cones primary axis
        dz_cone = np.sqrt(self.geo.dist**2 + dy_cone**2)
        # tilt is handled as a rotation but
        # has its travel distance (y) reset.
        comp_tilt = np.deg2rad(self.geo.tilt) * self.geo.dist
        # eccentricity of the resulting conic section
        ecc = np.round(np.cos(np.pi/2 - omega) / np.cos(theta), 10)
        # we need e**2-1 to calculate the width of the hyperbola
        e21 = ecc**2-1
        # y ('height') components/distances from central axis of the cone
        # intersecting the detector plane and the distance to
        # the y intersection of the conic section.
        y1 = dz_cone * np.sin(theta) / np.cos(omega + theta)
        y2 = dz_cone * np.sin(theta) / np.cos(omega - theta)

        # add x/y offsets
        # revert tilt rotation
        y0 = dy_cone - self.geo.yoff + comp_tilt 
        x0 = self.geo.xoff

        # add margin to slightly extend
        # conics outside of visible area
        _xdim = self.plo.xdim * 1.05
        _ydim = self.plo.ydim * 1.05
        # evaluate the eccentricity and parameterise
        # the resulting conic accordingly
        if abs(ecc) == 0:
            # circle
            h = (y1+y2)/2
            # check if the circle is visible
            if h - np.sqrt(y0**2 + x0**2) > np.sqrt(self.plo.ydim**2 + self.plo.xdim**2):
                return False, False, False
            t = np.linspace(0, 2*np.pi, 2*steps)
            x = x0 + h * np.sin(t)
            y = y0 + (y1-y2)/2 + h * np.cos(t)
        elif 0 < abs(ecc) < 1:
            # ellipse
            yd = (y1-y2)/2
            h = (y1+y2)/2
            w = dz_cone * np.sin(theta) * (y1+y2) / (2 * np.sqrt(y1*y2) * np.cos(theta))
            # ellipses that expand ouside the visible area:
            # add a margin to make sure the connecting line
            # of segmented ellipses is outside the visible area
            #
            # I hope this is faster than generating and handing
            # over a 'connect' array to the setData function
            # of the plotItem
            _xlim1 = (_xdim + x0) / w
            _xlim2 = (_xdim - x0) / w
            if _xlim1 < 1 and _xlim2 < 1:
                l = -np.arcsin(_xlim1)
                r =  np.arcsin(_xlim2)
                t = np.hstack([np.linspace(l, r, steps), np.linspace(-r+np.pi, -l+np.pi, steps)])
            else:
                t = np.linspace(0, 2*np.pi, 2*steps)
            x = x0 + w * np.sin(t)
            y = y0 + yd + h * np.cos(t)
        elif abs(ecc) == 1:
            # parabola
            yd = np.sign(ecc) * self.geo.dist * np.tan(abs(omega) - theta)
            a = np.sign(ecc) * self.geo.dist * np.tan(theta)
            l = -(_xdim + x0) / a
            r =  (_xdim - x0) / a
            t = np.linspace(l, r, steps)
            x = x0 + a*t
            y = y0 - dy_cone + yd + a/2 * t**2
        elif 1 < abs(ecc) < 100:
            # hyperbola
            h = np.sign(omega) * (y1+y2)/2
            w = h * np.sqrt(e21)
            l = -np.arcsinh((_xdim + x0) / w)
            r =  np.arcsinh((_xdim - x0) / w)
            t = np.linspace(l, r, steps)
            x = x0 + w * np.sinh(t)
            y = y0 + (y1-y2)/2 - h * np.cosh(t)
        elif abs(ecc) >= 100:
            # line
            t = np.linspace(-_xdim, _xdim, 2)
            a = -np.sign(ecc) * min(abs(y1), abs(y2))
            x = t
            y = y0 + np.ones(len(t)) * a

        # check if conic is visible
        cx = np.argwhere((x >= -_xdim) & (x <= _xdim))
        cy = np.argwhere((y >= -_ydim) & (y <= _ydim))
        if len(cy) == 0 or len(cx) == 0:
            # outside detector area
            return False, False, False
        
        # adjust the label position to maintain readibility
        # this works for most cases but is not the most optimal solution yet
        # OR: use the actual beam position to determine label position
        # beam_pos_y = -(self.geo.yoff + np.tan(np.deg2rad(self.geo.rota))*self.geo.dist)
        if omega <= 0:
            label_pos = max(y) if theta < np.pi/2 else min(y)
        else:
            label_pos = min(y) if theta < np.pi/2 else max(y)
        
        # return x, y and the label position
        return x, y, label_pos

    def show_tooltip(self, widget, event):
        if not widget.name or not self.plo.cont_ref_hkl:
            return False
        pos = QtCore.QPoint(*map(int, event.screenPos()))# - QtCore.QPoint(10,20)
        QtWidgets.QToolTip.showText(pos, str(widget.name))
        event.ignore()

    def update_screen(self, val=None):
        if val is not None:
            if self.sender().objectName() == 'dist':
                self.geo.dist = float(val)
            elif self.sender().objectName() == 'rota':
                self.geo.rota = float(val)
            elif self.sender().objectName() == 'tilt':
                self.geo.tilt = float(val)
            elif self.sender().objectName() == 'yoff':
                self.geo.yoff = float(val)
            elif self.sender().objectName() == 'xoff':
                self.geo.xoff = float(val)
            elif self.sender().objectName() == 'ener':
                self.geo.ener = float(val)

        # re-calculate cones and re-draw contours
        self.draw_contours()
        # draw reference contours
        if self.geo.reference != 'None':
            self.get_reference()
            self.draw_reference()

    def dragEnterEvent(self, event):
        # Drag-and-Drop cif-file
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        # Drag-and-Drop cif-file
        #  dropEvent()
        #  - check dropped file is a cif
        fpath = event.mimeData().urls()[0].toLocalFile()
        if os.path.splitext(fpath)[1] == '.cif':
            self.get_cif_reference(fpath)

    def init_par(self):
        # fetch the geometry, detector, plot specifications and limits
        self.get_defaults()
        # file name to store current settings
        # if file_dump doesn't exists, make a dump
        if not os.path.exists(self.file_dump):
            self.save_par()
        # if it exists load parameters
        else:
            self.load_par()
        # reset to default if requested
        if self.plo.reset_settings:
            self.get_defaults()
            self.save_par()
        # update with default if requested
        #  - e.g. to add missing entries
        if self.plo.update_settings:
            self.save_par()

    def save_par(self):
        # Writing geo as dict to file
        with open(self.file_dump, 'w') as wf:
            json.dump({'geo':self.geo.__dict__, 'plo':self.plo.__dict__, 'lmt':self.lmt.__dict__}, wf, indent=4)

    def load_par(self):
        # Opening JSON file as dict
        with open(self.file_dump, 'r') as of:
            pars = json.load(of)
        conv = {'geo':self.geo, 'plo':self.plo, 'lmt':self.lmt}
        for key, vals in pars.items():
            for p, x in vals.items():
                if p in conv[key].__dict__.keys():
                    setattr(conv[key], p, x)
                else:
                    print(f'WARNING: "{p}" is not a valid key!')

class container(object):
    pass

class SliderWidget(QtWidgets.QFrame):
    def __init__(self, parent, geo, plo, lmt):
        super().__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.geo = geo
        self.plo = plo
        self.lmt = lmt
        self.leaveEvent = self.toggle_panel
        self.enterEvent = self.toggle_panel
        frame = QtWidgets.QFrame()
        frame.setContentsMargins(0, 0, 0, 0)
        frame.setFixedHeight(self.plo.slider_margin)
        layout.addWidget(frame)
        frame.setStyleSheet(f'''
            QFrame {{
                border: {self.plo.slider_border_width}px solid {self.plo.slider_border_color};
                border-radius: {self.plo.slider_border_radius}px;
                background: {self.plo.slider_border_color};
            }}
        ''')
        self.box = QtWidgets.QFrame()
        self.box.setContentsMargins(0, 5, 0, 5)
        self.box.setStyleSheet(f'''
            QFrame {{
                border: {self.plo.slider_border_width}px solid {self.plo.slider_border_color};
                border-radius: {self.plo.slider_border_radius}px;
                background: {self.plo.slider_bg_color};
            }}
            QFrame:hover {{
                background: {self.plo.slider_bg_hover};
            }}
        ''')
        layout.addWidget(self.box)
        self.box.setHidden(True)
        self.box_toggle = False
        self.box_width_dynamic = 0
        self.box_height_show = int(np.ceil(parent.size().height()/3))
        self.box_height_hide = int(np.ceil(frame.size().height()))

        grid = QtWidgets.QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setRowStretch(1,10)
        self.box.setLayout(grid)
        
        _idx = 0
        if plo.enable_slider_ener:
            self.add_slider(grid, 'Energy\n[keV]', 'ener', _idx, self.geo.ener, lmt.ener_min, lmt.ener_max, lmt.ener_stp)
            self.box_width_dynamic += self.plo.slider_column_width
            _idx += 1
        if plo.enable_slider_dist:
            self.add_slider(grid, 'Distance\n[mm]', 'dist', _idx, self.geo.dist, lmt.dist_min, lmt.dist_max, lmt.dist_stp)
            self.box_width_dynamic += self.plo.slider_column_width
            _idx += 1
        if plo.enable_slider_yoff:
            self.add_slider(grid, 'Y offset\n[mm]', 'yoff', _idx, self.geo.yoff, lmt.yoff_min, lmt.yoff_max, lmt.yoff_stp)
            self.box_width_dynamic += self.plo.slider_column_width
            _idx += 1
        if plo.enable_slider_xoff:
            self.add_slider(grid, 'X offset\n[mm]', 'xoff', _idx, self.geo.xoff, lmt.xoff_min, lmt.xoff_max, lmt.xoff_stp)
            self.box_width_dynamic += self.plo.slider_column_width
            _idx += 1
        if plo.enable_slider_tilt:
            self.add_slider(grid, 'Tilt\n[˚]', 'tilt', _idx, self.geo.tilt, lmt.tilt_min, lmt.tilt_max, lmt.tilt_stp)
            self.box_width_dynamic += self.plo.slider_column_width
            _idx += 1
        if plo.enable_slider_rota:
            self.add_slider(grid, 'Rotation\n[˚]', 'rota', _idx, self.geo.rota, lmt.rota_min, lmt.rota_max, lmt.rota_stp)
            self.box_width_dynamic += self.plo.slider_column_width
            _idx += 1
        
        self.resize(self.box_width_dynamic, self.box_height_hide)
        self.center_frame()

    def center_frame(self):
        self.move(int((self.parent().size().width()-self.box_width_dynamic)/2), self.parent().offset_win32)

    def update_slider(self, label, value):
        label.setText(str(int(value)))

    def add_slider(self, layout, label, token, idx, lval, lmin, lmax, lstp):
        font = QtGui.QFont()
        font.setPixelSize(self.plo.slider_label_size)

        label_name = QtWidgets.QLabel(label)
        label_name.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        label_name.setFont(font)
        label_name.setStyleSheet(f'''
            QLabel {{
                color: {self.plo.slider_label_color};
                border: 0px solid transparent;
            }}
        ''')
        
        label_value = QtWidgets.QLabel()
        label_value.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        label_value.setFont(font)
        label_value.setStyleSheet(f'''
            QLabel {{
                color: {self.plo.slider_label_color};
                border: 0px solid transparent;
            }}
        ''')

        slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Vertical, objectName=token)
        slider.setValue(999)
        slider.valueChanged.connect(self.parent().update_screen)
        slider.valueChanged.connect(lambda value: self.update_slider(label_value, value))
        slider.setRange(int(lmin), int(lmax))
        slider.setSingleStep(int(lstp))
        slider.setPageStep(int(lstp))
        slider.setValue(int(lval))
        
        layout.addWidget(label_name, 0, idx, QtCore.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(slider, 1, idx, QtCore.Qt.AlignmentFlag.AlignHCenter)
        layout.addWidget(label_value, 2, idx, QtCore.Qt.AlignmentFlag.AlignCenter)

        return slider

    def toggle_panel(self, event):
        if type(event) == QtGui.QEnterEvent:
            #self.box.setHidden(not self.box.isHidden())
            self.box.setHidden(False)
            self.resize(self.box_width_dynamic, self.box_height_show)
        elif type(event) == QtCore.QEvent and not self.box_toggle:
            self.box.setHidden(True)
            self.resize(self.box_width_dynamic, self.box_height_hide)
        else:
            pass

    def mousePressEvent(self, event):
        if event.button() == QtCore.Qt.MouseButton.LeftButton:
            self.startPos = event.pos()
            self.box_toggle = not self.box_toggle

    def mouseMoveEvent(self, event):
        if event.buttons() == QtCore.Qt.MouseButton.LeftButton:
            # relative movement
            delta = event.pos() - self.startPos
            # window limits
            lim = QtCore.QPoint(self.parent().size().width() - self.size().width(),
                                self.parent().size().height() - self.size().height())
            # new temporary position
            _pos = self.pos() + delta
            # x small
            if _pos.x() < 0:
                delta.setX(-self.pos().x())
            # x large
            if (lim - _pos).x() < 0:
                delta.setX(lim.x() - self.pos().x())
            # y small
            if _pos.y() < 0:
                delta.setY(-self.pos().y())
            # y large
            if (lim - _pos).y() < 0:
                delta.setY(lim.y() - self.pos().y())
            # move window to new position
            self.move(self.pos() + delta)
            # keep the box open after dragging
            self.box_toggle = True

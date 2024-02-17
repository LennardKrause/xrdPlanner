import os
import sys
import json
import glob
import shutil
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
        self.path_home = os.path.dirname(__file__)

        # add an icon
        self.setWindowIcon(QtGui.QIcon(':/icons/xrdPlanner.png'))
        #self.setMouseTracking(True)

        # enable antialiasing
        pg.setConfigOptions(antialias=True, imageAxisOrder='row-major')

        # Drag-and-Drop cif-file
        #  dropEvent()
        #  - check dropped file is a cif
        #  calc_ref_from_cif()
        #  - use Dans_Diffraction to get d_spacings
        #  - dif.Crystal()
        #  - xtl.Scatter.powder()
        self.setAcceptDrops(True)

        # set path to settings folder
        self.path_settings = os.path.join(self.path_home, 'settings',)
        # set path to active settings token
        self.path_settings_token = os.path.join(self.path_settings, '.active',)
        # set path to default settings file
        self.path_settings_default = os.path.join(self.path_settings, 'default.json')
        # set path to current settings file
        self.active_settings = self.get_active_settings_file()
        self.path_settings_current = os.path.join(self.path_settings, self.active_settings)
        # set path to detector database
        self.path_detdb = os.path.join(self.path_home, 'detector_db.json')

        # move settings file from old location
        _old_settings_file = os.path.join(self.path_home, 'settings.json')
        if os.path.exists(_old_settings_file):
            shutil.move(_old_settings_file, self.path_settings_default)

        # delete active settings file token
        #  - token with currently active settings
        #    file will be placed on successful exit
        #  - if no token file is found -> crashed last time
        #  - load default settings
        self.delete_active_settings_file()

        # save/load parameters to/from file
        self.init_par()

        # Store available colormaps
        # PyQtGraph.colormap.listMaps(): Experimental, subject to change.
        self.colormaps = sorted(pg.colormap.listMaps())

        # get the translations to link
        # the settings parameter keys
        # to their description
        self.get_tooltips()

        # menubar is displayed within the main window on Windows
        # so we need to make space for it
        # no idea about other OS, if there are issues fix them here
        if sys.platform in ['win32', 'linux', 'linux2']:
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

        # Add current beamstop size to list
        if self.geo.bssz and not isinstance(self.geo.bssz, str) and self.geo.bssz not in self.geo.bs_list:
            self.geo.bs_list.append(self.geo.bssz)
            self.geo.bs_list.sort()

        # What standards should be available as reference
        # The d spacings will be imported from pyFAI
        self.ref_library = calibrant.names()
        # dict to store custom reference data
        self.ref_custom = {}
        self.ref_custom_hkl = {}

        # initialise all that depends on the settings
        # call this function to apply changes were made
        # to the settings file -> reload_settings()
        self.init_modifiables()
        
        # populate the menus with detectors, references and units
        self.init_menus()
        
        # initialize the screen
        self.init_screen()

        # disable/reset delta_d/d toggles
        #self.action_funct_deltadd_show.setEnabled(False)

        # add the slider frame
        # this calls draw_conics(), 
        # make sure that everything
        # that is needed is initialised
        self.sliderWidget = SliderWidget(self)

    def init_modifiables(self):
        # get colormap
        self.cont_cmap = pg.colormap.get(self.geo.colormap, skipCache=True)
        # reverse the colormap useful to increase visibility in darkmode
        if self.geo.darkmode:
            self.cont_cmap.reverse()

        # experimental darkmode?
        self.apply_theme(self.geo.darkmode)

        # set window color
        self.ax.setBackground(self.plot_bg_color)

        # generate contour levels
        self.cont_geom_num = np.linspace(self.plo.conic_tth_min, self.plo.conic_tth_max, self.plo.conic_tth_num)

        # translate unit for plot title
        self.unit_names = ['2\u03B8 [\u00B0]',
                           'd [\u212B]',
                           'q [\u212B\u207B\u00B9]',
                           'sin(\u03B8)/\u03BB [\u212B\u207B\u00B9]']
        if self.geo.unit >= len(self.unit_names):
            print(f'Error: Valid geo.unit range is from 0 to {len(self.unit_names)-1}, geo.unit={self.geo.unit}')
            raise SystemExit

        # get the detector specs
        # - update: overwrite existing file after load
        # - reset: overwrite existing file with defaults
        self.detector_db = self.get_det_library(update=self.plo.update_det_bank, reset=self.plo.reset_det_bank)

        # pick current detector
        self.det = self.get_specs_det()
        
        # remove unavailable detectors
        if self.geo.det_bank:
            # add chosen detector/size combination to the bank
            if self.geo.det_type not in self.geo.det_bank.keys():
                self.geo.det_bank[self.geo.det_type] = [self.geo.det_size]
            else:
                # if entry is str -> list
                if isinstance(self.geo.det_bank[self.geo.det_type], str):
                    self.geo.det_bank[self.geo.det_type] = [self.geo.det_bank[self.geo.det_type]]
                # if current not in list -> append
                if self.geo.det_size not in self.geo.det_bank[self.geo.det_type]:
                    self.geo.det_bank[self.geo.det_type].append(self.geo.det_size)
            # iterate over detector db
            for det in list(self.detector_db.keys()):
                if det not in self.geo.det_bank.keys():
                    self.detector_db.pop(det)
                else:
                    for size in list(self.detector_db[det]['size'].keys()):
                        if size not in self.geo.det_bank[det]:
                            self.detector_db[det]['size'].pop(size)

        # init the hkl tooltip
        font = QtGui.QFont()
        font.setPixelSize(self.plo.conic_hkl_label_size)
        font.setBold(True)
        QtWidgets.QToolTip.setFont(font)

    def init_par(self):
        # fetch the geometry, detector, plot specifications and limits
        self.get_defaults_all()
        # file name to store current settings
        # if file_dump doesn't exists, make a dump
        if not os.path.exists(self.path_settings_current):
            self.save_settings()
        # if it exists load parameters
        else:
            self.load_settings()
        # reset to default if requested
        if self.plo.reset_settings:
            self.get_defaults_all()
            self.save_settings()
        # update with default if requested
        #  - e.g. to add missing entries
        if self.plo.update_settings:
            self.save_settings()

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

        # define font to label conics
        font = QtGui.QFont()
        font.setPixelSize(self.plo.conic_label_size)
        font.setBold(True)

        # container for contour lines
        self.patches = {'beamcenter':None,
                        'beamstop':None,
                        'overlay':None,
                        #'deltadd':None,
                        'poni':None,
                        'bs_label':None,
                        'conic':[],
                        'reference':[],
                        'labels':[]}

        # add beam stop scatter plot
        # pxMode=False
        # Size is in scene coordinates and the spots scale with the
        # view. Otherwise, spots are always the same size regardless
        # of scaling, and size is given in px.
        self.patches['beamstop'] = pg.PlotDataItem(useCache=True,
                                                   pen=pg.mkPen(self.beamstop_edge_color,
                                                                width=self.plo.conic_linewidth),
                                                   brush=pg.mkBrush(self.beamstop_color),
                                                   fillOutline=True)
        self.ax.addItem(self.patches['beamstop'])

        # add empty plot per reference contour line
        for i in range(self.plo.conic_ref_num):
            ref = pg.PlotCurveItem(useCache=True)
            self.ax.addItem(ref)
            self.patches['reference'].append(ref)
            self.patches['reference'][i].setClickable(True, width=self.plo.conic_ref_linewidth)
            self.patches['reference'][i].sigClicked.connect(self.show_tooltip)
            self.patches['reference'][i].name = None

        # add isocurve for delta d/d
        #self.patches['deltadd'] = QtWidgets.QGraphicsEllipseItem()
        #_pen = pg.mkPen((255,255,255,255), width=self.plo.conic_linewidth, style=QtCore.Qt.PenStyle.DashLine)
        #_pen.setCapStyle(QtCore.Qt.PenCapStyle.FlatCap)
        #_pen.setDashPattern([4,8])
        #self.patches['deltadd'].setPen(_pen)
        #self.ax.addItem(self.patches['deltadd'])

        # add empty plot per contour line
        for i in range(self.plo.conic_tth_num):
            curve = pg.PlotCurveItem(useCache=True)
            self.ax.addItem(curve)
            self.patches['conic'].append(curve)
            temp_label = pg.TextItem(anchor=(0.5,0.5), fill=pg.mkBrush(self.conic_label_fill))
            temp_label.setFont(font)
            self.patches['labels'].append(temp_label)
            self.ax.addItem(temp_label)

        # label for beamstop contour
        bs_label = pg.TextItem(anchor=(0.5,0.5),
                               color=self.unit_label_color,
                               fill=self.conic_label_fill)
        bs_label.setFont(font)
        self.patches['bs_label'] = bs_label
        self.ax.addItem(bs_label)

        # add poni scatter plot
        self.patches['poni'] = pg.ScatterPlotItem(symbol = self.plo.poni_marker,
                                                  size = self.plo.poni_size,
                                                  brush = pg.mkBrush(self.cont_cmap.map(0, mode='qcolor')),
                                                  pen = pg.mkPen(None))
        self.ax.addItem(self.patches['poni'])

        # add beam center scatter plot
        self.patches['beamcenter'] = pg.ScatterPlotItem(symbol = self.plo.beamcenter_marker,
                                                        size = self.plo.beamcenter_size,
                                                        brush = pg.mkBrush(self.cont_cmap.map(0, mode='qcolor')),
                                                        pen = pg.mkPen(None))
        self.ax.addItem(self.patches['beamcenter'])

        self.patches['overlay'] = pg.ImageItem()
        self.patches['overlay'].hoverEvent = self.hoverEvent
        self.ax.addItem(self.patches['overlay'])
        # set overlay colors
        _high = QtGui.QColor(self.plot_bg_color)
        _high.setAlphaF(0.0)
        if self.plo.overlay_toggle_warn:
            cm = pg.ColorMap(pos=[0.0, 1.0], color=[self.overlay_threshold_color, _high])
        else:
            cm = pg.ColorMap(pos=[0.0, 1.0], color=[self.plot_bg_color, _high])
        cm.setMappingMode('clip')
        self.patches['overlay'].setColorMap(cm)

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
        self.init_unit_label()

        # add label for corrections
        self.init_correction_label()

        # create cones and draw contour lines
        self.update_screen()
        self.set_window_title()

    def init_menus(self, reset=False):

        if reset:
            self.menu_bar.clear()
        # if xrdPlanner is added as a widget to a
        # GUI use and append to the parent menuBar
        if self.parent():
            self.menu_bar = self.parent().menuBar()
        else:
            self.menu_bar = self.menuBar()

        # append 'menu' and 'value' (string!) as tuple to this list
        # update_menu_entries() will then reset the
        # checkmark to the active entry after reload of
        # settings.

        # self.update_menu_entries() will access
        # the menus and update the checkmarks upon
        # settings reload via self.reload_settings()
        menu_det = self.menu_bar.addMenu('Detector')
        group_det = QtGui.QActionGroup(self)
        group_det.setExclusive(True)

        #############
        # DETECTORS #
        #############
        # menu Detectors
        for d in sorted(self.detector_db):
            d_menu = QtWidgets.QMenu(d, self)
            menu_det.addMenu(d_menu)
            for s in self.detector_db[d]['size']:
                det_action = QtGui.QAction(s, self, checkable=True)
                self.set_menu_action(det_action, self.change_detector, d, s)
                d_menu.addAction(det_action)
                group_det.addAction(det_action)
                if d == self.geo.det_type and s == self.geo.det_size:
                    det_action.setChecked(True)
        
        #############
        # REFERENCE #
        #############
        # menu Reference
        menu_ref = self.menu_bar.addMenu('Reference')
        self.group_ref = QtGui.QActionGroup(self)
        self.group_ref.setExclusive(True)
        # menu Reference: add None
        ref_action = QtGui.QAction('None', self, checkable=True)
        self.set_menu_action(ref_action, self.change_reference, 'None')
        menu_ref.addAction(ref_action)
        self.group_ref.addAction(ref_action)
        if self.geo.reference.lower() == 'none':
            ref_action.setChecked(True)
        # menu Reference: add pyFAI library
        sub_menu_pyFAI = QtWidgets.QMenu('pyFAI', self)
        menu_ref.addMenu(sub_menu_pyFAI)
        for ref_name in sorted(self.ref_library):
            ref_action = QtGui.QAction(ref_name, self, checkable=True)
            self.set_menu_action(ref_action, self.change_reference, ref_name)
            sub_menu_pyFAI.addAction(ref_action)
            self.group_ref.addAction(ref_action)
            if ref_name == self.geo.reference:
                ref_action.setChecked(True)

        # menu Reference: add Custom
        self.sub_menu_custom = QtWidgets.QMenu('Custom', self)
        menu_ref.addMenu(self.sub_menu_custom)
        
        ############
        # BEAMSTOP #
        ############
        # menu Beamstop
        menu_bs = self.menu_bar.addMenu('Beamstop')
        group_bs = QtGui.QActionGroup(self)
        group_bs.setExclusive(True)
        # menu Beamstop: add None
        bs_action = QtGui.QAction('None', self, checkable=True)
        self.set_menu_action(bs_action, self.change_beamstop, 'None')
        menu_bs.addAction(bs_action)
        group_bs.addAction(bs_action)
        if isinstance(self.geo.bssz, str) and self.geo.bssz.lower() == 'none':
            bs_action.setChecked(True)
        # menu Beamstop: add sizes list
        sub_menu_bs = QtWidgets.QMenu('Sizes [mm]', self)
        menu_bs.addMenu(sub_menu_bs)
        for bs_size in sorted(self.geo.bs_list):
            bs_sub_action = QtGui.QAction(str(bs_size), self, checkable=True)
            self.set_menu_action(bs_sub_action, self.change_beamstop, bs_size)
            sub_menu_bs.addAction(bs_sub_action)
            group_bs.addAction(bs_sub_action)
            if bs_size == self.geo.bssz:
                bs_sub_action.setChecked(True)
        
        ###########
        #  UNITS  #
        ###########
        # menu Units
        # non-standard menu -> self.geo.unit is int
        self.menu_unit = self.menu_bar.addMenu('Unit')
        group_unit = QtGui.QActionGroup(self)
        group_unit.setExclusive(True)
        for unit_index, unit_name in enumerate(self.unit_names):
            unit_action = QtGui.QAction(unit_name, self, checkable=True)
            self.set_menu_action(unit_action, self.change_units, unit_index)
            self.menu_unit.addAction(unit_action)
            group_unit.addAction(unit_action)
            if unit_index == self.geo.unit:
                unit_action.setChecked(True)

        ############
        #   VIEW   #
        ############
        # menu View
        menu_view = self.menu_bar.addMenu('View')
        ################
        # VIEW - THEME #
        ################
        # non-standard menu -> self.geo.darkmode is bool
        menu_theme = menu_view.addMenu('Theme')
        group_theme = QtGui.QActionGroup(self)
        group_theme.setExclusive(True)
        for (theme, invert) in [('Light', False), ('Dark', True)]:
            theme_action = QtGui.QAction(theme, self, checkable=True)
            self.set_menu_action(theme_action, self.apply_theme, invert, True)
            group_theme.addAction(theme_action)
            menu_theme.addAction(theme_action)
            if invert == self.geo.darkmode:
                theme_action.setChecked(True)

        ###################
        # VIEW - COLORMAP #
        ###################
        self.menu_cmap = menu_view.addMenu('Colormap')
        group_cmap = QtGui.QActionGroup(self)
        group_cmap.setExclusive(True)
        for cmap_name in self.colormaps: # PyQtGraph.colormap.listMaps(): Experimental, subject to change.
            cmap_action = QtGui.QAction(cmap_name, self, checkable=True)
            self.set_menu_action(cmap_action, self.change_cmap, cmap_name)
            self.menu_cmap.addAction(cmap_action)
            group_cmap.addAction(cmap_action)
            if cmap_name == self.geo.colormap:
                cmap_action.setChecked(True)

        ##################
        # VIEW - OVERLAY #
        ##################
        menu_overlays = menu_view.addMenu('Overlay')
        # unit value hover toggle
        self.action_unit_hover = QtGui.QAction('Unit hover', self, checkable=True)
        self.set_menu_action(self.action_unit_hover, self.toggle_unit_hover)
        if self.plo.show_unit_hover:
            self.action_unit_hover.setChecked(True)
        else:
            self.action_unit_hover.setChecked(False)
        menu_overlays.addAction(self.action_unit_hover)
        menu_overlays.addSeparator()
        # polarisation map toggle
        self.action_show_pol = QtGui.QAction('Polarisation', self, checkable=True)
        self.set_menu_action(self.action_show_pol, self.toggle_overlay_polarisation)
        if self.plo.show_polarisation:
            self.action_show_pol.setChecked(True)
        else:
            self.action_show_pol.setChecked(False)
        menu_overlays.addAction(self.action_show_pol)
        # Solidangle map toggle
        self.action_show_ang = QtGui.QAction('Solid angle', self, checkable=True)
        self.set_menu_action(self.action_show_ang, self.toggle_overlay_solidangle)
        if self.plo.show_solidangle:
            self.action_show_ang.setChecked(True)
        else:
            self.action_show_ang.setChecked(False)
        menu_overlays.addAction(self.action_show_ang)
        # Overlay warn color toggle
        menu_overlays.addSeparator()
        self.action_overlay_warn = QtGui.QAction('Highlight', self, checkable=True)
        self.set_menu_action(self.action_overlay_warn, self.toggle_overlay_highlight)
        if self.plo.overlay_toggle_warn:
            self.action_overlay_warn.setChecked(True)
        else:
            self.action_overlay_warn.setChecked(False)
        menu_overlays.addAction(self.action_overlay_warn)

        ####################
        # VIEW - FUNCTIONS #
        ####################
        #menu_functions = menu_view.addMenu('Functions')
        # set deltadd parameters toggle
        #self.action_funct_deltadd_set = QtGui.QAction('Setup \u03B4d/d', self)
        #self.set_menu_action(self.action_funct_deltadd_set, self.window_function_deltadd)
        #menu_functions.addAction(self.action_funct_deltadd_set)
        # show deltadd toggle
        #self.action_funct_deltadd_show = QtGui.QAction('Show \u03B4d/d', self, checkable=True)
        #self.set_menu_action(self.action_funct_deltadd_show, self.toggle_function_deltadd)
        #if self.plo.show_deltadd:
        #    self.action_funct_deltadd_show.setChecked(True)
        #else:
        #    self.action_funct_deltadd_show.setChecked(False)
        #menu_functions.addAction(self.action_funct_deltadd_show)

        ############
        # SETTINGS #
        ############
        # menu Settings
        menu_Settings = self.menu_bar.addMenu('Settings')
        # submenu load settings files
        self.menu_custom_settings = menu_Settings.addMenu('Load')
        group_cset = QtGui.QActionGroup(self)
        group_cset.setExclusive(True)
        for cset_name in self.get_settings_files():
            cset_action = QtGui.QAction(cset_name, self, checkable=True)
            self.set_menu_action(cset_action, self.change_settings_file, cset_name)
            self.menu_custom_settings.addAction(cset_action)
            group_cset.addAction(cset_action)
            if cset_name == self.active_settings:
                cset_action.setChecked(True)
        
        ###################
        # SETTINGS - SAVE #
        ###################
        save_action = QtGui.QAction('Save', self)
        self.set_menu_action(save_action, self.save_current_settings)
        menu_Settings.addAction(save_action)
        #####################
        # SETTINGS - IMPORT #
        #####################
        import_action = QtGui.QAction('Import', self)
        self.set_menu_action(import_action, self.import_settings_file)
        menu_Settings.addAction(import_action)
        #####################
        # SETTINGS - EXPORT #
        #####################
        export_action = QtGui.QAction('Export', self)
        self.set_menu_action(export_action, self.show_export_window)
        menu_Settings.addAction(export_action)
        ###################
        # SETTINGS - EDIT #
        ###################
        menu_Settings.addSeparator()
        menu_Edit = menu_Settings.addMenu('Edit')
        # prepare platform dependent file reader
        if sys.platform == 'win32':
            tokens = [('Detector db', os.system, f'notepad {self.path_detdb}'),
                      ('Settings', self.edit_settings_file, f'notepad')]
        elif sys.platform == 'linux':
            tokens = [('Detector db', os.system, f'xdg-open {self.path_detdb}'),
                      ('Settings', self.edit_settings_file, f'xdg-open')]
        else:
            tokens = [('Detector db', os.system, f'open -t {self.path_detdb}'),
                      ('Settings', self.edit_settings_file, f'open -t')]
        for (name, funct, command) in tokens:
            edit_action = QtGui.QAction(name, self)
            self.set_menu_action(edit_action, funct, command)
            menu_Edit.addAction(edit_action)
        #####################
        # SETTINGS - DELETE #
        #####################
        menu_Settings.addSeparator()
        export_action = QtGui.QAction('Delete', self)
        self.set_menu_action(export_action, self.delete_settings_files)
        menu_Settings.addAction(export_action)

        # menu About
        menu_help = self.menu_bar.addMenu('Help')
        action_about = QtGui.QAction('xrdPlanner', self)
        self.set_menu_action(action_about, self.show_about_window)
        menu_help.addAction(action_about)

    def init_unit_label(self):
        font = QtGui.QFont()
        font.setPixelSize(self.plo.unit_label_size)
        self.unit_label = pg.TextItem(anchor=(0.0,0.0), color=self.unit_label_color, fill=self.unit_label_fill)
        self.unit_label.setText(self.unit_names[self.geo.unit])
        self.unit_label.setFont(font)
        self.ax.addItem(self.unit_label)
        self.unit_label.setPos(-self.xdim, self.ydim)

    def init_correction_label(self):
        font = QtGui.QFont()
        font.setPixelSize(self.plo.unit_label_size)
        self.cor_label = pg.TextItem(anchor=(1.0,1.0), color=self.unit_label_color, fill=self.unit_label_fill)
        self.cor_label.setText('')
        self.cor_label.setFont(font)
        self.ax.addItem(self.cor_label)
        self.cor_label.hide()
        self.cor_label.setToolTip('Correction factors:\nP: Polarisation\nS: Solid angle')
        self.cor_label.setPos(self.xdim, -self.ydim)

    def apply_theme(self, use_dark, redraw=False):
        # set darkmode
        self.geo.darkmode = use_dark
        _color_dark = QtGui.QColor(self.thm.color_dark)
        _color_light = QtGui.QColor(self.thm.color_light)
        # set highlight text color depending on the lightness of the colormap
        _highlight_color = self.cont_cmap.map(0.0, mode='qcolor')
        if _highlight_color.lightnessF() < 0.5:
            _highlight_text = _color_light.lighter(150)
        else:
            _highlight_text = _color_dark.darker(150)
        # define color palette
        if use_dark:
            # reference contour
            self.conic_label_fill = self.thm.dark_conic_label_fill
            self.conic_ref_color = self.thm.dark_conic_ref_color
            self.det_module_color = self.thm.dark_det_module_color
            self.det_module_fill = self.thm.dark_det_module_fill
            # general
            self.plot_bg_color = self.thm.dark_plot_bg_color
            self.beamstop_color = self.thm.dark_beamstop_color
            self.beamstop_edge_color = self.thm.dark_beamstop_edge_color
            self.unit_label_color = self.thm.dark_unit_label_color
            self.unit_label_fill = self.thm.dark_unit_label_fill
            self.overlay_threshold_color = self.thm.dark_overlay_threshold_color
            # slider
            self.slider_border_color = self.thm.dark_slider_border_color
            self.slider_bg_color = self.thm.dark_slider_bg_color
            self.slider_bg_hover = self.thm.dark_slider_bg_hover
            self.slider_label_color = self.thm.dark_slider_label_color
            # palette
            palette = QtGui.QPalette()
            palette.setColor(QtGui.QPalette.ColorRole.Window,          _color_dark)
            palette.setColor(QtGui.QPalette.ColorRole.WindowText,      _color_light)
            palette.setColor(QtGui.QPalette.ColorRole.Button,          _color_dark)
            palette.setColor(QtGui.QPalette.ColorRole.ButtonText,      _color_light)
            palette.setColor(QtGui.QPalette.ColorRole.Base,            _color_dark)
            palette.setColor(QtGui.QPalette.ColorRole.AlternateBase,   _color_dark.lighter(110))
            palette.setColor(QtGui.QPalette.ColorRole.Highlight,       _highlight_color)
            palette.setColor(QtGui.QPalette.ColorRole.HighlightedText, _highlight_text)
            palette.setColor(QtGui.QPalette.ColorRole.Text,            _color_light)
            palette.setColor(QtGui.QPalette.ColorRole.PlaceholderText, _color_light)
            palette.setColor(QtGui.QPalette.ColorRole.BrightText,      _color_light)
            palette.setColor(QtGui.QPalette.ColorRole.ToolTipBase,     _color_dark)
            palette.setColor(QtGui.QPalette.ColorRole.ToolTipText,     _color_light)
        else:
            # reference contour
            self.conic_label_fill = self.thm.light_conic_label_fill
            self.conic_ref_color = self.thm.light_conic_ref_color
            self.det_module_color = self.thm.light_det_module_color
            self.det_module_fill = self.thm.light_det_module_fill
            # general
            self.plot_bg_color = self.thm.light_plot_bg_color
            self.beamstop_color = self.thm.light_beamstop_color
            self.beamstop_edge_color = self.thm.light_beamstop_edge_color
            self.unit_label_color = self.thm.light_unit_label_color
            self.unit_label_fill = self.thm.light_unit_label_fill
            self.overlay_threshold_color = self.thm.light_overlay_threshold_color
            # slider
            self.slider_border_color = self.thm.light_slider_border_color
            self.slider_bg_color = self.thm.light_slider_bg_color
            self.slider_bg_hover = self.thm.light_slider_bg_hover
            self.slider_label_color = self.thm.light_slider_label_color

            # palette
            palette = QtGui.QPalette()
            palette.setColor(QtGui.QPalette.ColorRole.Window,          _color_light)
            palette.setColor(QtGui.QPalette.ColorRole.WindowText,      _color_dark)
            palette.setColor(QtGui.QPalette.ColorRole.Button,          _color_light)
            palette.setColor(QtGui.QPalette.ColorRole.ButtonText,      _color_dark)
            palette.setColor(QtGui.QPalette.ColorRole.Base,            _color_light)
            palette.setColor(QtGui.QPalette.ColorRole.AlternateBase,   _color_light.darker(110))
            palette.setColor(QtGui.QPalette.ColorRole.Highlight,       _highlight_color)
            palette.setColor(QtGui.QPalette.ColorRole.HighlightedText, _highlight_text)
            palette.setColor(QtGui.QPalette.ColorRole.Text,            _color_dark)
            palette.setColor(QtGui.QPalette.ColorRole.PlaceholderText, _color_dark)
            palette.setColor(QtGui.QPalette.ColorRole.BrightText,      _color_dark)
            palette.setColor(QtGui.QPalette.ColorRole.ToolTipBase,     _color_light)
            palette.setColor(QtGui.QPalette.ColorRole.ToolTipText,     _color_dark)

        # apply palette to app
        app = QtWidgets.QApplication.instance()
        app.setPalette(palette)
        app.setStyle('Fusion')
        
        # redraw the canvas
        if redraw:
            self.redraw_canvas()

    #############
    #  TOGGLE   #
    #############
    def toggle_unit_hover(self):
        self.plo.show_unit_hover = not self.plo.show_unit_hover
        if self.plo.show_unit_hover:
            self.action_unit_hover.setChecked(True)
        else:
            self.action_unit_hover.setChecked(False)
        self.redraw_canvas()

    def toggle_overlay_polarisation(self):
        self.plo.show_polarisation = not self.plo.show_polarisation
        if self.plo.show_polarisation:
            self.action_show_pol.setChecked(True)
        else:
            self.action_show_pol.setChecked(False)
        self.redraw_canvas()

    def toggle_overlay_solidangle(self):
        self.plo.show_solidangle = not self.plo.show_solidangle
        if self.plo.show_solidangle:
            self.action_show_ang.setChecked(True)
        else:
            self.action_show_ang.setChecked(False)
        self.redraw_canvas()

    def toggle_overlay_highlight(self):
        self.plo.overlay_toggle_warn = not self.plo.overlay_toggle_warn
        if self.plo.overlay_toggle_warn:
            self.action_overlay_warn.setChecked(True)
        else:
            self.action_overlay_warn.setChecked(False)
        self.redraw_canvas()

    #def toggle_function_deltadd(self):
    #    if not self.action_funct_deltadd_show.isEnabled():
    #        return
    #    self.plo.show_deltadd = not self.plo.show_deltadd
    #    if self.plo.show_deltadd:
    #        self.action_funct_deltadd_show.setChecked(True)
    #    else:
    #        self.action_funct_deltadd_show.setChecked(False)
    #    self.redraw_canvas()

    #############
    #   REDO    #
    #############
    def redraw_canvas(self):
        # save darkmode toggle
        #self.save_settings()
        self.init_modifiables()
        # clear the screen
        self.ax.clear()
        # re-initialise
        self.init_screen()
        # center the slider frame
        self.sliderWidget.apply_style()
        self.sliderWidget.init_sliders()
        self.sliderWidget.center_frame()

    def reload_settings(self):
        # load settings
        ###############################
        # clearing the det_bank makes #
        # sure that 'old' par files   #
        # without the det_bank entry  #
        # get access to all detectors #
        # -> backwards compatibility. #
        ###############################
        self.geo.det_bank = {}
        ###############################
        self.load_settings()
        # missing entries in the
        # settings file will be added
        self.save_settings()
        self.redraw_canvas()
        self.init_menus(reset=True)
        #self.update_menu_entries()

    #############
    #  CHANGE   #
    #############
    def change_beamstop(self, size):
        self.geo.bssz = size
        self.redraw_canvas()

    def change_cmap(self, cmap):
        self.geo.colormap = cmap
        for action in self.menu_cmap.actions():
            if action.text() == self.geo.colormap:
                action.setChecked(True)
        self.redraw_canvas()

    def change_detector(self, det_name, det_size):
        self.geo.det_type = det_name
        self.geo.det_size = det_size
        # get new detector specs
        self.det = self.get_specs_det()
        # clear the screen
        self.ax.clear()
        # re-initialise
        self.init_screen()
        # center the slider frame
        self.sliderWidget.center_frame()

    def change_units(self, unit_index):
        self.geo.unit = unit_index
        self.unit_label.setText(self.unit_names[unit_index])
        for num, action in enumerate(self.menu_unit.actions()):
            if num == self.geo.unit:
                action.setChecked(True)
        self.draw_conics()

    def change_reference(self, ref_name):
        self.geo.reference = ref_name
        self.get_reference()
        self.draw_reference()

    #############
    #    GET    #
    #############
    def get_reference(self):
        if self.geo.reference in self.ref_library:
            # get the d spacings for the calibrtant from pyFAI
            self.cont_ref_dsp = np.array(calibrant.get_calibrant(self.geo.reference).get_dSpacing()[:self.plo.conic_ref_num])
            self.cont_ref_hkl = None
        elif self.geo.reference in self.ref_custom:
            # get custom d spacings
            self.cont_ref_dsp = self.ref_custom[self.geo.reference]
            self.cont_ref_hkl = self.ref_custom_hkl[self.geo.reference]
        else:
            # set all d-spacings to -1
            self.cont_ref_dsp = np.zeros(self.plo.conic_ref_num)
            self.cont_ref_hkl = None
        # update window title
        self.set_window_title()

    def get_defaults_geo(self):
        ######################
        # Setup the geometry #
        ######################
        geo = container()
        geo.det_type = 'EIGER2'  # [str]  Pilatus3 / Eiger2
        geo.det_size = '4M'      # [str]  300K 1M 2M 6M / 1M 4M 9M 16M
        geo.ener = 25            # [keV]  Beam energy
        geo.dist = 100           # [mm]   Detector distance
        geo.voff = 0             # [mm]   Detector offset (vertical)
        geo.hoff = 0             # [mm]   Detector offset (horizontal)
        geo.rota = 0             # [deg]  Detector rotation
        geo.tilt = 0             # [deg]  Detector tilt
        geo.bssz = 3.0           # [mm]   Current beamstop size (or None)
        geo.bsdx = 40            # [mm]   Beamstop distance
        geo.unit = 1             # [0-3]  Contour legend
                                 #          0: 2-Theta
                                 #          1: d-spacing
                                 #          2: q-space
                                 #          3: sin(theta)/lambda
        geo.reference = 'None'   # [str]  Plot reference contours
                                 #          pick from pyFAI
        geo.darkmode = True      # [bool] Darkmode
        geo.colormap = 'viridis' # [cmap] Contour colormap
        geo.bs_list = [1.5,      # [list] Available beamstop sizes
                       2.0,
                       2.5,
                       3.0,
                       5.0]
        geo.det_bank = {}        # available detectors
                                 # 'key':['value'] or ['value1', 'value2']
                                 # 'key': Detector name
                                 # 'value': Detector model/size/type
                                 # e.g. {'PLATUS3:['1M', '2M'], 'EIGER2':['4M']}
                                 # empty dict enables all detectors
        return geo

    def get_defaults_plo(self):
        ################
        # Plot Details #
        ################
        plo = container()
        # - geometry contour section - 
        plo.conic_tth_min = 5               # [int]    Minimum 2-theta contour line
        plo.conic_tth_max = 100             # [int]    Maximum 2-theta contour line
        plo.conic_tth_num = 15              # [int]    Number of contour lines
        plo.beamcenter_marker = 'o'         # [marker] Beamcenter marker
        plo.beamcenter_size = 6             # [int]    Beamcenter size
        plo.poni_marker = 'x'               # [marker] Poni marker
        plo.poni_size = 8                   # [int]    Poni size
        plo.conic_linewidth = 2.0           # [float]  Contour linewidth (lw)
        plo.conic_label_size = 14           # [int]    Contour labelsize
        # - reference contour section - 
        plo.conic_ref_linewidth = 2.0       # [float]  Reference contour linewidth
        plo.conic_ref_num = 250             # [int]    Number of reference contours
        plo.conic_ref_cif_int = 0.01        # [float]  Minimum display intensity (cif)
        plo.conic_ref_cif_kev = 10.0        # [float]  Energy [keV] for intensity calculation
        plo.conic_ref_cif_irel = True       # [bool]   Linewidth relative to intensity
        plo.conic_ref_cif_lw_min = 0.1      # [float]  Minimum linewidth when using irel
        plo.conic_ref_cif_lw_mult = 3.0     # [float]  Linewidth multiplier when using irel
        plo.conic_hkl_show_int = False      # [bool]   Show intensity in hkl tooltip
        plo.conic_hkl_label_size = 14       # [int]    Font size of hkl tooltip
        # - module section - 
        plo.det_module_alpha = 0.20         # [float]  Detector module alpha
        plo.det_module_width = 1            # [int]    Detector module border width
        # - general section - 
        plo.conic_steps = 100               # [int]    Conic resolution
        plo.plot_size = 0                   # [int]    Plot size, px (0 for auto)
        plo.plot_size_fixed = True          # [bool]   Fix window size
        plo.unit_label_size = 16            # [int]    Label size, px
        plo.polarisation_fac = 0.99         # [float]  Horizontal polarisation factor
        plo.show_polarisation = True        # [bool]   Show polarisation overlay
        plo.show_solidangle = False         # [bool]   Show solid angle overlay
        plo.show_unit_hover = True          # [bool]   Show unit value on hover
        plo.overlay_resolution = 300        # [int]    Overlay resolution
        plo.overlay_toggle_warn = True      # [bool]   Overlay warn color threshold
        # - extra functions -
        #plo.show_deltadd = False            # [bool]   Show delta_d/d function
        #plo.sensor_thickness = 100e-6       # [float]  Detector sensor thickness [m]
        #plo.beam_divergence = 10e-6         # [float]  X-ray beam divergence [rad]
        #plo.scattering_diameter = 100e-6    # [float]  Scattering volume diameter [m]
        #plo.energy_resolution = 100e-6      # [float]  X-ray beam resolution [eV/keV]
        #plo.funct_deltadd_thresh = -4       # [float]  Log10 threshold for delta d/d contour
        # - slider section - 
        plo.slider_margin = 12              # [int]    Slider frame top margin
        plo.slider_border_width = 1         # [int]    Slider frame border width
        plo.slider_border_radius = 1        # [int]    Slider frame border radius (px)
        plo.slider_label_size = 14          # [int]    Slider frame label size
        plo.slider_column_width = 75        # [int]    Slider label column width
        plo.enable_slider_ener = True       # [bool]   Show energy slider
        plo.enable_slider_dist = True       # [bool]   Show distance slider
        plo.enable_slider_rota = True       # [bool]   Show rotation slider
        plo.enable_slider_voff = True       # [bool]   Show vertical offset slider
        plo.enable_slider_hoff = True       # [bool]   Show horizontal offset slider
        plo.enable_slider_tilt = True       # [bool]   Show tilt slider
        plo.enable_slider_bsdx = True       # [bool]   Show beamstop distance slider
        plo.slider_label_ener = 'Energy\n[keV]'            # [str] Label for energy slider
        plo.slider_label_dist = 'Distance\n[mm]'           # [str] Label for distance slider
        plo.slider_label_rota = 'Rotation\n[\u02da]'       # [str] Label for rotation slider
        plo.slider_label_voff = 'Vertical\noffset\n[mm]'   # [str] Label for vertical offset slider
        plo.slider_label_hoff = 'Horizontal\noffset\n[mm]' # [str] Label for horizontal offset slider
        plo.slider_label_tilt = 'Tilt\n[\u02da]'           # [str] Label for tilt slider
        plo.slider_label_bsdx = 'Beamstop\ndistance\n[mm]' # [str] Label for beamstop distance slider
        # - update/reset - 
        plo.update_settings = True          # [bool]   Update settings file after load
        plo.update_det_bank = True          # [bool]   Update detector bank after load
        plo.reset_settings = False          # [bool]   Reset settings file
        plo.reset_det_bank = False          # [bool]   Reset detector bank
        # - debug/testing -
        plo.set_debug = False               # [bool]   Debug mode

        return plo
    
    def get_defaults_thm(self):
        #################
        # Theme Details #
        #################
        thm = container()
        thm.color_dark = '#404040'                    # [color]  Global dark color
        thm.color_light = '#EEEEEE'                   # [color]  Global light color
        # light mode
        thm.light_conic_label_fill = '#FFFFFF'        # [color]  Contour label fill color
        thm.light_conic_ref_color = '#DCDCDC'         # [color]  Reference contour color
        thm.light_beamstop_color = '#FF000080'        # [color]  Beamstop color
        thm.light_beamstop_edge_color = '#FF0000'     # [color]  Beamstop edge color
        thm.light_det_module_color = '#404040'        # [color]  Detector module border color
        thm.light_det_module_fill = '#404040'         # [color]  Detector module background color
        thm.light_plot_bg_color = '#FFFFFF'           # [color]  Plot background color
        thm.light_unit_label_color = '#808080'        # [color]  Label color
        thm.light_unit_label_fill = '#FFFFFF'         # [color]  Label fill color
        thm.light_slider_border_color = '#808080'     # [color]  Slider frame border color
        thm.light_slider_bg_color = '#AAC0C0C0'       # [color]  Slider frame background color
        thm.light_slider_bg_hover = '#C0C0C0'         # [color]  Slider frame hover color
        thm.light_slider_label_color = '#000000'      # [color]  Slider frame label color
        thm.light_overlay_threshold_color = '#FF0000' # [color]  Map threshold color
        # dark mode
        thm.dark_conic_label_fill = '#000000'         # [color]  Contour label fill color
        thm.dark_conic_ref_color = '#303030'          # [color]  Reference contour color
        thm.dark_beamstop_color = '#FF0000AA'         # [color]  Beamstop color
        thm.dark_beamstop_edge_color = '#FF0000'      # [color]  Beamstop edge color
        thm.dark_det_module_color = '#EEEEEE'         # [color]  Detector module border color
        thm.dark_det_module_fill = '#EEEEEE'          # [color]  Detector module background color
        thm.dark_plot_bg_color = '#000000'            # [color]  Plot background color
        thm.dark_unit_label_color = '#C0C0C0'         # [color]  Label color
        thm.dark_unit_label_fill = '#000000'          # [color]  Label fill color
        thm.dark_slider_border_color = '#202020'      # [color]  Slider frame border color
        thm.dark_slider_bg_color = '#AA303030'        # [color]  Slider frame background color
        thm.dark_slider_bg_hover = '#303030'          # [color]  Slider frame hover color
        thm.dark_slider_label_color = '#C0C0C0'       # [color]  Slider frame label color
        thm.dark_overlay_threshold_color = '#FF0000'  # [color]  Map threshold color

        return thm
    
    def get_defaults_lmt(self):
        ##########
        # Limits #
        ##########
        lmt = container()
        lmt.ener_min =  5    # [int] Energy minimum [keV]
        lmt.ener_max =  100  # [int] Energy maximum [keV]
        lmt.ener_stp =  1    # [int] Energy step size [keV]

        lmt.dist_min =  40   # [int] Distance minimum [mm]
        lmt.dist_max =  1000 # [int] Distance maximum [mm]
        lmt.dist_stp =  1    # [int] Distance step size [mm]

        lmt.hoff_min = -150  # [int] Horizontal offset minimum [mm]
        lmt.hoff_max =  150  # [int] Horizontal offset maximum [mm]
        lmt.hoff_stp =  1    # [int] Horizontal offset step size [mm]

        lmt.voff_min = -250  # [int] Vertical offset minimum [mm]
        lmt.voff_max =  250  # [int] Vertical offset maximum [mm]
        lmt.voff_stp =  1    # [int] Vertical offset step size [mm]

        lmt.rota_min = -45   # [int] Rotation minimum [deg]
        lmt.rota_max =  45   # [int] Rotation maximum [deg]
        lmt.rota_stp =  1    # [int] Rotation step size [deg]

        lmt.tilt_min = -40   # [int] Tilt minimum [deg]
        lmt.tilt_max =  40   # [int] Tilt maximum [deg]
        lmt.tilt_stp =  1    # [int] Tilt step size [deg]

        lmt.bsdx_min =   5   # [int] Beamstop distance minimum [mm]
        lmt.bsdx_max = 1000  # [int] Beamstop distance maximum [mm]
        lmt.bsdx_stp =   1   # [int] Beamstop distance step size [mm]
        
        return lmt

    def get_defaults_all(self):
        # load the defaults
        # geo: geometry and detector specs
        self.geo = self.get_defaults_geo()
        self._geo = container()
        # plo: plot details
        self.plo = self.get_defaults_plo()
        # thm: theme details
        self.thm = self.get_defaults_thm()
        # lmt: geometry limits
        self.lmt = self.get_defaults_lmt()

    def get_specs_det(self):
        det_type = self.geo.det_type
        det_size = self.geo.det_size

        if det_type not in self.detector_db.keys():
            print(f'Unknown detector type: {det_type}')
            print(f'Current databank entries: {", ".join(self.detector_db.keys())}.')
            raise SystemExit
        
        if det_size not in self.detector_db[det_type]['size'].keys():
            print(f'Unknown detector type/size combination: {det_type}/{det_size}')
            print(f'Current {det_type} databank sizes: {", ".join(self.detector_db[det_type]["size"].keys())}.')
            raise SystemExit
        
        det = container()
        det.hms = self.detector_db[det_type]['hms']
        det.vms = self.detector_db[det_type]['vms']
        det.pxs = self.detector_db[det_type]['pxs']
        det.hgp = self.detector_db[det_type]['hgp']
        det.vgp = self.detector_db[det_type]['vgp']
        det.cbh = self.detector_db[det_type]['cbh']
        det.hmn, det.vmn = self.detector_db[det_type]['size'][det_size]
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
            'size' : {'300K':(1,3),
                        '1M':(2,5),
                        '2M':(3,8),
                        '6M':(5,12)},
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
            'size' : {'1M':(2,4),
                      '2M':(3,6),
                      '4M':(4,8)}
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
            'size' : { '1M':(1,2),
                       '4M':(2,4),
                       '9M':(3,6),
                      '16M':(4,8)},
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
            'size' : {'MX225-HS':(3,3),
                      'MX300-HS':(4,4)},
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
            'size' : { '7':(1,1),
                      '14':(1,2),
                      '28':(2,2)},
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
            'size' : { '7':(1,1),
                      '14':(1,2)},
            }
        
        ###################################
        # Specifications for Perkin-Elmer #
        ###################################
        detectors['Perkin-Elmer XRD'] = {
            'hms' : 204.8,  # [mm]  Module size (horizontal)
            'vms' : 204.8,  # [mm]  Module size (vertical)
            'pxs' : 100e-3, # [mm]  Pixel size
            'hgp' : 0,      # [pix] Gap between modules (horizontal)
            'vgp' : 0,      # [pix] Gap between modules (vertical)
            'cbh' : 0,      # [mm]  Central beam hole
            'size' : {'0822':(1,1),
                      '1611':(2,2),
                      '1620':(2,2),
                      '1621':(2,2),
                      '1622':(2,2),
                      '1642':(2,2)},
            }
        
        ############################
        # Specifications for Varex #
        ############################
        detectors['VAREX XRpad2'] = {
            'hms' : 428.8,  # [mm]  Module size (horizontal)
            'vms' : 428.8,  # [mm]  Module size (vertical)
            'pxs' : 100e-3, # [mm]  Pixel size
            'hgp' : 0,      # [pix] Gap between modules (horizontal)
            'vgp' : 0,      # [pix] Gap between modules (vertical)
            'cbh' : 0,      # [mm]  Central beam hole
            'size' : {'4343':(1,1)},
            }
        
        # make file dump
        if not os.path.exists(self.path_detdb) or reset:
            with open(self.path_detdb, 'w') as wf:
                json.dump(detectors, wf, indent=4)
        else:
            try:
                with open(self.path_detdb, 'r') as of:
                    for key, vals in json.load(of).items():
                        detectors[key] = vals
            except: # any error is critical here!
                print(f"Error parsing Detector db at: {self.path_detdb}")
                raise SystemExit
        
        if update:
            with open(self.path_detdb, 'w') as wf:
                json.dump(detectors, wf, indent=4)
        
        return detectors

    def get_tooltips(self):
        self.tooltips = dict()
        self.tooltips['geo'] = {
            'det_type':'[str] e.g.\n{}'.format('\n'.join(self.get_det_library(update=False).keys())),
            'det_size':'[str] e.g.\n{}'.format('\n'.join([f"{k}: {', '.join(v['size'])}" for k,v in self.get_det_library(update=False).items()])),
            'ener':'[keV] Beam energy',
            'dist':'[mm] Detector distance',
            'voff':'[mm] Detector offset (vertical)',
            'hoff':'[mm] Detector offset (horizontal)',
            'rota':'[deg] Detector rotation',
            'tilt':'[deg] Detector tilt',
            'bssz':'[mm] Current beamstop size (or None)',
            'bsdx':'[mm] Beamstop distance',
            'unit':'[0-3] Contour legend\n0: 2-Theta\n1: d-spacing\n2: q-space\n3: sin(theta)/lambda',
            'reference':'[str] Plot reference contours\npick from pyFAI or None:\n{}'.format(', '.join(calibrant.names())),
            'darkmode':'[bool] Darkmode',
            'colormap':'[cmap] Contour colormap:\n{}'.format(', '.join(self.colormaps)),
            'bs_list':'[list] Available beamstop sizes',
            'det_bank':'Available detector bank\nkey:[value] or [value1, value2]\nkey: Detector name\nvalue: Detector model/size/type\ne.g. {PLATUS3:[1M, 2M], EIGER2:[4M]}\nempty dict enables all detectors',
        }
        self.tooltips['plo'] = {
            'conic_tth_min':'[int] Minimum 2-theta contour line',
            'conic_tth_max':'[int] Maximum 2-theta contour line',
            'conic_tth_num':'[int] Number of contour lines',
            'beamcenter_marker':'[marker] Beamcenter marker',
            'beamcenter_size':'[int] Beamcenter size',
            'poni_marker':'[marker] Poni marker',
            'poni_size':'[int] Poni size',
            'conic_linewidth':'[float] Contour linewidth (lw)',
            'conic_label_size':'[int] Contour labelsize',
            'conic_ref_linewidth':'[float] Reference contour linewidth',
            'conic_ref_num':'[int] Number of reference contours',
            'conic_ref_cif_int':'[float] Minimum display intensity (cif)',
            'conic_ref_cif_kev':'[float] Energy [keV] for intensity calculation',
            'conic_ref_cif_irel':'[bool] Linewidth relative to intensity',
            'conic_ref_cif_lw_min':'[float] Minimum linewidth when using irel',
            'conic_ref_cif_lw_mult':'[float] Linewidth multiplier when using irel',
            'conic_hkl_show_int':'[bool] Show intensity in hkl tooltip',
            'conic_hkl_label_size':'[int] Font size of hkl tooltip',
            'det_module_alpha':'[float] Detector module alpha',
            'det_module_width':'[int] Detector module border width',
            'conic_steps':'[int] Conic resolution',
            'plot_size':'[int] Plot size, px (0 for auto)',
            'plot_size_fixed':'[bool] Fix window size',
            'unit_label_size':'[int] Label size, px',
            'polarisation_fac':'[float] Horizontal polarisation factor',
            'show_polarisation':'[bool] Show polarisation overlay',
            'show_solidangle':'[bool] Show solid angle overlay',
            'show_unit_hover':'[bool] Show unit value on hover',
            'overlay_resolution':'[int] Overlay resolution',
            'overlay_threshold':'[float] Overlay warn color threshold',
            'overlay_toggle_warn':'[bool] Toggle overlay highlight',
            'slider_margin':'[int] Slider frame top margin',
            'slider_border_width':'[int] Slider frame border width',
            'slider_border_radius':'[int] Slider frame border radius (px)',
            'slider_label_size':'[int] Slider frame label size',
            'slider_column_width':'[int] Slider label column width',
            'enable_slider_ener':'[bool] Show energy slider',
            'enable_slider_dist':'[bool] Show distance slider',
            'enable_slider_rota':'[bool] Show rotation slider',
            'enable_slider_voff':'[bool] Show vertical offset slider',
            'enable_slider_hoff':'[bool] Show horizontal offset slider',
            'enable_slider_tilt':'[bool] Show tilt slider',
            'enable_slider_bsdx':'[bool] Show beamstop distance slider',
            'slider_label_ener':'[str] Label for energy slider',
            'slider_label_dist':'[str] Label for distance slider',
            'slider_label_rota':'[str] Label for rotation slider',
            'slider_label_voff':'[str] Label for vertical offset slider',
            'slider_label_hoff':'[str] Label for horizontal offset slider',
            'slider_label_tilt':'[str] Label for tilt slider',
            'slider_label_bsdx':'[str] Label for beamstop distance slider',
            'update_settings':'[bool] Update settings file after load',
            'update_det_bank':'[bool] Update detector bank after load',
            'reset_settings':'[bool] Reset settings file',
            'reset_det_bank':'[bool] Reset detector bank',
            'set_debug':'[bool] Debug mode: Allow pan, zoom and use of toolbox menu'
        }
        self.tooltips['thm'] = {
            'color_dark':'[color] Global dark color',
            'color_light':'[color] Global light color',
            'light_conic_label_fill':'[color] Contour label fill color',
            'light_conic_ref_color':'[color] Reference contour color',
            'light_beamstop_color':'[color] Beamstop color',
            'light_beamstop_edge_color':'[color] Beamstop edge color',
            'light_det_module_color':'[color] Detector module border color',
            'light_det_module_fill':'[color] Detector module background color',
            'light_plot_bg_color':'[color] Plot background color',
            'light_unit_label_color':'[color] Label color',
            'light_unit_label_fill':'[color] Label fill color',
            'light_slider_border_color':'[color] Slider frame border color',
            'light_slider_bg_color':'[color] Slider frame background color',
            'light_slider_bg_hover':'[color] Slider frame hover color',
            'light_slider_label_color':'[color] Slider frame label color',
            'light_overlay_threshold_color':'[color] Overlay threshold color',
            'dark_conic_label_fill':'[color] Contour label fill color',
            'dark_conic_ref_color':'[color] Reference contour color',
            'dark_beamstop_color':'[color] Beamstop color',
            'dark_beamstop_edge_color':'[color] Beamstop edge color',
            'dark_det_module_color':'[color] Detector module border color',
            'dark_det_module_fill':'[color] Detector module background color',
            'dark_plot_bg_color':'[color] Plot background color',
            'dark_unit_label_color':'[color] Label color',
            'dark_unit_label_fill':'[color] Label fill color',
            'dark_slider_border_color':'[color] Slider frame border color',
            'dark_slider_bg_color':'[color] Slider frame background color',
            'dark_slider_bg_hover':'[color] Slider frame hover color',
            'dark_slider_label_color':'[color] Slider frame label color',
            'dark_overlay_threshold_color':'[color] Overlay threshold color',
        }
        self.tooltips['lmt'] = {
            'ener_min':'[int] Energy minimum [keV]',
            'ener_max':'[int] Energy maximum [keV]',
            'ener_stp':'[int] Energy step size [keV]',
            'dist_min':'[int] Distance minimum [mm]',
            'dist_max':'[int] Distance maximum [mm]',
            'dist_stp':'[int] Distance step size [mm]',
            'hoff_min':'[int] Horizontal offset minimum [mm]',
            'hoff_max':'[int] Horizontal offset maximum [mm]',
            'hoff_stp':'[int] Horizontal offset step size [mm]',
            'voff_min':'[int] Vertical offset minimum [mm]',
            'voff_max':'[int] Vertical offset maximum [mm]',
            'voff_stp':'[int] Vertical offset step size [mm]',
            'rota_min':'[int] Rotation minimum [deg]',
            'rota_max':'[int] Rotation maximum [deg]',
            'rota_stp':'[int] Rotation step size [deg]',
            'tilt_min':'[int] Tilt minimum [deg]',
            'tilt_max':'[int] Tilt maximum [deg]',
            'tilt_stp':'[int] Tilt step size [deg]',
            'bsdx_min':'[int] Beamstop distance minimum [mm]',
            'bsdx_max':'[int] Beamstop distance maximum [mm]',
            'bsdx_stp':'[int] Beamstop distance step size [mm]',
        }

        # DEBUG function:
        # check if all keys have a tooltip
        if self.plo.set_debug:
            print('DEBUG: Check tooltips for completeness')
            _d = {'geo':self.geo.__dict__,
                'plo':self.plo.__dict__,
                'thm':self.thm.__dict__,
                'lmt':self.lmt.__dict__}
            for k,v in _d.items():
                for i in v.keys():
                    if i not in self.tooltips[k].keys():
                        print(i, '-> missing')

    #############
    #   BUILD   #
    #############
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
                rect_item.setPen(pg.mkPen(color = self.det_module_color, width = self.plo.det_module_width))
                rect_item.setBrush(pg.mkBrush(color = self.det_module_fill))
                rect_item.setOpacity(self.plo.det_module_alpha)
                self.ax.addItem(rect_item)

    #############
    #  CONICS   #
    #############
    def draw_conics(self):
        # calculate the offset of the contours resulting from voff and rotation
        # shift the grid to draw the cones, to make sure the contours are drawn
        # within the visible area

        # convert theta in degrees to radians
        # for some reason I defined it negative some time ago
        # now there's no turning back!
        _omega = -np.deg2rad(self.geo.tilt + self.geo.rota)
        
        # beamcenter shift
        _comp_shift = -(self.geo.voff - self.geo.dist * np.tan(_omega) - np.deg2rad(self.geo.tilt) * self.geo.dist)

        # overlay
        if self.plo.show_polarisation or self.plo.show_solidangle or self.plo.show_unit_hover:# or self.plo.show_deltadd:
            _grd, self._tth, self._polcor, self._solang = self.calc_overlays(_omega, res=self.plo.overlay_resolution, pol=self.plo.polarisation_fac)
            self.patches['overlay'].setImage(_grd * self._polcor * self._solang,
                                             autoLevels=False,
                                             levels=[0.0,1.0],
                                             rect=(-self.xdim,
                                                   -self.ydim,
                                                    self.xdim * 2,
                                                    self.ydim * 2))
            #if self.plo.show_deltadd:
            #    _poni_shift = -(self.geo.voff - np.deg2rad(self.geo.tilt) * self.geo.dist)
            #    self.patches['isocurve'].setRect(self.geo.hoff - self._deltadd_width/2,
            #                                     _poni_shift - self._deltadd_width/2,
            #                                     self._deltadd_width,
            #                                     self._deltadd_width)
        else:
            self.patches['overlay'].setImage(None)
            self.cor_label.hide()
        
        # update beam center
        self.patches['beamcenter'].setData([self.geo.hoff],[_comp_shift])
        # update beam center
        self.patches['poni'].setData([self.geo.hoff],[-(self.geo.voff - np.deg2rad(self.geo.tilt)*self.geo.dist)])

        if self.geo.bssz and not (isinstance(self.geo.bssz, str) and self.geo.bssz.lower() == 'none'):
            # make sure it is a float (might be a string from the export window!)
            self.geo.bssz = float(self.geo.bssz)
            # update beam stop
            bs_theta = np.tan((self.geo.bssz/2) / self.geo.bsdx)

            # calculate the conic section corresponding to the theta angle
            # :returns False is conic is outside of visiblee area
            x, y, label_pos = self.calc_conic(_omega, bs_theta, steps=self.plo.conic_steps)
            if x is False:
                self.patches['beamstop'].setVisible(False)
                self.patches['bs_label'].setVisible(False)
            else:
                # plot the conic section
                self.patches['beamstop'].setData(x, y, fillLevel=y.max())
                self.patches['beamstop'].setVisible(True)

                _unit = self.calc_unit(bs_theta)
                self.patches['bs_label'].setPos(self.geo.hoff, label_pos)
                self.patches['bs_label'].setText(f'{_unit:.2f}')
                self.patches['bs_label'].setVisible(True)
        else:
            self.patches['beamstop'].setVisible(False)
            self.patches['bs_label'].setVisible(False)

        for _n, _ttd in enumerate(self.cont_geom_num):
            self.patches['conic'][_n].setVisible(False)
            self.patches['labels'][_n].setVisible(False)
            # current fraction for colormap
            _f = _n/len(self.cont_geom_num)

            # convert theta in degrees to radians
            theta = np.deg2rad(_ttd)

            # calculate the conic section corresponding to the theta angle
            # :returns False is conic is outside of visiblee area
            x, y, label_pos = self.calc_conic(_omega, theta, steps=self.plo.conic_steps)
            if x is False or len(x) == 0:
                continue

            # plot the conic section
            self.patches['conic'][_n].setData(x, y, pen=pg.mkPen(self.cont_cmap.map(_f, mode='qcolor'), width=self.plo.conic_linewidth))
            self.patches['conic'][_n].setVisible(True)

            _unit = self.calc_unit(theta)
            self.patches['labels'][_n].setPos(self.geo.hoff, label_pos)
            self.patches['labels'][_n].setText(f'{_unit:.2f}', color=self.cont_cmap.map(_f, mode='qcolor'))
            self.patches['labels'][_n].setVisible(True)

    def draw_reference(self):
        # plot reference contour lines
        # standard contour lines are to be drawn
        for _n in range(self.plo.conic_ref_num):
            self.patches['reference'][_n].setVisible(False)
            # number of d-spacings might be lower than the maximum number of allowed contours
            if _n < len(self.cont_ref_dsp):
                _d = self.cont_ref_dsp[_n]
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
                _omega = -np.deg2rad(self.geo.tilt + self.geo.rota)
                
                # calculate the conic section corresponding to the theta angle
                # :returns False is conic is outside of visiblee area
                x, y, _ = self.calc_conic(_omega, theta, steps=self.plo.conic_steps)
                if x is False:
                    continue

                # if hkl are available
                # put them in the proper container for the contour
                # so indexing gets it right
                irel = 1.0
                if self.cont_ref_hkl:
                    h, k, l, itot, irel = self.cont_ref_hkl[_n]
                    irel *= self.plo.conic_ref_cif_lw_mult
                    if self.plo.conic_hkl_show_int:
                        # alignment of the tooltip text is far from trivial
                        # detour via QTextEdit -> setAlignment and toHtml
                        # hence the <br> instead of \n
                        # 
                        # self.setStyleSheet('''QToolTip {... ) didn't work
                        self.patches['reference'][_n].name = f'({h: 2.0f} {k: 2.0f} {l: 2.0f})<br>{round(itot, 0):,.0f}'
                    else:
                        self.patches['reference'][_n].name = f'({h: 2.0f} {k: 2.0f} {l: 2.0f})'
                    if not self.plo.conic_ref_cif_irel:
                        irel = 1.0
                else:
                    self.patches['reference'][_n].name = None
                
                # plot the conic section
                self.patches['reference'][_n].setData(x, y, pen=pg.mkPen(self.conic_ref_color,
                                                                         width=max(self.plo.conic_ref_cif_lw_min, 
                                                                                   self.plo.conic_ref_linewidth * irel)))
                self.patches['reference'][_n].setVisible(True)
    
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
        # y ('height') components/distances from central axis of the cone
        # intersecting the detector plane and the distance to
        # the y intersection of the conic section.
        y1 = dz_cone * np.sin(theta) / np.cos(omega + theta)
        y2 = dz_cone * np.sin(theta) / np.cos(omega - theta)

        # add x/y offsets
        # revert tilt rotation
        y0 = dy_cone - self.geo.voff + comp_tilt 
        x0 = self.geo.hoff

        # add margin to slightly extend
        # conics outside of visible area
        _xdim = self.xdim * 1.05
        _ydim = self.ydim * 1.05
        # evaluate the eccentricity and parameterise
        # the resulting conic accordingly
        if abs(ecc) == 0:
            # circle
            h = (y1+y2)/2
            # check if the circle is visible
            if h - np.sqrt(y0**2 + x0**2) > np.sqrt(self.ydim**2 + self.xdim**2):
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
            # check if the ellipse is visible
            if _xlim1 < 0 and _xlim2 < 0:
                return False, False, False
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
            if h == 0:
                return False, False, False
            w = h * np.sqrt(ecc**2-1)
            l = -np.arcsinh((_xdim + x0) / w)
            r =  np.arcsinh((_xdim - x0) / w)
            t = np.linspace(l, r, steps)
            x = x0 + w * np.sinh(t)
            y = y0 + (y1-y2)/2 - h * np.cosh(t)
        elif abs(ecc) >= 100:
            # line
            t = np.linspace(-_xdim, _xdim, 2)
            x = t
            y = y0 + np.ones(len(t)) * y1

        # check if conic is visible
        cx = np.argwhere((x >= -_xdim) & (x <= _xdim))
        cy = np.argwhere((y >= -_ydim) & (y <= _ydim))
        if len(cx) == 0 or len(cy) == 0:
            # outside detector area
            return False, False, False
        
        # adjust the label position to maintain readibility
        # this works for most cases but is not the most optimal solution yet
        # OR: use the actual beam position to determine label position
        # beam_pos_y = -(self.geo.voff + np.tan(np.deg2rad(self.geo.rota))*self.geo.dist)
        if omega <= 0:
            label_pos = max(y) if theta < np.pi/2 else min(y)
        else:
            label_pos = min(y) if theta < np.pi/2 else max(y)
        
        # return x, y and the label position
        return x, y, label_pos

    def calc_ref_from_cif(self, fpath):
        # called when a cif is dropped onto the window
        xtl = dif.Crystal(fpath)
        # :return xval: arrray : x-axis of powder scan (units)
        # :return inten: array : intensity values at each point in x-axis
        # :return reflections: (h, k, l, xval, intensity) array of reflection positions, grouped by min_overlap
        xval, inten, reflections = xtl.Scatter.powder(scattering_type='xray', units='dspace', powder_average=True, min_overlap=0.02, energy_kev=self.plo.conic_ref_cif_kev)
        # reject low intensities: based on median or mean?
        # median is always around unity -> useless
        # mean rejects many, add adjustable multiplicator?
        used = reflections[reflections[:,4] > reflections[:,4].max() * self.plo.conic_ref_cif_int]
        # sort by intensity -> ascending -> flip
        ordered = used[used[:, 4].argsort()][::-1]
        # pick the strongest
        ordered = ordered[:self.plo.conic_ref_num]
        # assign dspacing
        self.cont_ref_dsp = ordered[:,3]
        # hkl -> integer
        # cast hkl array to list of tuples (for easy display)
        irel = ordered[:,4]/ordered[:,4].max()
        self.cont_ref_hkl = list(zip(ordered[:,0], ordered[:,1], ordered[:,2], ordered[:,4], irel))

        self.geo.reference = os.path.basename(fpath)
        self.ref_custom[self.geo.reference] = self.cont_ref_dsp
        self.ref_custom_hkl[self.geo.reference] = self.cont_ref_hkl

        ref_action = QtGui.QAction(self.geo.reference, self, checkable=True)
        self.set_menu_action(ref_action, self.change_reference, self.geo.reference)
        self.sub_menu_custom.addAction(ref_action)
        self.group_ref.addAction(ref_action)
        ref_action.setChecked(True)
        
        # update window title
        self.set_window_title()

        self.draw_reference()
    
    def show_tooltip(self, widget, event):
        if not widget.name or not self.cont_ref_hkl:
            event.ignore()
            return False
        text = QtWidgets.QTextEdit(str(widget.name))
        text.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        pos = QtCore.QPoint(*map(int, event.screenPos()))# - QtCore.QPoint(10,20)
        QtWidgets.QToolTip.showText(pos, text.toHtml())
        event.ignore()

    #############
    #  OVERLAY  #
    #############
    def calc_overlays(self, omega, res=150, pol=0.99):
        # scale overlay to detector dimensions
        _res_scale = self.ydim/self.xdim
        _res_v = int(round(_res_scale * res, 0))
        # neutral overlay grid
        grd = np.ones((_res_v, res))

        # make screen grid
        size_h = np.linspace(-self.xdim, self.xdim, res, endpoint=False)
        size_v = np.linspace(-self.ydim, self.ydim, _res_v, endpoint=False)
        _gx, _gy = np.meshgrid(size_h, size_v, sparse=True)
        # build vector -> 3 x n x m
        _vec = np.full((3, _res_v, res), self.geo.dist)
        _vec[0,:,:] = _gx - self.geo.hoff
        # Compensate for vertical offset and PONI offset caused by the tilt (sdd*tilt)
        _vec[1,:,:] = _gy + self.geo.voff - np.deg2rad(self.geo.tilt) * self.geo.dist

        # omega is the combination of rotation and tilt, in radians
        _rot = self.rot_100(omega)
        # reshape to allow matrix multiplication
        # _rot: 3 x 3 @ _norm: 3 x n*m -> _res: 3 x n x m
        _res = np.reshape(_rot @ np.reshape(_vec, (3,-1)), _vec.shape)

        # unit hover
        tth = 1.0
        if self.plo.show_unit_hover:
            # Distance POBI - pixel on grid
            R_a = np.sqrt(np.sum(_res[0:2]**2, axis=0)) * 1e-3 # m
            # POBI distance
            D_a = _res[2] * 1e-3 # m
            # 2theta - Angle between pixel, sample, and POBI
            tth = np.arctan(R_a/D_a)
            # remove very small values (tth < 0.057 deg) to avoid zero divide
            tth[tth < 1e-3] = np.nan

        # polarisation
        pc = 1.0
        if self.plo.show_polarisation:
            _mag = np.sqrt(np.sum(_res**2, axis=0))
            _norm = _res / _mag
            # add pol fractions (-> 1.0)
            # this notation is equivalent to cos(psi)**2,
            # psi angle of polarization direction to the observer
            # _res is direct cos(psi)
            pc_hor = _norm[0,:,:]**2 * pol
            pc_ver = _norm[1,:,:]**2 * (1-pol)
            pc = 1.0 - (pc_hor + pc_ver)

        # solid angle
        sa = 1.0
        if self.plo.show_solidangle:
            _mag = np.sqrt(np.sum(_vec**2, axis=0))
            sa = 1 / _mag**3
            sa = sa / np.max(sa)
        
        """# delta d / d
        dd_grid = 1.0
        dd_width = 0.0
        if self.plo.show_deltadd:
            
            c = self.plo.capillary_size
            t = self.plo.scattering_diameter
            p = self.det.pxs * 1e-3
            phi = self.plo.beam_divergence
            dE_E = self.plo.energy_resolution
            
            # use the "unrotated" vector coordinates to define
            # R and D in the detector plane
            # radial component of the unrotated detector
            # i.e. radius away from the PONI
            R_a = np.sqrt(np.sum(_vec[0:2]**2, axis=0)) * 1e-3 # m
            # PONI distance
            D_a = self.geo.dist * 1e-3 # m
            # 2theta-alpha - Angle between pixel, sample, and PONI
            _tth_a = np.arctan(R_a/D_a)
            # remove very small values (tth < 0.057 deg) to avoid zero divide
            _tth_a[_tth_a < 1e-3] = np.nan
            # estimate of the SDD spread (independent of rotation)
            dD = np.sqrt(1/4 * (c**2 + t**2))
            # estimate of the pixel radial spread in the detector plane
            dR = np.sqrt(1/4*(c**2/(np.cos(_tth_a)**2)          \
                            + p**2 + t**2*np.tan(_tth_a)**2     \
                            + phi**2*D_a**2/(np.cos(_tth_a)**4) \
                             )                                  \
                        )
            
            # delta d over d
            _dd = np.sqrt(1/4*np.tan(_tth_a/2)**(-2)                \
                          * np.sin(_tth_a)**2 * np.cos(_tth_a)**2   \
                          * ((dD/D_a)**2 + (dR/R_a)**2) + (dE_E)**2 \
                          )
            # H2, FWHM
            #_dd = 32*np.log(2)*( np.cos(_tth_a)**4/(16*D_a**2) \
            #             * ((np.tan(_tth_a)**2 \
            #             * (c**2 + 2*t**2) \
            #             + p**2 \
            #             + c**2/np.cos(_tth_a)**2 \
            #             + (D_a**2*phi**2)/np.cos(_tth_a)**4) ))
            dd_grid = np.log10(_dd * 180**2 / np.pi**2)
            _ra = R_a[dd_grid < self.plo.funct_deltadd_thresh]
            if any(_ra):
                dd_width = np.nanmin(_ra) * 2 * 1e3
        """

        return grd, tth, pc, sa#, dd_grid, dd_width

    def rot_100(self, a, cc=1):
        #Omega in radians
        ca = np.cos(a)
        sa = np.sin(a)
        if cc: sa = -sa
        return np.array([[1,   0,  0],
                         [0,  ca, sa],
                         [0, -sa, ca]])

    def window_function_deltadd(self):
        param_window = QtWidgets.QDialog()

        param_dict = {'Detector':[(0, 'Sensor thickness', self.plo.sensor_thickness, 1e-6, '\u00B5m')],
                      'Beam':[(1, 'Divergence', self.plo.beam_divergence, 1e-6, '\u00B5rad'),
                              (2, 'Energy resolution', self.plo.energy_resolution, 1e-6, '\u00B5eV/keV')],
                      'Sample':[(3, 'Scattering volume \u2300', self.plo.scattering_diameter, 1e-6, '\u00B5m')],
                      'Display':[(4, 'Threshold', self.plo.funct_deltadd_thresh, 1, 'log10(\u03B4d/d)')]}
        param_dict_change = {}

        layout = QtWidgets.QVBoxLayout()
        for title, entry in param_dict.items():
            box = QtWidgets.QGroupBox(title=title)
            box_layout = QtWidgets.QVBoxLayout()
            box_layout.setContentsMargins(0,0,0,0)
            for idx, label, value, div, unit in entry:
                entry_box = QtWidgets.QFrame()
                entry_layout = QtWidgets.QHBoxLayout()
                box_combobox = QtWidgets.QSpinBox(singleStep=1, minimum=-10, maximum=1000, value=int(value/div))
                entry_layout.addWidget(QtWidgets.QLabel(label))
                entry_layout.addWidget(box_combobox)
                entry_layout.addWidget(QtWidgets.QLabel(unit))
                entry_box.setLayout(entry_layout)
                box_layout.addWidget(entry_box)
                param_dict_change[idx] = [box_combobox, div]
            box.setLayout(box_layout)
            layout.addWidget(box)

        button_box = QtWidgets.QDialogButtonBox()
        button_box.addButton(QtWidgets.QDialogButtonBox.StandardButton.Apply)
        button_box.setCenterButtons(True)
        button_box.clicked.connect(lambda: self.window_function_deltadd_accept(param_dict_change, param_window))
        layout.addWidget(button_box)

        windowIcon = QtGui.QIcon()
        windowIcon.addPixmap(QtGui.QPixmap(':/icons/xrdPlanner.png'))
        
        param_window.setWindowIcon(windowIcon)
        param_window.setLayout(layout)
        param_window.setFixedSize(param_window.sizeHint())
        param_window.exec()
    
    def window_function_deltadd_accept(self, dict, window):
        self.plo.sensor_thickness = round(dict[0][0].value() * dict[0][1], 6)
        self.plo.beam_divergence = round(dict[1][0].value() * dict[1][1], 6)
        self.plo.energy_resolution = round(dict[2][0].value() * dict[2][1], 6)
        self.plo.capillary_size = round(dict[3][0].value() * dict[3][1], 6)
        self.plo.funct_deltadd_thresh = round(dict[4][0].value() * dict[4][1], 6)
        self.action_funct_deltadd_show.setEnabled(True)
        self.plo.show_deltadd = False
        self.toggle_function_deltadd()
        window.close()
    
    #############
    #  UPDATE   #
    #############
    def update_screen(self, val=None):
        if val is not None:
            if self.sender().objectName() == 'dist':
                self.geo.dist = float(val)
                if self.plo.enable_slider_bsdx:
                    # set beamstop max distance to current detector distance
                    # or the maximum allowed value, whichever is smaller
                    current_max = min(self.geo.dist, self.lmt.bsdx_max)
                    self.sliderWidget.update_slider_limits(self.sliderWidget.sl_bsdx,
                                                           self.lmt.bsdx_min,
                                                           current_max)
            elif self.sender().objectName() == 'rota':
                self.geo.rota = float(val)
            elif self.sender().objectName() == 'tilt':
                self.geo.tilt = float(val)
            elif self.sender().objectName() == 'voff':
                self.geo.voff = float(val)
            elif self.sender().objectName() == 'hoff':
                self.geo.hoff = float(val)
            elif self.sender().objectName() == 'ener':
                self.geo.ener = float(val)
            elif self.sender().objectName() == 'bsdx':
                self.geo.bsdx = float(val)

        # re-calculate cones and re-draw contours
        self.draw_conics()
        # draw reference contours
        if self.geo.reference != 'None':
            self.get_reference()
            self.draw_reference()

    #############
    #  EXPORT   #
    #############
    def show_export_window(self):
        # set flags=QtCore.Qt.WindowType.Tool for the window
        # to not loose focus when the FileDialog is closed
        self.export_window = QtWidgets.QMainWindow(parent=self, flags=QtCore.Qt.WindowType.Tool)
        self.export_window.setWindowTitle('Export current settings to file')
        layout_vbox = QtWidgets.QVBoxLayout()
        layout_hbox = QtWidgets.QHBoxLayout()
        layout_hbox.setContentsMargins(0,0,0,0)
        frame = QtWidgets.QFrame()
        frame.setLayout(layout_hbox)

        central_widget = QtWidgets.QWidget()
        central_widget.setLayout(layout_vbox)
        self.export_window.setCentralWidget(central_widget)

        # Fonts
        font_header = QtGui.QFont()
        font_header.setPixelSize(self.plo.slider_label_size)
        font_header.setBold(True)
        font_normal = QtGui.QFont()
        font_normal.setPixelSize(self.plo.slider_label_size)
        font_normal.setBold(False)

        # Detector tree box
        qbox_det = QtWidgets.QGroupBox()
        qbox_det.setTitle('Detector bank')
        qbox_det.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
        qbox_det.setFont(font_header)
        layout_qbox_det = QtWidgets.QVBoxLayout()
        qbox_det.setLayout(layout_qbox_det)

        # Detector tree
        self.tree_det = QtWidgets.QTreeWidget()
        self.tree_det.setToolTip('Add detectors and models to the detector bank.\nSelect none to make all detectors and models available.')
        self.tree_det.setColumnCount(1)
        self.tree_det.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.MultiSelection)
        self.tree_det.setHeaderLabels(['Detector type / size'])
        self.tree_det.header().setFont(font_header)
        self.tree_det.setAlternatingRowColors(True)
        _cur_det_type = None
        _cur_det_size = None
        for det, specs in self.get_det_library(update=False, reset=False).items():
            item = QtWidgets.QTreeWidgetItem([str(det)])
            item.setFont(0, font_header)
            if self.geo.det_type == str(det):
                _cur_det_type = item
            for val in specs['size'].keys():
                child = QtWidgets.QTreeWidgetItem([str(val)])
                child.setFont(0, font_normal)
                if self.geo.det_size == str(val) and self.geo.det_type == str(det):
                    _cur_det_size = child
                item.addChild(child)
            self.tree_det.addTopLevelItem(item)
        self.tree_det.expandAll()
        self.tree_det.itemClicked.connect(self.tree_ghost_select)

        # highlight current detector type/size
        if _cur_det_type is not None:
            _cur_det_type.setSelected(True)
        if _cur_det_size is not None:
            _cur_det_size.setSelected(True)

        # Beamstop list tree box
        qbox_bsb = QtWidgets.QGroupBox()
        qbox_bsb.setTitle('Beamstop bank')
        qbox_bsb.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
        qbox_bsb.setFont(font_header)
        layout_qbox_bsb = QtWidgets.QVBoxLayout()
        layout_qbox_bsb.setSpacing(0)
        qbox_bsb.setLayout(layout_qbox_bsb)

        # doubleSpinBox editor for the beamstop bank list widget
        class MyListDelegate(QtWidgets.QStyledItemDelegate):
            def displayText(self, value, locale):
                return f'{value:.1f}'
            
            def createEditor(self, parent, option, index):
                box = QtWidgets.QDoubleSpinBox(parent)
                box.setDecimals(1)
                box.setSingleStep(0.1)
                box.setMinimum(0)
                box.setMaximum(100)
                return box

        # Beamstop bank list widget
        self.tree_bsb = QtWidgets.QListWidget()
        self.tree_bsb.setAlternatingRowColors(True)
        self.tree_bsb.setItemDelegate(MyListDelegate())
        self.tree_bsb.setToolTip('Specify the available beamstop sizes.')
        self.tree_bsb.itemChanged.connect(self.tree_bsb.sortItems)
        for bs_size in self.geo.bs_list:
            item = QtWidgets.QListWidgetItem()
            item.setData(2, bs_size)
            item.setFont(font_normal)
            item.setFlags(item.flags()|QtCore.Qt.ItemFlag.ItemIsEditable)
            self.tree_bsb.addItem(item)

        # Add row to beamstop bank list widget
        def row_add(aQListWidget, font):
            vmax = 0
            for i in range(aQListWidget.count()):
                if aQListWidget.item(i).data(0) > vmax:
                    vmax = aQListWidget.item(i).data(0)
            item = QtWidgets.QListWidgetItem()
            item.setData(2, vmax + 1)
            item.setFont(font)
            item.setFlags(item.flags()|QtCore.Qt.ItemFlag.ItemIsEditable)
            aQListWidget.addItem(item)

        # Remove row to beamstop bank list widget
        def row_rem(aQListWidget):
            _ = aQListWidget.takeItem(aQListWidget.currentRow())

        # add/remove buttons for the beamstop bank list widget
        button_add = QtWidgets.QToolButton()
        button_add.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
        button_add.setText('+')
        button_add.setFont(font_header)
        button_add.setToolTip('Add a new entry to the list.')
        button_add.clicked.connect(lambda: row_add(self.tree_bsb, font_normal))
        # remove button
        button_rem = QtWidgets.QToolButton()
        button_rem.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
        button_rem.setText('-')
        button_rem.setFont(font_header)
        button_rem.setToolTip('Remove the highlighted entry from the list.')
        button_rem.clicked.connect(lambda: row_rem(self.tree_bsb))
        # button in layout
        layout_hbox_but = QtWidgets.QHBoxLayout()
        layout_hbox_but.addWidget(button_add)
        layout_hbox_but.addWidget(button_rem)
        qbox_buttons = QtWidgets.QGroupBox()
        qbox_buttons.setLayout(layout_hbox_but)

        # Parameter tree box
        qbox_par = QtWidgets.QGroupBox()
        qbox_par.setTitle('Review parameters')
        qbox_par.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
        qbox_par.setFont(font_header)
        layout_qbox_par = QtWidgets.QVBoxLayout()
        qbox_par.setLayout(layout_qbox_par)

        # Make only one column editable
        # Specify the editors for different data types
        class MyTreeDelegate(QtWidgets.QItemDelegate):
            def createEditor(self, parent, option, index):
                if index.column() == 1:
                    # bool needs to be evaluated before int as True and False
                    # will be interpreted as integers
                    if isinstance(index.data(0), bool):
                        return super(MyTreeDelegate, self).createEditor(parent, option, index)
                    elif isinstance(index.data(0), int):
                        box = QtWidgets.QSpinBox(parent)
                        box.setSingleStep(1)
                        box.setMinimum(int(-1e9))
                        box.setMaximum(int(1e9))
                        return box
                    elif isinstance(index.data(0), float):
                        box = QtWidgets.QDoubleSpinBox(parent)
                        box.setDecimals(2)
                        box.setSingleStep(0.01)
                        box.setMinimum(-1e9)
                        box.setMaximum(1e9)
                        return box
                    else:
                        return super(MyTreeDelegate, self).createEditor(parent, option, index)
                return None
            
        # Parameter tree widget
        self.tree_par = QtWidgets.QTreeWidget()
        self.tree_par.setColumnCount(2)
        self.tree_par.setHeaderLabels(['Parameter', 'Value'])
        self.tree_par.header().setFont(font_header)
        self.tree_par.setAlternatingRowColors(True)
        self.tree_par.setItemDelegate(MyTreeDelegate())
        # det_bank and bs_list need some special treatment
        # to facilitate their editing
        dont_show = ['det_bank', 'bs_list']
        for section, values in {'geo':self.geo.__dict__, 'plo':self.plo.__dict__, 'thm':self.thm.__dict__, 'lmt':self.lmt.__dict__}.items():
            item = QtWidgets.QTreeWidgetItem([str(section)])
            item.setFont(0, font_header)
            #item.setFlags(QtCore.Qt.ItemFlag.ItemIsSelectable|QtCore.Qt.ItemFlag.ItemIsEnabled)
            for par, val in values.items():
                if par in dont_show:
                    continue
                child = QtWidgets.QTreeWidgetItem()
                child.setData(0, 0, str(par))
                if section in self.tooltips and par in self.tooltips[section]:
                    child.setToolTip(0, self.tooltips[section][par])
                    child.setToolTip(1, self.tooltips[section][par])
                child.setData(1, 2, val)
                child.setFont(0, font_normal)
                child.setFont(1, font_normal)
                child.setFlags(item.flags()|QtCore.Qt.ItemFlag.ItemIsEditable)
                item.addChild(child)
            self.tree_par.addTopLevelItem(item)
        self.tree_par.expandAll()
        self.tree_par.resizeColumnToContents(0)
        self.tree_par.resizeColumnToContents(1)

        button_export = QtWidgets.QToolButton()
        button_export.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
        button_export.setText('Export settings to file')
        button_export.setFont(font_header)
        button_export.clicked.connect(self.export_to)

        qbox_btn = QtWidgets.QGroupBox()
        layout_qbox_btn = QtWidgets.QVBoxLayout()
        qbox_btn.setLayout(layout_qbox_btn)
        
        layout_vbox.addWidget(frame)
        layout_hbox.addWidget(qbox_det, stretch=2)
        layout_hbox.addWidget(qbox_bsb, stretch=1)
        layout_hbox.addWidget(qbox_par, stretch=3)
        layout_vbox.addWidget(qbox_btn)
        layout_qbox_det.addWidget(self.tree_det)
        layout_qbox_bsb.addWidget(self.tree_bsb)
        layout_qbox_bsb.addWidget(qbox_buttons)
        layout_qbox_par.addWidget(self.tree_par)
        layout_qbox_btn.addWidget(button_export)

        self.export_window.show()

    def tree_ghost_select(self, item):
        parent = item.parent()
        # toplevel -> multiselect
        if parent == None:
            if item.isSelected():
                [item.child(i).setSelected(True) for i in range(item.childCount())]
            else:
                [item.child(i).setSelected(False) for i in range(item.childCount())]
        # sublevel -> single + top select
        else:
            if any([parent.child(i).isSelected() for i in range(parent.childCount())]):
                parent.setSelected(True)
            else:
                parent.setSelected(False)

    def export_to(self):
        # make detector bank/dict from selection
        det_bank = {}
        for item in self.tree_det.selectedItems():
            if item.parent() == None:
                det_name = item.text(0)
                det_types = [item.child(i).text(0) for i in range(item.childCount()) if item.child(i).isSelected()]
                det_bank.update({det_name:det_types})

        # get beamstop list from listwidget
        bs_list = [round(self.tree_bsb.item(i).data(0), 5) for i in range(self.tree_bsb.count())]

        # make settings dict from treewidget
        settings = dict()
        for i in range(self.tree_par.topLevelItemCount()):
            top = self.tree_par.topLevelItem(i)
            key = top.data(0, 0)
            settings[key] = dict()
            for j in range(top.childCount()):
                k = top.child(j).data(0, 0)
                v = top.child(j).data(1, 0)
                if isinstance(v, float):
                    v = round(v, 5)
                settings[key][k] = v

        # save parameters to settings file
        target, filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Export settings file', self.path_settings_current, "Settings files (*.json)")
        if not target:
            return
        
        # add det_bank and bs_list to settings
        settings['geo']['det_bank'] = det_bank
        settings['geo']['bs_list'] = bs_list
        with open(target, 'w') as wf:
            json.dump(settings, wf, indent=4)
        
        # close the window
        if self.export_window:
            self.export_window.close()

        # activate exported settings file
        basename = os.path.basename(target)
        self.change_settings_file(basename)

        # add new settings file to menu and check it
        # add if not in list
        if basename not in [action.text() for action in self.menu_custom_settings.actions()]:
            cset_action = QtGui.QAction(basename, self, checkable=True)
            self.set_menu_action(cset_action, self.change_settings_file, basename)
            self.menu_custom_settings.addAction(cset_action)
            self.group_cset.addAction(cset_action)
            cset_action.setChecked(True)
        # check if in list
        else:
            for action in self.menu_custom_settings.actions():
                if action.text() == basename:
                    action.setChecked(True)

    #############
    #  UTILITY  #
    #############
    def show_about_window(self):
        msgBox = QtWidgets.QDialog()
        #msgBox.setWindowTitle('About')
        
        windowIcon = QtGui.QIcon()
        windowIcon.addPixmap(QtGui.QPixmap(':/icons/xrdPlanner.png'))
        
        font_title = QtGui.QFont()
        font_title.setPointSize(48)
        icon = QtWidgets.QLabel()
        icon.setPixmap(QtGui.QPixmap(':/icons/xrdPlanner.png'))
        title = QtWidgets.QLabel(f'<b>xrdPlanner</b>')
        title.setFont(font_title)
        suptitle = QtWidgets.QLabel(f'<b>Version {xrdPlanner.__version__}</b> (released {xrdPlanner.__date__})')
        github = QtWidgets.QLabel(f'<br>A tool to project X-ray diffraction cones on a detector screen at different \
                                    geometries (tilt, rotation, offset) and X-ray energies. For more information visit \
                                    us on <a href="https://github.com/LennardKrause/xrdPlanner">Github</a>.')
        github.setWordWrap(True)
        github.setOpenExternalLinks(True)
        published = QtWidgets.QLabel(f'<br>The article is published in<br><a href="https://doi.org/10.1107/S1600577523011086">\
                                       <i>J. Synchrotron Rad.</i> (2024). <b>31</b></a>')
        published.setOpenExternalLinks(True)
        authors = QtWidgets.QLabel(f'<br><b>Authors:</b><br>{"<br>".join(xrdPlanner.__authors__)}')
        email = QtWidgets.QLabel(f'<br>Feedback? <a href="mailto:{xrdPlanner.__email__}?subject=xrdPlanner feedback">{xrdPlanner.__email__}</a>')
        email.setOpenExternalLinks(True)

        box_layout = QtWidgets.QVBoxLayout()
        box_layout.setSpacing(0)
        for widget in [title, suptitle, github, published, authors, email]:
            box_layout.addWidget(widget)
        
        box = QtWidgets.QGroupBox()
        box.setFlat(True)
        box.setLayout(box_layout)

        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(icon)
        layout.addWidget(box)

        msgBox.setWindowIcon(windowIcon)
        msgBox.setLayout(layout)
        msgBox.setFixedSize(msgBox.sizeHint())
        msgBox.exec()

    def resize_window(self):
        # figure out proper plot dimensions
        self.xdim = (self.det.hms * self.det.hmn + self.det.pxs * self.det.hgp * self.det.hmn + self.det.cbh)/2
        self.ydim = (self.det.vms * self.det.vmn + self.det.pxs * self.det.vgp * self.det.vmn + self.det.cbh)/2
        
        # limit the axis x and y
        self.ax.setXRange(-self.xdim, self.xdim, padding=0, update=True)
        self.ax.setYRange(-self.ydim, self.ydim, padding=0, update=True)

        if self.plo.plot_size <= 0:
            _app = QtWidgets.QApplication.instance()
            _height = _app.primaryScreen().availableGeometry().height()
            self.plo.plot_size = int(np.ceil(_height*0.9))
        
        # get proper dimensions
        self.width = int(np.ceil(self.plo.plot_size * self.xdim / self.ydim))
        self.height = self.plo.plot_size + self.plo.slider_margin//2 + self.offset_win32

        # fix the window size
        if self.plo.plot_size_fixed:
            self.setMaximumHeight(self.height)
            self.setMinimumHeight(self.height)
            self.setMaximumWidth(self.width)
            self.setMinimumWidth(self.width)

        # resize the window
        self.resize(self.width, self.height)

    def set_window_title(self):
        if self.geo.reference.lower() == 'none':
            self.setWindowTitle(f'{self.det.name} - {self.active_settings}')
        else:
            self.setWindowTitle(f'{self.det.name} - {self.geo.reference} - {self.active_settings}')

    def set_menu_action(self, action, target, *args):
        action.triggered.connect(lambda: target(*args))

    def calc_unit(self, tth, return_dict=False):
        # calc_unit expects 2-Theta in radians

        # Conversion factor keV to Angstrom: 12.398
        # sin(t)/l: np.sin(Theta) / lambda -> (12.398/geo_energy)
        stl = np.sin(tth/2)/(12.398/self.geo.ener)
        # d-spacing: l = 2 d sin(t) -> 1/2(sin(t)/l)
        dsp = 1/(2*stl)
        units = {0:np.rad2deg(tth), 1:dsp, 2:stl*4*np.pi, 3:stl}
        if return_dict:
            return units
        else:
            return units[self.geo.unit]

    #############
    # SETTINGS  #
    #############
    def get_active_settings_file(self):
        # settings token file exists
        # get settings file
        if os.path.exists(self.path_settings_token):
            with open(self.path_settings_token, 'r') as rf:
                name = rf.read()
            if os.path.exists(os.path.join(self.path_settings, name)):
                return name
        # otherwise reset and return default
        # -> last exit ended with crash
        # get defaults
        self.get_defaults_all()
        # set active settings to default
        self.active_settings = os.path.basename(self.path_settings_default)
        # reset/save default file
        self.save_settings()
        return self.active_settings
    
    def delete_active_settings_file(self):
        if os.path.exists(self.path_settings_token):
            os.remove(self.path_settings_token)

    def set_active_settings_file(self):
        if not os.path.exists(self.path_settings):
            os.makedirs(self.path_settings)
        with open(self.path_settings_token, 'w') as wf:
            wf.write(self.active_settings)

    def delete_settings_files(self):
        # todo change button from 'open' to 'delete'
        fnames, filter = QtWidgets.QFileDialog.getOpenFileNames(self, 'Delete settings files', self.path_settings_current, "Settings files (*.json)")
        if fnames:
            for fname in fnames:
                os.remove(fname)
                for action in self.menu_custom_settings.actions():
                    if action.text() == os.path.basename(fname):
                        self.menu_custom_settings.removeAction(action)

    def save_current_settings(self):
        # self.geo is edited with the sliders
        # self._geo holds the initial values
        # usually I do not want to overwrite 
        # the startup values -> I write _geo
        # to the settings file.
        # unless this function is called!
        self._geo.__dict__.update(self.geo.__dict__)
        self.save_settings()

    def save_settings(self, target=None):
        if target == None:
            target = os.path.join(self.path_settings, self.active_settings)
        target_base = os.path.dirname(target)
        # create folder if not existing
        if not os.path.exists(target_base):
            os.makedirs(target_base)
        # Writing geo as dict to file
        with open(target, 'w') as wf:
            json.dump({'geo':self._geo.__dict__, 'plo':self.plo.__dict__, 'thm':self.thm.__dict__, 'lmt':self.lmt.__dict__}, wf, indent=4)

    def load_settings(self, skip=[]):
        # Some parameters need to be protected to save
        # the user some waiting time
        #
        # A check is performed if the value of the key
        # is within 'minimum' and 'maximum', if not it
        # is set to 'default'
        #
        # Add 'key':('minimum', 'maximum', 'default') to the dict
        _warn = {'conic_ref_cif_kev':( 5,   25,  12),
                     #'conic_tth_min':( 1,   10,   5),
                     #'conic_tth_max':(10,  180,  90),
                     #'conic_tth_num':( 1,  100,  20),
                     #'conic_ref_num':( 1,  500, 200),
                       'conic_steps':(10, 1000, 100),
                          'ener_stp':( 1,  100,   1),
                          'dist_stp':( 1,  100,   1),
                          'hoff_stp':( 1,  100,   1),
                          'voff_stp':( 1,  100,   1),
                          'rota_stp':( 1,  100,   1),
                          'tilt_stp':( 1,  100,   1),
                          'bsdx_stp':( 1,  100,   1),
                }
        # Add 'key':'val to the dict to enforce 'key' to be set to 'val'
        _force = {'show_deltadd':False}
        # Opening JSON file as dict
        try:
            with open(self.path_settings_current, 'r') as of:
                pars = json.load(of)
        except: # any error is critical here!
            print(f"Error parsing settings file at: {self.path_settings_current}")
            raise SystemExit
        conv = {'geo':self.geo, 'plo':self.plo, 'thm':self.thm, 'lmt':self.lmt}
        for key, vals in pars.items():
            if key in skip:
                continue
            for p, x in vals.items():
                if p in conv[key].__dict__.keys():
                    # make sure the strings are in order.
                    # this is bad!
                    if isinstance(x, str) and x.lower() == 'none':
                        x = 'None'
                    if p in _warn and x not in range(_warn[p][0], _warn[p][1]+1):
                            print(f'WARNING: {p} set to {x}!\nAllowed values are within {_warn[p][0], _warn[p][1]}, parameter set to {_warn[p][2]}.')
                            x = _warn[p][2]
                    if p in _force:
                        x = _force[p]
                    setattr(conv[key], p, x)
                else:
                    print(f'WARNING: "{p}" is not a valid key!')
        if 'geo' not in skip:
            # store the initial values of geo
            self._geo.__dict__.update(self.geo.__dict__)

    def edit_settings_file(self, command):
        os.system(f'{command} {self.path_settings_current}')

    def get_settings_files(self):
        return sorted(map(os.path.basename, glob.glob(os.path.join(self.path_settings, '*.json'))))

    def change_settings_file(self, name):
        self.active_settings = name
        self.path_settings_current = os.path.join(self.path_settings, name)
        self.reload_settings()
    
    def import_settings_file(self):
        fname, filter = QtWidgets.QFileDialog.getOpenFileName(self, 'Import settings file', '', "Settings files (*.json)")
        if fname:
            # copy to settings folder
            shutil.copy(fname, self.path_settings)
            bname = os.path.basename(fname)
            # Add to menu
            cset_action = QtGui.QAction(bname, self, checkable=True)
            self.set_menu_action(cset_action, self.change_settings_file, bname)
            self.menu_custom_settings.addAction(cset_action)
            self.group_cset.addAction(cset_action)
            cset_action.setChecked(True)
            # change settings and reload
            self.change_settings_file(bname)

    #############
    #   EVENT   #
    #############
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
            self.calc_ref_from_cif(fpath)

    def keyPressEvent(self, event):
        k = event.key()
        # bind 'c' to cycle colormaps
        if k == QtCore.Qt.Key.Key_C:
            if event.modifiers() == QtCore.Qt.KeyboardModifier.ShiftModifier:
                inc = -1
            else:
                inc = +1
            idx = self.colormaps.index(self.geo.colormap) + inc
            if idx >= len(self.colormaps):
                idx = 0
            elif idx < 0:
                idx = len(self.colormaps) - 1
            self.change_cmap(self.colormaps[idx])
        elif k == QtCore.Qt.Key.Key_P:
            self.toggle_overlay_polarisation()
        elif k == QtCore.Qt.Key.Key_A:
            self.toggle_overlay_solidangle()
        elif k == QtCore.Qt.Key.Key_H:
            self.toggle_overlay_highlight()
        elif k == QtCore.Qt.Key.Key_T:
            self.change_units(0)
        elif k == QtCore.Qt.Key.Key_D:
            self.change_units(1)
        elif k == QtCore.Qt.Key.Key_Q:
            self.change_units(2)
        elif k == QtCore.Qt.Key.Key_S:
            self.change_units(3)
        elif k == QtCore.Qt.Key.Key_U:
            self.toggle_unit_hover()
        elif k == QtCore.Qt.Key.Key_F1:
            self.show_about_window()
        #elif k == QtCore.Qt.Key.Key_R:
        #    self.toggle_function_deltadd()

    def closeEvent(self, event):
        # Save current settings file for
        # auto load on next startup.
        # If the program crashes no
        # sctive_settings token exists
        # and the defaul is loaded.
        self.set_active_settings_file()
        event.accept()

    def hoverEvent(self, event):
        """Hover event linked to cormap
        and should only be active while
        either or both maps are displayed """
        if event.isExit():
            self.cor_label.hide()
            self.unit_label.setText(f'{self.unit_names[self.geo.unit]}')
            return
        
        # hoverEvent is only called if
        # either or both maps are active
        # -> always show on isEnter event
        if event.isEnter():
            self.cor_label.show()

        # cormap displays the product of both corrections
        # but the individual values can be retrieves from their arrays
        # -> _polcor and _solang
        y, x = map(int, np.clip(event.pos(), [0,0], np.array(self.patches['overlay'].image.shape)[::-1]-1))

        # unit label value
        if not isinstance(self._tth, float):
            unit = self.calc_unit(self._tth[x,y])
            self.unit_label.setText(f'{self.unit_names[self.geo.unit]} {unit:.2f}')
        
        _text = []
        # calc_overlays returns either a np.array (if active)
        # or a float (inactive) for _polcor and _solang. 
        if not isinstance(self._polcor, float):
            _text.append(f'P: {self._polcor[x,y]:.2f}')
        if not isinstance(self._solang, float):
            _text.append(f'S: {self._solang[x,y]:.2f}')
        #if not isinstance(self._deltadd, float):
        #    _text.append(f'log10(\u03B4d/d): {self._deltadd[x,y]:.2f}')
        self.cor_label.setText('\n'.join(_text))

    """
    def update_menu_entries(self):
        # set checkmark: references none
        for action in self.menu_ref.actions():
            if action.text() == self.geo.reference:
                action.setChecked(True)
        # set checkmark: references pyFAI
        for action in self.sub_menu_pyFAI.actions():
            if action.text() == self.geo.reference :
                action.setChecked(True)
        # set checkmark: references custom
        for action in self.sub_menu_custom.actions():
            if action.text() == self.geo.reference :
                action.setChecked(True)
        # set checkmark: beamstop none
        for action in self.menu_bs.actions():
            if action.text() == str(self.geo.bssz):
                action.setChecked(True)
        # set checkmark: beamstop custom
        for action in self.sub_menu_bs.actions():
            if action.text() == str(self.geo.bssz):
                action.setChecked(True)
        # set checkmark: colormap
        for action in self.menu_cmap.actions():
            if action.text() == self.geo.colormap:
                action.setChecked(True)
        # set checkmark: active settings
        for action in self.menu_custom_settings.actions():
            if action.text() == self.active_settings:
                action.setChecked(True)

        # set checkmark: detectors
        # - move through submenus
        for menu in self.menu_det.actions():
            if menu.text() == self.geo.det_type:
                for action in menu.menu().actions():
                    if action.text() == self.geo.det_size:
                        action.setChecked(True)
        
        # set checkmark: units
        # - self.geo.unit is int
        for num, action in enumerate(self.menu_unit.actions()):
            if num == self.geo.unit:
                action.setChecked(True)

        # set checkmark: theme
        # - self.geo.darkmode is bool
        conv = {'Light': False, 'Dark': True}
        for action in self.menu_theme.actions():
            if conv[action.text()] == self.geo.darkmode:
                action.setChecked(True)
        
        # check show_polarisation
        if self.plo.show_polarisation:
            self.polmap_action.setChecked(True)
        else:
            self.polmap_action.setChecked(False)
        
        # check show_solidangle
        if self.plo.show_solidangle:
            self.solang_action.setChecked(True)
        else:
            self.solang_action.setChecked(False)

        # menu Beamstop: add sizes list
        if self.geo.bssz not in self.geo.bs_list and not isinstance(self.geo.bssz, str):
            self.geo.bs_list.append(self.geo.bssz)
            self.geo.bs_list.sort()
        
        # update menu Beamstop
        self.sub_menu_bs.clear()
        for bs_size in sorted(self.geo.bs_list):
            bs_sub_action = QtGui.QAction(str(bs_size), self, checkable=True)
            self.set_menu_action(bs_sub_action, self.change_beamstop, bs_size)
            self.sub_menu_bs.addAction(bs_sub_action)
            self.group_bs.addAction(bs_sub_action)
            if bs_size == self.geo.bssz:
                bs_sub_action.setChecked(True)
        
        # update menu Detectors
        self.menu_det.clear()
        for d in sorted(self.detector_db):
            d_menu = QtWidgets.QMenu(d, self)
            self.menu_det.addMenu(d_menu)
            for s in self.detector_db[d]['size']:
                det_action = QtGui.QAction(s, self, checkable=True)
                self.set_menu_action(det_action, self.change_detector, d, s)
                d_menu.addAction(det_action)
                self.group_det.addAction(det_action)
                if d == self.geo.det_type and s == self.geo.det_size:
                    det_action.setChecked(True)
    """

class container(object):
    pass

class SliderWidget(QtWidgets.QFrame):
    def __init__(self, parent):
        super().__init__(parent)
        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(0)
        self.leaveEvent = self.toggle_panel
        self.enterEvent = self.toggle_panel
        self.frame = QtWidgets.QFrame()
        self.frame.setContentsMargins(0, 0, 0, 0)
        self.frame.setFixedHeight(self.parent().plo.slider_margin)
        self.layout.addWidget(self.frame)
        self.box = QtWidgets.QFrame()
        self.box.setContentsMargins(0, 5, 0, 5)
        self.layout.addWidget(self.box)
        self.box.setHidden(True)
        self.box_toggle = False
        self.box_height_show = int(np.ceil(parent.size().height()/3))
        self.box_height_hide = int(np.ceil(self.frame.size().height()))

        # prevent moving the window
        # when clicking the sliders
        self.startPos = None

        self.grid = QtWidgets.QGridLayout()
        self.grid.setContentsMargins(0, 0, 0, 0)
        self.grid.setRowStretch(1,10)
        self.box.setLayout(self.grid)

        self.apply_style()
        self.init_sliders()
        self.center_frame()

    def apply_style(self):
        self.box.setStyleSheet(f'''
            QFrame {{
                border: {self.parent().plo.slider_border_width}px solid {self.parent().slider_border_color};
                border-radius: {self.parent().plo.slider_border_radius}px;
                background: {self.parent().slider_bg_color};
            }}
            QFrame:hover {{
                background: {self.parent().slider_bg_hover};
            }}
        ''')
        self.frame.setStyleSheet(f'''
            QFrame {{
                border: {self.parent().plo.slider_border_width}px solid {self.parent().slider_border_color};
                border-radius: {self.parent().plo.slider_border_radius}px;
                background: {self.parent().slider_border_color};
            }}
        ''')

    def init_sliders(self):
        # remove sliders and labels
        for i in reversed(range(self.grid.count())): 
            self.grid.itemAt(i).widget().deleteLater()
        # add new set of sliders with updated values
        _idx = 0
        self.box_width_dynamic = 0
        if self.parent().plo.enable_slider_ener:
            self.sl_ener = self.add_slider(self.grid,
                                           self.parent().plo.slider_label_ener,'ener', _idx,
                                           self.parent().geo.ener,
                                           self.parent().lmt.ener_min,
                                           self.parent().lmt.ener_max,
                                           self.parent().lmt.ener_stp)
            self.box_width_dynamic += self.parent().plo.slider_column_width
            _idx += 1
        if self.parent().plo.enable_slider_dist:
            self.sl_dist = self.add_slider(self.grid,
                                           self.parent().plo.slider_label_dist, 'dist', _idx,
                                           self.parent().geo.dist,
                                           self.parent().lmt.dist_min,
                                           self.parent().lmt.dist_max,
                                           self.parent().lmt.dist_stp)
            self.box_width_dynamic += self.parent().plo.slider_column_width
            _idx += 1
        if self.parent().plo.enable_slider_voff:
            self.sl_voff = self.add_slider(self.grid,
                                           self.parent().plo.slider_label_voff, 'voff', _idx,
                                           self.parent().geo.voff,
                                           self.parent().lmt.voff_min,
                                           self.parent().lmt.voff_max,
                                           self.parent().lmt.voff_stp)
            self.box_width_dynamic += self.parent().plo.slider_column_width
            _idx += 1
        if self.parent().plo.enable_slider_hoff:
            self.sl_hoff = self.add_slider(self.grid,
                                           self.parent().plo.slider_label_hoff, 'hoff', _idx,
                                           self.parent().geo.hoff,
                                           self.parent().lmt.hoff_min,
                                           self.parent().lmt.hoff_max,
                                           self.parent().lmt.hoff_stp)
            self.box_width_dynamic += self.parent().plo.slider_column_width
            _idx += 1
        if self.parent().plo.enable_slider_tilt:
            self.sl_tilt = self.add_slider(self.grid,
                                           self.parent().plo.slider_label_tilt, 'tilt', _idx,
                                           self.parent().geo.tilt,
                                           self.parent().lmt.tilt_min,
                                           self.parent().lmt.tilt_max,
                                           self.parent().lmt.tilt_stp)
            self.box_width_dynamic += self.parent().plo.slider_column_width
            _idx += 1
        if self.parent().plo.enable_slider_rota:
            self.sl_rota = self.add_slider(self.grid,
                                           self.parent().plo.slider_label_rota, 'rota', _idx,
                                           self.parent().geo.rota,
                                           self.parent().lmt.rota_min,
                                           self.parent().lmt.rota_max,
                                           self.parent().lmt.rota_stp)
            self.box_width_dynamic += self.parent().plo.slider_column_width
            _idx += 1
        if self.parent().plo.enable_slider_bsdx:
            self.sl_bsdx = self.add_slider(self.grid,
                                           self.parent().plo.slider_label_bsdx, 'bsdx', _idx,
                                           self.parent().geo.bsdx,
                                           self.parent().lmt.bsdx_min,
                                           self.parent().lmt.bsdx_max,
                                           self.parent().lmt.bsdx_stp)
            self.box_width_dynamic += self.parent().plo.slider_column_width
            # set beamstop distance slider max limit to detector distance
            self.update_slider_limits(self.sl_bsdx, self.parent().lmt.bsdx_min, self.parent().geo.dist)
            _idx += 1
        
        self.resize(self.box_width_dynamic, self.box_height_hide)

    def center_frame(self):
        self.move(int((self.parent().size().width()-self.box_width_dynamic)/2), self.parent().offset_win32)

    def update_slider_label(self, label, value):
        label.setText(str(int(value)))

    def update_slider_limits(self, slider, lmin, lmax):
        slider.setRange(int(lmin), int(lmax))

    def add_slider(self, layout, label, token, idx, lval, lmin, lmax, lstp):
        font = QtGui.QFont()
        font.setPixelSize(self.parent().plo.slider_label_size)

        label_name = QtWidgets.QLabel(label)
        label_name.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        label_name.setFont(font)
        label_name.setStyleSheet(f'''
            QLabel {{
                color: {self.parent().slider_label_color};
                border: 0px solid none;
                background: transparent;
            }}
        ''')
        
        label_value = QtWidgets.QLabel()
        label_value.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        label_value.setFont(font)
        label_value.setText(str(int(lval)))
        label_value.setStyleSheet(f'''
            QLabel {{
                color: {self.parent().slider_label_color};
                border: 0px solid transparent;
                background: transparent;
            }}
        ''')

        slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Vertical, objectName=token)
        slider.setRange(int(lmin), int(lmax))
        slider.setSingleStep(int(lstp))
        slider.setPageStep(int(lstp))
        slider.setValue(int(lval))
        slider.valueChanged.connect(self.parent().update_screen)
        slider.valueChanged.connect(lambda value: self.update_slider_label(label_value, value))

        layout.addWidget(label_name, 0, idx, QtCore.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(slider, 1, idx, QtCore.Qt.AlignmentFlag.AlignHCenter)
        layout.addWidget(label_value, 2, idx, QtCore.Qt.AlignmentFlag.AlignCenter)

        return slider#(slider, label_name, label_value)

    def toggle_panel(self, event):
        if isinstance(event, QtGui.QEnterEvent):
            #self.box.setHidden(not self.box.isHidden())
            self.box.setHidden(False)
            self.resize(self.box_width_dynamic, self.box_height_show)
        elif isinstance(event, QtCore.QEvent) and not self.box_toggle:
            self.box.setHidden(True)
            self.resize(self.box_width_dynamic, self.box_height_hide)
        else:
            pass

    def mousePressEvent(self, event):
        super().mousePressEvent(event)
        if event.button() == QtCore.Qt.MouseButton.LeftButton:
            self.startPos = event.pos()
            self.box_toggle = not self.box_toggle

    def mouseMoveEvent(self, event):
        super().mouseMoveEvent(event)
        if event.buttons() == QtCore.Qt.MouseButton.LeftButton:
            if not self.startPos:
                event.ignore()
                return
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

    def mouseReleaseEvent(self, event):
        super().mouseReleaseEvent(event)
        self.startPos = None
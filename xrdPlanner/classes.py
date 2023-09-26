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

        # enable antialiasing
        pg.setConfigOptions(antialias=True)

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

        # save/load parameters to/from file
        self.init_par()

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
        if self.geo.bssz and self.geo.bssz != 'None' and self.geo.bssz not in self.geo.bs_list:
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

        # add the slider frame
        # this calls draw_conics(), 
        # make sure that everything
        # that is needed is initialised
        self.sliderWidget = SliderWidget(self)

    def get_active_settings_file(self):
        if os.path.exists(self.path_settings_token):
            with open(self.path_settings_token, 'r') as rf:
                name = rf.read()
            if os.path.exists(os.path.join(self.path_settings, name)):
                return name
        # otherwise return default
        _default = os.path.basename(self.path_settings_default)
        self.set_active_settings_file(_default)
        return _default
    
    def set_active_settings_file(self, active):
        if not os.path.exists(self.path_settings):
            os.makedirs(self.path_settings)
        with open(self.path_settings_token, 'w') as wf:
            wf.write(active)

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

        # init the hkl tooltip
        font = QtGui.QFont()
        font.setPixelSize(self.plo.conic_hkl_label_size)
        font.setBold(True)
        QtWidgets.QToolTip.setFont(font)

    def apply_theme(self, use_dark, redraw=False):
        # set darkmode
        self.geo.darkmode = use_dark
        _color_dark = QtGui.QColor(self.thm.color_dark)
        _color_light = QtGui.QColor(self.thm.color_light)
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
            # slider
            self.slider_border_color = self.thm.dark_slider_border_color
            self.slider_bg_color = self.thm.dark_slider_bg_color
            self.slider_bg_hover = self.thm.dark_slider_bg_hover
            self.slider_label_color = self.thm.dark_slider_label_color
            # palette
            palette = QtGui.QPalette()
            palette.setColor(QtGui.QPalette.ColorRole.Window,          _color_dark)
            palette.setColor(QtGui.QPalette.ColorRole.ButtonText,      _color_light)
            palette.setColor(QtGui.QPalette.ColorRole.Base,            _color_dark)
            palette.setColor(QtGui.QPalette.ColorRole.Text,            _color_light)
            palette.setColor(QtGui.QPalette.ColorRole.HighlightedText, _color_dark)
            palette.setColor(QtGui.QPalette.ColorRole.WindowText,      _color_light)
            #palette.setColor(QtGui.QPalette.ColorRole.Button,          self.cont_cmap.map(1.0, mode='qcolor'))
            palette.setColor(QtGui.QPalette.ColorRole.Highlight,       self.cont_cmap.map(0.0, mode='qcolor'))
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
            # slider
            self.slider_border_color = self.thm.light_slider_border_color
            self.slider_bg_color = self.thm.light_slider_bg_color
            self.slider_bg_hover = self.thm.light_slider_bg_hover
            self.slider_label_color = self.thm.light_slider_label_color

            # palette
            palette = QtGui.QPalette()
            palette.setColor(QtGui.QPalette.ColorRole.Window,          _color_light)
            palette.setColor(QtGui.QPalette.ColorRole.ButtonText,      _color_dark)
            palette.setColor(QtGui.QPalette.ColorRole.Base,            _color_light)
            palette.setColor(QtGui.QPalette.ColorRole.Text,            _color_dark)
            palette.setColor(QtGui.QPalette.ColorRole.HighlightedText, _color_light)
            palette.setColor(QtGui.QPalette.ColorRole.WindowText,      _color_dark)
            #palette.setColor(QtGui.QPalette.ColorRole.Button,          self.cont_cmap.map(1.0, mode='qcolor'))
            palette.setColor(QtGui.QPalette.ColorRole.Highlight,       self.cont_cmap.map(0.0, mode='qcolor'))

        # apply palette to app
        app = QtWidgets.QApplication.instance()
        app.setPalette(palette)
        app.setStyle('Fusion')
        
        # redraw the canvas
        if redraw:
            self.redraw_canvas()
        
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

        # add empty plot per contour line
        for i in range(self.plo.conic_tth_num):
            curve = pg.PlotCurveItem(useCache=True)
            self.ax.addItem(curve)
            self.patches['conic'].append(curve)
            temp_label = pg.TextItem(anchor=(0.5,0.5), fill=pg.mkBrush(self.conic_label_fill))
            temp_label.setFont(font)
            self.patches['labels'].append(temp_label)
            self.ax.addItem(temp_label)

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

        # label for beamstop contour
        bs_label = pg.TextItem(anchor=(0.5,0.5),
                               color=self.unit_label_color,
                               fill=self.conic_label_fill)
        bs_label.setFont(font)
        self.patches['bs_label'] = bs_label
        self.ax.addItem(bs_label)

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
        self.menu_det = self.menu_bar.addMenu('Detector')
        group_det = QtGui.QActionGroup(self)
        group_det.setExclusive(True)

        # menu Detectors
        for d in sorted(self.detector_db):
            d_menu = QtWidgets.QMenu(d, self)
            self.menu_det.addMenu(d_menu)
            for s in self.detector_db[d]['size']:
                det_action = QtGui.QAction(s, self, checkable=True)
                self.set_menu_action(det_action, self.change_detector, d, s)
                d_menu.addAction(det_action)
                group_det.addAction(det_action)
                if d == self.geo.det_type and s == self.geo.det_size:
                    det_action.setChecked(True)
        
        # menu Reference
        self.menu_ref = self.menu_bar.addMenu('Reference')
        self.group_ref = QtGui.QActionGroup(self)
        self.group_ref.setExclusive(True)
        # menu Reference: add None
        ref_action = QtGui.QAction('None', self, checkable=True)
        self.set_menu_action(ref_action, self.change_reference, 'None')
        self.menu_ref.addAction(ref_action)
        self.group_ref.addAction(ref_action)
        if self.geo.reference == 'None':
            ref_action.setChecked(True)
        # menu Reference: add pyFAI library
        self.sub_menu_pyFAI = QtWidgets.QMenu('pyFAI', self)
        self.menu_ref.addMenu(self.sub_menu_pyFAI)
        for ref_name in sorted(self.ref_library):
            ref_action = QtGui.QAction(ref_name, self, checkable=True)
            self.set_menu_action(ref_action, self.change_reference, ref_name)
            self.sub_menu_pyFAI.addAction(ref_action)
            self.group_ref.addAction(ref_action)
            if ref_name == self.geo.reference:
                ref_action.setChecked(True)

        # menu Reference: add Custom
        self.sub_menu_custom = QtWidgets.QMenu('Custom', self)
        self.menu_ref.addMenu(self.sub_menu_custom)
        
        # menu Beamstop
        self.menu_bs = self.menu_bar.addMenu('Beamstop')
        self.group_bs = QtGui.QActionGroup(self)
        self.group_bs.setExclusive(True)
        # menu Beamstop: add None
        bs_action = QtGui.QAction('None', self, checkable=True)
        self.set_menu_action(bs_action, self.change_beamstop, 'None')
        self.menu_bs.addAction(bs_action)
        self.group_bs.addAction(bs_action)
        if self.geo.bssz == 'None':
            bs_action.setChecked(True)
        # menu Beamstop: add sizes list
        self.sub_menu_bs = QtWidgets.QMenu('Sizes [mm]', self)
        self.menu_bs.addMenu(self.sub_menu_bs)
        for bs_size in sorted(self.geo.bs_list):
            bs_sub_action = QtGui.QAction(str(bs_size), self, checkable=True)
            self.set_menu_action(bs_sub_action, self.change_beamstop, bs_size)
            self.sub_menu_bs.addAction(bs_sub_action)
            self.group_bs.addAction(bs_sub_action)
            if bs_size == self.geo.bssz:
                bs_sub_action.setChecked(True)
        
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

        # menu View
        menu_view = self.menu_bar.addMenu('View')
        # submenu Theme
        # non-standard menu -> self.geo.darkmode is bool
        self.menu_theme = menu_view.addMenu('Theme')
        group_theme = QtGui.QActionGroup(self)
        group_theme.setExclusive(True)
        for (theme, invert) in [('Light', False), ('Dark', True)]:
            theme_action = QtGui.QAction(theme, self, checkable=True)
            self.set_menu_action(theme_action, self.apply_theme, invert, True)
            group_theme.addAction(theme_action)
            self.menu_theme.addAction(theme_action)
            if invert == self.geo.darkmode:
                theme_action.setChecked(True)

        # submenu Colormap
        self.menu_cmap = menu_view.addMenu('Colormap')
        group_cmap = QtGui.QActionGroup(self)
        group_cmap.setExclusive(True)
        for cmap_name in sorted(pg.colormap.listMaps()):
            cmap_action = QtGui.QAction(cmap_name, self, checkable=True)
            self.set_menu_action(cmap_action, self.change_cmap, cmap_name)
            self.menu_cmap.addAction(cmap_action)
            group_cmap.addAction(cmap_action)
            if cmap_name == self.geo.colormap:
                cmap_action.setChecked(True)
        
        # menu Settings
        menu_Settings = self.menu_bar.addMenu('Settings')
        # import action
        import_action = QtGui.QAction('Import settings file', self)
        self.set_menu_action(import_action, self.import_settings_file)
        menu_Settings.addAction(import_action)
        # submenu Settings files
        self.menu_custom_settings = menu_Settings.addMenu('Load settings file')
        self.group_cset = QtGui.QActionGroup(self)
        self.group_cset.setExclusive(True)
        for cset_name in self.get_settings_files():
            cset_action = QtGui.QAction(cset_name, self, checkable=True)
            self.set_menu_action(cset_action, self.change_settings_file, cset_name)
            self.menu_custom_settings.addAction(cset_action)
            self.group_cset.addAction(cset_action)
            if cset_name == self.active_settings:
                cset_action.setChecked(True)
        # save action
        save_action = QtGui.QAction('Save current settings', self)
        self.set_menu_action(save_action, self.save_current_settings)
        menu_Settings.addAction(save_action)
        # submenu Edit files
        menu_Settings.addSeparator()
        # submenu Edit
        menu_Edit = menu_Settings.addMenu('Edit files')
        # prepare platform dependent file reader
        if sys.platform == 'win32':
            tokens = [('Edit detector db file', os.system, f'notepad {self.path_detdb}'),
                      ('Edit current settings file', self.edit_settings_file, f'notepad')]
        elif sys.platform == 'linux':
            tokens = [('Edit detector db file', os.system, f'xdg-open {self.path_detdb}'),
                      ('Edit current settings file', self.edit_settings_file, f'xdg-open')]
        else:
            tokens = [('Edit detector db file', os.system, f'open -t {self.path_detdb}'),
                      ('Edit current settings file', self.edit_settings_file, f'open -t')]
        for (name, funct, command) in tokens:
            edit_action = QtGui.QAction(name, self)
            self.set_menu_action(edit_action, funct, command)
            menu_Edit.addAction(edit_action)
        # reset action
        #reset_action = QtGui.QAction('Reset to default settings', self)
        #self.set_menu_action(reset_action, self.reset_to_default)
        #menu_default.addAction(reset_action)

        # menu About
        self.menu_help = self.menu_bar.addMenu('Help')
        self.action_about = QtGui.QAction('xrdPlanner', self)
        self.set_menu_action(self.action_about, self.about_window)
        self.menu_help.addAction(self.action_about)

    def about_window(self):
        msgBox = QtWidgets.QMessageBox()
        msgBox.setTextFormat(QtCore.Qt.TextFormat.RichText)
        msgBox.setStandardButtons(QtWidgets.QMessageBox.StandardButton.NoButton)
        msgBox.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.LinksAccessibleByMouse)
        windowIcon = QtGui.QIcon()
        windowIcon.addPixmap(QtGui.QPixmap(':/icons/xrdPlanner.png'))
        msgBox.setWindowIcon(windowIcon)
        msgBox.setIconPixmap(QtGui.QPixmap(':/icons/xrdPlanner.png'))
        msgBox.setWindowTitle('About')
        msgBox.setText('xrdPlanner 1.3.0')
        msgBox.setDetailedText("<a href='https://github.com/LennardKrause/xrdPlanner'>https://github.com/LennardKrause/xrdPlanner</a>")
        msgBox.setInformativeText('Author:\nLennard Krause, 2023\nlkrause@chem.au.dk')
        msgBox.exec()

    def edit_settings_file(self, command):
        os.system(f'{command} {self.path_settings_current}')

    def get_settings_files(self):
        return sorted(map(os.path.basename, glob.glob(os.path.join(self.path_settings, '*.json'))))

    def change_settings_file(self, name):
        self.active_settings = name
        self.path_settings_current = os.path.join(self.path_settings, name)
        self.set_active_settings_file(name)
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
        
        # menu Beamstop: add sizes list
        if self.geo.bssz not in self.geo.bs_list:
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

    def add_unit_label(self):
        font = QtGui.QFont()
        font.setPixelSize(self.plo.unit_label_size)
        self.unit_label = pg.TextItem(anchor=(0.0,0.0), color=self.unit_label_color, fill=self.unit_label_fill)
        self.unit_label.setText(self.unit_names[self.geo.unit])
        self.unit_label.setFont(font)
        self.ax.addItem(self.unit_label)
        self.unit_label.setPos(-self.xdim, self.ydim)

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
        width = int(np.ceil(self.plo.plot_size * self.xdim / self.ydim))
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
            self.setWindowTitle(f'{self.det.name} - {self.active_settings}')
        else:
            self.setWindowTitle(f'{self.det.name} - {self.geo.reference} - {self.active_settings}')

    def redraw_canvas(self):
        # save darkmode toggle
        self.save_par()
        self.init_modifiables()
        # clear the screen
        self.ax.clear()
        # re-initialise
        self.init_screen()
        # center the slider frame
        self.sliderWidget.apply_style()
        self.sliderWidget.init_sliders()
        self.sliderWidget.center_frame()

    def change_beamstop(self, size):
        self.geo.bssz = size
        self.redraw_canvas()

    def change_cmap(self, cmap):
        self.geo.colormap = cmap
        self.redraw_canvas()

    def reload_settings(self):
        # load settings
        self.load_par()
        self.redraw_canvas()
        self.update_menu_entries()

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
        self.draw_conics()

    def change_reference(self, ref_name):
        self.geo.reference = ref_name
        self.get_reference()
        self.draw_reference()

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
        geo.ener = 21            # [keV]  Beam energy
        geo.dist = 75            # [mm]   Detector distance
        geo.voff = 0             # [mm]   Detector offset (vertical)
        geo.hoff = 0             # [mm]   Detector offset (horizontal)
        geo.rota = 25            # [deg]  Detector rotation
        geo.tilt = 0             # [deg]  Detector tilt
        geo.bssz = 'None'        # [mm]   Current beamstop size (or 'None')
        geo.bsdx = 40            # [mm]   Beamstop distance
        geo.unit = 1             # [0-3]  Contour legend
                                 #          0: 2-Theta
                                 #          1: d-spacing
                                 #          2: q-space
                                 #          3: sin(theta)/lambda
        geo.reference = 'None'   # [str]  Plot reference contours
                                 #          pick from pyFAI
        geo.darkmode = False     # [bool] Darkmode
        geo.colormap = 'viridis' # [cmap] Contour colormap
        geo.bs_list = [1.5,      # [list] Available beamstop sizes
                       2.0,
                       2.5,
                       3.0,
                       5.0]
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
        plo.conic_ref_cif_irel = True       # [int]    Linewidth relative to intensity
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
        plo.plot_size_fixed = True          # [int]    Fix window size
        plo.unit_label_size = 16            # [int]    Label size, px
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
        thm.color_dark = '#404040'                # [color]  Global dark color
        thm.color_light = '#EEEEEE'               # [color]  Global light color
        # light mode
        thm.light_conic_label_fill = '#FFFFFF'    # [color]  Contour label fill color
        thm.light_conic_ref_color = '#404040'     # [color]  Reference contour color
        thm.light_beamstop_color = '#DCDCDC'      # [color]  Beamstop color
        thm.light_beamstop_edge_color = '#EEEEEE' # [color]  Beamstop edge color
        thm.light_det_module_color = '#404040'    # [color]  Detector module border color
        thm.light_det_module_fill = '#404040'     # [color]  Detector module background color
        thm.light_plot_bg_color = '#FFFFFF'       # [color]  Plot background color
        thm.light_unit_label_color = '#808080'    # [color]  Label color
        thm.light_unit_label_fill = '#FFFFFF'     # [color]  Label fill color
        thm.light_slider_border_color = '#808080' # [color]  Slider frame border color
        thm.light_slider_bg_color = '#AAC0C0C0'   # [color]  Slider frame background color
        thm.light_slider_bg_hover = '#C0C0C0'     # [color]  Slider frame hover color
        thm.light_slider_label_color = '#000000'  # [color]  Slider frame label color
        # dark mode
        thm.dark_conic_label_fill = '#000000'     # [color]  Contour label fill color
        thm.dark_conic_ref_color = '#202020'      # [color]  Reference contour color
        thm.dark_beamstop_color = '#202020'       # [color]  Beamstop color
        thm.dark_beamstop_edge_color = '#404040'  # [color]  Beamstop edge color
        thm.dark_det_module_color = '#EEEEEE'     # [color]  Detector module border color
        thm.dark_det_module_fill = '#EEEEEE'      # [color]  Detector module background color
        thm.dark_plot_bg_color = '#000000'        # [color]  Plot background color
        thm.dark_unit_label_color = '#C0C0C0'     # [color]  Label color
        thm.dark_unit_label_fill = '#000000'      # [color]  Label fill color
        thm.dark_slider_border_color = '#202020'  # [color]  Slider frame border color
        thm.dark_slider_bg_color = '#AA303030'    # [color]  Slider frame background color
        thm.dark_slider_bg_hover = '#303030'      # [color]  Slider frame hover color
        thm.dark_slider_label_color = '#C0C0C0'   # [color]  Slider frame label color

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

    def draw_conics(self):
        # calculate the offset of the contours resulting from voff and rotation
        # shift the grid to draw the cones, to make sure the contours are drawn
        # within the visible area

        # convert theta in degrees to radians
        # for some reason I defined it negative some time ago
        # now there's no turning back!
        omega = -np.deg2rad(self.geo.tilt + self.geo.rota)

        # beamcenter shift
        _comp_shift = -(self.geo.voff - self.geo.dist * np.tan(omega) - np.deg2rad(self.geo.tilt) * self.geo.dist)
        
        # update beam center
        self.patches['beamcenter'].setData([self.geo.hoff],[_comp_shift])
        # update beam center
        self.patches['poni'].setData([self.geo.hoff],[-(self.geo.voff - np.deg2rad(self.geo.tilt)*self.geo.dist)])

        if self.geo.bssz and self.geo.bssz != 'None':
            # update beam stop
            bs_theta = np.tan((self.geo.bssz/2) / self.geo.bsdx)

            # calculate the conic section corresponding to the theta angle
            # :returns False is conic is outside of visiblee area
            x, y, label_pos = self.calc_conic(omega, bs_theta, steps=self.plo.conic_steps)
            if x is False:
                self.patches['beamstop'].setVisible(False)
                self.patches['bs_label'].setVisible(False)
            else:
                # plot the conic section
                self.patches['beamstop'].setData(x, y, fillLevel=y.max())
                self.patches['beamstop'].setVisible(True)

                # Conversion factor keV to Angstrom: 12.398
                # sin(t)/l: np.sin(Theta) / lambda -> (12.398/geo_energy)
                stl = np.sin(bs_theta/2)/(12.398/self.geo.ener)

                # d-spacing: l = 2 d sin(t) -> 1/2(sin(t)/l)
                dsp = 1/(2*stl)

                # prepare the values in the different units / labels
                #_units = {0:np.rad2deg(theta), 1:_dsp, 2:_stl*4*np.pi, 3:_stl}
                _units = {0:np.rad2deg(bs_theta), 1:dsp, 2:stl*4*np.pi, 3:stl}
                self.patches['bs_label'].setPos(self.geo.hoff, label_pos)
                self.patches['bs_label'].setText(f'{_units[self.geo.unit]:.2f}')
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
            self.patches['labels'][_n].setPos(self.geo.hoff, label_pos)
            self.patches['labels'][_n].setText(f'{_units[self.geo.unit]:.2f}', color=self.cont_cmap.map(_f, mode='qcolor'))
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
                omega = -np.deg2rad(self.geo.tilt + self.geo.rota)
                
                # calculate the conic section corresponding to the theta angle
                # :returns False is conic is outside of visiblee area
                x, y, _ = self.calc_conic(omega, theta, steps=self.plo.conic_steps)
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

    def show_tooltip(self, widget, event):
        if not widget.name or not self.cont_ref_hkl:
            event.ignore()
            return False
        text = QtWidgets.QTextEdit(str(widget.name))
        text.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        pos = QtCore.QPoint(*map(int, event.screenPos()))# - QtCore.QPoint(10,20)
        QtWidgets.QToolTip.showText(pos, text.toHtml())
        event.ignore()

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

    def init_par(self):
        # fetch the geometry, detector, plot specifications and limits
        self.get_defaults_all()
        # file name to store current settings
        # if file_dump doesn't exists, make a dump
        if not os.path.exists(self.path_settings_current):
            self.save_par()
        # if it exists load parameters
        else:
            self.load_par()
        # reset to default if requested
        if self.plo.reset_settings:
            self.get_defaults_all()
            self.save_par()
        # update with default if requested
        #  - e.g. to add missing entries
        if self.plo.update_settings:
            self.save_par()

    def save_current_settings(self):
        # self.geo is edited with the sliders
        # self._geo holds the initial values
        # usually I do not want to overwrite 
        # the startup values -> I write _geo
        # to the settings file.
        # unless this function is called!
        self._geo.__dict__.update(self.geo.__dict__)
        self.save_par()

    def save_par(self):
        # Writing geo as dict to file
        with open(self.path_settings_current, 'w') as wf:
            json.dump({'geo':self._geo.__dict__, 'plo':self.plo.__dict__, 'thm':self.thm.__dict__, 'lmt':self.lmt.__dict__}, wf, indent=4)

    def load_par(self, skip=[]):
        # Some parameters need to be protected to save
        # the user some waiting time
        #
        # A check is performed if the value of the key
        # is within 'minimum' and 'maximum', if not it
        # is set to 'default'
        #
        # Add 'key':('minimum', 'maximum', 'default') to the dict
        _warn = {'conic_ref_cif_kev':( 5,   25,  12),
                     'conic_tth_min':( 1,   10,   5),
                     'conic_tth_max':(10,  180,  90),
                     'conic_tth_num':( 1,  100,  20),
                     'conic_ref_num':( 1,  500, 200),
                       'conic_steps':(10, 1000, 100),
                          'ener_stp':( 1,  100,   1),
                          'dist_stp':( 1,  100,   1),
                          'hoff_stp':( 1,  100,   1),
                          'voff_stp':( 1,  100,   1),
                          'rota_stp':( 1,  100,   1),
                          'tilt_stp':( 1,  100,   1),
                          'bsdx_stp':( 1,  100,   1),
                }
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
                    if type(x) == str and x.lower() == 'none':
                        x = 'None'
                    if p in _warn and x not in range(_warn[p][0], _warn[p][1]+1):
                            print(f'WARNING: {p} set to {x}!\nAllowed values are within {_warn[p][0], _warn[p][1]}, parameter set to {_warn[p][2]}.')
                            x = _warn[p][2]
                    setattr(conv[key], p, x)
                else:
                    print(f'WARNING: "{p}" is not a valid key!')
        if 'geo' not in skip:
            # store the initial values of geo
            self._geo.__dict__.update(self.geo.__dict__)

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
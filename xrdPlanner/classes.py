import os
import sys
import json
import glob
import shutil
import numpy as np
from scipy.optimize import curve_fit
import pyqtgraph as pg
import Dans_Diffraction as dif
from PyQt6 import QtWidgets, QtCore, QtGui
from pyFAI import calibrant
import xrdPlanner.resources

# Add the Absorption window and connect scattering diameter slider (from FWHM)
# change pxrd scatterplot highlight to use dedicated highlighter (scatterplot)
# check what windows need update on theme change
# out-class pxrd window

# work out how to update fwhm window only in visible

# debug fn to show who is calling what and when to find out why.
#import inspect
#def print_stack(context=True):
#    print('>----|')
#    if context:
#        print('\n'.join([f"{x.lineno:>5}| {x.function} > {''.join(x.code_context).strip()}" for x in inspect.stack()][1:][::-1]))
#    else:
#        print('\n'.join([f"{x.lineno:>5}| {x.function}" for x in inspect.stack()][1:][::-1]))

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        """
        Initializes the main application window and sets up various components and configurations.
        """
        super().__init__(*args, **kwargs)

        # set path home
        self.path_home = os.path.dirname(__file__)
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

        # list to keep track of currently highlighted contours
        self.highlight_timers = []
        # default unit cell parameters for custom cell window
        self.default_custom_cell = [6,6,6,90,90,90]
        # set path to settings folder
        self.path_settings = os.path.join(self.path_home, 'settings',)
        # set path to active settings token
        self.path_settings_token = os.path.join(self.path_settings, '.active',)
        # set path to default settings file
        self.path_settings_default = os.path.join(self.path_settings, 'default.json')
        # set path to current settings file
        self.active_settings = self.settings_get_active()
        self.path_settings_current = os.path.join(self.path_settings, self.active_settings)
        # set path to detector database
        self.path_detdb = os.path.join(self.path_home, 'detector_db.json')
        # set path to cif file paths
        self.path_cif_db = os.path.join(self.path_settings, 'cif_db.json')
        # initialize powder diffraction plot window
        self.pxrd_win = None
        # dicts to store custom reference data
        self.ref_cif = {}
        self.ref_cell = {}
        self.cont_ref_dsp = None
        self.cont_ref_hkl = None
        self.xtl = None
        # What standards should be available as reference
        # The d spacings will be imported from pyFAI
        self.ref_pyfai = calibrant.names()

        # move settings file from old location
        _old_settings_file = os.path.join(self.path_home, 'settings.json')
        if os.path.exists(_old_settings_file):
            shutil.move(_old_settings_file, self.path_settings_default)

        # delete active settings file token
        #  - token with currently active settings
        #    file will be placed on successful exit
        #  - if no token file is found -> crashed last time
        #  - load default settings
        self.settings_del_active()

        # save/load parameters to/from file
        self.params_init()

        # Store available colormaps
        # PyQtGraph.colormap.listMaps(): Experimental, subject to change.
        self.colormaps = sorted(pg.colormap.listMaps())

        # get the translations to link
        # the settings parameter keys
        # to their description
        self.get_tooltips()

        # get X-ray attenuation lengths
        # for the FWHM / detector sensor
        # thickness calculation
        self.get_att_lengths()
        
        # set to False to 'test' windows os behaviour
        self.menuBar().setNativeMenuBar(self.plo.use_native_menubar)

        # menubar is displayed within the main window on Windows
        # so we need to make space for it
        # no idea about other OS, if there are issues fix them here
        if sys.platform in ['win32', 'linux', 'linux2'] or not self.menuBar().isNativeMenuBar():
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
        #self.ax.scene().sigMouseMoved.connect(print)
        
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

        # initialise all that depends on the settings
        # call this function to apply changes were made
        # to the settings file -> settings_reload()
        self.modifiables_init()
        
        # initialize about window
        self.about_win = AboutWindow(parent=self,
                                     path_settings=self.path_settings,
                                     path_home=self.path_home,
                                     pixmap=self.pixmap,
                                     icon=self.icon)
        # initialize geometry window
        self.geometry_win = GeometryWindow(parent=self)
        # initialize hotkey window
        self.hotkeys_win = HotkeysWindow(parent=self)
        # initialize detdb window
        self.detdb_win = DetdbWindow(parent=self)
        # initialize export window
        self.export_win = ExportWindow(parent=self)
        # initialize fwhm parameter window
        self.fwhm_win = FwhmWindow(parent=self)
        # initialize unit cell window
        self.uc_win = UnitCellWindow(parent=self)
        # initialize absorption window
        #self.abs_win = AbsorptionWindow(parent=self)

        # get the hotkeys
        # self.hotkey_desc: list of tuples [(hotkey, description)]
        # self.hotkey_dict: dictionary of hotkeys (key, modifier): function
        self.get_hotkeys()
        self.hotkeys_win.add_hotkeys(self.hotkey_dict)

        # populate the menus with detectors, references and units
        self.menu_init()
        
        # load stored cif file links
        self.get_ref_db_from_file()

        # initialize the screen
        self.main_screen_init()

        # add an icon
        self.setWindowIcon(self.icon)

        # disable/reset delta_d/d toggles
        self.action_funct_fwhm_show.setEnabled(False)
        self.action_funct_fwhm_export.setEnabled(False)

        # add the slider frame
        # this calls draw_conics(), 
        # make sure that everything
        # that is needed is initialised
        self.sliderWidget = SliderWidget(self)

    def modifiables_init(self):
        """
        Initializes modifiable parameters for the plotting environment.
        This method performs the following tasks:
        - Retrieves and optionally reverses the colormap based on the dark mode setting.
        - Applies the theme based on the dark mode setting.
        - Sets the background color of the plotting window.
        - Generates contour levels for the plot.
        - Translates units for the plot title.
        - Validates the unit index and raises an error if it is out of range.
        - Retrieves and updates the detector specifications.
        - Selects the current detector parameters.
        - Removes unavailable detectors from the detector database based on the detector bank.
        - Updates the generic window settings.
        - Initializes the tooltip for hkl labels with specified font settings.
        Raises:
            SystemExit: If the unit index is out of the valid range.
        """
        # get colormap
        self.cont_cmap = pg.colormap.get(self.geo.colormap, skipCache=True)
        # reverse the colormap useful to increase visibility in darkmode
        if self.geo.darkmode:
            self.cont_cmap.reverse()

        # experimental darkmode?
        self.theme_apply(self.geo.darkmode)

        # set window color
        self.ax.setBackground(self.plot_bg_color)

        # generate contour levels
        self.cont_geom_num = np.linspace(self.plo.conic_tth_min, self.plo.conic_tth_max, self.plo.conic_tth_num)

        # translate unit for plot title
        self.unit_names = ['2\u03B8 [\u00B0]',
                           'd [\u212B]',
                           'Q [\u212B\u207B\u00B9]',
                           'sin(\u03B8)/\u03BB [\u212B\u207B\u00B9]']
        if self.geo.unit >= len(self.unit_names):
            print(f'Error: Valid geo.unit range is from 0 to {len(self.unit_names)-1}, geo.unit={self.geo.unit}')
            raise SystemExit

        # get the detector specs
        # - update: overwrite existing file after load
        # - reset: overwrite existing file with defaults
        self.detector_db = self.get_det_library(update=self.plo.update_det_bank, reset=self.plo.reset_det_bank)

        # pick current detector
        self.det = self.get_det_params()
        
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
        if self.ref_cif:
            for name in self.ref_cif.keys():
                self.add_cif_to_menu(name)

        self.update_win_generic()

        # init the hkl tooltip
        font = QtGui.QFont()
        font.setPixelSize(self.plo.conic_hkl_label_size)
        font.setBold(True)
        QtWidgets.QToolTip.setFont(font)

    def params_init(self):
        """
        Initialize parameters for the xrdPlanner class.

        This method performs the following steps:
        1. Fetches the default geometry, detector, plot specifications, and limits.
        2. Checks if the settings file exists:
           - If it does not exist, it saves the current settings to a file.
           - If it exists, it loads the parameters from the file.
        3. Resets to default settings if requested.
        4. Updates the settings file with default values if requested, 
           which can be useful for adding missing entries.

        Attributes:
            self.path_settings_current (str): Path to the current settings file.
            self.plo.reset_settings (bool): Flag to reset settings to default.
            self.plo.update_settings (bool): Flag to update settings with default values.
        """
        # fetch the geometry, detector, plot specifications and limits
        self.get_defaults_all()
        # file name to store current settings
        # if file_dump doesn't exists, make a dump
        if not os.path.exists(self.path_settings_current):
            self.settings_save_to_file()
        # if it exists load parameters
        else:
            self.settings_load_from_file()
        # reset to default if requested
        if self.plo.reset_settings:
            self.get_defaults_all()
            self.settings_save_to_file()
        # update with default if requested
        #  - e.g. to add missing entries
        if self.plo.update_settings:
            self.settings_save_to_file()

    def main_screen_init(self):
        """
        Initializes the main screen for the application.

        This method sets up the plot for contours and beam center, configures various
        plot settings, and initializes containers and items for contour lines, beam stop,
        reference lines, and other plot elements.

        Steps performed:
        1. Locks the aspect ratio of the plot.
        2. Hides the bottom and left axes.
        3. Disables mouse pan/zoom and right-click context menu.
        4. Hides the autoscale button.
        5. Enables debug mode functions if debug mode is set.
        6. Defines the font for labeling conics.
        7. Initializes containers for various plot elements.
        8. Adds a scatter plot for the beam stop.
        9. Adds empty plots for reference contour lines and makes them clickable.
        10. Adds empty plots for contour lines and labels them.
        11. Adds a label for the beam stop contour.
        12. Adds scatter plots for the poni and beam center.
        13. Adds an overlay image item and sets its color map.
        14. Builds detector modules.
        15. Resizes the window and plot to fit proper dimensions.
        16. Initializes unit and correction labels.
        17. Creates cones and draws contour lines.
        18. Sets the window title.
        """
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
            #self.ax.setMenuEnabled(True)

        # define font to label conics
        font = QtGui.QFont()
        font.setPixelSize(self.plo.conic_label_size)
        font.setBold(True)

        # container for contour lines
        self.patches = {'beamcenter':None,
                        'beamstop':None,
                        'overlay':None,
                        #'fwhm':None,
                        #'isocurve':None,
                        'poni':None,
                        'bs_label':None,
                        'conic':[],
                        'reference':[],
                        'ref_hl_label':None,
                        'ref_hl_curve':None,
                        'labels':[],
                        'polar_grid':[]}

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
            ref = HoverableCurveItem(self, useCache=True, mouseWidth=1)
            self.ax.addItem(ref)
            self.patches['reference'].append(ref)
            self.patches['reference'][i].name = None
            self.patches['reference'][i].index = None

        # add ref highlight curve
        ref_hkl_curve = pg.PlotCurveItem(useCache=True,
                                         pen=pg.mkPen(color=self.conic_highlight,
                                                      width=self.plo.conic_ref_linewidth*2))
        ref_hkl_curve.setVisible(False)
        ref_hkl_curve.index = None
        self.patches['ref_hl_curve'] = ref_hkl_curve
        self.ax.addItem(ref_hkl_curve)

        # add hkl annotation label
        ref_hkl_label = pg.TextItem(anchor=(0.5,1.0), color=self.unit_label_color)
        ref_hkl_label.setFont(font)
        ref_hkl_label.setVisible(False)
        ref_hkl_label.setZValue(1) # place label on top of all other items
        self.patches['ref_hl_label'] = ref_hkl_label
        self.ax.addItem(ref_hkl_label)

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
        self.patches['poni'] = pg.ScatterPlotItem(symbol=self.plo.poni_marker,
                                                  size=self.plo.poni_size,
                                                  brush=pg.mkBrush(self.cont_cmap.map(0, mode='qcolor')),
                                                  pen=pg.mkPen(None))
        self.ax.addItem(self.patches['poni'])

        # add beam center scatter plot
        self.patches['beamcenter'] = pg.ScatterPlotItem(symbol=self.plo.beamcenter_marker,
                                                        size=self.plo.beamcenter_size,
                                                        brush=pg.mkBrush(self.cont_cmap.map(0, mode='qcolor')),
                                                        pen=pg.mkPen(None))
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
        self.resize_win()

        # add unit label
        self.label_unit_init()

        # add label for corrections
        self.label_corr_init()

        # add empty lines for azimuthal grid
        self.polar_grid_init()

        # create cones and draw contour lines
        self.update_screen()
        self.set_win_title()

    def polar_grid_init(self):
        """initialize a polar grid that can be easily hidden or shown"""
        # lines for azimuthal grid
        for i in range(self.plo.azimuth_num):
            line = pg.PlotCurveItem(useCache=True,
                                    pen=pg.mkPen(self.grid_color, width=1))
            self.ax.addItem(line)
            self.patches['polar_grid'].append(line)

    ############
    #  LABELS  #
    ############
    def label_unit_init(self):
        """
        Initializes and sets up the unit label for the plot.

        This method creates a text item for the unit label using the specified font,
        color, and fill properties. It then sets the text of the label to the name
        of the current unit and adds it to the plot axis at a specified position.

        Attributes:
            self.unit_label (pg.TextItem): The text item for the unit label.
            self.unit_label_color (QtGui.QColor): The color of the unit label text.
            self.unit_label_fill (QtGui.QBrush): The fill color of the unit label background.
            self.unit_names (dict): A dictionary mapping unit identifiers to their names.
            self.geo.unit (str): The current unit identifier.
            self.plo.unit_label_size (int): The pixel size of the unit label font.
            self.xdim (float): The x-coordinate for positioning the unit label.
            self.ydim (float): The y-coordinate for positioning the unit label.
            self.ax (pg.PlotItem): The plot axis to which the unit label is added.
        """
        font = QtGui.QFont()
        font.setPixelSize(self.plo.unit_label_size)
        self.unit_label = pg.TextItem(anchor=(0.0,0.0), color=self.unit_label_color, fill=self.unit_label_fill)
        self.unit_label.setText(self.unit_names[self.geo.unit])
        self.unit_label.setFont(font)
        self.ax.addItem(self.unit_label)
        self.unit_label.setPos(-self.xdim, self.ydim)

    def label_corr_init(self):
        """
        Initializes and configures the correction label for the plot.

        This method sets up a text label with specific font size, color, and fill properties.
        The label is initially hidden and positioned at the specified coordinates.
        A tooltip is also added to the label to provide additional information.

        Attributes:
            cor_label (pg.TextItem): The text item used as the correction label.
            plo.unit_label_size (int): The pixel size for the label font.
            unit_label_color (str or QColor): The color of the label text.
            unit_label_fill (str or QColor): The fill color of the label background.
            xdim (float): The x-coordinate for positioning the label.
            ydim (float): The y-coordinate for positioning the label.
        """
        font = QtGui.QFont()
        font.setPixelSize(self.plo.unit_label_size)
        self.cor_label = pg.TextItem(anchor=(1.0,1.0), color=self.unit_label_color, fill=self.unit_label_fill)
        self.cor_label.setText('')
        self.cor_label.setFont(font)
        self.ax.addItem(self.cor_label)
        self.cor_label.hide()
        self.cor_label.setToolTip('P: Polarisation\nS: Solid angle\nF: FWHM [\u00B0]')
        self.cor_label.setPos(self.xdim, -self.ydim)

    def label_conic_pos_auto(self, x, y):
        """
        Adjusts the label position to maintain readability and ensures the label
        stays within visible boundaries of the detector dimensions. This method 
        tries to place labels as close as possible to the horizontal center.

        Parameters:
        x (numpy.ndarray): Array of x coordinates of the conic section.
        y (numpy.ndarray): Array of y coordinates of the conic section.

        Returns:
        list: A list containing the x and y coordinates of the label position 
              if the conic section is visible, otherwise returns False.
        """
        # This tries to place labels as close as possible
        #  to the horizontal center.
        # Adjust the label position to maintain readibility
        #  XX% of detector dimensions make sure the label 
        #  stays within visible boundaries.
        # Available label positions depend on the stepsize
        #  used to calculate the conic section. E.g. finer
        #  stepping makes labels follow the conics smoother.
        _condition = np.nonzero((x >= -self.xdim*0.96)&
                                (x <=  self.xdim*0.96)&
                                (y >= -self.ydim*0.98)&
                                (y <=  self.ydim*0.98))
        # regions in x and y matching the condition
        _visible_x = x[_condition]
        _visible_y = y[_condition]
        # check if conic is visible
        if _visible_x.size == 0 and _visible_y.size == 0:
            # outside visible area
            return False
        # if conic is not cut-off, use horizontal offset as x value.
        # this makes the display is cleaner, as the positions otherwise
        # rely on the sampling of the conic e.g. rough stepping.
        if _visible_y.max() == y.max():
            label_x = self.geo.hoff
            label_y = _visible_y.max()
        elif _visible_y.min() == y.min():
            label_x = self.geo.hoff
            label_y = _visible_y.min()
        else:
            # place label as close as possible to the horizontal center
            # get the two smallest indices (absolute -> positive and negative).
            # we don't know if the positive or the negative value is closer to
            # the target so we get all possible four (2x, 2y) and find the best
            # combination in the next step.
            _sorted = abs(_visible_x - self.geo.hoff).argsort()[:4]
            # get the index of the largest sum (-> upper right)
            _visible = np.argmax(_visible_x[_sorted]+_visible_y[_sorted])
            # get the index of the target value
            _idx = _sorted[_visible]
            # depending on the curvature use either maximum or
            # minimum value as label y position.
            label_x = _visible_x[_idx]
            label_y = _visible_y[_idx]
        return [label_x, label_y]
        
    def label_conic_pos_static(self, x, y, xdim, ydim, omega, theta):
        """
        Determine the label position for a conic section on a detector.

        This method checks if the conic section defined by the coordinates (x, y) is within the detector area
        defined by dimensions (xdim, ydim). If the conic section is within the detector area, it calculates
        an appropriate label position based on the given angles omega and theta to maintain readability.

        Parameters:
        x (numpy.ndarray): Array of x-coordinates of the conic section.
        y (numpy.ndarray): Array of y-coordinates of the conic section.
        xdim (float): Half-width of the detector area.
        ydim (float): Half-height of the detector area.
        omega (float): Angle in degrees, used to determine label position.
        theta (float): Angle in radians, used to determine label position.

        Returns:
        list or bool: A list containing the label position [hoff, y] if the conic section is within the detector area,
                  otherwise False.
        """
        # check if conic is visible
        cx = np.argwhere((x >= -xdim) & (x <= xdim))
        cy = np.argwhere((y >= -ydim) & (y <= ydim))
        if len(cx) == 0 or len(cy) == 0:
            # outside detector area
            return False
        # adjust the label position to maintain readibility
        # this works for most cases but is not the most optimal solution yet
        # OR: use the actual beam position to determine label position
        # beam_pos_y = -(self.geo.voff + np.tan(np.deg2rad(self.geo.rota))*self.geo.dist)
        if omega <= 0:
            label_pos = [self.geo.hoff, max(y)] if theta < np.pi/2 else [self.geo.hoff, min(y)]
        else:
            label_pos = [self.geo.hoff, min(y)] if theta < np.pi/2 else [self.geo.hoff, max(y)]
        return label_pos

    def label_set_position(self, pos):
        if pos.lower() in ['u', 'up', 'upper', 't', 'top']:
            self.unit_label.setPos(-self.xdim, self.ydim)
            self.unit_label.setAnchor((0.0, 0.0))
        elif pos.lower() in ['d', 'down', 'b', 'bottom']:
            self.unit_label.setPos(-self.xdim, -self.ydim)
            self.unit_label.setAnchor((0.0, 1.0))
        else:
            return

    ###########
    #  THEME  #
    ###########
    def theme_apply(self, use_dark, redraw=False):
        """
        Apply the theme to the application.
        This method sets the application's theme to either dark mode or light mode 
        based on the `use_dark` parameter. It updates various UI elements such as 
        icons, colors, and palettes to match the selected theme. Optionally, it can 
        also redraw the canvas if the `redraw` parameter is set to True.
        Parameters:
        -----------
        use_dark : bool
            If True, apply the dark theme. If False, apply the light theme.
        redraw : bool, optional
            If True, redraw the canvas after applying the theme. Default is False.
        Returns:
        --------
        None
        """
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
            # icon
            self.pixmap = QtGui.QPixmap(':/icons/xrdPlanner_dark')
            self.icon = QtGui.QIcon(':/icons/xrdPlanner_dark')
            # reference contour
            self.conic_label_fill = QtGui.QColor(self.thm.dark_conic_label_fill)
            self.conic_ref_color = QtGui.QColor(self.thm.dark_conic_ref_color)
            self.conic_highlight = QtGui.QColor(self.thm.dark_conic_highlight)
            self.det_module_color = QtGui.QColor(self.thm.dark_det_module_color)
            self.det_module_fill = QtGui.QColor(self.thm.dark_det_module_fill)
            # general
            self.plot_bg_color = QtGui.QColor(self.thm.dark_plot_bg_color)
            self.beamstop_color = QtGui.QColor(self.thm.dark_beamstop_color)
            self.beamstop_edge_color = QtGui.QColor(self.thm.dark_beamstop_edge_color)
            self.unit_label_color = QtGui.QColor(self.thm.dark_unit_label_color)
            self.unit_label_fill = QtGui.QColor(self.thm.dark_unit_label_fill)
            self.overlay_threshold_color = QtGui.QColor(self.thm.dark_overlay_threshold_color)
            self.grid_color = QtGui.QColor(self.thm.dark_grid_color)
            # slider
            self.slider_border_color = QtGui.QColor(self.thm.dark_slider_border_color)
            self.slider_bg_color = QtGui.QColor(self.thm.dark_slider_bg_color)
            self.slider_bg_hover = QtGui.QColor(self.thm.dark_slider_bg_hover)
            self.slider_label_color = QtGui.QColor(self.thm.dark_slider_label_color)
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
            # icon
            self.pixmap = QtGui.QPixmap(':/icons/xrdPlanner')
            self.icon = QtGui.QIcon(':/icons/xrdPlanner')
            # reference contour
            self.conic_label_fill = QtGui.QColor(self.thm.light_conic_label_fill)
            self.conic_ref_color = QtGui.QColor(self.thm.light_conic_ref_color)
            self.conic_highlight = QtGui.QColor(self.thm.light_conic_highlight)
            self.det_module_color = QtGui.QColor(self.thm.light_det_module_color)
            self.det_module_fill = QtGui.QColor(self.thm.light_det_module_fill)
            # general
            self.plot_bg_color = QtGui.QColor(self.thm.light_plot_bg_color)
            self.beamstop_color = QtGui.QColor(self.thm.light_beamstop_color)
            self.beamstop_edge_color = QtGui.QColor(self.thm.light_beamstop_edge_color)
            self.unit_label_color = QtGui.QColor(self.thm.light_unit_label_color)
            self.unit_label_fill = QtGui.QColor(self.thm.light_unit_label_fill)
            self.overlay_threshold_color = QtGui.QColor(self.thm.light_overlay_threshold_color)
            self.grid_color = QtGui.QColor(self.thm.light_grid_color)
            # slider
            self.slider_border_color = QtGui.QColor(self.thm.light_slider_border_color)
            self.slider_bg_color = QtGui.QColor(self.thm.light_slider_bg_color)
            self.slider_bg_hover = QtGui.QColor(self.thm.light_slider_bg_hover)
            self.slider_label_color = QtGui.QColor(self.thm.light_slider_label_color)
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

        # change palette for all child windows
        self.change_palette_recursive(self, app.palette())

        # redraw the canvas
        if redraw:
            self.redraw_canvas()
            self.about_win.update_logo(self.pixmap)
    
    def change_palette_recursive(self, root, palette):
        """
        Recursively changes the palette of a given root widget and all its child widgets.

        Args:
            root (QtWidgets.QWidget): The root widget whose palette will be changed.
            palette (QtGui.QPalette): The new palette to be applied to the root widget and its children.

        Returns:
            None
        """
        root.setPalette(palette)
        for child in root.children():
            if isinstance(child, QtWidgets.QWidget):
                self.change_palette_recursive(child, palette)
    
    ##########
    #  MENU  #
    ##########

    def menu_init(self, reset=False):
        """
        Initializes the menu bar for the application.
        Parameters:
        reset (bool): If True, clears the existing menu bar before initializing.
        This method sets up various menus and submenus including:
        - Detector: Allows selection of different detectors and their sizes.
        - Reference: Allows selection of reference types including 'None', 'From pyFAI', 'From cif', and 'From cell'.
        - Beamstop: Allows selection of beamstop sizes and 'None'.
        - Unit: Allows selection of different units.
        - View: Contains submenus for Theme, Colormap, Overlay, and Functions.
        - Settings: Contains options to load, save, import, export settings, and edit files.
        - Help: Provides information about the application, geometry conventions, and hotkeys.
        Each menu item is associated with specific actions that are triggered when the item is selected.
        """

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
        # settings reload via self.settings_reload()
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
                self.menu_set_action(det_action, self.change_detector, d, s)
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
        self.menu_set_action(ref_action, self.change_reference, 'None')
        menu_ref.addAction(ref_action)
        self.group_ref.addAction(ref_action)
        if self.geo.reference.lower() == 'none':
            ref_action.setChecked(True)
        # menu Reference: add pyFAI library
        sub_menu_pyFAI = QtWidgets.QMenu('From pyFAI', self)
        menu_ref.addMenu(sub_menu_pyFAI)
        for ref_name in sorted(self.ref_pyfai):
            ref_action = QtGui.QAction(ref_name, self, checkable=True)
            self.menu_set_action(ref_action, self.change_reference, ref_name)
            sub_menu_pyFAI.addAction(ref_action)
            self.group_ref.addAction(ref_action)
            if ref_name == self.geo.reference:
                ref_action.setChecked(True)

        # menu Reference: add cif
        self.sub_menu_cif = QtWidgets.QMenu('From cif', self)
        menu_ref.addMenu(self.sub_menu_cif)
        if self.ref_cif:
            for name in self.ref_cif.keys():
                self.add_cif_to_menu(name)
        self.sub_menu_cif.setDisabled(self.sub_menu_cif.isEmpty())

        # menu Reference: add cell
        self.sub_menu_cell = QtWidgets.QMenu('From cell', self)
        menu_ref.addMenu(self.sub_menu_cell)
        if self.ref_cell:
            for name in self.ref_cell.keys():
                ref_action = QtGui.QAction(name, self, checkable=True)
                self.menu_set_action(ref_action, self.change_reference, name)
                self.sub_menu_cell.addAction(ref_action)
                self.group_ref.addAction(ref_action)
                if name == self.geo.reference:
                    ref_action.setChecked(True)
        self.sub_menu_cell.setDisabled(self.sub_menu_cell.isEmpty())

        menu_ref.addSeparator()
        # menu Reference: add None
        cell_action = QtGui.QAction('Calculate from Cell', self)
        cell_action.triggered.connect(self.uc_win.show)
        menu_ref.addAction(cell_action)
        
        ############
        # BEAMSTOP #
        ############
        # menu Beamstop
        menu_bs = self.menu_bar.addMenu('Beamstop')
        group_bs = QtGui.QActionGroup(self)
        group_bs.setExclusive(True)
        # menu Beamstop: add None
        bs_action = QtGui.QAction('None', self, checkable=True)
        self.menu_set_action(bs_action, self.change_beamstop, 'None')
        menu_bs.addAction(bs_action)
        group_bs.addAction(bs_action)
        if isinstance(self.geo.bssz, str) and self.geo.bssz.lower() == 'none':
            bs_action.setChecked(True)
        # menu Beamstop: add sizes list
        sub_menu_bs = QtWidgets.QMenu('Sizes [mm]', self)
        menu_bs.addMenu(sub_menu_bs)
        for bs_size in sorted(self.geo.bs_list):
            bs_sub_action = QtGui.QAction(str(bs_size), self, checkable=True)
            self.menu_set_action(bs_sub_action, self.change_beamstop, bs_size)
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
            self.menu_set_action(unit_action, self.change_units, unit_index)
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
            self.menu_set_action(theme_action, self.theme_apply, invert, True)
            group_theme.addAction(theme_action)
            menu_theme.addAction(theme_action)
            if invert == self.geo.darkmode:
                theme_action.setChecked(True)

        # Colored reference
        self.action_colored_reference = QtGui.QAction('Colored reference', self, checkable=True)
        self.menu_set_action(self.action_colored_reference, self.toggle_colored_reference)
        if self.plo.colored_reference:
            self.action_colored_reference.setChecked(True)
        else:
            self.action_colored_reference.setChecked(False)
        menu_view.addAction(self.action_colored_reference)

        ###################
        # VIEW - COLORMAP #
        ###################
        # indicate the hot key by underlining the letter with '&'
        self.menu_cmap = menu_view.addMenu('&Colormap')
        group_cmap = QtGui.QActionGroup(self)
        group_cmap.setExclusive(True)
        for cmap_name in self.colormaps: # PyQtGraph.colormap.listMaps(): Experimental, subject to change.
            cmap_action = QtGui.QAction(cmap_name, self, checkable=True)
            self.menu_set_action(cmap_action, self.change_cmap, cmap_name)
            self.menu_cmap.addAction(cmap_action)
            group_cmap.addAction(cmap_action)
            if cmap_name == self.geo.colormap:
                cmap_action.setChecked(True)

        ##################
        # VIEW - OVERLAY #
        ##################
        menu_overlays = menu_view.addMenu('Overlay')
        # unit value hover toggle
        self.action_unit_hover = QtGui.QAction('&Unit hover', self, checkable=True)
        self.menu_set_action(self.action_unit_hover, self.toggle_unit_hover)
        if self.plo.show_unit_hover:
            self.action_unit_hover.setChecked(True)
        else:
            self.action_unit_hover.setChecked(False)
        menu_overlays.addAction(self.action_unit_hover)
        # azimuthal grid toggle
        self.action_grid = QtGui.QAction('&Grid', self, checkable=True)
        self.menu_set_action(self.action_grid, self.toggle_grid)
        if self.plo.show_grid:
            self.action_grid.setChecked(True)
        else:
            self.action_grid.setChecked(False)
        menu_overlays.addAction(self.action_grid)
        # separator
        menu_overlays.addSeparator()
        # polarisation map toggle
        self.action_show_pol = QtGui.QAction('&Polarisation', self, checkable=True)
        self.menu_set_action(self.action_show_pol, self.toggle_overlay_polarisation)
        if self.plo.show_polarisation:
            self.action_show_pol.setChecked(True)
        else:
            self.action_show_pol.setChecked(False)
        menu_overlays.addAction(self.action_show_pol)
        # Solidangle map toggle
        self.action_show_ang = QtGui.QAction('Solid &angle', self, checkable=True)
        self.menu_set_action(self.action_show_ang, self.toggle_overlay_solidangle)
        if self.plo.show_solidangle:
            self.action_show_ang.setChecked(True)
        else:
            self.action_show_ang.setChecked(False)
        menu_overlays.addAction(self.action_show_ang)
        # Overlay warn color toggle
        menu_overlays.addSeparator()
        self.action_overlay_warn = QtGui.QAction('&Highlight', self, checkable=True)
        self.menu_set_action(self.action_overlay_warn, self.toggle_overlay_highlight)
        if self.plo.overlay_toggle_warn:
            self.action_overlay_warn.setChecked(True)
        else:
            self.action_overlay_warn.setChecked(False)
        menu_overlays.addAction(self.action_overlay_warn)

        ####################
        # VIEW - FUNCTIONS #
        ####################
        menu_functions = menu_view.addMenu('Functions')
        #set fwhm parameters toggle
        self.action_funct_fwhm_set = QtGui.QAction('Setup &FWHM', self)
        self.menu_set_action(self.action_funct_fwhm_set, self.fwhm_win.show)
        menu_functions.addAction(self.action_funct_fwhm_set)
        #show fwhm toggle
        self.action_funct_fwhm_show = QtGui.QAction('Show FWHM', self, checkable=True)
        self.menu_set_action(self.action_funct_fwhm_show, self.toggle_fwhm)
        if self.plo.show_fwhm:
           self.action_funct_fwhm_show.setChecked(True)
        else:
           self.action_funct_fwhm_show.setChecked(False)
        menu_functions.addAction(self.action_funct_fwhm_show)
        #export fwhm toggle
        self.action_funct_fwhm_export = QtGui.QAction('Export FWHM', self)
        self.menu_set_action(self.action_funct_fwhm_export, self.fwhm_win.export_grid)
        menu_functions.addAction(self.action_funct_fwhm_export)
        
        # PXRD pattern
        self.action_pxrd_pattern = QtGui.QAction('P&XRD pattern', self)
        self.menu_set_action(self.action_pxrd_pattern, self.win_pxrd_plot)
        menu_view.addAction(self.action_pxrd_pattern)

        # Absorption window
        #self.action_abs_win = QtGui.QAction('A&bsorption', self)
        #self.menu_set_action(self.action_abs_win, self.abs_win.show)
        #menu_view.addAction(self.action_abs_win)

        ############
        # SETTINGS #
        ############
        # menu Settings
        menu_Settings = self.menu_bar.addMenu('Settings')
        # submenu load settings files
        self.menu_custom_settings = menu_Settings.addMenu('Load')
        self.group_cset = QtGui.QActionGroup(self)
        self.group_cset.setExclusive(True)
        for cset_name in self.settings_get_files():
            cset_action = QtGui.QAction(cset_name, self, checkable=True)
            self.menu_set_action(cset_action, self.settings_change_file, cset_name)
            self.menu_custom_settings.addAction(cset_action)
            self.group_cset.addAction(cset_action)
            if cset_name == self.active_settings:
                cset_action.setChecked(True)
        
        ###################
        # SETTINGS - SAVE #
        ###################
        save_action = QtGui.QAction('Save', self)
        self.menu_set_action(save_action, self.settings_save_current)
        menu_Settings.addAction(save_action)
        #####################
        # SETTINGS - IMPORT #
        #####################
        import_action = QtGui.QAction('Import', self)
        self.menu_set_action(import_action, self.settings_import_win)
        menu_Settings.addAction(import_action)
        #####################
        # SETTINGS - EXPORT #
        #####################
        export_action = QtGui.QAction('Export editor', self)
        self.menu_set_action(export_action, self.export_win.show)
        menu_Settings.addAction(export_action)
        
        menu_Settings.addSeparator()
        #####################
        # SETTINGS - DETDB  #
        #####################
        detdb_action = QtGui.QAction('Detector db editor', self)
        #self.menu_set_action(detdb_action, self.win_detdb_show)
        self.menu_set_action(detdb_action, self.detdb_win.show)
        menu_Settings.addAction(detdb_action)
        ###################
        # SETTINGS - EDIT #
        ###################
        menu_Edit = menu_Settings.addMenu('Edit files')
        ##########################################
        # prepare platform dependent file reader #
        ##########################################
        if sys.platform == 'win32':
            tokens = [('Detector db', os.system, f'notepad {self.path_detdb}'),
                      ('Settings', self.settings_edit_file, f'notepad')]
        elif sys.platform == 'linux':
            tokens = [('Detector db', os.system, f'xdg-open {self.path_detdb}'),
                      ('Settings', self.settings_edit_file, f'xdg-open')]
        else:
            tokens = [('Detector db', os.system, f'open -t {self.path_detdb}'),
                      ('Settings', self.settings_edit_file, f'open -t')]
        for (name, funct, command) in tokens:
            edit_action = QtGui.QAction(name, self)
            self.menu_set_action(edit_action, funct, command)
            menu_Edit.addAction(edit_action)
        #####################
        # SETTINGS - RESET  #
        #####################
        reset_ddb_action = QtGui.QAction('Reset detector db', self)
        self.menu_set_action(reset_ddb_action, self.reset_detector_db)
        menu_Settings.addAction(reset_ddb_action)
        #####################
        # SETTINGS - DELETE #
        #####################
        menu_Settings.addSeparator()
        export_action = QtGui.QAction('Delete', self)
        self.menu_set_action(export_action, self.settings_delete_file)
        menu_Settings.addAction(export_action)

        # menu About
        menu_help = self.menu_bar.addMenu('Help')
        action_about = QtGui.QAction('xrdPlanner', self)
        #self.menu_set_action(action_about, self.show_about_win)
        self.menu_set_action(action_about, self.about_win.show)
        menu_help.addAction(action_about)
        action_geometry = QtGui.QAction('Geometry conventions', self)
        self.menu_set_action(action_geometry, self.geometry_win.show)
        menu_help.addAction(action_geometry)
        action_hotkeys = QtGui.QAction('Hotkeys', self)
        self.menu_set_action(action_hotkeys, self.hotkeys_win.show)
        menu_help.addAction(action_hotkeys)

    def menu_set_action(self, action, target, *args):
        """
        Connects a given action's triggered signal to a target function with optional arguments.

        Parameters:
        action (QAction): The action whose triggered signal will be connected.
        target (callable): The function to be called when the action is triggered.
        *args: Additional arguments to be passed to the target function when called.
        """
        action.triggered.connect(lambda: target(*args))

    ############
    #  TOGGLE  #
    ############
    def toggle_unit_hover(self):
        """
        Toggles the visibility of unit hover in the plot.

        This method switches the state of `show_unit_hover` in the plot object (`plo`).
        If `show_unit_hover` is enabled, it sets the corresponding action (`action_unit_hover`)
        to checked; otherwise, it unchecks it. Finally, it triggers a redraw of the canvas.
        """
        self.plo.show_unit_hover = not self.plo.show_unit_hover
        if self.plo.show_unit_hover:
            self.action_unit_hover.setChecked(True)
        else:
            self.action_unit_hover.setChecked(False)
        self.redraw_canvas()

    def toggle_overlay_polarisation(self):
        """
        Toggles the polarisation overlay on or off.

        This method switches the state of the polarisation overlay by toggling
        the `show_polarisation` attribute of the `plo` object. It also updates
        the checked state of the `action_show_pol` action to reflect the new
        state and triggers a redraw of the canvas.
        """
        self.plo.show_polarisation = not self.plo.show_polarisation
        if self.plo.show_polarisation:
            self.action_show_pol.setChecked(True)
        else:
            self.action_show_pol.setChecked(False)
        self.redraw_canvas()

    def toggle_overlay_solidangle(self):
        """
        Toggles the visibility of the solid angle overlay on the plot.

        This method switches the state of the `show_solidangle` attribute of the `plo` object.
        It also updates the checked state of the `action_show_ang` action accordingly and 
        triggers a redraw of the canvas to reflect the changes.
        """
        self.plo.show_solidangle = not self.plo.show_solidangle
        if self.plo.show_solidangle:
            self.action_show_ang.setChecked(True)
        else:
            self.action_show_ang.setChecked(False)
        self.redraw_canvas()

    def toggle_overlay_highlight(self):
        """
        Toggles the overlay highlight warning state.

        This method switches the state of the overlay highlight warning between
        enabled and disabled. It updates the corresponding action's checked state
        and triggers a canvas redraw to reflect the change.
        """
        self.plo.overlay_toggle_warn = not self.plo.overlay_toggle_warn
        if self.plo.overlay_toggle_warn:
            self.action_overlay_warn.setChecked(True)
        else:
            self.action_overlay_warn.setChecked(False)
        self.redraw_canvas()

    def toggle_fwhm(self):
        """
        Toggles the visibility of the Full Width at Half Maximum (FWHM) on the plot.

        This method checks if the FWHM action is enabled. If it is, it toggles the 
        `show_fwhm` attribute of the plot object (`plo`). It also updates the 
        checked state of the FWHM action in the UI and triggers a redraw of the canvas.
        """
        if not self.action_funct_fwhm_show.isEnabled():
            return
        self.plo.show_fwhm = not self.plo.show_fwhm
        if self.plo.show_fwhm:
            self.action_funct_fwhm_show.setChecked(True)
        else:
            self.action_funct_fwhm_show.setChecked(False)
        self.redraw_canvas()

    def toggle_colored_reference(self):
        self.plo.colored_reference = not self.plo.colored_reference
        if self.plo.colored_reference:
            self.action_colored_reference.setChecked(True)
        else:
            self.action_colored_reference.setChecked(False)
        self.redraw_canvas()

    def toggle_grid(self):
        self.plo.show_grid = not self.plo.show_grid
        if self.plo.show_grid:
            self.action_grid.setChecked(True)
        else:
            self.action_grid.setChecked(False)
        self.redraw_canvas()

    ############
    #   REDO   #
    ############
    def redraw_canvas(self):
        """
        Redraws the canvas by performing the following steps:
        1. Initializes modifiable elements.
        2. Clears the current screen.
        3. Re-initializes the main screen.
        4. Applies styles to the slider widget.
        5. Initializes sliders in the slider widget.
        6. Centers the slider frame.
        """
        # save darkmode toggle
        #self.settings_save_to_file()
        self.modifiables_init()
        # clear the screen
        self.ax.clear()
        # re-initialise
        self.main_screen_init()
        # center the slider frame
        self.sliderWidget.apply_style()
        self.sliderWidget.init_sliders()
        self.sliderWidget.center_frame()

    def reset_detector_db(self):
        """
        Resets the detector database and updates the current detector settings.

        This method performs the following actions:
        1. Resets the `detector_db` by calling `get_det_library` with `update=False` and `reset=True`.
        2. Checks if the current detector type (`self.geo.det_type`) is in the updated `detector_db`.
           - If not, sets `self.geo.det_type` to the first available detector type.
        3. Checks if the current detector size (`self.geo.det_size`) is available for the current detector type.
           - If not, sets `self.geo.det_size` to the first available size for the current detector type.
        4. Calls `change_detector` to apply the changes.

        """
        self.detector_db = self.get_det_library(update=False, reset=True)
        # if current detector is not in the available list of detectors
        # pick the first valid one to continue
        if self.geo.det_type not in self.detector_db:
            self.geo.det_type = next(iter(self.detector_db.keys()))
        # if the current size is not available pick a valid size
        if self.geo.det_size not in self.detector_db[self.geo.det_type]['size']:
            self.geo.det_size = next(iter(self.detector_db[self.geo.det_type]['size'].keys()))
        self.change_detector()

    #############
    #  UPDATE   #
    #############
    def update_screen(self, val=None):
        """
        Updates the screen based on the provided value and the sender's object name.

        Args:
            val (float, optional): The new value to update. Defaults to None.

        Updates:
            - Updates the geometry attributes (`dist`, `rota`, `tilt`, `voff`, `hoff`, `ener`, `bsdx`)
              based on the sender's object name.
            - Adjusts the beamstop slider limits if the sender is 'dist' and the slider is enabled.
            - Re-calculates and re-draws cones and contours.
            - Updates child windows.
            - Draws reference contours if a reference is set.
        """
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
                #if self.abs_win.isVisible():
                #    self.abs_win.update_keV(self.geo.ener)
            elif self.sender().objectName() == 'bsdx':
                self.geo.bsdx = float(val)

        # re-calculate cones and re-draw contours
        self.draw_conics()
        # update child windows
        self.update_win_generic()
        # draw reference contours
        if self.geo.reference != 'None':
            self.get_reference()
            self.draw_reference()

    def update_win_generic(self):
        """
        Updates the window with generic settings.

        This method performs updates to the PXRD window and the FWHM window
        by calling the respective update methods.
        """
        if hasattr(self, 'pxrd_win') and self.pxrd_win is not None and self.pxrd_win.isVisible():
            self.win_pxrd_update()
        if hasattr(self, 'fwhm_win'):
            self.fwhm_win.update()

    ##################
    #  DRAW CONICS   #
    ##################
    def draw_conics(self):
        """
        Draws conic sections on the grid based on the current geometry and plotting settings.
        This method calculates the offset of the contours resulting from vertical offset (voff) 
        and rotation, shifts the grid to draw the cones, and ensures the contours are drawn 
        within the visible area. It also handles the overlay of additional information such as 
        polarisation, solid angle, unit hover, and full width at half maximum (FWHM).
        The method performs the following steps:
        1. Converts theta from degrees to radians.
        2. Calculates the beam center shift.
        3. Handles the overlay of additional information if required.
        4. Updates the beam center and PONI (Point of Normal Incidence).
        5. Hides the beamstop contour and label initially.
        6. Checks if the beamstop is on screen and updates its visibility and position.
        7. Calculates the maximum resolution for the given geometry.
        8. Plots conic sections at given 2theta values.
        The method uses several helper functions to calculate overlays, conic sections, 
        label positions, and units.
        Attributes:
            self.geo: Geometry settings including tilt, rotation, vertical offset, distance, 
                      horizontal offset, beamstop size, and beamstop distance.
            self.plo: Plotting settings including overlay resolution, polarisation factor, 
                      conic steps, conic label auto-placement, and conic line width.
            self.patches: Dictionary of plot elements including overlay, beam center, PONI, 
                          beamstop, and labels.
            self.cont_geom_num: Array of 2theta values for plotting conic sections.
            self.cont_cmap: Colormap for conic sections.
        """
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
        if self.plo.show_polarisation or self.plo.show_solidangle or self.plo.show_unit_hover or self.plo.show_fwhm:
            _grd, self._tth, self._azi, self._polcor, self._solang, self._fwhm = self.calc_overlays(_omega, res=self.plo.overlay_resolution, pol=self.plo.polarisation_fac)
            self.patches['overlay'].setImage(_grd * self._polcor * self._solang,
                                             autoLevels=False,
                                             levels=[0.0,1.0],
                                             rect=(-self.xdim,
                                                   -self.ydim,
                                                    self.xdim * 2,
                                                    self.ydim * 2))
        else:
            self._tth = None
            self._azi = None
            self.patches['overlay'].setImage(None)
            self.cor_label.hide()
        
        # update beam center
        self.patches['beamcenter'].setData([self.geo.hoff],[_comp_shift])
        # update beam center
        self.patches['poni'].setData([self.geo.hoff],[-(self.geo.voff - np.deg2rad(self.geo.tilt)*self.geo.dist)])

        # hide beamstop contour and label
        self.patches['beamstop'].setVisible(False)
        self.patches['bs_label'].setVisible(False)
        # check if beamstop is on screen, change visibility
        if self.geo.bssz and not (isinstance(self.geo.bssz, str) and self.geo.bssz.lower() == 'none'):
            # make sure it is a float (might be a string from the export window!)
            self.geo.bssz = float(self.geo.bssz)
            # update beam stop
            self.bs_theta = np.tan((self.geo.bssz/2) / self.geo.bsdx)

            # calculate the conic section corresponding to the theta angle
            # :returns False is conic is outside of visiblee area
            x, y = self.calc_conic(_omega, self.bs_theta, steps=self.plo.conic_steps)
            if x is not False:
                # figure out the label positions
                if self.plo.conic_label_auto:
                    label_pos = self.label_conic_pos_auto(x, y)
                else:
                    label_pos = self.label_conic_pos_static(x, y, self.xdim, self.ydim, _omega, self.bs_theta)
                
                # continue if label can be placed
                if label_pos is not False:
                    # plot the conic section
                    self.patches['beamstop'].setData(x, y, fillLevel=y.max())
                    self.patches['beamstop'].setVisible(True)

                    _unit = self.calc_unit(self.bs_theta)
                    self.patches['bs_label'].setPos(*label_pos)
                    self.patches['bs_label'].setText(f'{_unit:.2f}')
                    self.patches['bs_label'].setVisible(True)
        
        # calculate the maximum resolution for the given geometry
        if self.plo.conic_tth_auto:
            theta_max = np.rad2deg(self.calc_tth_max(scale=0.90))
            # make new array of 2theta values
            self.cont_geom_num = np.linspace(theta_max/self.plo.conic_tth_num, theta_max, self.plo.conic_tth_num)
        
        # plot conics at given 2theta values
        for _n, _ttd in enumerate(self.cont_geom_num):
            self.patches['conic'][_n].setVisible(False)
            self.patches['labels'][_n].setVisible(False)
            # current fraction for colormap
            _f = _n/len(self.cont_geom_num)

            # convert theta in degrees to radians
            theta = np.deg2rad(_ttd)

            # calculate the conic section corresponding to the theta angle
            # :returns False is conic is outside of visiblee area
            x, y = self.calc_conic(_omega, theta, steps=self.plo.conic_steps)
            if x is False or x.size == 0:
                continue

            # figure out the label positions
            if self.plo.conic_label_auto:
                label_pos = self.label_conic_pos_auto(x, y)
            else:
                label_pos = self.label_conic_pos_static(x, y, self.xdim, self.ydim, _omega, theta)
            # continue if label can be placed
            if label_pos is False:
                continue

            # plot the conic section
            _unit = self.calc_unit(theta)
            self.patches['labels'][_n].setPos(*label_pos)
            if self.plo.colored_reference:
                self.patches['conic'][_n].setData(x, y, pen=pg.mkPen(self.conic_ref_color, width=self.plo.conic_linewidth))
                self.patches['labels'][_n].setText(f'{_unit:.2f}', color=self.conic_ref_color)
            else:
                self.patches['conic'][_n].setData(x, y, pen=pg.mkPen(self.cont_cmap.map(_f, mode='qcolor'), width=self.plo.conic_linewidth))
                self.patches['labels'][_n].setText(f'{_unit:.2f}', color=self.cont_cmap.map(_f, mode='qcolor'))
            self.patches['conic'][_n].setVisible(True)
            self.patches['labels'][_n].setVisible(True)
            
            ## plot azimuthal grid lines
            if self.plo.show_grid:
                # calculate the azimuthal grid points
                grid_vectors = self.calc_azi_grid(_omega)
                for i,[a, v] in enumerate(grid_vectors.items()):
                    # plot the azimuthal grid point
                    # Currently the angle a is not used
                    # but might be useful for future reference
                    self.patches['polar_grid'][i].setData(v[:,0],
                                                        v[:,1],)
                    self.patches['polar_grid'][i].setVisible(True)
            else:
                for i in range(self.plo.azimuth_num):
                    self.patches['polar_grid'][i].setVisible(False)

    def draw_reference(self):
        """
        Draws reference contour lines on the plot.
        This method plots standard contour lines based on the number of reference 
        conic sections specified. It handles the visibility and data of each contour 
        line, calculates the necessary angles and positions, and updates the plot 
        accordingly.
        The method performs the following steps:
        1. Iterates over the number of reference conic sections.
        2. Sets the visibility of each contour line to False initially.
        3. Checks if the current index is within the range of available d-spacings.
        4. Calculates the wavelength and theta angle for each d-spacing.
        5. Converts theta from degrees to radians and adjusts for tilt and rotation.
        6. Calculates the conic section corresponding to the theta angle.
        7. Updates the contour line data and visibility if the conic section is within the visible area.
        8. If available, updates the contour line with hkl values and intensity information.
        9. Plots the conic section with the appropriate line width and color.
        Notes:
        - The method assumes that `self.geo.ener` is the energy value used for wavelength calculation.
        - The method uses `np.arcsin` for angle calculation and `np.deg2rad` for degree to radian conversion.
        - The method relies on `self.calc_conic` to compute the conic section coordinates.
        - The method uses `pg.mkPen` for setting the pen properties of the contour lines.
        Attributes:
        - self.plo.conic_ref_num: Number of reference conic sections.
        - self.cont_ref_dsp: List of d-spacings for the reference contours.
        - self.geo.ener: Energy value for wavelength calculation.
        - self.geo.tilt: Tilt angle for the geometry.
        - self.geo.rota: Rotation angle for the geometry.
        - self.plo.conic_steps: Number of steps for conic section calculation.
        - self.cont_ref_hkl: List of hkl values and intensities for the reference contours.
        - self.plo.conic_ref_cif_irel: Flag to show relative intensity.
        - self.plo.conic_ref_cif_lw_mult: Line width multiplier for relative intensity.
        - self.plo.conic_hkl_show_int: Flag to show intensity in the tooltip.
        - self.plo.conic_ref_cif_lw_min: Minimum line width for the reference contours.
        - self.plo.conic_ref_linewidth: Base line width for the reference contours.
        - self.conic_ref_color: Color for the reference contour lines.
        - self.patches['reference']: Dictionary to store the reference contour line patches.
        """
        if len(self.patches['reference']) == 0:
            return
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
                x, y = self.calc_conic(_omega, theta, steps=self.plo.conic_steps)
                if x is False:
                    continue

                # if hkl are available
                # put them in the proper container for the contour
                # so indexing gets it right
                width = 1.0
                if self.cont_ref_hkl is not None:
                    h, k, l, itot, irel = self.cont_ref_hkl[_n]
                    if self.plo.conic_ref_cif_irel and irel > 0:
                        width = irel * self.plo.conic_ref_cif_lw_mult
                    if self.plo.conic_hkl_show_int:
                        # alignment of the tooltip text is far from trivial
                        # detour via QTextEdit -> setAlignment and toHtml
                        # hence the <br> instead of \n
                        # 
                        # self.setStyleSheet('''QToolTip {... ) didn't work
                        self.patches['reference'][_n].name = f'{h: 2.0f} {k: 2.0f} {l: 2.0f}\n{round(itot, 0):,.0f}'
                        self.patches['reference'][_n].index = _n
                    else:
                        self.patches['reference'][_n].name = f'{h: 2.0f} {k: 2.0f} {l: 2.0f}'
                        self.patches['reference'][_n].index = _n
                else:
                    self.patches['reference'][_n].name = None
                    self.patches['reference'][_n].index = None
                
                # current fraction for colormap
                if self.plo.colored_reference:
                    # current fraction for colormap
                    _f = _n/len(self.cont_ref_dsp)
                    # plot the conic section
                    self.patches['reference'][_n].setData(x, y, pen=pg.mkPen(self.cont_cmap.map(_f, mode='qcolor'),
                                                                             width=max(self.plo.conic_ref_cif_lw_min, 
                                                                                       self.plo.conic_ref_linewidth * width)))
                else:
                    self.patches['reference'][_n].setData(x, y, pen=pg.mkPen(self.conic_ref_color,
                                                                             width=max(self.plo.conic_ref_cif_lw_min, 
                                                                                       self.plo.conic_ref_linewidth * width)))
                self.patches['reference'][_n].setVisible(True)
        
        # update highlighted curve position and hkl label
        _index_hl = self.patches['ref_hl_curve'].index
        if _index_hl is None:
            pass
        elif self.patches['reference'][_index_hl].isVisible():
            self.patches['reference'][_index_hl].highlight()
        else:
            self.patches['reference'][_index_hl].lowlight()
    
    ############
    #  CHANGE  #
    ############
    def change_beamstop(self, size):
        """
        Change the size of the beamstop and redraw the canvas.

        Parameters:
        size (float): The new size of the beamstop.
        """
        self.geo.bssz = size
        self.redraw_canvas()

    def change_cmap(self, token):
        """
        Change the colormap of the geo object based on the provided token.

        Parameters:
        token (str or int): If a string is provided, it sets the colormap directly to the string value.
                    If an integer is provided, it changes the colormap by the index offset.
                    Positive integers move forward in the colormap list, and negative integers move backward.

        Returns:
        None
        """
        # cmap name or index
        if isinstance(token, str):
            self.geo.colormap = token
        elif isinstance(token, int):
            idx = self.colormaps.index(self.geo.colormap) + token
            if idx >= len(self.colormaps):
                idx = 0
            elif idx < 0:
                idx = len(self.colormaps) - 1
            self.geo.colormap = self.colormaps[idx]
        else:
            return
        # change cmap
        for action in self.menu_cmap.actions():
            if action.text() == self.geo.colormap:
                action.setChecked(True)
        self.redraw_canvas()
        #self.update_win_generic()

    def change_detector(self, det_name=None, det_size=None):
        """
        Change the detector settings and update the display.

        This method allows changing the detector type and size. It then updates
        the detector parameters, clears the current display, re-initializes the 
        main screen, centers the slider frame, and updates the window.

        Args:
            det_name (str, optional): The name/type of the detector. Defaults to None.
            det_size (tuple, optional): The size of the detector. Defaults to None.
        """
        if det_name is not None:
            self.geo.det_type = det_name
        if det_size is not None:
            self.geo.det_size = det_size
        # get new detector specs
        self.det = self.get_det_params()
        # clear the screen
        self.ax.clear()
        # re-initialise
        self.main_screen_init()
        # center the slider frame
        self.sliderWidget.center_frame()
        #self.update_win_generic()

    def change_units(self, unit_index):
        """
        Change the units of measurement for the geometry and update the UI accordingly.

        Args:
            unit_index (int): The index of the new unit to be set.

        Updates:
            - Sets the geometry unit to the specified unit index.
            - Updates the unit label text to reflect the new unit.
            - Checks the corresponding action in the unit menu.
            - Redraws the conics.
            - Updates the generic window elements.
        """
        self.geo.unit = unit_index
        self.unit_label.setText(self.unit_names[unit_index])
        for num, action in enumerate(self.menu_unit.actions()):
            if num == self.geo.unit:
                action.setChecked(True)
        self.draw_conics()
        self.update_win_generic()

    def change_reference(self, ref_name):
        """
        Change the reference geometry to the specified reference name.

        This method updates the reference geometry of the object to the given
        reference name and subsequently calls methods to get, draw, and update
        the reference.

        Args:
            ref_name (str): The name of the new reference geometry.
        """
        self.geo.reference = ref_name
        self.get_reference()
        self.draw_reference()
        self.update_win_generic()

    #############
    #    GET    #
    #############
    def get_reference(self):
        """
        Retrieves and sets the reference d-spacings and HKL values based on the provided reference type.
        The method checks the type of reference (library, CIF, or cell) and retrieves the corresponding d-spacings 
        and HKL values. If the reference is not found, it sets the d-spacings to zero. Additionally, it handles 
        closing the PXRD window if the reference is not from a CIF file and updates the window title.
        Attributes:
            geo.reference (str): The reference identifier.
            ref_library (dict): A dictionary containing reference data from the library.
            ref_cif (dict): A dictionary containing reference data from CIF files.
            ref_cell (dict): A dictionary containing reference data from cell files.
            plo.conic_ref_num (int): The number of conic references.
            cont_ref_dsp (np.array): The array of d-spacings for the reference.
            cont_ref_hkl (list or None): The list of HKL values for the reference.
            pxrd_win (object or None): The PXRD window object.
        Raises:
            KeyError: If the reference is not found in any of the reference dictionaries.
        """
        # reset the Dans_Diffraction object
        self.xtl = None
        # check what type of reference is selected
        if self.geo.reference in self.ref_pyfai:
            # get the d spacings for the calibrtant from pyFAI
            self.cont_ref_dsp = np.array(calibrant.get_calibrant(self.geo.reference).get_dSpacing()[:self.plo.conic_ref_num])
            self.cont_ref_hkl = None
        elif self.geo.reference in self.ref_cif:
            if not self.ref_cif[self.geo.reference].has_hkl or not self.ref_cif[self.geo.reference].has_dsp:
                self.calc_ref_from_cif(self.ref_cif[self.geo.reference].cif, open_pxrd=True)
            else:
                # get custom d spacings
                self.cont_ref_dsp = self.ref_cif[self.geo.reference].dsp
                self.cont_ref_hkl = self.ref_cif[self.geo.reference].hkl
        elif self.geo.reference in self.ref_cell:
                self.cont_ref_dsp = self.ref_cell[self.geo.reference].dsp
                self.cont_ref_hkl = self.ref_cell[self.geo.reference].hkl
        else:
            # set all d-spacings to 0
            self.cont_ref_dsp = np.zeros(self.plo.conic_ref_num)
            self.cont_ref_hkl = None
        
        # close pxrd window if reference is not from cif
        if self.pxrd_win is not None and self.geo.reference not in self.ref_cif:
            if self.pxrd_win.isVisible():
                self.pxrd_win.close()
        
        # update window title
        self.set_win_title()

    def get_defaults_geo(self):
        """
        Sets up and returns the default geometry configuration for the detector.

        Returns:
            Container: An object containing the default geometry settings
        """
        ######################
        # Setup the geometry #
        ######################
        geo = Container()
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
        """
        Get default plot settings.

        Returns:
            Container: An object containing default plot settings.

        Default Plot Settings:
            - geometry contour section
            - reference contour section
            - module section
            - general section
            - pxrd plot
            - extra functions
            - slider section
            - update/reset
            - debug/testing
        """
        ################
        # Plot Details #
        ################
        plo = Container()
        # - geometry contour section - 
        plo.conic_tth_min = 5               # [int]    Minimum 2-theta contour line
        plo.conic_tth_max = 100             # [int]    Maximum 2-theta contour line
        plo.conic_tth_num = 12              # [int]    Number of contour lines
        plo.conic_tth_auto = True           # [bool]   Dynamic contour levels
        plo.beamcenter_marker = 'o'         # [marker] Beamcenter marker
        plo.beamcenter_size = 6             # [int]    Beamcenter size
        plo.poni_marker = 'x'               # [marker] Poni marker
        plo.poni_size = 8                   # [int]    Poni size
        plo.conic_linewidth = 2.0           # [float]  Contour linewidth (lw)
        plo.conic_label_size = 14           # [int]    Contour labelsize
        plo.conic_label_auto = True         # [bool]   Dynamic label positions
        plo.colored_reference = False       # [bool]   Color reference cones instead
        # - reference contour section - 
        plo.conic_ref_linewidth = 2.0       # [float]  Reference contour linewidth
        plo.conic_ref_timeout = 500         # [int]    highlight contour on click (msec)
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
        plo.plot_padding = 0                # [int]    Padding of the detector screen
        plo.unit_label_size = 16            # [int]    Label size, px
        plo.polarisation_fac = 0.99         # [float]  Horizontal polarisation factor
        plo.show_polarisation = True        # [bool]   Show polarisation overlay
        plo.show_solidangle = False         # [bool]   Show solid angle overlay
        plo.show_unit_hover = True          # [bool]   Show unit value on hover
        plo.show_grid = False               # [bool]   Show azimuthal grid
        plo.azimuth_num = 13                # [int]    Number of azimuthal grid lines
        plo.overlay_resolution = 300        # [int]    Overlay resolution
        plo.overlay_toggle_warn = True      # [bool]   Overlay warn color threshold
        # - pxrd plot -
        plo.pxrd_marker_symbol = 'arrow_up' # [marker] Symbol to mark peaks
        plo.pxrd_marker_offset = 0.05       # [float]  offset of marker from x-axis
        # - extra functions -
        plo.show_fwhm = False               # [bool]   Show delta_d/d function
        plo.sensor_thickness = 1000e-6      # [float]  Detector sensor thickness [m]
        plo.sensor_material = 'CdTe'        # [str]    Detector sensor material
        plo.beam_divergence = 10e-6         # [float]  X-ray beam divergence [rad]
        plo.scattering_diameter = 100e-6    # [float]  Scattering volume diameter [m]
        plo.energy_resolution = 1.4e-4      # [float]  X-ray beam resolution [eV/keV]
        #plo.funct_fwhm_thresh = 1e-3       # 
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
        plo.use_native_menubar = True       # [bool]   Use native menubar
        plo.set_debug = False               # [bool]   Debug mode

        return plo
    
    def get_defaults_thm(self):
        """
        Returns a Container object with default theme settings for both light and dark modes.

        The theme settings include various color configurations for different UI elements such as:
        - Global colors
        - Contour labels
        - Reference contours
        - Beamstop colors and edges
        - Detector module borders and backgrounds
        - Plot backgrounds
        - Label colors and fills
        - Slider frame borders, backgrounds, and hover colors
        - Slider frame label colors
        - Map threshold colors

        Returns:
            Container: An object containing the default theme settings.
        """
        #################
        # Theme Details #
        #################
        thm = Container()
        thm.color_dark = '#404040'                    # [color]  Global dark color
        thm.color_light = '#EEEEEE'                   # [color]  Global light color
        # light mode
        thm.light_conic_label_fill = '#FFFFFF'        # [color]  Contour label fill color
        thm.light_conic_ref_color = '#DCDCDC'         # [color]  Reference contour color
        thm.light_conic_highlight = '#FF0000'         # [color]  Reference contour highlight color
        thm.light_beamstop_color = '#80FF0000'        # [color]  Beamstop color
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
        thm.light_grid_color = '#AAAAAA'              # [color]  Grid color
        # dark mode
        thm.dark_conic_label_fill = '#000000'         # [color]  Contour label fill color
        thm.dark_conic_ref_color = '#303030'          # [color]  Reference contour color
        thm.dark_conic_highlight = '#FF0000'          # [color]  Reference contour highlight color
        thm.dark_beamstop_color = '#AAFF0000'         # [color]  Beamstop color
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
        thm.dark_grid_color = '#606060'               # [color]  Grid color

        return thm
    
    def get_defaults_lmt(self):
        """
        Get default limits for various parameters.
        Returns:
            Container: An object containing default limits
        """
        ##########
        # Limits #
        ##########
        lmt = Container()
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
        """
        Load and set all default configurations.

        This method initializes the following default configurations:
        - Geometry and detector specifications (`geo` and `_geo`)
        - Plot details (`plo`)
        - Theme details (`thm`)
        - Geometry limits (`lmt`)

        Each configuration is loaded using its respective `get_defaults_*` method.
        """
        # load the defaults
        # geo: geometry and detector specs
        self.geo = self.get_defaults_geo()
        self._geo = Container()
        # plo: plot details
        self.plo = self.get_defaults_plo()
        # thm: theme details
        self.thm = self.get_defaults_thm()
        # lmt: geometry limits
        self.lmt = self.get_defaults_lmt()

    def get_det_params(self):
        """
        Retrieves the detector parameters based on the detector type and size from the detector database.
        This method checks if the detector type and size exist in the detector database. If either the type or size
        is not found, it prints an error message and exits the program. If the type and size are found, it creates
        a Container object and populates it with the corresponding detector parameters from the database.
        Returns:
            Container: An object containing the detector parameters.
        Raises:
            SystemExit: If the detector type or size is not found in the database, or if there is an error reading
                the detector database.
        """
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
        try:
            det = Container()
            det.hmp = self.detector_db[det_type]['hmp']
            det.vmp = self.detector_db[det_type]['vmp']
            det.pxs = self.detector_db[det_type]['pxs']
            det.hgp = self.detector_db[det_type]['hgp']
            det.vgp = self.detector_db[det_type]['vgp']
            det.cbh = self.detector_db[det_type]['cbh']
            det.hmn, det.vmn = self.detector_db[det_type]['size'][det_size]
            det.name = f'{det_type} {det_size}'
        except:
            print(f'Error reading detector db -> resetting file.')
            print(f' - Please restart xrdPlanner.')
            self.get_det_library(update=False, reset=True)
            raise SystemExit

        return det

    def get_det_library(self, update=True, reset=False):
        """
        Retrieves the detector specifications library.
        This method constructs a dictionary containing specifications for various detectors.
        It can either update the existing detector database file or reset it based on the parameters provided.
        Parameters:
        update (bool): If True, updates the detector database file with the current specifications.
                       Default is True.
        reset (bool): If True, resets the detector database file to the current specifications.
                      Default is False.
        Returns:
        dict: A dictionary containing the specifications for various detectors.
        Detector Specifications:
        - PILATUS3
        - PILATUS4
        - EIGER2
        - MPCCD
        - RAYONIX
        - PHOTON-II
        - PHOTON-III
        - Perkin-Elmer XRD
        - VAREX XRpad2
        - CITIUS
        Each detector specification includes:
        - pxs: Pixel size in mm
        - hmp: Module size (horizontal) in pixels
        - vmp: Module size (vertical) in pixels
        - hgp: Module gap (horizontal) in pixels
        - vgp: Module gap (vertical) in pixels
        - cbh: Central beam hole in pixels
        - size: Dictionary mapping detector model to its dimensions (rows, columns)
        Raises:
        SystemExit: If there is an error parsing the detector database file.
        """
        ###########################
        # Detector Specifications #
        ###########################
        detectors = dict()
        ###############################
        # Specifications for Pilatus3 #
        ###############################
        detectors['PILATUS3'] = {
            'pxs' : 172e-3, # [mm] Pixel size
            'hmp' : 487,    # [px] Module size (horizontal
            'vmp' : 195,    # [px] Module size (vertical)
            'hgp' : 7,      # [px] Module gap (horizontal)
            'vgp' : 17,     # [px] Module gap (vertical)
            'cbh' : 0,      # [px] Central beam hole
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
            'pxs' : 150e-3, # [mm] Pixel size
            'hmp' : 513,    # [px] Module size (horizontal
            'vmp' : 255,    # [px] Module size (vertical)
            'hgp' : 7,      # [px] Module gap (horizontal)
            'vgp' : 20,     # [px] Module gap (vertical)
            'cbh' : 0,      # [px] Central beam hole
            'size' : {'1M':(2,4),
                      '2M':(3,6),
                      '4M':(4,8)}
            }
        
        #############################
        # Specifications for Eiger2 #
        #############################
        detectors['EIGER2'] = {
            'pxs' : 75e-3,  # [mm] Pixel size
            'hmp' : 1028,   # [px] Module size (horizontal
            'vmp' : 512,    # [px] Module size (vertical)
            'hgp' : 12,     # [px] Module gap (horizontal)
            'vgp' : 38,     # [px] Module gap (vertical)
            'cbh' : 0,      # [px] Central beam hole
            'size' : { '1M':(1,2),
                       '4M':(2,4),
                       '9M':(3,6),
                      '16M':(4,8)},
            }
        
        #############################
        # Specifications for MPCCD #
        #############################
        detectors['MPCCD'] = {
            'pxs' : 50e-3,   # [mm] Pixel size
            'hmp' : 1024,    # [px] Module size (horizontal
            'vmp' : 512,     # [px] Module size (vertical)
            'hgp' : 18,      # [px] Module gap (horizontal)
            'vgp' : 27,      # [px] Module gap (vertical)
            'cbh' : 60,      # [px] Central beam hole
            'size' : {'4M':(2,4)},
            }

        ##############################
        # Specifications for RAYONIX #
        ##############################
        detectors['RAYONIX'] = {
            'pxs' : 39e-3,  # [mm] Pixel size
            'hmp' : 1920,   # [px] Module size (horizontal
            'vmp' : 1920,   # [px] Module size (vertical)
            'hgp' : 0,      # [px] Module gap (horizontal)
            'vgp' : 0,      # [px] Module gap (vertical)
            'cbh' : 0,      # [px] Central beam hole
            'size' : {'MX225-HS':(3,3),
                      'MX300-HS':(4,4)},
            }
        
        #############################
        # Specifications for PHOTON #
        #############################
        detectors['PHOTON-II'] = {
            'pxs' : 135e-3, # [mm] Pixel size
            'hmp' : 768,    # [px] Module size (horizontal
            'vmp' : 512,    # [px] Module size (vertical)
            'hgp' : 0,      # [px] Module gap (horizontal)
            'vgp' : 0,      # [px] Module gap (vertical)
            'cbh' : 0,      # [px] Central beam hole
            'size' : { '7':(1,1),
                      '14':(1,2)},
            }
        
        #############################
        # Specifications for PHOTON #
        #############################
        detectors['PHOTON-III'] = {
            'pxs' : 135e-3, # [mm] Pixel size
            'hmp' : 768,    # [px] Module size (horizontal
            'vmp' : 512,    # [px] Module size (vertical)
            'hgp' : 0,      # [px] Module gap (horizontal)
            'vgp' : 0,      # [px] Module gap (vertical)
            'cbh' : 0,      # [px] Central beam hole
            'size' : { '7':(1,1),
                      '14':(1,2),
                      '28':(2,2)},
            }
        
        ###################################
        # Specifications for Perkin-Elmer #
        ###################################
        detectors['Perkin-Elmer XRD'] = {
            'pxs' : 100e-3, # [mm] Pixel size
            'hmp' : 2048,   # [px] Module size (horizontal
            'vmp' : 2048,   # [px] Module size (vertical)
            'hgp' : 0,      # [px] Module gap (horizontal)
            'vgp' : 0,      # [px] Module gap (vertical)
            'cbh' : 0,      # [px] Central beam hole
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
            'pxs' : 100e-3, # [mm] Pixel size
            'hmp' : 4288,   # [px] Module size (horizontal)
            'vmp' : 4288,   # [px] Module size (vertical)
            'hgp' : 0,      # [px] Module gap (horizontal)
            'vgp' : 0,      # [px] Module gap (vertical)
            'cbh' : 0,      # [px] Central beam hole
            'size' : {'4343':(1,1)},
            }
        
        #############################
        # Specifications for CITIUS #
        #############################
        detectors['CITIUS'] = {
            'pxs' : 72.6e-3, # [mm] Pixel size
            'hmp' : 728,     # [px] Module size (horizontal)
            'vmp' : 384,     # [px] Module size (vertical)
            'hgp' : 22,      # [px] Module gap (horizontal)
            'vgp' : 36,      # [px] Module gap (vertical)
            'cbh' : 0,       # [px] Central beam hole
            'size' : {'20.2':(6,12)},
            }
        
        # make file dump
        if not os.path.exists(self.path_detdb) or reset:
            with open(self.path_detdb, 'w') as wf:
                json.dump(detectors, wf, indent=4)
        # read from detector db file
        else:
            try:
                temp = {}
                with open(self.path_detdb, 'r') as of:
                    for key, vals in json.load(of).items():
                        temp[key] = vals
                # 'hms' key indicated old mm detector module format
                # -> backup and update!
                if 'hms' in next(iter(temp.values())):
                    # this is due to the old [mm] 'hms' detector module format
                    # after an update to version 2.0.0 make backup of old
                    # detector entries. All changes to detector_db.json are lost!
                    with open(f'{self.path_detdb}.bak', 'w') as wf:
                        json.dump(temp, wf, indent=4)
                    # save db in new detector format
                    update = True
                else:
                    detectors.update(temp)
            except: # any error is critical here!
                print(f"Error parsing Detector db at: {self.path_detdb}")
                raise SystemExit
        
        if update:
            with open(self.path_detdb, 'w') as wf:
                json.dump(detectors, wf, indent=4)
        
        return detectors

    def get_tooltips(self):
        """
        Generate tooltips for various configuration parameters.

        This method initializes a dictionary of tooltips for different categories
        of configuration parameters, including 'geo', 'plo', 'thm', and 'lmt'.
        Each category contains key-value pairs where the key is the parameter name
        and the value is a string describing the parameter and its expected input.

        Categories:
        - 'geo': Geometric parameters related to the detector setup.
        - 'plo': Plotting parameters for contour lines and markers.
        - 'thm': Theme-related parameters for color settings.
        - 'lmt': Limit parameters for various settings.

        The method also includes a debug function to check if all keys have a tooltip.

        Example:
            self.get_tooltips()

        Raises:
            KeyError: If a parameter key is missing a tooltip in debug mode.
        """
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

    def get_hotkeys(self):
        """
        Initializes a dictionary mapping hotkeys to their respective functions and descriptions.
        The dictionary structure is as follows:
        - Keys are tuples of the form ('key', 'modifier'), where 'key' is a valid QtKey string
          and 'modifier' is an optional key modifier (e.g., 'SHIFT', 'CTRL', 'ALT').
        - Values are dictionaries with the following keys:
          - 'fn': A tuple containing the function to be executed and its argument.
          - 'dc': A string description of the hotkey's function.
        The dictionary is divided into sections, each indicated by a key with a value of None.
        Available sections and their hotkeys:
        - 'Windows':
          - ('F1', None): Show about window
          - ('F2', None): Show geometry window
          - ('F3', None): Show hotkey window
          - ('X', None): Show PXRD window
          - ('F', None): Show FWHM window
        - 'Unit':
          - ('T', None): Change to 2-theta
          - ('D', None): Change to d-spacing
          - ('Q', None): Change to Q-space
          - ('S', None): Change to sin(θ)/λ
        - 'Label position':
          - ('Up', None): Top left
          - ('Down', None): Bottom left
        - 'Toggles':
          - ('P', None): Toggle polarisation overlay
          - ('A', None): Toggle solid angle overlay
          - ('H', None): Toggle highlight / transparency
          - ('U', None): Toggle unit hover display
          - ('G', None): Toggle azimuthal grid
        - 'Colormaps':
          - ('C', None): Next colormap
          - ('C', 'SHIFT'): Previous colormap
          - ('R', None): Toggle colored reference
        Returns:
            None
        """
        # (str):None indicated a new section / tile using str
        # and tuple('key', 'modifier'):{'fn':(function, arg), 'dc':'description'}
        # to map hotkeys to functions and descriptions
        #
        # 'key' is the key to be pressed, must be a valid QtKey string and mappable
        # to the QtGui.QKeySequence(key).toString() function used in the keyReleaseEvent
        #
        # the following keymodifier mappings are available:
        # QtCore.Qt.KeyboardModifier.AltModifier:'ALT'
        # QtCore.Qt.KeyboardModifier.ShiftModifier:'SHIFT'
        # QtCore.Qt.KeyboardModifier.ControlModifier:'CTRL'
        self.hotkey_dict = {('Windows'):None,
                            ('F1', None):  {'fn':(self.about_win.show, None),
                                            'dc':'Show about window'},
                            ('F2', None):  {'fn':(self.geometry_win.show, None),
                                            'dc':'Show geometry window'},
                            ('F3', None):  {'fn':(self.hotkeys_win.show, None),
                                            'dc':'Show hotkey window'},
                            ('X', None):   {'fn':(self.win_pxrd_plot, None),
                                            'dc':'Show PXRD window'},
                            ('F', None):   {'fn':(self.fwhm_win.show, None),
                                            'dc':'Show FWHM window'},
                            #('B', None):   {'fn':(self.abs_win.show, None),
                            #                'dc':'Show Absorption window'},
                            ('O', None):   {'fn':(self.detdb_win.show, None),
                                            'dc':'Show Detector DB'},
                            ('E', None):   {'fn':(self.export_win.show, None),
                                            'dc':'Show Export window'},
                            
                            ('Unit'):None,
                            ('T', None):   {'fn':(self.change_units, 0),
                                            'dc':'Change to 2-theta'},
                            ('D', None):   {'fn':(self.change_units, 1),
                                            'dc':'Change to d-spacing'},
                            ('Q', None):   {'fn':(self.change_units, 2),
                                            'dc':'Change to Q-space'},
                            ('S', None):   {'fn':(self.change_units, 3),
                                            'dc':'Change to sin(\u03B8)/\u03BB'},
                            
                            ('Label position'):None,
                            ('Up', None):  {'fn':(self.label_set_position, 'top'),
                                            'dc':'Top left'},
                            ('Down', None):{'fn':(self.label_set_position, 'bottom'),
                                            'dc':'Bottom left'},
                            
                            ('Toggles'):None,
                            ('P', None):   {'fn':(self.toggle_overlay_polarisation, None),
                                            'dc':'Toggle polarisation overlay'},
                            ('A', None):   {'fn':(self.toggle_overlay_solidangle, None),
                                            'dc':'Toggle solid angle overlay'},
                            ('H', None):   {'fn':(self.toggle_overlay_highlight, None),
                                            'dc':'Toggle highlight / transparency'},
                            ('U', None):   {'fn':(self.toggle_unit_hover, None),
                                            'dc':'Toggle unit hover display'},
                            ('G', None):   {'fn':(self.toggle_grid, None),
                                            'dc':'Toggle azimuthal grid'},
                            
                            ('Colormaps'):None,
                            ('C', None):   {'fn':(self.change_cmap, 1),
                                            'dc':'Next colormap'},
                            ('C', 'SHIFT'):{'fn':(self.change_cmap, -1),
                                            'dc':'Previous colormap'},
                            ('R', None):   {'fn':(self.toggle_colored_reference, None),
                                            'dc':'Toggle colored reference'},
                            }

    #############
    #   BUILD   #
    #############
    def build_detector(self):
        """
        Builds the detector by placing modules in a grid pattern based on the detector's configuration.

        The method calculates the positions of the detector modules in millimeters, taking into account
        the pixel size, module gaps, and central beam hole size. It then creates and adds graphical 
        rectangle items representing the modules to the scene.

        The placement logic ensures that the beam position is either between the modules (for even 
        numbers of modules) or at the center module (for odd numbers of modules). Additionally, it 
        handles specific configurations for MPCCD detectors used at SACLA/SPring-8.

        The graphical representation of each module is customized with specified colors, pen width, 
        fill, and opacity.

        Attributes:
            self.det.hmp (float): Horizontal module size in pixels.
            self.det.vmp (float): Vertical module size in pixels.
            self.det.hgp (float): Horizontal gap between modules in pixels.
            self.det.vgp (float): Vertical gap between modules in pixels.
            self.det.cbh (float): Central beam hole size in pixels.
            self.det.hmn (int): Number of horizontal modules.
            self.det.vmn (int): Number of vertical modules.
            self.det.pxs (float): Pixel size in millimeters.
            self.det_module_color (QColor): Color of the module outline.
            self.plo.det_module_width (int): Width of the module outline.
            self.det_module_fill (QColor): Fill color of the module.
            self.plo.det_module_alpha (float): Opacity of the module.
            self.ax (QGraphicsScene): The scene to which the modules are added.

        Returns:
            None
        """
        # build detector modules
        # pixel -> mm
        _hms = self.det.hmp * self.det.pxs
        _vms = self.det.vmp * self.det.pxs
        _hgs = self.det.hgp * self.det.pxs
        _vgs = self.det.vgp * self.det.pxs
        _cbh = self.det.cbh * self.det.pxs
        # beam position is between the modules (even) or at the center module (odd)
        # determined by the "+det.hmn%2" part
        for i in range(-self.det.hmn//2+self.det.hmn%2, self.det.hmn-self.det.hmn//2):
            for j in range(-self.det.vmn//2+self.det.vmn%2, self.det.vmn-self.det.vmn//2):
                # - place modules along x (i) and y (j) keeping the gaps in mind ( + (det.hgp*det.pxs)/2)
                # - the " - ((det.hmp+det.hgp*det.pxs)/2)" positions the origin (the beam) at the center of a module
                #   and "det.hmn%2" makes sure this is only active for detectors with an odd number of modules
                # - define sets of panels that collectively move to realize a central hole offset for MPCCD detectors
                #   that are used at SACLA/SPring-8:
                #   x = (...) + (det.cbh/2)*(2*(j&det.vmn)//det.vmn-1)
                #   y = (...) + (det.cbh/2)*(1-2*(i&det.hmn)//det.hmn)
                # - negative values of det.cbh for 'clockwise' offset order
                origin_x = i * (_hms + _hgs) \
                             - ((_hms + _hgs)/2) * (self.det.hmn % 2) \
                             + (_hgs)/2 \
                             + (_cbh/2) * (2*(j & self.det.vmn) // self.det.vmn-1)
                origin_y = j * (_vms + _vgs) \
                             - ((_vms + _vgs)/2) * (self.det.vmn%2) \
                             + (_vgs/2) \
                             + (_cbh/2) * (1-2*(i & self.det.hmn) // self.det.hmn)
                # add the module
                rect_item = QtWidgets.QGraphicsRectItem(origin_x, origin_y,  _hms, _vms)
                rect_item.setPen(pg.mkPen(color = self.det_module_color, width = self.plo.det_module_width))
                rect_item.setBrush(pg.mkBrush(color = self.det_module_fill))
                rect_item.setOpacity(self.plo.det_module_alpha)
                self.ax.addItem(rect_item)

    ##########
    #  PXRD  #
    ##########
    def win_pxrd_plot(self, keep=False):
        """
        Displays or updates the Powder X-ray Diffraction (PXRD) plot window.
        If the PXRD window is already visible and `keep` is False, the window will be closed.
        Otherwise, the window will be updated or created if it does not exist.
        Parameters:
        -----------
        keep : bool, optional
            If True, keeps the PXRD window open even if it is already visible. Default is False.
        Notes:
        ------
        - The PXRD window includes a plot of the PXRD data, with options to interact with the plot.
        - A disclaimer box is included in the window, indicating that the feature is in the test phase.
        - A button to set up Full Width at Half Maximum (FWHM) is also included in the disclaimer box.
        """
        if self.pxrd_win is not None:
            self.win_pxrd_update()
            self.pxrd_win.adjustSize()
            self.pxrd_win.show(keep=keep)
            return

        if self.cont_ref_dsp is None or self.cont_ref_hkl is None:
            return
        
        font = QtGui.QFont()
        font.setPixelSize(self.plo.conic_hkl_label_size)
        font.setBold(True)

        self.pxrd_win = HotkeyDialog(parent=self, flags=QtCore.Qt.WindowType.Tool)
        self.pxrd_win.setStyleSheet('QGroupBox { font-weight: bold; }')

        layout = QtWidgets.QVBoxLayout()
        
        powder_box = QtWidgets.QGroupBox('')
        powder_box_layout = QtWidgets.QVBoxLayout()
        powder_box_layout.setContentsMargins(6,6,6,6)
        powder_box.setLayout(powder_box_layout)
        
        self.pxrd_curve = pg.PlotCurveItem(mousewidth=20)
        self.pxrd_curve.setPen(pg.mkPen(color=self.palette().text().color()))
        self.pxrd_curve.setClickable(True)
        self.pxrd_curve.sigClicked.connect(self.win_pxrd_curve_clicked)
        self.pxrd_scatt = ClickableScatterPlotItem(size=20, tip=None)
        self.pxrd_scatt.setPen(pg.mkPen(color=self.palette().highlight().color()))
        self.pxrd_scatt.setBrush(pg.mkBrush(color=self.palette().highlight().color()))
        self.pxrd_scatt.setSymbol(self.plo.pxrd_marker_symbol)
        self.pxrd_scatt.sigClicked.connect(self.win_pxrd_hkl_clicked)
        self.pxrd_scatt_highlighted = None

        self.pxrd_label = pg.TextItem(anchor=(0.5,0.9))
        self.pxrd_label.setFont(font)
        self.pxrd_regio = pg.LinearRegionItem(pen=pg.mkPen(None))
        self.pxrd_regio.setMovable(False)
        self.pxrd_plot = pg.PlotWidget(background=self.palette().base().color())
        self.pxrd_ghost_legend = self.pxrd_plot.addLegend(offset=(0,1), verSpacing=0, labelTextSize='12pt')
        self.pxrd_plot.viewport().setAttribute(QtCore.Qt.WidgetAttribute.WA_AcceptTouchEvents, False)
        self.pxrd_plot.setMenuEnabled(False)
        self.pxrd_plot.setLabel(axis='left', text='Intensity [arb. units]', color=self.palette().text().color())
        self.pxrd_plot.showGrid(x=True, y=True, alpha=0.5)
        self.pxrd_plot.getPlotItem().getViewBox().setDefaultPadding(0.0)
        self.pxrd_plot.addItem(self.pxrd_curve, zlevel=1)
        # if points are not visible, e.g. outside of currentrange, the plot size will
        # not be adjusted properly and stays larger than expected.
        # ignoreBound allows the plot to autorange disregarding the 'invisble' points
        # however we need to increase the min y padding to make the scatterplot visible!
        self.pxrd_plot.addItem(self.pxrd_scatt, ignoreBounds=True)
        self.pxrd_plot.addItem(self.pxrd_label)
        # same as above, allow the plot to disregard the beamstop if out of bounds
        self.pxrd_plot.addItem(self.pxrd_regio, ignoreBounds=True)
        # this invisible line will do the padding for us -> thanks Frederik!
        self.pxrd_line = pg.InfiniteLine(pos=0, angle=0, pen=pg.mkPen(None))
        self.pxrd_plot.addItem(self.pxrd_line)
        

        self.pxrd_ghost_curves = [] # list to store ghost lines
        self.pxrd_ghost_bank = [('blinky', '#fb0009'), ('pinky', '#fc60ff'), ('inky', '#22ffff'),
                                ('clyde', '#f27a10'), ('sue', '#9208ff'), ('dinky', '#e1e1e1'),
                                ('miru', '#80ff14'), ('tim', '#dd8308'), ('funky', '#18b804'),
                                ('kinky', '#ffdc14'), ('orson', '#96a382')]

        powder_box_layout.addWidget(self.pxrd_plot)
        layout.addWidget(powder_box)

        # Disclimer box info
        citation_box = QtWidgets.QGroupBox()
        citation_box.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        citation_box_layout = QtWidgets.QHBoxLayout()
        citation_box_layout.setContentsMargins(6,6,6,6)
        citation_box.setLayout(citation_box_layout)
        # add ghost plot
        button_ghost = QtWidgets.QPushButton('Ghosts')
        button_ghost.setToolTip('Add or remove ghost lines to the PXRD plot.\nSet visibility in the legend.\nRemove all invisible ghosts or last added.')
        button_ghost.clicked.connect(self.win_pxrd_add_ghost)
        menu = QtWidgets.QMenu()
        menu.addAction('Add', self.win_pxrd_add_ghost)
        menu.addAction('Remove', self.win_pxrd_rem_ghost)
        button_ghost.setMenu(menu)
        citation_box_layout.addWidget(button_ghost, alignment=QtCore.Qt.AlignmentFlag.AlignLeft)
        # add citation
        citation = QtWidgets.QLabel('This feature is currently in <b>test phase</b>, feedback is very welcome!')
        citation.setOpenExternalLinks(True)
        citation_box_layout.addWidget(citation, alignment=QtCore.Qt.AlignmentFlag.AlignLeft)
        # Add the fwhm button
        button_fwhm = QtWidgets.QPushButton('Setup FHWM')
        button_fwhm.clicked.connect(self.fwhm_win.show)
        citation_box_layout.addWidget(button_fwhm, alignment=QtCore.Qt.AlignmentFlag.AlignRight)
        layout.addWidget(citation_box)
        
        self.win_pxrd_update()
        self.pxrd_win.setWindowIcon(self.icon)
        self.pxrd_win.setLayout(layout)
        self.pxrd_win.show()

    def win_pxrd_update(self, update_dict=None):
        """
        Updates the PXRD (Powder X-ray Diffraction) window with new data and settings.
        Parameters:
        update_dict (dict, optional): A dictionary containing updated parameters for the PXRD calculation. 
                                      If None, default parameters from the object are used. The dictionary 
                                      should have the following keys:
                                      - '0': Sensor thickness
                                      - '1': Sensor material
                                      - '2': Beam divergence
                                      - '3': Energy resolution
                                      - '4': Scattering diameter
        Returns:
        None
        """
        if self.pxrd_win is None:
            return
        
        # change color if theme changed
        self.pxrd_curve.setPen(self.palette().text().color())
        self.pxrd_plot.setBackground(self.palette().base().color())
        self.pxrd_scatt.setPen(self.palette().highlight().color())
        self.pxrd_scatt.setBrush(self.palette().highlight().color())

        if self.geo.reference.lower() == 'none':
            self.pxrd_win.setWindowTitle(f'Load a cif to show PXRD pattern')
            self.pxrd_curve.setData(x=None, y=None, data=None)
            self.pxrd_scatt.setData(x=None, y=None, data=None)
            return
        
        if self.geo.reference in self.ref_cif and self.ref_cif[self.geo.reference].is_complete:
            dsp = self.ref_cif[self.geo.reference].dsp # np.array
            peak_ttr, idx = self.dsp2tth(dsp)
            inten = np.array(self.ref_cif[self.geo.reference].hkl)[:,4][idx] # list[ (h k l I Irel) ]
            hkl = list(map(tuple, np.array(self.ref_cif[self.geo.reference].hkl)[:,:3][idx].astype(int)))
        else:
            self.pxrd_win.setWindowTitle(f'Load a cif to show PXRD pattern')
            self.pxrd_curve.setData(x=None, y=None, data=None)
            self.pxrd_scatt.setData(x=None, y=None, data=None)
            return
        
        if update_dict is None:
            thk = self.plo.sensor_thickness
            mat = self.plo.sensor_material
            div = self.plo.beam_divergence
            dEE = self.plo.energy_resolution
            dia = self.plo.scattering_diameter
        else:
            thk = update_dict['0']
            mat = update_dict['1']
            div = update_dict['2']
            dEE = update_dict['3']
            dia = update_dict['4']

        self.pxrd_win.setWindowTitle(f'PXRD pattern of {self.geo.reference}')
        self.pxrd_plot.setTitle(f'Estimated diffration pattern for a <b>perpendicular\
                                 & flat</b> detector geometry (F\u00B2 & d-spacing @\
                                 {self.plo.conic_ref_cif_kev} keV).', color=self.palette().text().color())
        
        fwhm = self.calc_FWHM(dis=self.geo.dist * 1e-3,
                              dia=dia,
                              thk=thk,
                              mat=mat,
                              pxs=self.det.pxs * 1e-3,
                              tth=peak_ttr,
                              nrg=self.geo.ener,
                              div=div,
                              dEE=dEE,
                              deg=False)
        # calc visible extent
        if self._tth is not None and not isinstance(self._tth, float):
            max_res_r = np.nanmax(self._tth)
            min_res_r = max(np.nanmin(self._tth), fwhm.mean()/10)
        else:
            max_res_r = self.calc_tth_max()
            min_res_r = fwhm.mean()/10

        # add beamstop shadow
        if self.geo.bssz and not (isinstance(self.geo.bssz, str) and self.geo.bssz.lower() == 'none'):
            if min_res_r < self.bs_theta:
                regio_min = self.calc_unit(min_res_r)
                regio_max = self.calc_unit(self.bs_theta)
                self.pxrd_regio.setRegion((regio_min, regio_max))
                regio_color = self.beamstop_color
                regio_color.setAlphaF(0.10)
                self.pxrd_regio.setBrush(regio_color)
                self.pxrd_regio.setVisible(True)
        else:
            self.pxrd_regio.setVisible(False)
        
        ttr = np.arange(min_res_r, max_res_r, fwhm.mean()/10)

        gauss = np.sum(self.gaussian(ttr, np.atleast_2d(peak_ttr).T, np.atleast_2d(fwhm).T) * np.atleast_2d(inten).T, axis=0)
        xval = self.calc_unit(ttr)
        self.pxrd_curve.setData(x=xval, y=gauss)
        
        peak_xval = self.calc_unit(peak_ttr)
        self.pxrd_scatt_offset = -(gauss.max() + gauss.min()) * self.plo.pxrd_marker_offset
        self.pxrd_scatt.setData(x=peak_xval,
                                y=np.zeros(len(peak_xval))+self.pxrd_scatt_offset,
                                data=hkl)
        # remove peaks that are ouside of the detector screen
        cond = (peak_xval <= xval.max()) & (peak_xval >= xval.min())
        self.pxrd_scatt.setPointsVisible(cond)
        self.pxrd_line.setPos(self.pxrd_scatt_offset*2)

        # add axis label
        self.pxrd_plot.setLabel(axis='bottom', text=self.unit_names[self.geo.unit], color=self.palette().text().color())
        
        # flip axis for d-spacing
        if self.geo.unit == 1:
            self.pxrd_plot.getPlotItem().getViewBox().invertX(True)
            #self.pxrd_plot.setXRange(20, xval.min(), padding=0)
        else:
            self.pxrd_plot.getPlotItem().getViewBox().invertX(False)

        self.pxrd_plot.getPlotItem().getViewBox().setLimits(xMin=xval.min(),
                                                            xMax=xval.max(),
                                                            yMin=self.pxrd_scatt_offset*2,
                                                            yMax=gauss.max())

        if self.pxrd_scatt_highlighted is not None:
            self.win_pxrd_highlight(self.pxrd_scatt_highlighted.index())

    def win_pxrd_hkl_clicked(self, widget, points, event):
        if not widget.name:
            return
        elif len(points) > 0:
            point = self.get_closest_point_x(points, event.pos())
            if event.buttons() == QtCore.Qt.MouseButton.RightButton:
                self.patches['reference'][point.index()].lowlight()
                return
            # highlight in main window
            # calls the highlight method of the HoverableCurveItem
            # calls win_pxrd_highlight
            self.patches['reference'][point.index()].highlight()

    def win_pxrd_curve_clicked(self, widget, event):
        """
        Handles the event when a PXRD curve is clicked in the window.

        This method is triggered when the user clicks on the PXRD curve in the window.
        It processes the click event, determines the position of the click, and identifies
        the corresponding points on the PXRD scatter plot. It then calls another method
        to handle the click event on the identified points.

        Args:
            widget: The widget that received the click event.
            event: The event object containing details about the click event.

        Returns:
            None
        """
        pos = event.pos()
        pos.setY(self.pxrd_scatt_offset)
        points = self.pxrd_scatt.pointsAt(pos)
        self.win_pxrd_hkl_clicked(self.pxrd_scatt, points, event)

    def win_pxrd_set_hkl_label(self, point):
        self.pxrd_label.setText(' '.join(map(str, point.data())))
        self.pxrd_label.setPos(point.pos())
        self.pxrd_label.setVisible(True)
    
    def win_pxrd_highlight(self, index):
        # highlight from main window
        if self.pxrd_win is not None:
            self.win_pxrd_lowlight()
            point = self.pxrd_scatt.points()[index]
            point.setPen(pg.mkPen(self.conic_highlight))
            point.setBrush(pg.mkBrush(self.conic_highlight))
            self.pxrd_scatt_highlighted = point
            self.win_pxrd_set_hkl_label(point)

    def win_pxrd_lowlight(self):
        # called by HoverableCurveItem:lowlight
        if self.pxrd_win is not None:
            if self.pxrd_scatt_highlighted is not None:
                _pen = pg.mkPen(color=self.palette().highlight().color())
                _brush = pg.mkBrush(color=self.palette().highlight().color())
                self.pxrd_scatt_highlighted.setPen(_pen)
                self.pxrd_scatt_highlighted.setBrush(_brush)
                self.pxrd_scatt_highlighted = None
            self.pxrd_label.setVisible(False)

    def win_pxrd_add_ghost(self):
        """
        Adds a ghost curve to the PXRD plot from the ghost bank.

        This method pops the first ghost curve from the `pxrd_ghost_bank`, sets its data and color,
        and then appends it to the `pxrd_ghost_curves` list. The ghost curve is then added to the 
        PXRD plot with a z-level based on the number of ghost curves.

        Attributes:
            pxrd_ghost_bank (list): A list of tuples containing the name and color of ghost curves.
            pxrd_curve (object): An object containing the x and y data for the PXRD curve.
            pxrd_ghost_curves (list): A list to store the ghost curves added to the plot.
            pxrd_plot (object): The plot object where the ghost curves are added.

        Returns:
            None
        """
        if len(self.pxrd_ghost_bank) > 0:
            name, color = self.pxrd_ghost_bank.pop(0)
            pxrd_ghost = pg.PlotCurveItem(name=name)
            pxrd_ghost.setData(x=self.pxrd_curve.xData, y=self.pxrd_curve.yData)
            pxrd_ghost.setPen(color)
            self.pxrd_ghost_curves.append(pxrd_ghost)
            self.pxrd_plot.addItem(pxrd_ghost, zlevel=-len(self.pxrd_ghost_curves), ignoreBounds=True)

    def win_pxrd_rem_ghost(self):
        """
        Removes ghost curves from the PXRD plot that are not visible.

        This method iterates through the list of current PXRD ghost curves and 
        removes those that are not visible from the plot. The removed ghost curves 
        are then added to the ghost bank with their name and pen color. If no ghost 
        curves are removed and there are still ghost curves present, the last ghost 
        curve in the list is removed from the plot.

        Attributes:
            pxrd_ghost_curves (list): List of current PXRD ghost curves.
            pxrd_plot (object): The PXRD plot object.
            pxrd_ghost_bank (list): List to store removed ghost curves with their 
                                    name and pen color.
        """
        remove_later = []
        for ghost in self.pxrd_ghost_curves:
            if not ghost.isVisible():
                remove_later.append(ghost)
        for ghost in remove_later:
            self.pxrd_plot.getPlotItem().removeItem(ghost)
            self.pxrd_ghost_curves.remove(ghost)
            self.pxrd_ghost_bank.insert(0, (ghost.name(), ghost.opts['pen'].color().name()))
        if len(remove_later) == 0 and len(self.pxrd_ghost_curves) > 0:
                ghost = self.pxrd_ghost_curves[-1]
                self.pxrd_plot.getPlotItem().removeItem(ghost)
                self.pxrd_ghost_curves.remove(ghost)

    #######
    # CIF #
    #######
    def get_ref_db_from_file(self):
        """
        Reads a reference database from a file and updates the internal reference dictionary.

        This method checks if the file specified by `self.path_cif_db` exists. If it does, it opens the file,
        reads its contents as JSON, and iterates through the items. For each item, it creates a `Ref` object
        and adds it to the internal reference dictionary `self.ref_cif`. Additionally, it updates the reference
        menu with the new reference names.

        Raises:
            FileNotFoundError: If the file specified by `self.path_cif_db` does not exist.
        """
        if os.path.exists(self.path_cif_db):
            with open(self.path_cif_db, 'r') as of:
                temp = json.load(of)
            for name, item in list(temp.items()):
                self.ref_cif[name] = Ref(name=name, cif=item['cif'])
                self.add_cif_to_menu(name)

    def save_ref_db_to_file(self):
        """
        Saves the reference CIF database to a file.

        This method iterates through the `ref_cif` dictionary, checks if each entry has a CIF file,
        and if so, adds it to a temporary dictionary. The temporary dictionary is then written to
        the file specified by `self.path_cif_db` in JSON format with an indentation of 4 spaces.

        Attributes:
            self.path_cif_db (str): The file path where the CIF database will be saved.
            self.ref_cif (dict): A dictionary containing reference CIF entries.

        Raises:
            IOError: If there is an issue writing to the file.
        """
        temp = {}
        with open(self.path_cif_db, 'w') as wf:
            for k in self.ref_cif:
                if self.ref_cif[k].has_cif:
                    temp[k] = {}
                    temp[k]['cif'] = self.ref_cif[k].cif
            json.dump(temp, wf, indent=4)

    def remove_ref_from_db(self, name):
        """
        Removes a reference from the database and updates the UI accordingly.
        Args:
            name (str): The name of the reference to be removed.
        This method performs the following actions:
        1. Removes the reference from the internal reference dictionary.
        2. Saves the updated reference database to a file.
        3. Closes the reference window if it is currently displayed.
        4. Updates the geometry reference to 'None' if the removed reference was the current reference.
        5. Changes the current reference in the UI.
        6. Activates the main window to ensure the menu is updated correctly.
        7. Clears and repopulates the reference menu with the remaining references.
        """
        self.ref_cif.pop(name)
        self.save_ref_db_to_file()

        # close window if currently displayed
        if self.geo.reference == name:
            if self.pxrd_win is not None:
                if self.pxrd_win.isVisible():
                    self.pxrd_win.close()
            self.geo.reference = 'None'
            self.change_reference(self.geo.reference)
        
        # focus main window to update the menu correctly
        self.activateWindow()
        # reset the menu
        self.sub_menu_cif.clear()
        # populate the menu
        for name in self.ref_cif.keys():
            self.add_cif_to_menu(name)

        # toggle enable/disable if empty
        self.sub_menu_cif.setDisabled(self.sub_menu_cif.isEmpty())

    def add_cif_to_menu(self, name):
        """
        Adds a reference to the submenu if it does not already exist.
        This method checks if a submenu with the given name already exists. If it does not,
        it creates a new submenu with the specified name and adds actions to open and delete
        the reference. The 'Open' action is checkable and is linked to the `change_reference`
        method. The 'Delete' action is linked to the `remove_ref_from_db` method.
        Args:
            name (str): The name of the reference to be added to the submenu.
        """
        
        # add only if no menu with same name exists
        if name in [item.text() for item in self.sub_menu_cif.actions()]:
            return
        
        ref_menu = QtWidgets.QMenu(name, self)
        self.sub_menu_cif.addMenu(ref_menu)
        self.sub_menu_cif.setEnabled(True)
        cif = self.ref_cif[name].cif
        if os.path.exists(cif):
            ref_action_open = QtGui.QAction('Open', self, checkable=True)
            self.menu_set_action(ref_action_open, self.change_reference, name)
            ref_menu.addAction(ref_action_open)
            self.group_ref.addAction(ref_action_open)
            if self.geo.reference == name:
                ref_action_open.setChecked(True)
        ref_action_del = QtGui.QAction('Delete', self)
        self.menu_set_action(ref_action_del, self.remove_ref_from_db, name)
        ref_menu.addAction(ref_action_del)

        # disable/enable menu if empty
        self.sub_menu_cif.setDisabled(self.sub_menu_cif.isEmpty())

    def calc_ref_from_cif(self, fpath, open_pxrd=True):
        """
        Calculate reference data from a CIF (Crystallographic Information File).

        This method processes a CIF file to extract crystallographic data and 
        generate powder diffraction patterns. The results are stored in the 
        object's attributes and can be optionally displayed in a plot.

        Args:
            fpath (str): The file path to the CIF file.
            open_pxrd (bool, optional): If True, opens the powder X-ray diffraction 
                plot. Defaults to True.

        Returns:
            None

        Attributes:
            cont_ref_dsp (numpy.ndarray): Array of d-spacing values for the strongest reflections.
            cont_ref_hkl (list of tuples): List of tuples containing (h, k, l, intensity, relative intensity) 
                for the strongest reflections.
            geo.reference (str): The basename of the CIF file.
            ref_cif (dict): Dictionary containing reference data with the CIF file basename as the key.
        """
        # called when a cif is dropped onto the window
        self.xtl = dif.Crystal(fpath)
        # :return xval: arrray : x-axis of powder scan (units)
        # :return inten: array : intensity values at each point in x-axis
        # :return reflections: (h, k, l, xval, intensity) array of reflection positions, grouped by min_overlap
        xval, inten, reflections = self.xtl.Scatter.powder(scattering_type='xray', units='dspace', powder_average=True, min_overlap=0.02, energy_kev=self.plo.conic_ref_cif_kev)
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
        # update entry, might be incomplete (load on startup only has cif path)
        # and we already did the heavy lifting
        self.ref_cif[self.geo.reference] = Ref(name=self.geo.reference, dsp=self.cont_ref_dsp, hkl=self.cont_ref_hkl, cif=fpath)
        # add to menu
        self.add_cif_to_menu(self.geo.reference)

        # update window title
        self.set_win_title()

        if open_pxrd:
            self.win_pxrd_plot(keep=True)

        self.draw_reference()

        self.save_ref_db_to_file()

    #############
    #  UTILITY  #
    #############
    def resize_win(self):
        """
        Adjusts the dimensions of the plotting window and sets the appropriate ranges for the axes.
        This method calculates the proper dimensions for the plot based on the detector's properties
        and the plot's padding. It then limits the x and y axes ranges accordingly. If the plot size
        is not set, it determines an appropriate size based on the screen's available height. Finally,
        it sets the window's width and height, and if the plot size is fixed, it ensures the window
        cannot be resized beyond these dimensions.
        Attributes:
            xdim (float): The calculated width dimension for the plot.
            ydim (float): The calculated height dimension for the plot.
            win_width (int): The final width of the window.
            win_height (int): The final height of the window.
        """
        # figure out proper plot dimensions
        # this adds a padding of hmn x vmn
        # to remove:  + self.det.hgp * (self.det.hmn -1)
        #             + self.det.vgp * (self.det.hmn -1)
        self.xdim = (self.det.hmp * self.det.hmn + self.det.hgp * (self.det.hmn-1) + self.det.cbh)/2 * self.det.pxs + self.plo.plot_padding
        self.ydim = (self.det.vmp * self.det.vmn + self.det.vgp * (self.det.vmn-1) + self.det.cbh)/2 * self.det.pxs + self.plo.plot_padding
        
        # limit the axis x and y
        self.ax.setXRange(-self.xdim, self.xdim, padding=0, update=True)
        self.ax.setYRange(-self.ydim, self.ydim, padding=0, update=True)

        if self.plo.plot_size <= 0:
            _app = QtWidgets.QApplication.instance()
            _height = _app.primaryScreen().availableGeometry().height()
            self.plo.plot_size = int(np.ceil(_height*0.9))
        
        # get proper dimensions
        self.win_width = int(np.ceil(self.plo.plot_size * self.xdim / self.ydim))
        self.win_height = self.plo.plot_size + self.plo.slider_margin + self.offset_win32

        # fix the window size
        if self.plo.plot_size_fixed:
            self.setMaximumHeight(self.win_height)
            self.setMinimumHeight(self.win_height)
            self.setMaximumWidth(self.win_width)
            self.setMinimumWidth(self.win_width)

        # resize the window
        self.resize(self.win_width, self.win_height)

    def set_win_title(self):
        """
        Sets the window title based on the detector name, geometry reference, and active settings.

        If the geometry reference is 'none' (case insensitive), the window title will be set to 
        "<detector name> - <active settings>". Otherwise, it will be set to 
        "<detector name> - <geometry reference> - <active settings>".

        Attributes:
            geo (object): An object containing the geometry reference.
            det (object): An object containing the detector name.
            active_settings (str): The current active settings.
        """
        if self.geo.reference.lower() == 'none':
            self.setWindowTitle(f'{self.det.name} - {self.active_settings}')
        else:
            self.setWindowTitle(f'{self.det.name} - {self.geo.reference} - {self.active_settings}')

    def calc_unit(self, tth):
        """
        Calculate the unit based on the given 2-Theta value.

        This function converts the given 2-Theta value (in radians) to various units
        based on the geometry energy and the selected unit type.

        Parameters:
        tth (float): The 2-Theta value in radians.

        Returns:
        float: The calculated unit value based on the selected unit type.
               The unit type is determined by `self.geo.unit`:
               - 0: 2-Theta in degrees
               - 1: d-spacing in Angstroms
               - 2: sin(Theta)/lambda multiplied by 4π
               - 3: sin(Theta)/lambda
        """
        # calc_unit expects 2-Theta in radians

        # Conversion factor keV to Angstrom: 12.398
        # sin(t)/l: np.sin(Theta) / lambda -> (12.398/geo_energy)
        stl = np.sin(tth/2)/(12.398/self.geo.ener)
        # d-spacing: l = 2 d sin(t) -> 1/2(sin(t)/l)
        dsp = 1/(2*stl)
        units = {0:np.rad2deg(tth), 1:dsp, 2:stl*4*np.pi, 3:stl}
        return units[self.geo.unit]

    def calc_tth_max(self, scale=1.0):
        """
        Calculate the maximum 2theta angle for the given geometry
         - used to generate 2theta values to draw conics
         - m is a multiplier used to shrink the 2theta angle
           to keep the maximum resolution conic visible
        
        Returns
        -------
        float, maximum 2-theta in radians
        """
        # make screen grid
        size_h = np.array([-self.xdim, self.xdim])*scale
        size_v = np.array([-self.ydim, self.ydim])*scale
        _gx, _gy = np.meshgrid(size_h, size_v, sparse=True)
        # build vector -> 3 x n x m
        _vec = np.full((3, 2, 2), self.geo.dist)
        _vec[0,:,:] = _gx - self.geo.hoff
        # Compensate for vertical offset and PONI offset caused by the tilt (sdd*tilt)
        _vec[1,:,:] = _gy + self.geo.voff - np.deg2rad(self.geo.tilt) * self.geo.dist
        # apply combined rotation and tilt
        _rot = self.rot_100(-np.deg2rad(self.geo.tilt + self.geo.rota))
        # reshape to allow matrix multiplication
        # _rot: 3 x 3 @ _norm: 3 x n*m -> _res: 3 x n x m
        _res = np.reshape(_rot @ np.reshape(_vec, (3,-1)), _vec.shape)
        # Distance POBI - pixel on grid
        R_a = np.sqrt(np.sum(_res[0:2]**2, axis=0))
        # POBI distance
        D_a = _res[2]
        # 2theta - Angle between pixel, sample, and POBI
        tth = np.arctan2(R_a, D_a)
        # 2theta max
        return tth.max()

    def calc_conic(self, omega, theta, steps=100):
        """
        Calculate the conic section formed by the intersection of a plane and a cone.
        This method computes the coordinates of the conic section (circle, ellipse, parabola, hyperbola, or line)
        formed by the intersection of a plane and a cone based on the given parameters.
        Parameters:
        omega (float): The angle of the cone's axis relative to the plane.
        theta (float): The angle of the intersecting plane.
        steps (int, optional): The number of steps for parameterization of the conic section. Default is 100.
        Returns:
        tuple: A tuple containing two arrays (x, y) representing the coordinates of the conic section.
               If the conic section is not visible, returns (False, False).
        Notes:
        - The method skips drawing smaller/larger ±90 degree contours and rejects overlap of the 'backscattering'.
        - The eccentricity of the resulting conic section is evaluated to parameterize the conic accordingly.
        - The method handles circles, ellipses, parabolas, hyperbolas, and lines based on the eccentricity value.
        - For ellipses and hyperbolas, the method checks if the conic section is visible within the given dimensions.
        References:
        - https://math.stackexchange.com/questions/4079720/on-the-equation-of-the-ellipse-formed-by-intersecting-a-plane-and-cone
        - https://www.geogebra.org/
        """
        # Apart from textbooks on conic sections,
        # https://math.stackexchange.com/questions/4079720/on-the-equation-of-the-ellipse-formed-by-intersecting-a-plane-and-cone
        # here is a great discussion on how to get
        # to the equation for the resulting ellipse
        # formed by the intersection of a plane and
        # a cone. We get the circle for free and just
        # adapt for faster calculation.
        # From there we need to work out the hyperbola
        # and parabola and here only 3d-drawing the
        # problem in geogebra made me understand what's
        # going on: https://www.geogebra.org/
        # 
        # skip drawing smaller/larger +-90 deg contours
        # reject overlap of the 'backscattering'
        # -> limitation of the current implementation
        if theta > np.pi/2 + abs(omega):
            return False, False

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
        #_ydim = self.ydim * 1.05
        # evaluate the eccentricity and parameterise
        # the resulting conic accordingly
        if abs(ecc) == 0:
            # circle
            h = (y1+y2)/2
            # check if the circle is visible
            if h - np.sqrt(y0**2 + x0**2) > np.sqrt(self.ydim**2 + self.xdim**2):
                return False, False
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
                return False, False
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
                return False, False
            w = h * np.sqrt(ecc**2-1)
            l = -np.arcsinh((_xdim + x0) / w)
            r =  np.arcsinh((_xdim - x0) / w)
            t = np.linspace(l, r, steps)
            x = x0 + w * np.sinh(t)
            y = y0 + (y1-y2)/2 - h * np.cosh(t)
        elif abs(ecc) >= 100:
            # line
            t = np.linspace(-_xdim, _xdim, steps)
            x = t
            y = y0 + np.ones(len(t)) * y1
        
        return x, y
    
    def calc_overlays(self, omega, res=150, pol=0.99):
        """
        Calculate overlays for the detector grid.
        Parameters:
        -----------
        omega : float
            The combination of rotation and tilt, in radians.
        res : int, optional
            The resolution of the grid, default is 150.
        pol : float, optional
            The polarization factor, default is 0.99.
        Returns:
        --------
        tuple
            A tuple containing:
            - grd : ndarray
                The neutral overlay grid.
            - tth : ndarray
                The 2theta angle between pixel, sample, and POBI.
            - pc : ndarray
                The polarization correction factor.
            - sa : ndarray
                The solid angle correction factor.
            - fwhm : ndarray
                The full width at half maximum (FWHM) for the detector.
        """
        # scale overlay to detector dimensions
        if res is not None and res > 0:
            _res_scale = self.ydim/self.xdim
            _res_v = int(round(_res_scale * res, 0))
        else:
            res    = (self.det.hmp * self.det.hmn + self.det.hgp * (self.det.hmn-1) + self.det.cbh)
            _res_v = (self.det.vmp * self.det.vmn + self.det.vgp * (self.det.vmn-1) + self.det.cbh)

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
        azi = None
        if self.plo.show_unit_hover:
            # Distance POBI - pixel on grid
            R_a = np.sqrt(np.sum(_res[0:2]**2, axis=0)) * 1e-3 # m
            # POBI distance
            D_a = _res[2] * 1e-3 # m
            # 2theta - Angle between pixel, sample, and POBI
            tth = np.arctan2(R_a, D_a)
            # remove very small values (tth < 0.057 deg) to avoid zero divide
            tth[tth < 1e-3] = np.nan
            if self.plo.show_grid:
                # calculate the azimuthal angle eta
                azi = -np.arctan2(_res[0], _res[1])

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
        
        # delta d / d
        fwhm = 1.0
        #fwhm_width = 0.0
        if self.plo.show_fwhm:
            # dia: sample scattering diameter
            dia = self.plo.scattering_diameter
            # thk: detector sensor thickness
            thk = self.plo.sensor_thickness
            # mat: detector sensor material
            mat = self.plo.sensor_material
            # nrg: X-ray energy
            nrg = self.geo.ener
            # pxs: detector pixel size
            pxs = self.det.pxs * 1e-3
            # div: X-ray beam divergence
            div = self.plo.beam_divergence
            # dEE: X-ray energy resolution
            dEE = self.plo.energy_resolution
            # use the "unrotated" vector coordinates to define
            # R and D in the detector plane
            # radial component of the unrotated detector
            # i.e. radius away from the PONI
            r = np.sqrt(np.sum(_vec[0:2]**2, axis=0)) * 1e-3 # m
            # dis: PONI distance
            dis = self.geo.dist * 1e-3 # m
            # 2theta-alpha - Angle between pixel, sample, and PONI
            tth_a = np.arctan2(r, dis)
            # remove very small values (tth < 0.057 deg) to avoid zero divide
            tth_a[tth_a < 1e-3] = np.nan
            # H2, FWHM
            fwhm = self.calc_FWHM(dis, dia, thk, mat, pxs, tth_a, nrg, div, dEE)

        return grd, tth, azi, pc, sa, fwhm

    def calc_azi_grid(self, omega):
        """calculate the azimuthal grid points and return a dictionary with the vectors"""

        # First define vectors from the sample to the azimuthal
        # grid points for a vertical detector geometry. Then rotate
        # the vectors by the detector tilt angle. The azimuthal grid
        # points are then defined by the intersection of the rotated
        # vectors with the detector plane.

        _comp_shift = -(self.geo.voff - self.geo.dist * np.tan(omega) - np.deg2rad(self.geo.tilt) * self.geo.dist)
        # calculate the azimuthal grid points
        _azi = np.linspace(-np.pi, np.pi, self.plo.azimuth_num) # radians
        # calculate unit vectors to tth=5 deg
        _azi_vec = np.array([-np.sin(np.pi/36)*np.sin(_azi), # x
                                np.sin(np.pi/36)*np.cos(_azi), # y
                                np.cos(np.pi/36)*np.ones(_azi.shape[0]), # z
                                ]) 
        
        # scale the vectors such that z = sdd for all azimuthal grid points
        _azi_vec *=  self.geo.dist/_azi_vec[2,:]
        
        # rotate the vectors by the detector tilt+rotation angle
        _rot = self.rot_100(omega)
        _azi_vec = np.dot(_rot.T, _azi_vec)
        
        # rescale the vectors such that z = sdd for all azimuthal grid points
        # to make the points intersect the detector plane
        _azi_vec *=  self.geo.dist/np.abs(_azi_vec[2,:])
        # make sure the vectors are pointing in the right direction
        _azi_vec[2] = np.abs(_azi_vec[2])

        # account for the horizontal offset and the beam center shift
        _azi_vec[0] += self.geo.hoff
        _azi_vec[1] -= self.geo.voff

        grid_vectors = {}
        v0 = np.array([self.geo.hoff, _comp_shift])
        for i,a in enumerate(_azi):
            # find the vector in the detector plane from the beam center
            # to the azimuthal grid point
            v = np.array([_azi_vec[0,i], _azi_vec[1,i]])-np.array([self.geo.hoff, -(self.geo.voff - self.geo.dist * np.tan(omega))])
            # extend the vector to the edge of the detector including offsets
            scale = np.sqrt(self.xdim**2+self.ydim**2)+np.sqrt(self.geo.hoff**2+_comp_shift**2)
            v = v/np.linalg.norm(v)*scale
            # add the beam center shift
            v = np.array([self.geo.hoff, _comp_shift]) + v
            # store the vector
            grid_vectors[np.round(a*180/np.pi, 2)] = np.stack([v0, v])
        return grid_vectors

    def dsp2tth(self, dsp):
        """
        Converts d-spacing to 2-theta
         uses arcsin -> restricted to
         the interval -pi/2 - pi/2
        
        Returns
        -------
          np.arr, np.arr: tth values, valid indices
        """
        lambda_2d = (12.398/self.geo.ener) / (2*np.atleast_1d(dsp))
        idx = np.nonzero(lambda_2d < 1)
        tth = 2 * np.arcsin(lambda_2d[idx])
        return tth, idx

    #################
    #    UTILITY    #
    #  INDEPENDENT  #
    #################
    def get_att_lengths(self):
        # X-ray attenuation lengths z for Si and CdTe in meter [m] where z = ln(1/e)/mu
        # calculated in 1 keV steps from 1-150 keV
        # table values from Chantler (2000) https://doi.org/10.1063/1.1321055
        # calculated using the xraydb python module https://xraypy.github.io/XrayDB/
        self.att_lengths = {'Si':[2.92336338e-06, 1.52940940e-06, 4.41817449e-06, 9.68094927e-06,
                                  1.81408951e-05, 3.02711037e-05, 4.66249577e-05, 6.77473199e-05,
                                  9.52638566e-05, 1.30549774e-04, 1.73447018e-04, 2.24595659e-04,
                                  2.84605387e-04, 3.54031165e-04, 4.33380209e-04, 5.23140541e-04,
                                  6.23712068e-04, 7.35455743e-04, 8.58680607e-04, 9.93585970e-04,
                                  1.14037836e-03, 1.29916870e-03, 1.46997950e-03, 1.65282867e-03,
                                  1.84755219e-03, 2.05413078e-03, 2.27435657e-03, 2.51166980e-03,
                                  2.76076718e-03, 3.02154479e-03, 3.29339416e-03, 3.57595496e-03,
                                  3.86839569e-03, 4.17070154e-03, 4.48147168e-03, 4.80065791e-03,
                                  5.12767462e-03, 5.46137164e-03, 5.80102998e-03, 6.15160747e-03,
                                  6.50892811e-03, 6.87076056e-03, 7.23575505e-03, 7.60377581e-03,
                                  7.97432820e-03, 8.34584473e-03, 8.71928990e-03, 9.09192353e-03,
                                  9.46569452e-03, 9.83772142e-03, 1.02072666e-02, 1.05776192e-02,
                                  1.09433726e-02, 1.13091543e-02, 1.16704952e-02, 1.20265003e-02,
                                  1.23811631e-02, 1.27324096e-02, 1.30769437e-02, 1.34176491e-02,
                                  1.37548957e-02, 1.40888666e-02, 1.44155703e-02, 1.47396457e-02,
                                  1.50576116e-02, 1.53688649e-02, 1.56783348e-02, 1.59792892e-02,
                                  1.62758741e-02, 1.65681895e-02, 1.68560238e-02, 1.71398554e-02,
                                  1.74122248e-02, 1.76875177e-02, 1.79512540e-02, 1.82145339e-02,
                                  1.84671438e-02, 1.87236452e-02, 1.89677703e-02, 1.92079643e-02,
                                  1.94519568e-02, 1.96857927e-02, 1.99143713e-02, 2.01393629e-02,
                                  2.03608585e-02, 2.05788754e-02, 2.07833706e-02, 2.09941098e-02,
                                  2.12017351e-02, 2.14014048e-02, 2.15968609e-02, 2.17951736e-02,
                                  2.19821326e-02, 2.21717840e-02, 2.23511828e-02, 2.25366433e-02,
                                  2.27091166e-02, 2.28903416e-02, 2.30570644e-02, 2.32331077e-02,
                                  2.33946370e-02, 2.35541564e-02, 2.37235475e-02, 2.38780490e-02,
                                  2.40302885e-02, 2.41803970e-02, 2.43399625e-02, 2.44882422e-02,
                                  2.46323882e-02, 2.47746094e-02, 2.49149003e-02, 2.50533840e-02,
                                  2.51901681e-02, 2.53252522e-02, 2.54586507e-02, 2.55903974e-02,
                                  2.57204103e-02, 2.58490340e-02, 2.59762024e-02, 2.61018936e-02,
                                  2.62261777e-02, 2.63492097e-02, 2.64707961e-02, 2.65910939e-02,
                                  2.67101422e-02, 2.68280390e-02, 2.69448830e-02, 2.70602934e-02,
                                  2.71637201e-02, 2.72706288e-02, 2.73828591e-02, 2.74939145e-02,
                                  2.76039513e-02, 2.77130685e-02, 2.78212863e-02, 2.79285825e-02,
                                  2.80224089e-02, 2.81218802e-02, 2.82264203e-02, 2.83300579e-02,
                                  2.84220811e-02, 2.85159760e-02, 2.86172018e-02, 2.87175212e-02,
                                  2.88171887e-02, 2.89147878e-02, 2.89995540e-02, 2.90922933e-02,
                                  2.91891453e-02, 2.92852862e-02],
                          'CdTe':[8.86140064e-08, 4.41597213e-07, 1.15951587e-06, 8.36374589e-07,
                                  8.36034310e-07, 1.32724815e-06, 1.97362567e-06, 2.79221097e-06,
                                  3.79736962e-06, 5.07172559e-06, 6.58646052e-06, 8.34536221e-06,
                                  1.03633319e-05, 1.26445947e-05, 1.51992203e-05, 1.80585164e-05,
                                  2.12380712e-05, 2.47559111e-05, 2.86220163e-05, 3.28393095e-05,
                                  3.74211141e-05, 4.23847110e-05, 4.77396937e-05, 5.34604046e-05,
                                  5.95357580e-05, 6.58583493e-05, 2.02635396e-05, 2.23180704e-05,
                                  2.45073128e-05, 2.68096054e-05, 2.92163420e-05, 1.99269432e-05,
                                  2.15332477e-05, 2.32726190e-05, 2.51068287e-05, 2.70499442e-05,
                                  2.90850735e-05, 3.12143579e-05, 3.34390888e-05, 3.57606328e-05,
                                  3.81810110e-05, 4.07109020e-05, 4.33566740e-05, 4.61081874e-05,
                                  4.89663832e-05, 5.19342068e-05, 5.50092594e-05, 5.81952784e-05,
                                  6.14950533e-05, 6.49069814e-05, 6.84326118e-05, 7.20745366e-05,
                                  7.58320671e-05, 7.97392476e-05, 8.37930712e-05, 8.79696403e-05,
                                  9.22722203e-05, 9.66991799e-05, 1.01254296e-04, 1.05935695e-04,
                                  1.10745003e-04, 1.15686534e-04, 1.20756232e-04, 1.25977939e-04,
                                  1.31365089e-04, 1.36890012e-04, 1.42551567e-04, 1.48351972e-04,
                                  1.54290824e-04, 1.60366939e-04, 1.66584684e-04, 1.72943591e-04,
                                  1.79440083e-04, 1.86083525e-04, 1.92866175e-04, 1.99789869e-04,
                                  2.06859516e-04, 2.14074597e-04, 2.21432043e-04, 2.28971821e-04,
                                  2.36785318e-04, 2.44789163e-04, 2.52956288e-04, 2.61270529e-04,
                                  2.69748772e-04, 2.78374465e-04, 2.87160014e-04, 2.96098601e-04,
                                  3.05197594e-04, 3.14457023e-04, 3.23869498e-04, 3.33440320e-04,
                                  3.43178015e-04, 3.53053462e-04, 3.63091886e-04, 3.73283283e-04,
                                  3.83651351e-04, 3.94156444e-04, 4.04836771e-04, 4.15653703e-04,
                                  4.26625960e-04, 4.37759844e-04, 4.49053794e-04, 4.60495597e-04,
                                  4.72077099e-04, 4.83837118e-04, 4.95744109e-04, 5.07781340e-04,
                                  5.19964672e-04, 5.32316813e-04, 5.44823489e-04, 5.57483665e-04,
                                  5.70297686e-04, 5.83249863e-04, 5.96338131e-04, 6.09594120e-04,
                                  6.23013435e-04, 6.36594830e-04, 6.50297249e-04, 6.64142222e-04,
                                  6.78147980e-04, 6.92282395e-04, 7.06561344e-04, 7.20987946e-04,
                                  7.35548731e-04, 7.50255560e-04, 7.65122247e-04, 7.80063756e-04,
                                  7.95147555e-04, 8.10396744e-04, 8.25793330e-04, 8.41276246e-04,
                                  8.56872987e-04, 8.72612687e-04, 8.88496254e-04, 9.04519016e-04,
                                  9.20619462e-04, 9.36881518e-04, 9.53292628e-04, 9.69782092e-04,
                                  9.86397315e-04, 1.00312754e-03, 1.01996070e-03, 1.03689757e-03,
                                  1.05395590e-03, 1.07113920e-03, 1.08845123e-03, 1.10583974e-03,
                                  1.12333211e-03, 1.14094222e-03]}

    def rot_100(self, a, cc=1):
        """
        Generate a rotation matrix for a rotation around the [100] axis.

        Parameters:
        a (float): The angle of rotation in radians.
        cc (int, optional): Clockwise rotation if 1 (default), counterclockwise if 0.

        Returns:
        numpy.ndarray: A 3x3 rotation matrix.
        """
        #Omega in radians
        ca = np.cos(a)
        sa = np.sin(a)
        if cc: sa = -sa
        return np.array([[1,   0,  0],
                         [0,  ca, sa],
                         [0, -sa, ca]])

    def calc_hkld(self, ucp, res=0.2e-10, dec=4, cen='P'):
        """
        Generate hkls for unit cell parameters
        Calculate resolution in d-spacing (dsp)
        Round dsp to decimals, only use unique numbers

        :param ucp: Unit cell parameters, list, [a, b, c, alpha, beta, gamma]
        :param res: Maximum resolution
        :param dec: d-spacing sampling rate, remove multiplicity
        :param cen: Remove systematic absences according to centring
        
        :return: array of [[h k l dsp]]
        """
        def cart_from_cell(cell):
            """
            Convert lattice constants to Cartesian coordinates.

            Parameters:
            cell (numpy.ndarray): A 1D array with 6 elements representing the lattice constants.
                          The first three elements are the lengths of the cell edges (a, b, c)
                          in angstroms, and the last three elements are the angles (alpha, beta, gamma)
                          in degrees.

            Returns:
            tuple: A tuple containing three numpy arrays (av, bv, cv) representing the Cartesian coordinates
                   of the cell vectors.

            Raises:
            ValueError: If the input array does not have exactly 6 elements.
            """
            if cell.shape != (6,):
                raise ValueError('Lattice constants must be 1d array with 6 elements')
            a, b, c = cell[:3]*1E-10
            alpha, beta, gamma = np.radians(cell[3:])
            av = np.array([a, 0, 0], dtype=float)
            bv = np.array([b * np.cos(gamma), b * np.sin(gamma), 0], dtype=float)
            # calculate vector c
            x = np.cos(beta)
            y = (np.cos(alpha) - x * np.cos(gamma)) / np.sin(gamma)
            z = np.sqrt(1. - x**2. - y**2.)
            cv = np.array([x, y, z], dtype=float)
            cv /= np.linalg.norm(cv)
            cv *= c
            return av, bv, cv
        
        def matrix_from_cell(cell):
            """
            Generate a transformation matrix from a given cell.

            Parameters:
            cell (array-like): A 3x3 array representing the cell parameters.

            Returns:
            numpy.ndarray: A 3x3 transformation matrix rounded to 6 decimal places.

            Notes:
            - The function first converts the input cell to a numpy array.
            - It then calculates the Cartesian vectors from the cell parameters.
            - The reciprocal lattice vectors (a*, b*, c*) are computed using cross products and dot products.
            - These vectors are assembled into a 3x3 transformation matrix.
            """
            cell = np.array(cell)
            av, bv, cv = cart_from_cell(cell)
            a_star = (np.cross(bv, cv)) / (np.cross(bv, cv).dot(av))
            b_star = (np.cross(cv, av)) / (np.cross(cv, av).dot(bv))
            c_star = (np.cross(av, bv)) / (np.cross(av, bv).dot(cv))
            A = np.zeros((3, 3), dtype='float')  # transform matrix
            A[:, 0] = a_star
            A[:, 1] = b_star
            A[:, 2] = c_star
            return np.round(A,6)

        def applyExtinctionRules(hkl, centering='P'):
            """
            Apply extinction rules based on the centering type to filter out certain Miller indices (hkl).
            Parameters:
            hkl (array-like): An array of Miller indices to be filtered.
            centering (str, optional): The centering type of the crystal lattice. 
                           Default is 'P'. Options are:
                           - 'P': Primitive
                           - 'I': Body-centered
                           - 'A': A-face centered
                           - 'B': B-face centered
                           - 'C': C-face centered
                           - 'F': Face-centered
            Returns:
            numpy.ndarray: Filtered array of Miller indices after applying the extinction rules.
            """
            hkl = np.atleast_2d(hkl)

            if centering == 'P':
                pass
        
            elif centering == 'I':
                # h+k+l = even
                hkl = hkl[np.sum(hkl, axis=1)%2 == 0]
        
            elif centering == 'A':
                # k + l = even
                hkl = hkl[np.sum(hkl[:,1:], axis=1)%2 == 0]
        
            elif centering == 'B':
                # h + l = even
                hkl = hkl[np.sum(hkl[:,[0,2]], axis=1)%2 == 0]
        
            elif centering == 'C':
                # h + k = even
                hkl = hkl[np.sum(hkl[:,:2], axis=1)%2 == 0]
        
            elif centering == 'F':
                # h, k, l all odd or all even
                hkl = hkl[np.sum(hkl%2 == 0, axis=1)%3 == 0]
        
            return hkl
        
        A = matrix_from_cell(ucp)
        q_cutoff = 1. / res
        max_h = min(int(np.ceil(q_cutoff / np.linalg.norm(A[:,0]))), 127)
        max_k = min(int(np.ceil(q_cutoff / np.linalg.norm(A[:,1]))), 127)
        max_l = min(int(np.ceil(q_cutoff / np.linalg.norm(A[:,2]))), 127)
        # hkl grid
        #hh = np.arange(max_h, -max_h-1, -1)
        hh = np.arange(    0, -max_h-1, -1)
        kk = np.arange(max_k, -max_k-1, -1)
        ll = np.arange(max_l, -max_l-1, -1)
        # this determines (for no obvious reason)
        # the order of the array
        ks, hs, ls = np.meshgrid(kk, hh, ll)
        hkl = np.ones((hs.size, 3), dtype=np.int8)
        hkl[:,0] = hs.reshape(-1)
        hkl[:,1] = ks.reshape(-1)
        hkl[:,2] = ls.reshape(-1)
        # remove 0 0 0 reflection
        #hkl = np.delete(hkl, len(hkl)//2, 0)
        rem = max_k * (2*max_l+1) + max_l
        hkl = np.delete(hkl, rem, axis=0)
        # remove high resolution hkls
        # too expensive to be useful
        #hkl = hkl[(np.linalg.norm(A.dot(hkl.T).T, axis=1) <= q_cutoff)]
        # apply systematic absences from centring
        hkl = applyExtinctionRules(hkl, centering=cen)
        # calculate the d-spacing
        # go from meters to Angstrom
        # cast to int to speed up the next step
        #  -> np.unique sorting
        dsp = ((1/np.linalg.norm(A.dot(hkl.T).T, axis=1))*10**(10+dec)).astype(np.uint32)
        # get reduced indices
        dsp, idx = np.unique(dsp, return_index=True)
        # stack the hkl and dsp
        out = np.hstack([hkl[idx], (dsp*10**(-dec)).reshape(-1,1)])[::-1]
        return out

    def calc_FWHM(self, dis, dia, thk, mat, pxs, tth, nrg, div, dEE, deg=True):
            """
            Calculate FWHM

            Parameters
            ----------
            dis: poni distance
            dia: sample scattering diameter
            thk: detector sensor thickness
            mat: detector sensor material
            pix: detector pixel size
            tth: 2-theta angle
            nrg: X-ray energy
            div: X-ray beam divergence
            dEE: X-ray energy resolution
            thk: detector sensor thickness
            deg: return degrees (True) or radians (False)
            
            Returns
            -------
            array: fwhm
            """
            if len(self.att_lengths[mat]) > int(nrg):
                thk = min(thk, self.att_lengths[mat][int(nrg)])
            A = 2*np.log(2) / dis**2 * (pxs**2-2*thk**2-dia**2)
            B = 2*np.log(2) / dis**2 * (2*thk**2 + 2*dia**2)
            C = 2*np.log(2) * div**2
            M = (4*np.sqrt(2*np.log(2)) * dEE)**2 * ((1-np.cos(tth))/(1+np.cos(tth)))
            X = np.cos(tth)
            H2 = A*X**4 + B*X**2 + C + M
            fwhm = np.sqrt(H2)
            if deg is True:
                return fwhm * 180 / np.pi
            else:
                return fwhm

    def gaussian(self, x, m, s):
        """
        Calculate the value of a Gaussian function.

        Parameters:
        x (float or array-like): The input value(s) where the Gaussian function is evaluated.
        m (float): The mean (center) of the Gaussian distribution.
        s (float): The standard deviation (spread or width) of the Gaussian distribution.

        Returns:
        float or array-like: The value(s) of the Gaussian function at the given input value(s).
        """
        return 1/(np.sqrt(2*np.pi)*s)*np.exp(-np.square((x - m)/s)/2)

    def get_closest_point_x(self, points, pos):
        """
        Find the point in a list of points that is closest to a given position.

        Args:
            points (list): A list of objects, each having a method `pos()` that returns an object with method `x()`.
            pos (object): An object with a method `x()` representing the position to compare against.

        Returns:
            object: The point from the list that is closest to the given position.
        """
        dist = [np.linalg.norm(np.array(p.pos().x()) - np.array(pos.x())) for p in points]
        return points[np.argmin(dist)]

    #############
    # SETTINGS  #
    #############
    def settings_get_active(self):
        """
        Retrieves the active settings file name.

        This method checks if a settings token file exists. If it does, it reads the
        token file to get the name of the active settings file and verifies its existence.
        If the active settings file is found, its name is returned.

        If the settings token file does not exist or the active settings file is not found,
        the method resets to default settings, saves the default settings to a file, and
        returns the name of the default settings file.

        Returns:
            str: The name of the active settings file or the default settings file.
        """
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
        self.settings_save_to_file()
        return self.active_settings
    
    def settings_del_active(self):
        """
        Deletes the settings token file if it exists.

        This method checks if the settings token file, specified by the 
        attribute `self.path_settings_token`, exists. If the file exists, 
        it removes the file from the filesystem.

        Returns:
            None
        """
        if os.path.exists(self.path_settings_token):
            os.remove(self.path_settings_token)

    def settings_set_active(self):
        """
        Sets the active settings by creating the necessary directories and writing the active settings to a token file.

        This method performs the following steps:
        1. Checks if the directory specified by `self.path_settings` exists. If not, it creates the directory.
        2. Opens the file specified by `self.path_settings_token` in write mode.
        3. Writes the active settings (stored in `self.active_settings`) to the file.

        Raises:
            OSError: If there is an issue creating the directory or writing to the file.
        """
        if not os.path.exists(self.path_settings):
            os.makedirs(self.path_settings)
        with open(self.path_settings_token, 'w') as wf:
            wf.write(self.active_settings)

    def settings_delete_file(self):
        """
        Opens a file dialog to select and delete settings files.

        This method changes the file dialog button from 'open' to 'delete' and allows the user to select
        multiple settings files (with a .json extension) to delete. If files are selected, they are removed
        from the filesystem and corresponding actions are removed from the custom settings menu.

        Returns:
            None
        """
        # todo change button from 'open' to 'delete'
        fnames, filter = QtWidgets.QFileDialog.getOpenFileNames(self, 'Delete settings files', self.path_settings_current, "Settings files (*.json)")
        if fnames:
            for fname in fnames:
                os.remove(fname)
                for action in self.menu_custom_settings.actions():
                    if action.text() == os.path.basename(fname):
                        self.menu_custom_settings.removeAction(action)

    def settings_save_current(self):
        """
        Save the current geometry settings.

        This method updates the initial geometry settings (`_geo`) with the 
        current geometry settings (`geo`). It is typically used when the 
        current settings need to be saved as the new startup values. After 
        updating `_geo`, the settings are saved to a file by calling 
        `settings_save_to_file`.
        """
        # self.geo is edited with the sliders
        # self._geo holds the initial values
        # usually I do not want to overwrite 
        # the startup values -> I write _geo
        # to the settings file.
        # unless this function is called!
        self._geo.__dict__.update(self.geo.__dict__)
        self.settings_save_to_file()

    def settings_save_to_file(self, target=None):
        """
        Save the current settings to a file.

        If no target file is specified, the settings will be saved to the default location.

        Args:
            target (str, optional): The file path where the settings should be saved. 
                        If None, the settings will be saved to the default location.

        Returns:
            None

        Notes:
            - If the target directory does not exist, it will be created.
            - If the target file is not writable, the function will return without saving.
            - The settings are saved in JSON format, including 'geo', 'plo', 'thm', and 'lmt' attributes.
        """
        if target == None:
            target = os.path.join(self.path_settings, self.active_settings)
        target_base = os.path.dirname(target)
        # create folder if not existing
        if not os.path.exists(target_base):
            os.makedirs(target_base)
        # Don't overwrite beamlime standard settings
        if not os.access(target, os.W_OK):
            #print(f'{os.path.basename(target)} is protected.')
            return
        # Writing geo as dict to file
        with open(target, 'w') as wf:
            json.dump({'geo':self._geo.__dict__, 'plo':self.plo.__dict__, 'thm':self.thm.__dict__, 'lmt':self.lmt.__dict__}, wf, indent=4)

    def settings_load_from_file(self, skip=[]):
        """
        Load settings from a JSON file and apply them to the object's attributes.

        This method reads a JSON file specified by `self.path_settings_current` and updates
        the object's attributes based on the contents of the file. It performs checks to ensure
        that certain parameters are within specified ranges and enforces specific values for
        some parameters.

        Parameters:
        -----------
        skip : list, optional
            A list of keys to skip when loading settings. Default is an empty list.

        Raises:
        -------
        SystemExit
            If there is an error parsing the settings file.

        Notes:
        ------
        - The method uses two dictionaries, `_warn` and `_force`, to manage parameter constraints:
            - `_warn` contains parameters with their allowed minimum, maximum, and default values.
            - `_force` contains parameters that should be set to specific values regardless of the file content.
        - If a parameter value is outside the allowed range specified in `_warn`, it is set to the default value,
          and a warning message is printed.
        - If a parameter is in `_force`, it is set to the enforced value.
        - The method updates the attributes of `self.geo`, `self.plo`, `self.thm`, and `self.lmt` based on the
          contents of the JSON file.
        - If 'geo' is not in the `skip` list, the initial values of `self.geo` are stored in `self._geo`.
        """
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
        _force = {'show_fwhm':False}
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

    def settings_edit_file(self, command):
        """
        Executes a system command to edit the current settings file.

        Args:
            command (str): The command to execute, typically an editor command.
        """
        os.system(f'{command} {self.path_settings_current}')

    def settings_get_files(self):
        """
        Retrieve a sorted list of JSON filenames from the settings directory.

        This method searches the directory specified by `self.path_settings` for 
        files with a `.json` extension, extracts their base names, and returns 
        them in a sorted list.

        Returns:
            list: A sorted list of JSON filenames found in the settings directory.
        """
        return sorted(map(os.path.basename, glob.glob(os.path.join(self.path_settings, '*.json'))))

    def settings_change_file(self, name):
        """
        Change the active settings file and reload settings.

        This method updates the active settings file to the specified name,
        constructs the full path to the new settings file, and reloads the settings.

        Args:
            name (str): The name of the new settings file to activate.
        """
        self.active_settings = name
        self.path_settings_current = os.path.join(self.path_settings, name)
        self.settings_reload()
    
    def settings_import_win(self):
        """
        Opens a file dialog to import a settings file in JSON format. If a file is selected,
        it copies the file to the settings folder, adds it to the custom settings menu, 
        and changes the current settings to the newly imported file.

        Steps:
        1. Opens a file dialog to select a JSON settings file.
        2. Copies the selected file to the settings folder.
        3. Adds the file to the custom settings menu as a checkable action.
        4. Changes the current settings to the newly imported file and reloads the settings.

        Returns:
            None
        """
        fname, filter = QtWidgets.QFileDialog.getOpenFileName(self, 'Import settings file', '', "Settings files (*.json)")
        if fname:
            # copy to settings folder
            shutil.copy(fname, self.path_settings)
            bname = os.path.basename(fname)
            # Add to menu
            cset_action = QtGui.QAction(bname, self, checkable=True)
            self.menu_set_action(cset_action, self.settings_change_file, bname)
            self.menu_custom_settings.addAction(cset_action)
            self.group_cset.addAction(cset_action)
            cset_action.setChecked(True)
            # change settings and reload
            self.settings_change_file(bname)

    def settings_reload(self):
        """
        Reloads the settings for the application.

        This method performs the following actions:
        1. Clears the `det_bank` to ensure backward compatibility with older parameter files 
           that may not have the `det_bank` entry, allowing access to all detectors.
        2. Loads settings from a file.
        3. Saves the settings back to the file, adding any missing entries.
        4. Redraws the canvas to reflect the updated settings.
        5. Reinitializes the menu with a reset.

        Note: The `update_menu_entries` method call is commented out.
        """
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
        self.settings_load_from_file()
        # missing entries in the
        # settings file will be added
        self.settings_save_to_file()
        self.redraw_canvas()
        self.menu_init(reset=True)

    #############
    #   EVENT   #
    #############
    def dragEnterEvent(self, event):
        """
        Handles the drag enter event for the widget.

        This method is triggered when a drag-and-drop action enters the widget.
        It checks if the dragged data contains URLs (e.g., files) and accepts
        the event if it does. Otherwise, it ignores the event.

        Args:
            event (QDragEnterEvent): The drag enter event containing information
                                     about the dragged data.
        """
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        """
        Handle the drop event for a drag-and-drop operation.

        This method processes the dropped file, checks if it is a CIF file,
        and if so, calls the method to calculate reference data from the CIF file.

        Args:
            event (QDropEvent): The drop event containing the dropped file.

        Returns:
            None
        """
        fpath = event.mimeData().urls()[0].toLocalFile()
        if os.path.splitext(fpath)[1] == '.cif':
            self.calc_ref_from_cif(fpath, open_pxrd=True)

    def keyReleaseEvent(self, event):
        """
        Handles key press events and executes corresponding functions based on 
        key and modifier combinations defined in the hotkey_dict.
        Args:
            event (QKeyEvent): The key event that triggered this method.
        The method checks the key and modifier combination of the event and 
        looks it up in the hotkey_dict. If a matching combination is found, 
        the associated function is executed with or without arguments as defined 
        in the hotkey_dict.
        The modifiers are mapped as follows:
            - NoModifier or KeypadModifier: None
            - AltModifier: 'ALT'
            - ShiftModifier: 'SHIFT'
            - ControlModifier: 'CTRL'
        """
        k = event.key()
        m = event.modifiers()
        
        # set m to None if no modifier or only keypad modifier is pressed
        if m == QtCore.Qt.KeyboardModifier.NoModifier or m == QtCore.Qt.KeyboardModifier.KeypadModifier:
            m = None
        elif m == QtCore.Qt.KeyboardModifier.AltModifier:
            m = 'ALT'
        elif m == QtCore.Qt.KeyboardModifier.ShiftModifier:
            m = 'SHIFT'
        elif m == QtCore.Qt.KeyboardModifier.ControlModifier:
            m = 'CTRL'
        
        # check if key and modifier combination exists in hotkey_dict
        k = QtGui.QKeySequence(k).toString()
        if (k,m) in self.hotkey_dict:
            fn, arg = self.hotkey_dict.get((k,m))['fn']
            if arg is not None:
                fn(arg)
            else:
                fn()

    def closeEvent(self, event):
        """
        Handles the close event for the application window.

        This method is called when the application window is closed. It saves the 
        current settings to ensure they are automatically loaded on the next startup. 
        If the program crashes, no active settings token exists, and the default 
        settings are loaded.

        Args:
            event (QCloseEvent): The close event triggered when the window is closed.
        """
        self.settings_set_active()
        event.accept()

    def hoverEvent(self, event):
        """
        Handles hover events linked to the cormap and updates labels accordingly.
        This method is triggered when a hover event occurs and updates the 
        `cor_label` and `unit_label` based on the event's position and the 
        current state of the map overlays.
        Parameters:
        event (QHoverEvent): The hover event containing information about the 
                             cursor's position and state.
        Behavior:
        - Hides `cor_label` and updates `unit_label` when the cursor exits the area.
        - Shows `cor_label` when the cursor enters the area and any of the map 
          overlays are active.
        - Updates `unit_label` with the calculated unit value if `_tth` is not a float.
        - Updates `cor_label` with the correction values from `_polcor`, `_solang`, 
          and `_fwhm` if they are not floats.
        """
        # Hover event linked to cormap and should only be active while
        # either or both maps are displayed
        if event.isExit():
            self.cor_label.hide()
            self.unit_label.setText(f'{self.unit_names[self.geo.unit]}')
            return
        
        # hoverEvent is only called if
        # either or both maps are active
        # -> always show on isEnter event
        if event.isEnter() and (self.action_funct_fwhm_show.isChecked() or self.action_show_ang.isChecked() or self.action_show_pol.isChecked()):
            self.cor_label.show()

        # cormap displays the product of both corrections
        # but the individual values can be retrieves from their arrays
        # -> _polcor and _solang
        y, x = map(int, np.clip(event.pos(), [0,0], np.array(self.patches['overlay'].image.shape)[::-1]-1))

        # unit label value
        if not isinstance(self._tth, float):
            unit = self.calc_unit(self._tth[x,y])
            if self.plo.show_grid:
                self.unit_label.setText(f'{self.unit_names[self.geo.unit]} {unit:.2f}\nazi [\u00B0] {self._azi[x,y]*180/np.pi:.0f}')
            else:
                self.unit_label.setText(f'{self.unit_names[self.geo.unit]} {unit:.2f}')
        
        _text = []
        # calc_overlays returns either a np.array (if active)
        # or a float (inactive) for _polcor and _solang. 
        if not isinstance(self._polcor, float):
            _text.append(f'P: {self._polcor[x,y]:.2f}')
        if not isinstance(self._solang, float):
            _text.append(f'S: {self._solang[x,y]:.2f}')
        if not isinstance(self._fwhm, float):
            _text.append(f'H: {self._fwhm[x,y]:.4f}\u00B0')
            #_text.append(f'H: {self.calc_unit(self._fwhm[x,y]):.4f}')
        self.cor_label.setText('\n'.join(_text))

############
#   MISC   #
############
class Container(object):
    """
    A class used to represent a Container.
    """
    pass

class Ref(object):
    """
    Ref is a class that represents a reference object with attributes for name, dsp, hkl, and cif.
    """
    def __init__(self, name, dsp=None, hkl=None, cif=None):
        """
        Initialize a new instance of the class.

        Parameters:
        name (str): The name of the instance.
        dsp (optional): The dsp value. Default is None.
        hkl (optional): The hkl value. Default is None.
        cif (optional): The cif value. Default is None.
        """
        self.name = name
        self.dsp = dsp
        self.hkl = hkl
        self.cif = cif
        self.has_dsp = False
        self.has_hkl = False
        self.has_cif = False
        # set falgs
        if dsp is not None:
            self.has_dsp = True
        if hkl is not None:
            self.has_hkl = True
        if cif is not None:
            self.has_cif = True
        
    def add_cif(self, cif):
        """
        Adds a CIF (Crystallographic Information File) to the instance.

        Parameters:
        cif (str): The CIF data to be added.

        Sets the instance's `cif` attribute to the provided CIF data and 
        marks `has_cif` as True.
        """
        self.cif = cif
        self.has_cif = True

    def add_dsp(self, dsp):
        """
        Adds a DSP (Digital Signal Processor) to the instance.

        Parameters:
        dsp (object): The DSP object to be added.

        Sets:
        self.dsp: The provided DSP object.
        self.has_dsp (bool): Flag indicating that a DSP has been added, set to True.
        """
        self.dsp = dsp
        self.has_dsp = True

    def add_hkl(self, hkl):
        """
        Adds the Miller indices (hkl) to the instance and sets the has_hkl flag to True.

        Parameters:
        hkl (tuple): A tuple representing the Miller indices (h, k, l).
        """
        self.hkl = hkl
        self.has_hkl = True
    
    def rem_cif(self):
        """
        Remove the CIF (Crystallographic Information File) from the object.

        This method sets the `cif` attribute to None and the `has_cif` attribute to False,
        indicating that the object no longer has an associated CIF file.
        """
        self.cif = None
        self.has_cif = False

    def rem_hkl(self):
        """
        Remove the Miller indices (hkl) from the object.

        This method sets the `hkl` attribute to None and the `has_hkl` attribute to False,
        indicating that the object no longer has Miller indices assigned.
        """
        self.hkl = None
        self.has_hkl = False

    def rem_dsp(self):
        """
        Remove the dsp attribute from the instance.

        This method sets the `dsp` attribute to None and the `has_dsp` attribute to False.
        """
        self.dsp = None
        self.has_dsp = False
    
    def validate_cif(self):
        """
        Validates the existence of the CIF (Crystallographic Information File).

        This method checks if the CIF file specified by the `self.cif` attribute exists.
        If the file does not exist, it calls the `rem_cif` method to handle the missing file.

        Raises:
            FileNotFoundError: If the CIF file does not exist.
        """
        if not os.path.exists(self.cif):
            self.rem_cif()
    
    def is_complete(self):
        """
        Check if the object is complete.

        This method checks if all required attributes (`has_cif`, `has_dsp`, `has_hkl`) 
        are present and returns True if they are all present, otherwise False.

        Returns:
            bool: True if all required attributes are present, False otherwise.
        """
        return np.all([self.has_cif, self.has_dsp, self.has_hkl])

###############
#   WIDGETS   #
###############
class SliderWidget(QtWidgets.QFrame):
    """
    SliderWidget is a custom Qt widget that provides a frame with a set of vertical sliders,
    each with a label and value display. The widget can be toggled to show or hide the sliders
    and can be dragged within its parent window.
    """
    def __init__(self, parent):
        """
        Initializes the custom widget with the given parent.

        Args:
            parent (QWidget): The parent widget.

        Attributes:
            layout (QVBoxLayout): The main layout for the widget.
            leaveEvent (function): Event handler for when the mouse leaves the widget.
            enterEvent (function): Event handler for when the mouse enters the widget.
            frame (QFrame): A frame widget used within the layout.
            box (QFrame): A frame widget that can be toggled to show or hide.
            box_toggle (bool): A flag indicating the toggle state of the box.
            box_height_show (int): The height of the box when shown.
            box_height_hide (int): The height of the box when hidden.
            startPos (QPoint or None): The starting position for preventing window movement.
            grid (QGridLayout): The grid layout for the box frame.
        """
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
        """
        Applies the style to the box and frame widgets using the parent's style properties.

        This method sets the style sheets for the `self.box` and `self.frame` widgets based on the
        parent's properties such as `slider_border_width`, `slider_border_color`, `slider_border_radius`,
        `slider_bg_color`, and `slider_bg_hover`.

        The style includes:
        - Border width and color
        - Border radius
        - Background color
        - Hover background color for the box

        The styles are applied using the `setStyleSheet` method with formatted strings.

        Raises:
            AttributeError: If any of the required parent properties are missing.
        """
        self.box.setStyleSheet(f'''
            QFrame {{
                border: {self.parent().plo.slider_border_width}px solid {self.parent().slider_border_color.name(format=QtGui.QColor.NameFormat.HexArgb)};
                border-radius: {self.parent().plo.slider_border_radius}px;
                background: {self.parent().slider_bg_color.name(format=QtGui.QColor.NameFormat.HexArgb)};
            }}
            QFrame:hover {{
                background: {self.parent().slider_bg_hover.name(format=QtGui.QColor.NameFormat.HexArgb)};
            }}
        ''')
        self.frame.setStyleSheet(f'''
            QFrame {{
                border: {self.parent().plo.slider_border_width}px solid {self.parent().slider_border_color.name(format=QtGui.QColor.NameFormat.HexArgb)};
                border-radius: {self.parent().plo.slider_border_radius}px;
                background: {self.parent().slider_border_color.name(format=QtGui.QColor.NameFormat.HexArgb)};
            }}
        ''')

    def init_sliders(self):
        """
        Initializes and updates the sliders in the UI grid layout.
        This method performs the following steps:
        1. Removes all existing sliders and labels from the grid layout.
        2. Adds new sliders based on the enabled slider options in the parent object.
        3. Updates the dynamic width of the box based on the number of sliders added.
        4. Sets the maximum limit for the beamstop distance slider to the detector distance.
        Sliders are added for the following parameters if enabled:
        - Energy ('ener')
        - Distance ('dist')
        - Vertical offset ('voff')
        - Horizontal offset ('hoff')
        - Tilt ('tilt')
        - Rotation ('rota')
        - Beamstop distance ('bsdx')
        The method also updates the slider limits for the beamstop distance slider.
        Attributes:
        - self.grid: The grid layout where sliders are added.
        - self.box_width_dynamic: The dynamic width of the box based on the number of sliders.
        - self.box_height_hide: The height of the box when sliders are hidden.
        - self.parent(): The parent object containing slider configurations and limits.
        """
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
        """
        Centers the frame within its parent window.

        This method calculates the horizontal position required to center the frame
        within its parent window and moves the frame to that position. The vertical
        position is set to the parent's offset for Win32 systems.

        Returns:
            None
        """
        self.move(int((self.parent().size().width()-self.box_width_dynamic)/2), self.parent().offset_win32)

    def update_slider_label(self, label, value):
        """
        Updates the text of a given label to reflect the provided slider value.

        Args:
            label (QLabel): The label whose text needs to be updated.
            value (float): The value from the slider to set as the label's text.

        Returns:
            None
        """
        label.setText(str(int(value)))

    def update_slider_limits(self, slider, lmin, lmax):
        """
        Update the limits of a given slider.

        Parameters:
        slider (QSlider): The slider object whose limits are to be updated.
        lmin (float): The minimum limit to set for the slider.
        lmax (float): The maximum limit to set for the slider.
        """
        slider.setRange(int(lmin), int(lmax))

    def add_slider(self, layout, label, token, idx, lval, lmin, lmax, lstp):
        """
        Adds a vertical slider with a label and value display to the specified layout.
        Args:
            layout (QtWidgets.QGridLayout): The layout to which the slider and labels will be added.
            label (str): The text for the slider's label.
            token (str): The object name for the slider.
            idx (int): The column index in the layout where the slider and labels will be placed.
            lval (float): The initial value of the slider.
            lmin (float): The minimum value of the slider.
            lmax (float): The maximum value of the slider.
            lstp (float): The step size for the slider.
        Returns:
            QtWidgets.QSlider: The created slider widget.
        """
        font = QtGui.QFont()
        font.setPixelSize(self.parent().plo.slider_label_size)

        label_name = QtWidgets.QLabel(label)
        label_name.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        label_name.setFont(font)
        label_name.setStyleSheet(f'''
            QLabel {{
                color: {self.parent().slider_label_color.name(format=QtGui.QColor.NameFormat.HexArgb)};
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
                color: {self.parent().slider_label_color.name(format=QtGui.QColor.NameFormat.HexArgb)};
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
        """
        Toggles the visibility of a panel based on the type of event received.

        Args:
            event (QEvent): The event that triggers the toggle. If the event is of type
                            QtGui.QEnterEvent, the panel will be shown. If the event is of type
                            QtCore.QEvent and the panel is not already toggled, the panel will be hidden.

        Behavior:
            - If the event is a QtGui.QEnterEvent, the panel is shown and resized to `self.box_width_dynamic` 
              and `self.box_height_show`.
            - If the event is a QtCore.QEvent and `self.box_toggle` is False, the panel is hidden and resized 
              to `self.box_width_dynamic` and `self.box_height_hide`.
            - For other events, no action is taken.
        """
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
        """
        Handles the mouse press event.

        This method is called when a mouse button is pressed. It first calls the 
        parent class's mousePressEvent method to ensure any default behavior is 
        executed. If the left mouse button is pressed, it records the position of 
        the mouse press and toggles the state of the box_toggle attribute.

        Args:
            event (QMouseEvent): The event object containing details about the 
                                 mouse press event.
        """
        super().mousePressEvent(event)
        if event.button() == QtCore.Qt.MouseButton.LeftButton:
            self.startPos = event.pos()
            self.box_toggle = not self.box_toggle

    def mouseMoveEvent(self, event):
        """
        Handles the mouse move event for the widget.

        This method is called when the mouse is moved over the widget. It processes
        the event to enable dragging of the widget when the left mouse button is held
        down. The widget's position is updated based on the mouse movement, while
        ensuring it stays within the bounds of its parent widget.

        Args:
            event (QtGui.QMouseEvent): The mouse event containing information about
                                       the mouse movement.

        Returns:
            None
        """
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
        """
        Handles the mouse release event.

        This method is called when the mouse button is released. It calls the
        parent class's mouseReleaseEvent method and then resets the startPos
        attribute to None.

        Args:
            event (QMouseEvent): The mouse event that triggered this method.
        """
        super().mouseReleaseEvent(event)
        self.startPos = None

###############
#   WINDOWS   #
###############
class HotkeyDialog(QtWidgets.QDialog):
    """
    HotkeyDialog is a custom dialog class that extends QtWidgets.QDialog to handle
    key press events, specifically the Escape key to close the dialog.

    Methods:
        __init__(self, parent, *args, **kwargs):

        keyReleaseEvent(self, event):
            Handles key press events for the widget. Calls the parent's keyReleaseEvent
            method and checks if the pressed key is the Escape key. If the Escape key
            is pressed, the widget will be closed.
    """
    def __init__(self, parent, *args, **kwargs):
        """
        Initialize a new instance of the class.

        Args:
            parent: The parent object.
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.
        """
        super().__init__(parent)

        self.setWindowFlags(QtCore.Qt.WindowType.Dialog)
        #self.setWindowFlags(QtCore.Qt.WindowType.Tool)
        #self.setWindowModality(QtCore.Qt.WindowModality.NonModal)

    def show(self, keep=False):
        """
        Shows the dialog and sets the focus to the dialog window.
        """
        # Move to the last position if available
        #if self.has_pos is not None:
        if self.pos():
            self.move(self.pos())
        # Show/hide the dialog
        if self.isVisible() and not keep:
            self.close()
        else:
            super().show()
            self.setFocus()

    def update(self):
        pass
        
    def keyReleaseEvent(self, event):
        """
        Handles key press events for the widget.

        Parameters:
        event (QKeyEvent): The key event that triggered this method.

        This method calls the parent's keyReleaseEvent method and checks if the
        pressed key is the Escape key. If the Escape key is pressed, the widget
        will be closed.
        """
        self.parent().keyReleaseEvent(event)
        k = event.key()
        if k == QtCore.Qt.Key.Key_Escape:
            self.close()

class AboutWindow(HotkeyDialog):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.path_settings = kwargs.get('path_settings', None)
        self.path_home = kwargs.get('path_home', None)
        self.pixmap = kwargs.get('pixmap', QtGui.QPixmap())
        self.icon = kwargs.get('icon', QtGui.QIcon())
        self.add_content()
        self.setFixedSize(self.sizeHint())

    def add_content(self):
        # title font
        font_title = QtGui.QFont()
        font_title.setPointSize(48)
        title = QtWidgets.QLabel(f'<b>xrdPlanner</b>')
        title.setFont(font_title)
        suptitle = QtWidgets.QLabel(f'<b>Version {xrdPlanner.__version__}</b> (released {xrdPlanner.__date__})')
        # add information
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
        path = QtWidgets.QLabel(f'<br>Open settings file <a href=file:///{self.path_settings}>location</a>')
        path.setOpenExternalLinks(True)
        detdb = QtWidgets.QLabel(f'<br>Find <u>old</u> detector db backup file (detector_db.json.bak) <a href=file:///{self.path_home}>here</a>')
        detdb.setOpenExternalLinks(True)
        # add widgets to layout
        box_layout = QtWidgets.QVBoxLayout()
        box_layout.setSpacing(0)
        for widget in [title, suptitle, github, published, authors, email, path, detdb]:
            box_layout.addWidget(widget)
        # box containing info
        self.box = QtWidgets.QGroupBox()
        self.box.setFlat(True)
        self.box.setLayout(box_layout)
        # label containing the logo
        self.logo = QtWidgets.QLabel()
        self.logo.setPixmap(self.pixmap.scaledToHeight(self.box.sizeHint().height(), mode=QtCore.Qt.TransformationMode.SmoothTransformation))
        # add both to the window layout
        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.logo)
        layout.addWidget(self.box)
        # finish the window and make it unresizable
        self.setWindowIcon(self.icon)
        self.setLayout(layout)

    def update_logo(self, pixmap):
        self.pixmap = pixmap
        self.logo.setPixmap(self.pixmap.scaledToHeight(self.box.sizeHint().height(), mode=QtCore.Qt.TransformationMode.SmoothTransformation))

class GeometryWindow(HotkeyDialog):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.setWindowTitle('Geometry conventions')
        self.icon = kwargs.get('icon', QtGui.QIcon())
        self.add_content()
        self.setFixedSize(self.sizeHint())
    
    def add_content(self):
        font_title = QtGui.QFont()
        font_title.setPointSize(28)
        title = QtWidgets.QLabel(f'<b>Geometry conventions</b>')
        title.setFont(font_title)
        description = QtWidgets.QLabel('<b>Translations</b> without any <b>tilt</b>/<b>rotation</b> are the horizontal and vertical distances '
                                       'between the centre of the detector and the point of normal incidence PONI (<i>top</i>). The SDD is the distance from '
                                       'the sample to the PONI. A <b>rotation</b> moves the detector along the goniometer circle '
                                       '(constant SDD), keeping the PONI at the same position relative to the detector '
                                       'surface, here the detector centre (<i>lower left</i>). A <b>tilt</b> rolls the detector surface on '
                                       'the goniometer circle, hence the SDD is fixed, but the PONI shifts along the detector '
                                       'face (<i>lower right</i>).')
        description.setAlignment(QtCore.Qt.AlignmentFlag.AlignJustify)
        description.setWordWrap(True)
        pmap = QtGui.QPixmap(':/icons/xrdPlanner_geom').scaled(512, 512, aspectRatioMode=QtCore.Qt.AspectRatioMode.KeepAspectRatio, transformMode=QtCore.Qt.TransformationMode.SmoothTransformation)
        icon = QtWidgets.QLabel()
        icon.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        icon.setPixmap(pmap)
        
        published = QtWidgets.QLabel('For more information see:<a href="https://doi.org/10.1107/S1600577523011086"> '
                                     '<i>J. Synchrotron Rad.</i> (2024). <b>31</b></a> or <a href="https://github.com/LennardKrause/xrdPlanner">Github</a>.')
        published.setOpenExternalLinks(True)

        box_layout = QtWidgets.QVBoxLayout()
        box_layout.setSpacing(6)
        box_layout.addWidget(title)
        box_layout.addWidget(description)
        box_layout.addWidget(icon)
        box_layout.addWidget(published)

        box = QtWidgets.QGroupBox()
        box.setFlat(True)
        box.setLayout(box_layout)
        box.setContentsMargins(0,0,0,0)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(box)
        
        self.setWindowIcon(self.icon)
        self.setLayout(layout)

class HotkeysWindow(HotkeyDialog):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.setWindowTitle('Hotkeys')
        self.icon = kwargs.get('icon', QtGui.QIcon())

        self.add_content()

        if 'hotkey_dict' in kwargs:
            self.add_hotkeys(kwargs['hotkey_dict'])

    def add_hotkeys(self, hotkey_dict):
        font_bold = QtGui.QFont()
        font_bold.setBold(True)

        self.table.setRowCount(0)
        for i,(entry, value) in enumerate(hotkey_dict.items()):
            self.table.insertRow(i)
            # handle title / sections
            if isinstance(entry, str) and not value:
                self.table.setItem(i, 0, QtWidgets.QTableWidgetItem(entry))
                self.table.setSpan(i, 0, 1, 2)
                #table.item(i, 0).setBackground(self.palette().highlight().color())
                self.table.item(i, 0).setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                self.table.item(i, 0).setFont(font_bold)
            else:
                key, mod = entry
                desc = value['dc']
                # handle special keys
                if key.lower() == 'left':
                    key = '\u2190'
                elif key.lower() == 'up':
                    key = '\u2191'
                elif key.lower() == 'right':
                    key = '\u2192'
                elif key.lower() == 'down':
                    key = '\u2193'
                else:
                    key = QtGui.QKeySequence(key).toString()
                # add modification keys
                if mod:
                    if mod == 'SHIFT':
                        mod = '\u21E7'
                    key = f'{mod}+{key}'
                # add to table
                self.table.setItem(i, 0, QtWidgets.QTableWidgetItem(key))
                self.table.setItem(i, 1, QtWidgets.QTableWidgetItem(desc))

        self.table.resizeColumnsToContents()
        self.resize(self.sizeHint())

    def add_content(self):
        self.table = QtWidgets.QTableWidget()
        self.table.setColumnCount(2)
        self.table.setSortingEnabled(False)
        self.table.setHorizontalHeaderLabels(['Key','Action'])
        self.table.setAlternatingRowColors(False)
        self.table.verticalHeader().hide()
        self.table.setWordWrap(False)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setEditTriggers(QtWidgets.QTableWidget.EditTrigger.NoEditTriggers)
        self.table.setSelectionMode(QtWidgets.QTableWidget.SelectionMode.NoSelection)
        self.table.setSizeAdjustPolicy(QtWidgets.QTableWidget.SizeAdjustPolicy.AdjustToContents)
        self.table.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.table)

        self.setWindowIcon(self.icon)
        self.setLayout(layout)

class AbsorptionWindow(HotkeyDialog):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.setWindowTitle(f'Transmission of {self.parent().geo.reference}')
        self.icon = kwargs.get('icon', QtGui.QIcon())
        self.threshold = kwargs.get('threshold', 2)
        self.abs_range = kwargs.get('abs_range', None)
        self.abs_vals = kwargs.get('abs_vals', None)
        self.stepsize = kwargs.get('stepsize', 0.1)
        self.scat_diam = self.parent().plo.scattering_diameter
        self.current_keV = None

        keV_lower = self.parent().lmt.ener_min
        keV_upper = self.parent().lmt.ener_max
        if abs(keV_upper-keV_lower) < 1:
            self.keV_range = (np.floor(keV_lower*0.49), np.ceil(keV_lower*1.51))
        else:
            self.keV_range = (keV_lower, keV_upper)
        self.add_content()
        self.setFixedSize(self.sizeHint())
    
    def add_content(self):
        layout = QtWidgets.QVBoxLayout()

        self.abs_plot = pg.PlotWidget(background=self.palette().base().color())
        self.abs_plot.setTitle(f'', color=self.palette().text().color())
        self.abs_plot.viewport().setAttribute(QtCore.Qt.WidgetAttribute.WA_AcceptTouchEvents, False)
        self.abs_plot.setMenuEnabled(False)
        self.abs_plot.setLabel('left', '-log(I/I\u2080)')
        self.abs_plot.setLabel('bottom', 'Energy', units='keV')
        self.abs_plot.showGrid(x=True, y=True, alpha=0.5)
        self.abs_plot.getPlotItem().getViewBox().setDefaultPadding(0.01)
        layout.addWidget(self.abs_plot)
        self.abs_scat = pg.ScatterPlotItem(size=5,
                                           symbol='o',
                                           pen=pg.mkPen('g', width=1),
                                           brush='g')
        self.curve_label = pg.TextItem(anchor=(0.0,1.0), color=self.palette().text().color())

        line_threshold = pg.InfiniteLine(angle=0, movable=False)
        line_threshold.setPen(pg.mkPen(self.parent().conic_highlight, width=1, style=QtCore.Qt.PenStyle.DashLine))
        line_threshold.setPos(self.threshold)
        self.line_current_kev = pg.InfiniteLine(angle=90, movable=False, pen=pg.mkPen('g', width=1))
        self.abs_curve = pg.PlotCurveItem(useCache=True)
        # plot item order
        self.abs_plot.addItem(self.abs_curve)
        self.abs_plot.addItem(line_threshold)
        self.abs_plot.addItem(self.line_current_kev)
        self.abs_plot.addItem(self.abs_scat)
        self.abs_plot.addItem(self.curve_label)

        self.parm_box = QtWidgets.QGroupBox('')
        self.parm_layout = QtWidgets.QHBoxLayout()
        self.parm_box.setLayout(self.parm_layout)

        self.pack_desc = QtWidgets.QLabel('Packing density')
        self.pack_slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.pack_slider.setRange(1, 100)
        self.pack_slider.setSingleStep(1)
        self.pack_slider.setPageStep(10)
        self.pack_slider.setValue(100)
        self.pack_slider.valueChanged.connect(self.update)
        self.pack_label = QtWidgets.QLabel()
        self.pack_label.setText(f'{self.pack_slider.value()} %')
        self.parm_layout.addWidget(self.pack_desc)
        self.parm_layout.addWidget(self.pack_slider)
        self.parm_layout.addWidget(self.pack_label)

        self.diam_desc = QtWidgets.QLabel(f'Scattering volume \u2300: {int(self.parent().plo.scattering_diameter*1e6)} \u00B5m')
        self.parm_layout.addWidget(self.diam_desc)

        layout.addWidget(self.parm_box)

        # add description
        description = QtWidgets.QLabel('This feature is currently in <b>test phase</b>, feedback is very welcome!<br>\
                                        The estimated FWHM (H, in degrees) is shown in the bottom right corner.')
        description.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(description)

        self.setLayout(layout)
    
    def calc_absorption(self, stepsize=0.1):
        self.stepsize = stepsize
        self.abs_range = np.arange(*self.keV_range, self.stepsize)
        self.abs_vals = self.parent().xtl.Properties.absorption(energy_kev=self.abs_range)*1e3 # mm-1

    def update(self):
        if self.parent().xtl is None:
            self.close()
            return
        # called by: mainWindow slider, update_screen
        pack_density = self.pack_slider.value()*1e-2
        self.pack_label.setText(f'{pack_density*100:.0f}%')
        if self.abs_range is None or self.abs_vals is None:
            self.calc_absorption()
        self.abs_curve.setData(self.abs_range, self.abs_vals*pack_density*self.scat_diam*1e3)
        if self.current_keV != self.parent().geo.ener:
            self.update_keV(self.parent().geo.ener)

    def update_keV(self, keV):
        self.line_current_kev.setPos(keV)
        abs_val = self.abs_curve.getData()[1][np.argmin(np.abs(self.abs_range-keV))]
        self.curve_label.setPos(keV, abs_val)
        self.abs_scat.setData([keV], [abs_val])
        self.curve_label.setText(f'{np.exp(-abs_val)*100:.0f}%')
        self.current_keV = keV

    def update_scat_diameter(self, value):
        self.diam_desc.setText(f'Scattering volume \u2300: {value*1e6} \u00B5m')
        self.scat_diam = value
        #self.update()

    def show(self, keep=False):
        self.update()
        return super().show(keep)

# references to parent
class DetdbWindow(HotkeyDialog):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.setWindowTitle('Edit detector databank')
        self.add_content()
    
    def add_content(self):
        # set flags=QtCore.Qt.WindowType.Tool for the window
        # to not loose focus when the FileDialog is closed
        self.finished.connect(self.win_detdb_on_close)
        layout_vbox = QtWidgets.QVBoxLayout()
        layout_hbox = QtWidgets.QHBoxLayout()
        layout_hbox.setContentsMargins(0,0,0,0)
        frame = QtWidgets.QFrame()
        frame.setLayout(layout_hbox)
        self.setLayout(layout_vbox)

        self.setStyleSheet('QGroupBox { font-weight: bold; } ')

        # Fonts
        font_header = QtGui.QFont()
        font_header.setPixelSize(self.parent().plo.slider_label_size)
        font_header.setBold(True)
        font_normal = QtGui.QFont()
        font_normal.setPixelSize(self.parent().plo.slider_label_size)
        font_normal.setBold(False)

        # Detector tree box
        qbox_detbank = QtWidgets.QGroupBox()
        qbox_detbank.setTitle('Detector bank')
        qbox_detbank.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
        layout_qbox_detbank = QtWidgets.QVBoxLayout()
        qbox_detbank.setLayout(layout_qbox_detbank)

        self._db_dict = self.parent().detector_db.copy()

        tooltips = {'hmp' : 'Module size (horizontal) [px]',
                    'vmp' : 'Module size (vertical) [px]',
                    'pxs' : 'Pixel size [mm]',
                    'hgp' : 'Module gap (horizontal) [px]',
                    'vgp' : 'Module gap (vertical) [px]',
                    'cbh' : 'Central beam hole [px]',
                    'name': 'Detector size name [str]',
                    'hmn' : 'Module number (horizontal) [int]',
                    'vmn' : 'Module number (vertical) [int]',
                    }

        # detector bank list widget
        self.list_db = QtWidgets.QListWidget()
        self.list_db.setAlternatingRowColors(True)
        self.list_db.currentItemChanged.connect(self.win_detdb_par_update)
        self.list_db.itemChanged.connect(self.win_detdb_dict_update)
        layout_qbox_detbank.addWidget(self.list_db)
        for det in self._db_dict.keys():
            item = QtWidgets.QListWidgetItem()
            item.setData(0, det)
            self.list_db.addItem(item)

        # Detector size box
        qbox_detsize = QtWidgets.QGroupBox()
        qbox_detsize.setTitle('Size')
        qbox_detsize.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
        layout_qbox_detsize = QtWidgets.QVBoxLayout()
        qbox_detsize.setLayout(layout_qbox_detsize)

        # doubleSpinBox editor for the beamstop bank list widget
        class itemDelegateInt(QtWidgets.QStyledItemDelegate):
            def createEditor(self, parent, option, index):
                if isinstance(index.data(0), int):
                    box = QtWidgets.QSpinBox(parent)
                    box.setSingleStep(1)
                    box.setMinimum(1)
                    box.setMaximum(100)
                    return box
                else:
                    return super(itemDelegateInt, self).createEditor(parent, option, index)
        
        # detector size table widget
        self.table_size = QtWidgets.QTableWidget(columnCount=3)
        self.table_size.setHorizontalHeaderLabels(['name','hmn', 'vmn'])
        self.table_size.horizontalHeaderItem(0).setToolTip(tooltips['name'])
        self.table_size.horizontalHeaderItem(1).setToolTip(tooltips['hmn'])
        self.table_size.horizontalHeaderItem(2).setToolTip(tooltips['vmn'])
        self.table_size.setAlternatingRowColors(True)
        self.table_size.verticalHeader().hide()
        self.table_size.setSelectionBehavior(QtWidgets.QTableWidget.SelectionBehavior.SelectRows)
        self.table_size.setSelectionMode(QtWidgets.QTableWidget.SelectionMode.SingleSelection)
        self.table_size.setItemDelegate(itemDelegateInt())
        self.table_size.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.table_size.cellChanged.connect(self.win_detdb_size_update)
        self.table_size.itemSelectionChanged.connect(self.win_detdb_dim_update)
        layout_qbox_detsize.addWidget(self.table_size)

        # Detector parameter box
        qbox_detpar = QtWidgets.QGroupBox()
        qbox_detpar.setTitle('Parameter')
        qbox_detpar.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
        layout_qbox_detpar = QtWidgets.QGridLayout()
        qbox_detpar.setLayout(layout_qbox_detpar)

        # detector par table widget
        self.detpar_items = {}
        det = next(iter(self._db_dict.keys()))
        for idx, (par, val) in enumerate(self._db_dict[det].items()):
            if par == 'size':
                pass
            elif par == 'pxs':
                widget_label = QtWidgets.QLabel(par)
                widget_label.setToolTip(tooltips[par])
                widget_value = QtWidgets.QDoubleSpinBox()
                widget_value.setValue(val)
                widget_value.setDecimals(4)
                widget_value.setMinimum(0.0001)
                widget_value.setMaximum(10)
                widget_value.setSingleStep(0.001)
                widget_value.setSuffix(' mm')
                widget_value.setObjectName(par)
                widget_value.valueChanged.connect(self.win_detdb_par_change)
                layout_qbox_detpar.addWidget(widget_label, idx, 0)
                layout_qbox_detpar.addWidget(widget_value, idx, 1)
                self.detpar_items[par] = widget_value
            else:
                widget_label = QtWidgets.QLabel(par)
                widget_label.setToolTip(tooltips[par])
                widget_value = QtWidgets.QSpinBox()
                widget_value.setValue(val)
                widget_value.setMinimum(0)
                widget_value.setMaximum(100000)
                widget_value.setSuffix(' px')
                widget_value.setObjectName(par)
                widget_value.valueChanged.connect(self.win_detdb_par_change)
                layout_qbox_detpar.addWidget(widget_label, idx, 0)
                layout_qbox_detpar.addWidget(widget_value, idx, 1)
                self.detpar_items[par] = widget_value

        layout_qbox_detpar.setRowStretch(idx, 1)

        # Detector dimensions box
        qbox_detdim = QtWidgets.QGroupBox()
        qbox_detdim.setTitle('Detector area (h \u00D7 v)')
        qbox_detdim.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
        layout_qbox_detdim = QtWidgets.QVBoxLayout()
        qbox_detdim.setLayout(layout_qbox_detdim)
        self._dim_mm = QtWidgets.QLabel()
        self._dim_px = QtWidgets.QLabel()
        layout_qbox_detdim.addWidget(self._dim_mm)
        layout_qbox_detdim.addWidget(self._dim_px)
        layout_qbox_detpar.addWidget(qbox_detdim, layout_qbox_detpar.rowCount(), 0, QtCore.Qt.AlignmentFlag.AlignCenter, 2)

        # Add row to detector list widget
        def row_add(aQListWidget):
            new = f'CUSTOM{aQListWidget.count():>02}'
            item = QtWidgets.QListWidgetItem()
            item.setData(2, new)
            item.setFlags(item.flags()|QtCore.Qt.ItemFlag.ItemIsEditable)

            self._db_dict[new] = {
                'pxs' : 100e-3, # [mm] Pixel size
                'hmp' : 100,    # [px] Module size (horizontal)
                'vmp' : 100,    # [px] Module size (vertical)
                'hgp' : 0,      # [px] Gap between modules (horizontal)
                'vgp' : 0,      # [px] Gap between modules (vertical)
                'cbh' : 0,      # [px] Central beam hole
                'size' : {'ONE':(1,1)},
                }
            
            aQListWidget.addItem(item)
            aQListWidget.setCurrentItem(item)

        # Remove row to detector list widget
        def rem_item(aQListWidget):
            det = aQListWidget.takeItem(aQListWidget.currentRow())
            self._db_dict.pop(det.text())

        # Remove row to detector list widget
        def rem_row(aQTableWidget):
            row = aQTableWidget.currentRow()
            det = self.list_db.currentItem().text()
            size = self.table_size.item(row, 0).text()
            aQTableWidget.removeRow(row)
            self._db_dict[det]['size'].pop(size)

        # Detector button box
        qbox_detbutton = QtWidgets.QGroupBox()
        qbox_detbutton.setTitle('')
        layout_qbox_detbutton = QtWidgets.QHBoxLayout()
        qbox_detbutton.setLayout(layout_qbox_detbutton)
        # add/remove buttons detector list widget
        button_add = QtWidgets.QToolButton()
        button_add.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
        button_add.setText('+')
        button_add.setToolTip('Add a new entry to the list.')
        button_add.clicked.connect(lambda: row_add(self.list_db))
        # remove button
        button_rem = QtWidgets.QToolButton()
        button_rem.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
        button_rem.setText('-')
        button_rem.setToolTip('Remove the highlighted entry from the list.')
        button_rem.clicked.connect(lambda: rem_item(self.list_db))
        layout_qbox_detbutton.addWidget(button_add)
        layout_qbox_detbutton.addWidget(button_rem)
        layout_qbox_detbank.addWidget(qbox_detbutton)

        # Add row to detector size list widget
        def row_add_size(aQTableWidget, currentItem):
            new = f'C{aQTableWidget.rowCount():>02}'
            self._db_dict[currentItem.text()]['size'][new] = (1,1)
            self.win_detdb_par_update(currentItem)

        # Detector size button box
        qbox_sizebutton = QtWidgets.QGroupBox()
        qbox_sizebutton.setTitle('')
        layout_qbox_sizebutton = QtWidgets.QHBoxLayout()
        qbox_sizebutton.setLayout(layout_qbox_sizebutton)
        # add/remove buttons for the beamstop bank list widget
        button_size_add = QtWidgets.QToolButton()
        button_size_add.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
        button_size_add.setText('+')
        button_size_add.setToolTip('Add a new entry to the list.')
        button_size_add.clicked.connect(lambda: row_add_size(self.table_size, self.list_db.currentItem()))
        # remove button
        button_size_rem = QtWidgets.QToolButton()
        button_size_rem.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
        button_size_rem.setText('-')
        button_size_rem.setToolTip('Remove the highlighted entry from the list.')
        button_size_rem.clicked.connect(lambda: rem_row(self.table_size))
        layout_qbox_sizebutton.addWidget(button_size_add)
        layout_qbox_sizebutton.addWidget(button_size_rem)
        layout_qbox_detsize.addWidget(qbox_sizebutton)

        # Close button box
        qbox_accept = QtWidgets.QGroupBox()
        #qbox_accept.setTitle('')
        #qbox_accept.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
        layout_qbox_accept = QtWidgets.QHBoxLayout()
        qbox_accept.setLayout(layout_qbox_accept)

        # store button
        button_store = QtWidgets.QToolButton()
        button_store.setSizePolicy(QtWidgets.QSizePolicy.Policy.Preferred, QtWidgets.QSizePolicy.Policy.Fixed)
        button_store.setText('Save')
        button_store.setToolTip('This will overwrite the detector_db.json. Changes are saved and available upon restart of xrdPlanner.')
        button_store.clicked.connect(self.win_detdb_overwrite)

        # preview button
        button_preview = QtWidgets.QToolButton()
        button_preview.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
        button_preview.setText('Preview')
        button_preview.setToolTip('Show a preview of the current detector. All changes are available during this session, the data bank will not be updated.')
        button_preview.clicked.connect(self.win_detdb_preview)

        # add widgets
        layout_hbox.addWidget(qbox_detbank)
        layout_hbox.addWidget(qbox_detpar)
        layout_hbox.addWidget(qbox_detsize)
        layout_vbox.addWidget(frame)
        layout_vbox.addWidget(qbox_accept)
        layout_qbox_accept.addWidget(button_store)
        layout_qbox_accept.addWidget(button_preview)
        
        # select detector
        self.list_db.setCurrentRow(0)
        # calculate detector size for given settings
        self.win_detdb_dim_update()
        self.list_db.setFocus()
    
    def win_detdb_dim_update(self):
        """
        Updates the detector database dimensions displayed in the UI.

        This method retrieves the current detector and its dimensions from the UI elements,
        calculates the dimensions in both millimeters and pixels, and updates the corresponding
        UI text fields with these values.

        The dimensions are calculated based on the horizontal and vertical multipliers and gaps,
        as well as the pixel size from the detector database.

        Raises:
            ValueError: If the text from the table items cannot be converted to integers.
        """
        det = self.list_db.currentItem().text()
        size_row = self.table_size.currentRow()
        hmn = int(self.table_size.item(size_row, 1).text())
        vmn = int(self.table_size.item(size_row, 2).text())
        pxs = self._db_dict[det]['pxs']
        xdim = self._db_dict[det]['hmp'] * hmn + self._db_dict[det]['hgp'] * (hmn-1) + self._db_dict[det]['cbh']
        ydim = self._db_dict[det]['vmp'] * vmn + self._db_dict[det]['vgp'] * (vmn-1) + self._db_dict[det]['cbh']
        self._dim_mm.setText(f'{xdim * pxs:.1f} \u00D7 {ydim * pxs:.1f} mm\u00B2')
        self._dim_px.setText(f'{xdim:.0f} \u00D7 {ydim:.0f} px')
    
    def win_detdb_par_change(self, val):
        """
        Updates the detector database parameter with the given value and refreshes the display.

        This method retrieves the current detector from the list, identifies the parameter to be updated
        based on the sender's object name, and updates the corresponding value in the internal database
        dictionary. After updating the value, it calls the method to update the display dimensions.

        Args:
            val: The new value to set for the specified parameter.
        """
        det = self.list_db.currentItem().text()
        par = self.sender().objectName()
        self._db_dict[det][par] = val
        self.win_detdb_dim_update()

    def win_detdb_size_update(self):
        """
        Updates the 'size' entry in the detector database dictionary based on the current items in the table.

        This method retrieves the currently selected detector from the list and updates its 'size' entry in the 
        internal database dictionary (`_db_dict`). The 'size' entry is replaced with a new dictionary where each 
        key is a size name and each value is a tuple containing horizontal and vertical measurements.

        The method also calls `win_detdb_dim_update` to update the dimensions.

        Note:
            The old 'size' entry is completely replaced, and any previous data is lost.

        Raises:
            AttributeError: If the current item in the list or any item in the table is None.
            ValueError: If the horizontal or vertical measurements cannot be converted to integers.
        """
        # The table is replacing the dict[det]['size'] entry on change
        # I don't know how to properly handle the renaming of a 'size'
        # the 'old' name, that is to be replaced be 'new', is lost
        det = self.list_db.currentItem().text()
        self._db_dict[det]['size'] = {}
        for row in range(self.table_size.rowCount()):
            size = self.table_size.item(row, 0).text()
            hmn = int(self.table_size.item(row, 1).text())
            vmn = int(self.table_size.item(row, 2).text())
            self._db_dict[det]['size'][size] = (hmn, vmn)
        self.win_detdb_dim_update()

    def win_detdb_dict_update(self, QListWidgetItem):
        """
        Updates the detector database dictionary with the selected detector.

        Args:
            QListWidgetItem (QListWidgetItem): The item selected from the QListWidget, representing a detector.

        Updates:
            - Moves the current detector entry to the new detector key in the database dictionary.
            - Calls the method to update the detector dimensions in the window.
        """
        det = QListWidgetItem.text()
        self._db_dict[det] = self._db_dict.pop(self.current_detector)
        self.win_detdb_dim_update()

    def win_detdb_par_update(self, item):
        """
        Updates the detector database parameters in the UI based on the selected item.

        Args:
            item (QTableWidgetItem): The selected item from the detector database table.

        This method performs the following actions:
        - Sets the current detector to the text of the selected item.
        - Iterates through the parameters and values of the selected detector in the database.
        - If the parameter is 'size', it updates the size table in the UI:
            - Blocks signals from the size table to prevent unwanted signal emissions.
            - Clears the contents of the size table.
            - Sets the row count of the size table to 0.
            - Iterates through the size values and populates the size table with the new values.
            - Sets the current cell of the size table to the first cell of each row.
            - Unblocks signals from the size table.
        - For other parameters, it updates the corresponding UI elements:
            - Blocks signals from the UI element to prevent unwanted signal emissions.
            - Sets the value of the UI element to the parameter value.
            - Unblocks signals from the UI element.
        - Calls the `win_detdb_dim_update` method to update the dimensions in the UI.
        """
        self.current_detector = item.text()
        for par, val in self._db_dict[item.text()].items():
            if par == 'size':
                self.table_size.blockSignals(True)
                self.table_size.clearContents()
                self.table_size.setRowCount(0)
                for i, (name, (hmn,vmn)) in enumerate(val.items()):
                    self.table_size.insertRow(i)
                    self.table_size.setItem(i, 0, QtWidgets.QTableWidgetItem(name))
                    item_hmn = QtWidgets.QTableWidgetItem(hmn)
                    item_hmn.setData(QtCore.Qt.ItemDataRole.EditRole, int(hmn))
                    item_vmn = QtWidgets.QTableWidgetItem(vmn)
                    item_vmn.setData(QtCore.Qt.ItemDataRole.EditRole, int(vmn))
                    self.table_size.setItem(i, 1, item_hmn)
                    self.table_size.setItem(i, 2, item_vmn)
                    self.table_size.setCurrentCell(i, 0)
                self.table_size.blockSignals(False)
            else:
                self.detpar_items[par].blockSignals(True)
                self.detpar_items[par].setValue(val)
                self.detpar_items[par].blockSignals(False)
        self.win_detdb_dim_update()
    
    def win_detdb_on_close(self):
        """
        Handles the event when the detector database window is closed.

        This method performs the following actions:
        1. Copies the current state of the detector database to a backup dictionary.
        2. Reinitializes the menu with a reset option.

        Returns:
            None
        """
        self.parent().detector_db = self._db_dict.copy()
        self.parent().menu_init(reset=True)

    def win_detdb_preview(self):
        """
        Updates the detector database and changes the detector.

        This method copies the current detector database to `self.detector_db`, 
        retrieves the selected detector and its size from the UI, and then 
        updates the detector configuration using the `change_detector` method.

        Steps:
        1. Copies the current detector database to `self.detector_db`.
        2. Retrieves the selected detector from `self.list_db`.
        3. Retrieves the size of the selected detector from `self.table_size`.
        4. Calls `self.change_detector` with the selected detector and size.

        Note:
            This method assumes that `self.list_db` and `self.table_size` are 
            properly initialized and populated with the relevant data.

        """
        # add new detectors to
        # -> self.geo.det_bank
        # if it is not empty!
        self.parent().detector_db = self._db_dict.copy()
        detector = self.list_db.currentItem().text()
        size = self.table_size.item(self.table_size.currentRow(), 0).text()
        self.parent().change_detector(detector, size)

    def win_detdb_overwrite(self):
        """
        Overwrites the detector database file with the current detector database.

        This method first previews the detector database changes by calling 
        `win_detdb_preview`. Then, it writes the current state of `self.detector_db` 
        to the file specified by `self.path_detdb` in JSON format with an indentation 
        of 4 spaces. Finally, it closes the detector database window.

        Returns:
            None
        """
        self.win_detdb_preview()
        with open(self.parent().path_detdb, 'w') as wf:
            json.dump(self.parent().detector_db, wf, indent=4)
        self.close()

# references to parent
class ExportWindow(HotkeyDialog):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.setWindowTitle('Export current settings to file')
        self.add_content()
    
    def add_content(self):
        # set flags=QtCore.Qt.WindowType.Tool for the window
        # to not loose focus when the FileDialog is closed
        layout_vbox = QtWidgets.QVBoxLayout()
        layout_hbox = QtWidgets.QHBoxLayout()
        layout_hbox.setContentsMargins(0,0,0,0)
        frame = QtWidgets.QFrame()
        frame.setLayout(layout_hbox)

        self.setLayout(layout_vbox)

        # Fonts
        font_header = QtGui.QFont()
        font_header.setPixelSize(self.parent().plo.slider_label_size)
        font_header.setBold(True)
        font_normal = QtGui.QFont()
        font_normal.setPixelSize(self.parent().plo.slider_label_size)
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
        for det, specs in self.parent().get_det_library(update=False, reset=False).items():
            item = QtWidgets.QTreeWidgetItem([str(det)])
            item.setFont(0, font_header)
            if self.parent().geo.det_type == str(det):
                _cur_det_type = item
            for val in specs['size'].keys():
                child = QtWidgets.QTreeWidgetItem([str(val)])
                child.setFont(0, font_normal)
                if self.parent().geo.det_size == str(val) and self.parent().geo.det_type == str(det):
                    _cur_det_size = child
                item.addChild(child)
            self.tree_det.addTopLevelItem(item)
        self.tree_det.expandAll()
        self.tree_det.itemClicked.connect(self.win_export_ghost_select)

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
        class itemDelegateFloat(QtWidgets.QStyledItemDelegate):
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
        self.tree_bsb.setItemDelegate(itemDelegateFloat())
        self.tree_bsb.setToolTip('Specify the available beamstop sizes.')
        self.tree_bsb.itemChanged.connect(self.tree_bsb.sortItems)
        for bs_size in self.parent().geo.bs_list:
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
        class genericItemDelegate(QtWidgets.QItemDelegate):
            def createEditor(self, parent, option, index):
                if index.column() == 1:
                    # bool needs to be evaluated before int as True and False
                    # will be interpreted as integers
                    if isinstance(index.data(0), bool):
                        return super(genericItemDelegate, self).createEditor(parent, option, index)
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
                    elif isinstance(index.data(0), str) and index.data(0).startswith('#'):
                        current_color = QtGui.QColor(index.data(0))
                        box = QtWidgets.QColorDialog(current_color, parent)
                        box.setOption(QtWidgets.QColorDialog.ColorDialogOption.ShowAlphaChannel, on=True)
                        return box
                    else:
                        return super(genericItemDelegate, self).createEditor(parent, option, index)
                return None
            
            def setModelData(self, editor, model, index):
                if isinstance(editor, QtWidgets.QColorDialog):
                    if editor.result():
                        model.setData(index, editor.selectedColor().name(format=QtGui.QColor.NameFormat.HexArgb))
                        model.setData(index.siblingAtColumn(2), QtGui.QBrush(editor.selectedColor()), QtCore.Qt.ItemDataRole.BackgroundRole)
                else:
                    return super(genericItemDelegate, self).setModelData(editor, model, index)

        # Parameter tree widget
        self.tree_par = QtWidgets.QTreeWidget()
        self.tree_par.setColumnCount(3)
        #self.tree_par.setHeaderLabels(['Parameter', 'Value', '#'])
        #self.tree_par.header().setFont(font_header)
        self.tree_par.setHeaderHidden(True)
        self.tree_par.setAlternatingRowColors(True)
        self.tree_par.setItemDelegate(genericItemDelegate())

        # det_bank and bs_list need some special treatment
        # to facilitate their editing
        dont_show = ['det_bank', 'bs_list']
        # detailed header
        translation = {'geo':'Geometry',
                       'plo':'Plot layout',
                       'thm':'Theme',
                       'lmt':'Limits'}
        for section, values in {'geo':self.parent().geo.__dict__,
                                'plo':self.parent().plo.__dict__,
                                'thm':self.parent().thm.__dict__,
                                'lmt':self.parent().lmt.__dict__}.items():
            item = QtWidgets.QTreeWidgetItem([translation[section]])
            item.setFont(0, font_header)
            item.setText(1, section)
            #item.setFlags(QtCore.Qt.ItemFlag.ItemIsSelectable|QtCore.Qt.ItemFlag.ItemIsEnabled)
            for par, val in values.items():
                if par in dont_show:
                    continue
                child = QtWidgets.QTreeWidgetItem()
                child.setData(0, 0, str(par))
                if section in self.parent().tooltips and par in self.parent().tooltips[section]:
                    child.setToolTip(0, self.parent().tooltips[section][par])
                    child.setToolTip(1, self.parent().tooltips[section][par])
                child.setData(1, 2, val)
                if isinstance(val, str) and val.startswith('#'):
                    child.setBackground(2, QtGui.QBrush(QtGui.QColor(val)))
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
        button_export.clicked.connect(self.win_export_to)

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

    def win_export_ghost_select(self, item):
        """
        Handles the selection logic for a tree-like structure in a GUI.

        Parameters:
        item (QTreeWidgetItem): The item whose selection state has changed.

        Behavior:
        - If the item is a top-level item (i.e., it has no parent), selecting it will select all its children,
          and deselecting it will deselect all its children.
        - If the item is a child item (i.e., it has a parent), selecting any child will select the parent,
          and deselecting all children will deselect the parent.
        """
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

    def win_export_to(self):
        """
        Exports the current detector and beamstop settings to a JSON file.
        This method performs the following steps:
        1. Collects selected detector items from the tree widget and organizes them into a dictionary.
        2. Retrieves beamstop values from a list widget and rounds them to 5 decimal places.
        3. Constructs a settings dictionary from the tree widget, rounding float values to 5 decimal places.
        4. Prompts the user to select a file location to save the settings as a JSON file.
        5. Adds the detector bank and beamstop list to the settings dictionary.
        6. Writes the settings dictionary to the selected JSON file.
        7. Closes the export window if it is open.
        8. Activates the newly exported settings file.
        9. Adds the new settings file to the custom settings menu and checks it if it is not already in the list.
        Returns:
            None
        """
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
            key = top.data(1, 0)
            settings[key] = dict()
            for j in range(top.childCount()):
                k = top.child(j).data(0, 0)
                v = top.child(j).data(1, 0)
                if isinstance(v, float):
                    v = round(v, 5)
                settings[key][k] = v

        # save parameters to settings file
        target, filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Export settings file', self.parent().path_settings_current, "Settings files (*.json)")
        if not target:
            return
        
        # add det_bank and bs_list to settings
        settings['geo']['det_bank'] = det_bank
        settings['geo']['bs_list'] = bs_list
        with open(target, 'w') as wf:
            json.dump(settings, wf, indent=4)
        
        # close the window
        if self.isVisible():
            self.close()

        # activate exported settings file
        basename = os.path.basename(target)
        self.parent().settings_change_file(basename)

        # add new settings file to menu and check it
        # add if not in list
        if basename not in [action.text() for action in self.parent().menu_custom_settings.actions()]:
            cset_action = QtGui.QAction(basename, self, checkable=True)
            self.parent().menu_set_action(cset_action, self.parent().settings_change_file, basename)
            self.parent().menu_custom_settings.addAction(cset_action)
            self.parent().group_cset.addAction(cset_action)
            cset_action.setChecked(True)
        # check if in list
        else:
            for action in self.parent().menu_custom_settings.actions():
                if action.text() == basename:
                    action.setChecked(True)

    ###############
    #  UNIT CELL  #
    ###############

# references to parent
class UnitCellWindow(HotkeyDialog):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.setWindowTitle('Custom unit cell')
        self.add_content()
    
    def add_content(self):
        layout = QtWidgets.QVBoxLayout()
        layout.setSpacing(0)

        box = QtWidgets.QGroupBox(title='Enter unit cell parameters')
        box.setStyleSheet('QGroupBox { font-weight: bold; }')
        box_layout = QtWidgets.QVBoxLayout()
        box_layout.setSpacing(0)
        box_layout.setContentsMargins(0,0,0,0)
        self.uc_dict_change = {}
        self.uc_check_boxes = {}
        self.linked_axes = True
        #              label,        unit,  link, min, max
        uc_widgets = [('a',     ' \u212b', [1,2],   1,  99, 2),
                      ('b',     ' \u212b',  True,   1,  99, 2),
                      ('c',     ' \u212b',  True,   1,  99, 2),
                      ('\u03b1', '\u00b0', False,  60, 150, 1),
                      ('\u03b2', '\u00b0', False,  60, 150, 1),
                      ('\u03b3', '\u00b0', False,  60, 150, 1),
                      ('Centring',  False, False,   0,   0, 0),
                      ('Sampling',   True, False,   0,   0, 0),
                      ('Name',       None, False,   0,   0, 0)]
        for idx, (label, unit, link, minval, maxval, decimals) in enumerate(uc_widgets):
            entry_box = QtWidgets.QFrame()
            entry_layout = QtWidgets.QHBoxLayout()
            entry_layout.setSpacing(6)
            entry_layout.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
            entry_label = QtWidgets.QLabel(label)
            entry_layout.addWidget(entry_label)
            if unit is None:
                box_widget = QtWidgets.QLineEdit(text=f'Custom {len(self.parent().ref_cell)+1:>02}')
                entry_layout.addWidget(box_widget)
            elif unit is False:
                box_widget = QtWidgets.QComboBox()
                box_widget.addItems(['P','F','I','C','A','B'])
                entry_layout.addWidget(box_widget)
            elif unit is True:
                box_widget = QtWidgets.QComboBox()
                box_widget.addItems(['1','2','3','4','5'])
                box_widget.setCurrentIndex(2)
                entry_layout.addWidget(box_widget)
            else:
                box_widget = QtWidgets.QDoubleSpinBox(decimals=decimals, singleStep=10**(-decimals), minimum=minval, maximum=maxval, value=self.parent().default_custom_cell[idx], suffix=unit)
                box_widget.setProperty('linked', link)
                box_widget.setProperty('index', idx)
                #box_widget.setStepType(QtWidgets.QDoubleSpinBox.StepType.AdaptiveDecimalStepType)
                entry_layout.addWidget(box_widget)
                box_check = QtWidgets.QCheckBox()
                box_check.setProperty('index', idx)
                self.uc_check_boxes[idx] = box_check
                if link is True:
                    box_check.setChecked(True)
                else:
                    box_check.setChecked(False)
                    box_check.setEnabled(False)
                entry_layout.addWidget(box_check)
            
            entry_box.setLayout(entry_layout)
            box_layout.addWidget(entry_box)
            self.uc_dict_change[idx] = box_widget
        
        for idx in range(6):
            self.uc_dict_change[idx].valueChanged.connect(self.win_uc_set_link_sbox)
            self.uc_check_boxes[idx].stateChanged.connect(self.win_uc_set_link_cbox)

        box.setLayout(box_layout)
        layout.addWidget(box)

        button_box = QtWidgets.QFrame()
        button_layout = QtWidgets.QHBoxLayout()
        # Add the apply button
        button_apply = QtWidgets.QDialogButtonBox()
        button_apply.addButton(QtWidgets.QDialogButtonBox.StandardButton.Apply)
        button_apply.setCenterButtons(True)
        button_apply.clicked.connect(lambda: self.win_uc_apply(self.uc_dict_change))
        button_layout.addWidget(button_apply)
        # Add the OK button
        button_ok = QtWidgets.QDialogButtonBox()
        button_ok.addButton(QtWidgets.QDialogButtonBox.StandardButton.Ok)
        button_ok.setCenterButtons(True)
        button_ok.clicked.connect(self.win_uc_accept)
        button_layout.addWidget(button_ok)
        # Set the button box layout
        button_box.setLayout(button_layout)
        layout.addWidget(button_box)

        self.setWindowIcon(self.parent().icon)
        self.setWindowTitle('Add custom unit cell')
        self.setLayout(layout)
        self.setFixedSize(self.sizeHint())
        self.accepted.connect(self.win_uc_accept)
        self.rejected.connect(self.close)
    
    def win_uc_apply(self, dict):
        """
        Applies the unit cell parameters from the given dictionary to the current object.
        Parameters:
        dict (dict): A dictionary containing unit cell parameters and other settings. 
                 Expected keys and their corresponding values:
                 - 0 to 5: Objects with a `value()` method returning unit cell parameters.
                 - 6: Object with a `currentText()` method returning the center.
                 - 7: Object with a `currentText()` method returning the decimal precision.
                 - 8: Object with a `text()` method returning the reference name.
        This method performs the following actions:
        1. Extracts unit cell parameters from the dictionary and stores them in `self.default_custom_cell`.
        2. Calculates the hkl values and stores them in `self.cont_ref_dsp` and `self.cont_ref_hkl`.
        3. Updates the reference name in `self.geo.reference`.
        4. Creates a `Ref` object with the reference name, dsp values, and hkl values, and stores it in `self.ref_cell`.
        5. Updates the window title and draws the reference.
        """
        uc = []
        for i in range(6):
            uc.append(dict[i].value())
        hkld = self.parent().calc_hkld(uc, dec=int(dict[7].currentText()), cen=dict[6].currentText())
        self.parent().default_custom_cell = uc

        self.parent().cont_ref_dsp = hkld[:,3]
        # hkl -> integer
        # cast hkl array to list of tuples (for easy display)
        _zeros = np.zeros(hkld.shape[0])
        self.parent().cont_ref_hkl = list(zip(hkld[:,0], hkld[:,1], hkld[:,2], _zeros, _zeros))

        self.parent().geo.reference = dict[8].text()
        reference = Ref(name=self.parent().geo.reference, dsp=self.parent().cont_ref_dsp, hkl=self.parent().cont_ref_hkl)
        self.parent().ref_cell[self.parent().geo.reference] = reference
        
        # update window title
        self.parent().set_win_title()
        self.parent().draw_reference()

    def win_uc_accept(self):
        """
        Applies changes to the unit cell, updates the reference, and manages the UI actions accordingly.

        This method performs the following steps:
        1. Applies changes to the unit cell dictionary.
        2. Updates the reference geometry.
        3. Creates a new QAction for the updated reference.
        4. Adds the new QAction to the relevant menu and action group.
        5. Sets the new QAction as checked.
        6. Closes the unit cell window.
        """
        self.win_uc_apply(self.uc_dict_change)
        self.parent().change_reference(self.parent().geo.reference)
        self.win_uc_add_to_menu(self.parent().geo.reference)
        self.close()
    
    def win_uc_set_link_sbox(self):
        """
        Handles the state change of a widget's 'linked' property and updates associated checkboxes and child widgets.

        If the 'linked' property of the sender widget is True, it sets the property to False and unchecks the associated checkbox.
        If the 'linked' property is a list, it iterates through the list, and for each index, if the corresponding checkbox is checked,
        it updates the value of the associated child widget without emitting signals.

        The sender widget is expected to have the following properties:
        - 'linked': A boolean or a list of indices.
        - 'index': An index to identify the associated checkbox.

        The method interacts with the following attributes of the class:
        - self.uc_check_boxes: A list of checkboxes.
        - self.uc_dict_change: A dictionary of child widgets.

        Returns:
            None
        """
        widget = self.sender()
        linked = widget.property('linked')
        if linked is True:
            widget.setProperty('linked', False)
            self.uc_check_boxes[widget.property('index')].setChecked(False)
        elif isinstance(linked, list):
            for idx in linked:
                checkbox = self.uc_check_boxes[idx]
                child = self.uc_dict_change[idx]
                if checkbox.isChecked():
                    child.blockSignals(True)
                    child.setValue(widget.value())
                    child.blockSignals(False)
    
    def win_uc_set_link_cbox(self):
        """
        Updates the 'linked' property of an item in the uc_dict_change dictionary based on the state of a checkbox.

        This method retrieves the sender of the signal (assumed to be a checkbox), gets its 'index' property, and updates
        the corresponding item in the uc_dict_change dictionary by setting its 'linked' property to the checked state of the checkbox.

        Attributes:
            self.uc_dict_change (dict): A dictionary where the keys are indices and the values are objects with a 'linked' property.
        """
        check = self.sender()
        self.uc_dict_change[check.property('index')].setProperty('linked', check.isChecked())

    def win_uc_add_to_menu(self, name):
        # Doesn't work on windows, menus won't show up!
        ref_action = QtGui.QAction(name, self, checkable=True)
        self.parent().menu_set_action(ref_action, self.parent().change_reference, name)
        self.parent().sub_menu_cell.addAction(ref_action)
        self.parent().group_ref.addAction(ref_action)
        ref_action.setChecked(True)
        self.parent().sub_menu_cell.setDisabled(self.parent().sub_menu_cell.isEmpty())

    ##########
    #  FWHM  #
    ##########

# references to parent
class FwhmWindow(HotkeyDialog):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.setWindowTitle('Setup calculation of the instrumental broadening')
        self.add_content()
    
    def add_content(self):
        # dict: index: needed to identify the vars upon reassignment (accept) and to link it to the tooltips
        #       name: Display string
        #       variable: the variable to use
        #       exponent: exponent to display the value more intuitive
        #       decimals: how many decimals to display
        #       stepsize: single step size
        #       unit: unit to display
        self.param_dict = {'Detector':[(0, 'Sensor thickness', self.parent().plo.sensor_thickness, 1e-6, 0, 1, '\u00B5m'),
                                       (1, 'Sensor material', self.parent().plo.sensor_material, 0, 0, 1, '')],
                               'Beam':[(2, 'Divergence', self.parent().plo.beam_divergence, 1e-6, 0, 1, '\u00B5rad'),
                                       (3, 'Energy resolution \u0394E/E', self.parent().plo.energy_resolution, 1e-4, 2, 0.1, '\u00B710\u207b\u2074')],
                             'Sample':[(4, 'Scattering volume \u2300', self.parent().plo.scattering_diameter, 1e-6, 0, 5, '\u00B5m')],}
                           #'Display':[(5, 'FWHM threshold', self.plo.funct_fwhm_thresh, 1e-3, 'FWHM [m\u00B0]')]}
        self.param_dict_change = {}

        # dict needed to update the preview plot in the setup window
        self.param_dict_tr = {'0':self.parent().plo.sensor_thickness,
                              '1':self.parent().plo.sensor_material,
                              '2':self.parent().plo.beam_divergence,
                              '3':self.parent().plo.energy_resolution,
                              '4':self.parent().plo.scattering_diameter}
        self.setStyleSheet('QGroupBox { font-weight: bold; }')

        # tooltip dictionary
        tooltips = {0:'Estimated effective thickness of the detector absorption layer depending on sensor material and incident wavelength.',
                    1:'Detector sensor layer material.',
                    2:'Estimated beam divergence.',
                    3:'Expected bandwidth of the incident X-ray beam.',
                    4:'Intersection of the sample size and the beam size as the diameter.'}
        # add user input spinboxes for the variables
        #  - it uses doublespinboxes for all floats
        #    and a combobox to choose the sensor material.
        #    Distinction is made by setting the divisor to 0.
        layout = QtWidgets.QVBoxLayout()
        for title, entry in self.param_dict.items():
            box = QtWidgets.QGroupBox(title=title)
            box_layout = QtWidgets.QVBoxLayout()
            box_layout.setContentsMargins(0,0,0,0)
            for idx, label, value, div, decimals, step, unit in entry:
                entry_box = QtWidgets.QFrame()
                entry_layout = QtWidgets.QHBoxLayout()
                entry_layout.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                if div > 0:
                    box_combobox = QtWidgets.QDoubleSpinBox(decimals=decimals, singleStep=step, minimum=step, maximum=100000, value=value/div)
                    box_combobox.valueChanged.connect(self.update)
                    box_combobox.setObjectName(f'{idx} {div}')
                else:
                    box_combobox = QtWidgets.QComboBox()
                    box_combobox.addItems(self.parent().att_lengths.keys())
                    box_combobox.setCurrentIndex(box_combobox.findText(self.parent().plo.sensor_material))
                    box_combobox.currentTextChanged.connect(self.update)
                    box_combobox.setObjectName(f'{idx} {div}')
                box_combobox.setToolTip(tooltips[idx])
                entry_label = QtWidgets.QLabel(label)
                entry_label.setToolTip(tooltips[idx])
                entry_layout.addWidget(entry_label)
                entry_layout.addWidget(box_combobox)
                entry_unit = QtWidgets.QLabel(unit)
                entry_unit.setToolTip(tooltips[idx])
                entry_layout.addWidget(entry_unit)
                entry_box.setLayout(entry_layout)
                box_layout.addWidget(entry_box)
                self.param_dict_change[idx] = [box_combobox, div]
            box.setLayout(box_layout)
            layout.addWidget(box)

        # preview plot
        preview_box = QtWidgets.QGroupBox('Preview')
        preview_box_layout = QtWidgets.QVBoxLayout()
        preview_box_layout.setContentsMargins(6,6,6,6)
        preview_box.setLayout(preview_box_layout)
        self.fwhm_curve = pg.PlotCurveItem()
        self.fwhm_curve.setPen(pg.mkPen(color=self.palette().text().color(), width=3))
        self.fwhm_plot = pg.PlotWidget(background=self.palette().base().color())
        self.fwhm_plot.viewport().setAttribute(QtCore.Qt.WidgetAttribute.WA_AcceptTouchEvents, False)
        self.fwhm_plot.setMenuEnabled(False)
        self.fwhm_plot.setTitle('FWHM = A\u00D7X\u2074 + B\u00D7X\u00B2 + C + M', color=self.palette().text().color())
        self.fwhm_plot.setLabel(axis='bottom', text='2\u03B8 [\u00B0]', color=self.palette().text().color())
        self.fwhm_plot.setLabel(axis='left', text='FWHM [\u00B0]', color=self.palette().text().color())
        self.fwhm_plot.scale(sx=2, sy=3)
        self.fwhm_plot.addItem(self.fwhm_curve)
        # Thomas-Cox-Hastings curve
        self.tch_curve = pg.PlotCurveItem()
        self.tch_curve.setPen(pg.mkPen(color=self.palette().highlight().color(), width=6))
        self.fwhm_plot.addItem(self.tch_curve)
        self.fwhm_line = pg.InfiniteLine(pos=0, angle=90, pen=pg.mkPen(None))
        self.fwhm_plot.addItem(self.fwhm_line)
        preview_box_layout.addWidget(self.fwhm_plot)
        layout.addWidget(preview_box)

        # Thomas-Cox-Hastings
        tch_box = QtWidgets.QGroupBox('Estimate Thomas-Cox-Hastings parameters')
        tch_box_layout = QtWidgets.QHBoxLayout()
        tch_box.setLayout(tch_box_layout)
        tch_btn_est = QtWidgets.QPushButton('Estimate')
        tch_btn_est.clicked.connect(self.estimate_tch)
        tch_box_layout.addWidget(tch_btn_est)
        tch_box_layout.addSpacing(20)
        tch_label_U = QtWidgets.QLabel('U:', alignment=QtCore.Qt.AlignmentFlag.AlignRight|QtCore.Qt.AlignmentFlag.AlignVCenter)
        tch_box_layout.addWidget(tch_label_U)
        self.tch_value_U = QtWidgets.QLabel()
        tch_box_layout.addWidget(self.tch_value_U)
        tch_label_V = QtWidgets.QLabel('V:', alignment=QtCore.Qt.AlignmentFlag.AlignRight|QtCore.Qt.AlignmentFlag.AlignVCenter)
        tch_box_layout.addWidget(tch_label_V)
        self.tch_value_V = QtWidgets.QLabel()
        tch_box_layout.addWidget(self.tch_value_V)
        tch_label_W = QtWidgets.QLabel('W:', alignment=QtCore.Qt.AlignmentFlag.AlignRight|QtCore.Qt.AlignmentFlag.AlignVCenter)
        tch_box_layout.addWidget(tch_label_W)
        self.tch_value_W = QtWidgets.QLabel()
        tch_box_layout.addWidget(self.tch_value_W)
        layout.addWidget(tch_box)

        # add description
        description = QtWidgets.QLabel('This feature is currently in <b>test phase</b>, feedback is very welcome!<br>\
                                        The estimated FWHM (H, in degrees) is shown in the bottom right corner.')
        description.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(description)
        # Add the apply & reset buttons
        button_box = QtWidgets.QGroupBox()
        button_box_layout = QtWidgets.QHBoxLayout()
        button_box.setLayout(button_box_layout)
        button_apply = QtWidgets.QPushButton('Apply')
        button_apply.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        button_apply.clicked.connect(lambda: self.accept(self.param_dict_change))
        button_box_layout.addWidget(button_apply)
        layout.addWidget(button_box)
        # Disclimer box info
        citation_box = QtWidgets.QGroupBox()
        citation_box.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        citation_box_layout = QtWidgets.QVBoxLayout()
        citation_box_layout.setContentsMargins(6,6,6,6)
        citation_box.setLayout(citation_box_layout)
        citation = QtWidgets.QLabel('The interested user is referred to the article:<br>\
                                      <b>On the resolution function for powder diffraction with area detectors</b><br>\
                                      <a href="https://doi.org/10.1107/S2053273321007506">\
                                      D. Chernyshov <i> et al., Acta Cryst.</i> (2021). <b>A77</b>, 497-505</a>')
        citation.setOpenExternalLinks(True)
        citation.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        citation_box_layout.addWidget(citation)
        layout.addWidget(citation_box)

        self.setWindowIcon(self.parent().icon)
        self.setLayout(layout)
        self.setFixedSize(self.sizeHint())
        #self.update()
    
    def update(self, val=None):
        """
        Updates the Full Width at Half Maximum (FWHM) plot based on the current parameters and theme settings.

        Parameters:
        val (optional): A value to update the parameter dictionary. If provided, it should be a numerical value.

        Description:
        - Updates the pen color and background of the FWHM plot based on the current theme.
        - If `val` is provided, updates the parameter dictionary `param_dict_tr` with the new value.
        - Calculates the FWHM using the current parameters:
            - Sensor thickness
            - Sensor material
            - Beam divergence
            - Energy resolution
            - Scattering diameter
        - Updates the FWHM plot with the calculated values.
        - Calls `win_pxrd_update` to refresh the plot with the updated parameters.
        """
        # change color if theme was changed
        self.fwhm_curve.setPen(pg.mkPen(color=self.palette().text().color(), width=3))
        self.fwhm_plot.setBackground(self.palette().base().color())

        # Crude calculation of the FWHM
        # 
        # get updated values from param_dict_tr
        # 0, self.plo.sensor_thickness
        # 1, self.plo.sensor_material
        # 2, self.plo.beam_divergence
        # 3, self.plo.energy_resolution
        # 4, self.plo.scattering_diameter
        if val is not None:
            idx, div = self.sender().objectName().split()
            div = float(div)
            if div > 0:
                self.param_dict_tr[idx] = val * div
            else:
                self.param_dict_tr[idx] = val
            #if idx == '4' and hasattr(self.parent(), 'abs_win'):
            #    self.parent().abs_win.update_scat_diameter(val*div)

        if len(self.parent().att_lengths[self.param_dict_tr['1']]) > int(self.parent().geo.ener):
            t = min(self.param_dict_tr['0'], self.parent().att_lengths[self.param_dict_tr['1']][int(self.parent().geo.ener)])
        else:
            t = self.param_dict_tr['0']
        p = self.parent().det.pxs * 1e-3
        phi = self.param_dict_tr['2']
        dE_E = self.param_dict_tr['3']
        c = self.param_dict_tr['4']

        dist = self.parent().geo.dist * 1e-3
        tth = np.linspace(0, np.pi/2, 180)
        # or use:
        # _current_tth_max = self.calc_tth_max()
        # to get max tth for the current setup

        # H2, FWHM
        A = 2*np.log(2) / dist**2 * (p**2-2*t**2-c**2)
        B = 2*np.log(2) / dist**2 * (2*t**2 + 2*c**2)
        C = 2*np.log(2) * phi**2
        M = (4*np.sqrt(2*np.log(2)) * dE_E)**2 * ((1-np.cos(tth))/(1+np.cos(tth)))
        X = np.cos(tth)
        H2 = A*X**4 + B*X**2 + C + M
        fwhm = np.sqrt(H2) * 180 / np.pi
        self.fwhm_curve.setData(np.rad2deg(tth), fwhm)
        # update pxrd plot if exists and visible
        if hasattr(self.parent(), 'pxrd_win') and self.parent().pxrd_win is not None and self.parent().pxrd_win.isVisible():
            self.parent().win_pxrd_update(update_dict=self.param_dict_tr)

    def accept(self, udict):
        """
        Updates the properties of the `plo` object based on the values from the provided dictionary
        and closes the given window.

        Parameters:
        udict (dict): A dictionary where each key corresponds to an index and each value is a list
                      containing a widget and a divisor. The widget's value is used to update the 
                      properties of the `plo` object after being multiplied by the divisor.
        win (QWidget): The window that will be closed after updating the properties.

        Updates:
        - self.plo.sensor_thickness: Rounded value from udict[0].
        - self.plo.sensor_material: Text from udict[1].
        - self.plo.beam_divergence: Rounded value from udict[2].
        - self.plo.energy_resolution: Rounded value from udict[3].
        - self.plo.scattering_diameter: Rounded value from udict[4].

        Additional Actions:
        - Enables the action to show FWHM.
        - Sets `self.plo.show_fwhm` to False.
        - Calls `self.toggle_fwhm()` to toggle the FWHM display.
        - Closes the provided window.
        """
        # assign combo/spinbox values
        # Index:[combo/spinbox widget, divisor]
        self.parent().plo.sensor_thickness = round(udict[0][0].value() * udict[0][1], 6)
        self.parent().plo.sensor_material = str(udict[1][0].currentText())
        self.parent().plo.beam_divergence = round(udict[2][0].value() * udict[2][1], 6)
        self.parent().plo.energy_resolution = round(udict[3][0].value() * udict[3][1], 6)
        self.parent().plo.scattering_diameter = round(udict[4][0].value() * udict[4][1], 6)
        #self.parent().plo.funct_fwhm_thresh = round(udict[5][0].value() * udict[5][1], 6)
        self.parent().action_funct_fwhm_show.setEnabled(True)
        self.parent().action_funct_fwhm_export.setEnabled(True)
        self.parent().plo.show_fwhm = False
        self.parent().toggle_fwhm()
        self.close()
    
    def export_grid(self):
        # calculate the FWHM grid
        # res = None -> use the current detector pixel dimensions
        _fwhm_grid = self.parent().calc_overlays(omega=0, res=None, pol=self.parent().plo.polarisation_fac)[-1]
        # find target file
        default_path = os.path.join(os.path.expanduser('~'), f"{self.parent().det.name.replace(' ', '_')}_FWHM")
        target, filter = QtWidgets.QFileDialog.getSaveFileName(self, 'Export FWHM grid', default_path, "Compressed numpy array (*.npz)")
        if not target:
            return
        # save compressed array to target
        np.savez_compressed(target, fwhm=np.flipud(_fwhm_grid))
    
    def estimate_tch(self):
        tth_max_rad = self.parent().calc_tth_max()
        tth_min = 0.0
        tth_max = np.rad2deg(tth_max_rad)
        tth, fwhm = self.fwhm_curve.getData()
        tch_cond = (tth > tth_min) & (tth < tth_max)
        tch_tth = tth[tch_cond]
        instr_fwhm = fwhm[tch_cond]

        def TCH_pv_fit(tth,*UVW):
            """parser for the TCH pseudo-Voigt function"""
            fwhm,_ = self.calc_TCH_pV(tth,*UVW,0,0)
            return fwhm
        
        UVW = 0, 0, instr_fwhm[0]**2
        [u,v,w], _ = curve_fit(TCH_pv_fit, tch_tth, instr_fwhm, p0=UVW, maxfev=10000)
        tch_fwhm,_ = self.calc_TCH_pV(tch_tth, u, v, w, 0, 0)
        self.tch_curve.setData(tch_tth, tch_fwhm)
        self.tch_value_U.setText(f'{u:.6f}')
        self.tch_value_V.setText(f'{v:.6f}')
        self.tch_value_W.setText(f'{w:.6f}')

    def calc_TCH_pV(self, tth, U, V, W, X, Y):
        """
        Thompson-Cox-Hastings pseudo-Voigt
        Gaussian FWHM: U*tan(theta)^2 + V*tan(theta) + W
        Lorentzian FWHM: X*tan(theta) + Y/cos(theta)
        return pseudo Voigt FWHM
        """
        theta = tth*np.pi/360
        tt = np.tan(theta)
        ct = np.cos(theta)
        G = np.sqrt(U*tt**2 + V*tt + W)
        L = X*tt + Y/ct

        H = (G**5 + 2.69269*G**4*L + 2.42843*G**3*L**2 \
                + 4.47163*G**2*L**3 + 0.07842*G*L**4 + L**5)**(1/5)
        # notice that the eta is defined differently in the FullProf manual
        # such that eta = 1-eta_FullProf
        eta = 1-(1.36603*L/H - 0.47719*(L/H)**2 + 0.11116*(L/H)**3)
        return H, eta

    def highlight(self, index):
        # called by HoverableCurveItem:highlight
        self.fwhm_line.setPen(pg.mkPen(self.parent().conic_highlight))
        self.fwhm_line.setPos(float(np.rad2deg(self.parent().dsp2tth(self.parent().cont_ref_dsp[index])[0])))

    def lowlight(self):
        # called by HoverableCurveItem:lowlight
        self.fwhm_line.setPen(pg.mkPen(None))

##################
#   PLOT ITEMS   #
##################
class HoverableCurveItem(pg.PlotCurveItem):
    def __init__(self, parent=None, hoverable=True, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.parent = parent
        self.basePen = self.opts['pen']
        self.hoverable = hoverable
        self.setAcceptHoverEvents(hoverable)

    def setData(self, *args, **kargs):
        # added the pen re-assignment to the function.
        if 'pen' in kargs:
            self.basePen = kargs['pen']
        return super().setData(*args, **kargs)

    def highlight(self, pos=None):
        self.parent.patches['ref_hl_curve'].setData(self.xData, self.yData)
        self.parent.patches['ref_hl_curve'].setVisible(True)
        self.parent.patches['ref_hl_curve'].index = self.index
        self.parent.patches['ref_hl_label'].setText(str(self.name))
        self.parent.patches['ref_hl_label'].setVisible(True)
        if pos is None:
            self.parent.patches['ref_hl_label'].setPos(self.xData[0], self.yData[0])
        else:
            self.parent.patches['ref_hl_label'].setPos(pos)

        self.parent.win_pxrd_highlight(self.index)
        self.parent.fwhm_win.highlight(self.index)

    def lowlight(self):
        self.setPen(self.basePen)
        self.parent.patches['ref_hl_label'].setVisible(False)
        self.parent.patches['ref_hl_curve'].setVisible(False)
        self.parent.win_pxrd_lowlight()
        self.parent.fwhm_win.lowlight()

    def hoverEvent(self, ev):
        if self.hoverable and not ev.isExit() and ev.buttons() == QtCore.Qt.MouseButton.LeftButton:
            if self.mouseShape().contains(ev.pos()):
                self.highlight(ev.pos())
    
    def mouseClickEvent(self, ev):
        if ev.buttons() == QtCore.Qt.MouseButton.LeftButton:
            if self.mouseShape().contains(ev.pos()):
                self.highlight(ev.pos())
        else:
            self.lowlight()
        return super().mouseClickEvent(ev)

class ClickableScatterPlotItem(pg.ScatterPlotItem):
    """
    Reimplementation of the ScatterPlotItem class with added click event handling.
    The original class is limited to left mouse button click events.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.offset_y = 0
    
    def setData(self, *args, **kargs):
        super().setData(*args, **kargs)
        if len(self.data['y']) > 0:
            self.offset_y = self.data['y'][0]

    def mouseClickEvent(self, ev):
        if not ev.button() == QtCore.Qt.MouseButton.NoButton:
            pts = self.pointsAt(ev.pos())
            if len(pts) > 0:
                self.ptsClicked = pts
                ev.accept()
                self.sigClicked.emit(self, self.ptsClicked, ev)
            else:
                #print "no spots"
                ev.ignore()
        else:
            ev.ignore()

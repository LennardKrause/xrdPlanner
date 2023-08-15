# xrdPlanner
#### A tool to project X-ray diffraction cones on a detector screen at different geometries (tilt, rotation, offset) and X-ray energies
 - Main application is to visualize the maximum achievable resolution at a given geometry.
 - Is able to project diffraction cones for standard samples or directly from cif files.
 - Can turn out valuable when planning beamtimes at synchrotron facilities (e.g. [DanMAX](https://www.maxiv.lu.se/beamlines-accelerators/beamlines/danmax/)).
 - Helps in deciding on the geometry of an experiment.
 - The math used is not meant to bring people to the moon but to provide a quick and simple preview.
 - The module building code is designed for [Dectris](https://www.dectris.com) [PILATUS3](https://www.dectris.com/detectors/x-ray-detectors/pilatus3/) / [EIGER2](https://www.dectris.com/detectors/x-ray-detectors/eiger2/) or [SACLA](https://sacla.xfel.jp/?lang=en) MPCCD Detectors (central hole geometry) but one-module systems like the [Bruker](https://www.bruker.com/en.html) [Photon II](https://www.bruker.com/en/products-and-solutions/diffractometers-and-scattering-systems/single-crystal-x-ray-diffractometers/sc-xrd-components/detectors.html) and [Rayonix](https://www.rayonix.com/) [MX-HS](https://www.rayonix.com/rayonix-mx-hs-series/) are possible as well.
 - It uses [python3](https://www.python.org), [numpy](https://numpy.org), [pyqt6](https://www.riverbankcomputing.com/software/pyqt/), [pyqtgraph](https://pyqtgraph.readthedocs.io/en/latest/), [pyFAI](https://pyfai.readthedocs.io/en/v2023.1/) and [Dans_Diffraction](https://github.com/DanPorter/Dans_Diffraction).

## Short how-to:
 - pip install xrdPlanner.
 - Type _xrdPlanner_ in a terminal and hit enter.
 - Choose a detector and a model from the _Detector_ menu.
 - Pick a reference from the _Reference_ menu to plot its contours ([pyFAI](https://pyfai.readthedocs.io/en/v2023.1/)).
 - Drop a .cif file onto the window to draw its contours ([Dans_Diffraction](https://github.com/DanPorter/Dans_Diffraction)), click a contour to get a hkl tooltip.
 - Use the units from the _Units_ menu you are the most comfortable with.
 - Hover over the grey line at the top to show the sliders. Click it to make it stay open.
 - Drag the sliders to change energy and geometry.

## Customisation:
  - Edit the _settings.json_ file and the _detector_db.json_ files.
  - Use _Settings_ -> _Edit Settings_ or _Edit Detector db_.
  - A click on _Apply Changes_ lets you see the difference.
  - _geo_ determines the startup defaults.
  - _plo_ customises the general layout and visuals.
  - _lmt_ sets the limiting values of the geometry/energy sliders.
  - Add all the missing detectors to the _detector_db.json_.
  - Check the _Settings file documentation_ at the bottom of this page for details.

## Known Bugs:
  - On Windows: Switching Dark/Light mode requires restart to change the window frame color.

## After the Update:
   Sometimes I might change the name of a parameter and you will get a warning message upon startup looking something like this: _WARNING: "conic_ref_min_int" is not a valid key_! Either that key is no longer in use or its name got changed and is now reset to the default value. The settings file is updated and the warning should no longer appear after restart. Apart from this, your edited settings file will not be altered after updating.
#### Changed the name of the following parameters:
  - conic_ref_min_int -> conic_ref_cif_int
  - conic_ref_irel_lw_min -> conic_ref_cif_lw_min
  - conic_ref_use_irel -> conic_ref_cif_irel
#### Added a new key:
  - conic_ref_cif_kev: this key sets the energy at which Dans_Dffraction calculates the intensities from a cif, increasing the value allows for higher resolution reference conics. However, the calculation will get slower.

## Latest updates:
  - 2023-08-15 Bugfix: Fixed several bugs with regard to the save/load of the settings file (Again).
  - 2023-08-15 Update: Changed the way themes / styles and customisation works internally.
  - 2023-07-14 Update: Added a key _plo.conic_ref_cif_kev_ to edit the energy for the cif intensity calculation.
  - 2023-07-14 Bugfix: Fixed a bug in the calculation of the conics, sections close to 90 deg. would sometimes not be drawn.
  - 2023-06-30 Update: Reference hkl intensity determines linewidth (irel).
  - 2023-06-30 Bugfix: Reference lines stay after settings reload.
  - 2023-06-23 Bugfix: Fixed several bugs with regard to the reloading of the settings file.
  - 2023-06-21 Update: Settings files accessible from menu, changes can be applied on the fly.
  - 2023-06-14 Update: Big speed update.
  - 2023-06-01 Update: countourpy was dropped, the conics are now calculated directly instead of being evaluated on a grid.
  - 2023-05-25 Update: Dans_Diffraction is used in favour of gemmi as it allows the direct calculation of intensities from the cif.
  - 2023-04-26 Update: A hkl tooltip is shown on click on a contour (only for cif contours).
  - 2023-04-25 Bugfix: Segmented contours are now drawn properly.
  - 2023-04-20 Bugfix: Confined slider window mobility to main window area.
  - 2023-04-10 Bugfix: Main window aspect ratio on Windows (menu bar within window).
  - 2023-04-10 Bugfix: Label size could not be adjusted.
  - 2023-04-10 Bugfix: Large angle (2-Theta > 90) contour label positions.
  - 2023-04-09 Update: Drop a cif file onto the window to draw its contours (uses [Dans_Diffraction](https://github.com/DanPorter/Dans_Diffraction)).
  - 2023-04-05 Update: Uses pyqt6, pyqtgraph and contourpy, dropped matplotlib backend.

## Older updates
  - 2023-03-23 Update: Settings are saved to (if doesn't exist) or loaded from (if exists) a 'settings.json' file.
  - 2023-03-23 Update: Added horizontal offset support and slider.
  - 2022-06-07 Update: Added functionality to plot Standard (LaB6, CeO2, ...) contours (needs [pyFAI](https://pyfai.readthedocs.io/en/master/)).
  - 2022-04-28 Update: Changed contour line generation to accept a list of 2-theta values as input.
  - 2022-04-27 Update: Added support for [SACLA](https://sacla.xfel.jp/?lang=en) MPCCD Detectors (central hole geometry).
  - 2022-04-25 Bugfix: Calculation of the beamcenter (rotation and tilt).

## Examples
#### A PILATUS3 2M detector and a Silicon sample.
![Preview](https://github.com/LennardKrause/xrdPlanner/blob/main/examples/PILATUS3_2M_Si.png)

#### A rotated EIGER2 4M detector and a Aluminium sample (darkmode).
![Preview](https://github.com/LennardKrause/xrdPlanner/blob/main/examples/EIGER2_4M_Al.png)


## Settings file documentation

#### geo - startup defaults
    det_type = 'EIGER2'  # [str]  Pilatus3 / Eiger2 / etc.
                         #        -> Detector menu entry
    det_size = '4M'      # [str]  300K 1M 2M 6M / 1M 4M 9M 16M
                         #        -> Detector submenu entry
    ener = 21.0          # [keV]  Beam energy
    dist = 75.0          # [mm]   Detector distance
    yoff = 0.0           # [mm]   Detector offset (vertical)
    xoff = 0.0           # [mm]   Detector offset (horizontal)
    rota = 25.0          # [deg]  Detector rotation
    tilt = 0.0           # [deg]  Detector tilt
    unit = 1             # [0-3]  Contour legend
                         #         0: 2-Theta
                         #         1: d-spacing
                         #         2: q-space
                         #         3: sin(theta)/lambda
    reference = 'None'   # [str]  Plot reference contours
                         #          pick from pyFAI
    darkmode = False     # [bool] Darkmode
    colormap = 'viridis' # [cmap] Contour colormap

#### plo - plot settings
    # - geometry contour section - 
    conic_tth_min = 5               # [int]    Minimum 2-theta contour line
    conic_tth_max = 150             # [int]    Maximum 2-theta contour line
    conic_tth_num = 30              # [int]    Number of contour lines
    beamcenter_marker = 'o'         # [marker] Beam center marker
    beamcenter_size = 6             # [int]    Beam center size
    conic_linewidth = 4.0           # [float]  Contour linewidth
    conic_label_size = 14           # [int]    Contour label size
    
    # - reference contour section - 
    conic_ref_linewidth = 10.0      # [float]  Reference contour linewidth
    conic_ref_num = 100             # [int]    Number of reference contours
    conic_ref_cif_int = 0.01        # [float]  Minimum display intensity (cif)
    conic_ref_cif_kev = 10.0        # [float]  Energy [keV] for intensity calculation
    conic_ref_cif_irel = True       # [int]    Linewidth relative to intensity
    conic_ref_cif_lw_min = 2.0      # [float]  Minimum linewidth when using irel
    conic_hkl_show_int = False      # [bool]   Show intensity in hkl tooltip
    conic_hkl_label_size = 14       # [int]    Font size of hkl tooltip
    
    # - module section - 
    det_module_alpha = 0.20         # [float]  Detector module alpha
    det_module_width = 1            # [int]    Detector module border width
    
    # - general section - 
    conic_steps = 100               # [int]    Conic resolution
    plot_size = 768                 # [int]    Plot size, px
    plot_size_fixed = True          # [int]    Fix window size
    unit_label_size = 16            # [int]    Label size, px
    
    # - slider section - 
    slider_margin = 12              # [int]    Slider frame top margin
    slider_border_width = 1         # [int]    Slider frame border width
    slider_border_radius = 1        # [int]    Slider frame border radius (px)
    slider_label_size = 14          # [int]    Slider frame label size
    slider_column_width = 75        # [int]    Slider label column width
    enable_slider_ener = True       # [bool]   Show energy slider
    enable_slider_dist = True       # [bool]   Show distance slider
    enable_slider_rota = True       # [bool]   Show rotation slider
    enable_slider_yoff = True       # [bool]   Show vertical offset slider
    enable_slider_xoff = True       # [bool]   Show horizontal offset slider
    enable_slider_tilt = True       # [bool]   Show tilt slider
    
    # - update/reset - 
    update_settings = True          # [bool]   Update settings file after load
    update_det_bank = True          # [bool]   Update detector bank after load
    reset_settings = False          # [bool]   Reset settings file
    reset_det_bank = False          # [bool]   Reset detector bank
    
    # - debug/testing -
    set_debug = False               # [bool]   Debug mode

#### thm - theme
    color_dark = '#404040'                # [color]  Global dark color
    color_light = '#EEEEEE'               # [color]  Global light color
    
    # light mode
    light_conic_label_fill = '#FFFFFF'    # [str]    Contour label fill color
    light_conic_ref_color = '#DCDCDC'     # [color]  Reference contour color
    light_det_module_color = '#404040'    # [color]  Detector module border color
    light_det_module_fill = '#404040'     # [color]  Detector module background color
    light_plot_bg_color = '#FFFFFF'       # [str]    Plot background color
    light_unit_label_color = '#808080'    # [str]    Label color
    light_unit_label_fill = '#FFFFFF'     # [str]    Label fill color
    light_slider_border_color = '#808080' # [str]    Slider frame border color
    light_slider_bg_color = '#AAC0C0C0'   # [str]    Slider frame background color
    light_slider_bg_hover = '#C0C0C0'     # [str]    Slider frame hover color
    light_slider_label_color = '#000000'  # [str]    Slider frame label color
    
    # dark mode
    dark_conic_label_fill = '#000000'     # [str]    Contour label fill color
    dark_conic_ref_color = '#202020'      # [color]  Reference contour color
    dark_det_module_color = '#EEEEEE'     # [color]  Detector module border color
    dark_det_module_fill = '#EEEEEE'      # [color]  Detector module background color
    dark_plot_bg_color = '#000000'        # [str]    Plot background color
    dark_unit_label_color = '#C0C0C0'     # [str]    Label color
    dark_unit_label_fill = '#000000'      # [str]    Label fill color
    dark_slider_border_color = '#202020'  # [str]    Slider frame border color
    dark_slider_bg_color = '#AA303030'    # [str]    Slider frame background color
    dark_slider_bg_hover = '#303030'      # [str]    Slider frame hover color
    dark_slider_label_color = '#C0C0C0'   # [str]    Slider frame label color

#### lmt - limits

    # energy [keV]
    ener_min =  5.0    # [float] Energy minimum [keV]
    ener_max =  100.0  # [float] Energy maximum [keV]
    ener_stp =  1.0    # [float] Energy step size [keV]
    
    # distance [mm]
    dist_min =  40.0   # [float] Distance minimum [mm]
    dist_max =  1000.0 # [float] Distance maximum [mm]
    dist_stp =  1.0    # [float] Distance step size [mm]
    
    # X-offset [mm]
    xoff_min = -150.0  # [float] Horizontal offset minimum [mm]
    xoff_max =  150.0  # [float] Horizontal offset maximum [mm]
    xoff_stp =  1.0    # [float] Horizontal offset step size [mm]
    
    # Y-offset [mm]
    yoff_min = -250.0  # [float] Vertical offset minimum [mm]
    yoff_max =  250.0  # [float] Vertical offset maximum [mm]
    yoff_stp =  1.0    # [float] Vertical offset step size [mm]
    
    # rotation [deg]
    rota_min = -60.0   # [float] Rotation minimum [deg]
    rota_max =  60.0   # [float] Rotation maximum [deg]
    rota_stp =  1.0    # [float] Rotation step size [deg]
    
    # tilt [deg]
    tilt_min = -25.0   # [float] Tilt minimum [deg]
    tilt_max =  25.0   # [float] Tilt maximum [deg]
    tilt_stp =  1.0    # [float] Tilt step size [deg]

## Detector db entries

The detector_db.json is stored as a dictionary
- key: detector name e.g. PILATUS3
- value: dictionary {Entry:Value,}, see table below
 
 | Entry |  Value  | Hint |
 |-------|---------|------|
 |  hms  | 83.8    | [mm]   Module size (horizontal)
 |  vms  | 33.5    | [mm]   Module size (vertical)
 |  pxs  | 172e-3  | [mm]   Pixel size
 |  hgp  | 7       | [pix]  Gap between modules (horizontal)
 |  vgp  | 17      | [pix]  Gap between modules (vertical)
 |  cbh  | 0       | [mm]   Central beam hole

The size Entry is a dictionary {key:value,}
- key: detector size / type, e.g. 300K
- value: list [hmn, vmn]
  - hmn: [int]  Number of modules (horizontal)
  - vmn: [int]  Number of modules (vertical)
 
#### Example of the PILATUS3 entry
    "PILATUS3": {
        "hms": 83.8,
        "vms": 33.5,
        "pxs": 0.172,
        "hgp": 7,
        "vgp": 17,
        "cbh": 0,
        "size": {
            "300K": [1, 3],
              "1M": [2, 5],
              "2M": [3, 8],
              "6M": [5,12]
        }
    },


#### Example code for adding xrdPlanner as a widget into an existing GUI
###### xrdPlanner uses its own menu bar, setting the GUI as the parent for xrdPlanner makes it add its menus to the parents menu bar, and likely more in the future.

    import sys
    from PyQt6 import QtWidgets
    from xrdPlanner.classes import MainWindow as xrdPlanner
    
    class MainWindow(QtWidgets.QMainWindow):
        def __init__(self):
            super().__init__()
            # make layout and widget
            layout = QtWidgets.QGridLayout()
            central_widget = QtWidgets.QWidget()
            central_widget.setLayout(layout)
            self.setCentralWidget(central_widget)
            # add xrdPlanner to layout
            xrdPlanner_as_widget = xrdPlanner(parent=self)
            layout.addWidget(xrdPlanner_as_widget)
    
    if __name__ == '__main__':
        app = QtWidgets.QApplication(sys.argv)
        window = MainWindow()
        window.show()
        app.exec()

#### I hope this turns out to be useful for someone!

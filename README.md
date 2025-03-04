# xrdPlanner
## Version 2.4.0 (released 04.03.2025)
#### A tool to project X-ray diffraction cones on a detector screen at different geometries (tilt, rotation, offset) and X-ray energies
 - Article published in [_J. Synchrotron Rad._ (2024). **31**](https://doi.org/10.1107/S1600577523011086).
 - Main application is to visualize the maximum achievable resolution at a given geometry.
 - Is able to project diffraction cones for standard samples or directly from cif files.
 - Can turn out valuable when planning beamtimes at synchrotron facilities (e.g. [DanMAX](https://www.maxiv.lu.se/beamlines-accelerators/beamlines/danmax/)).
 - Helps in deciding on the geometry of an experiment.
 - Has a 1D [PXRD plot window](#pxrd-pattern-plot) for cif files to show the diffraction pattern and [intrumental peak broadening](#instrumental-broadening).
 - Meant to run standalone but is readily [insertable](#example-code) as a widget into an existing PyQt6 GUI.
 - The math used is not meant to bring people to the moon but to provide a quick and simple preview.
 - This is not meant to accurately simulate a diffraction experiment, the step sizes are integer values in mm or degrees.
 - The module building code is designed for [Dectris](https://www.dectris.com) [PILATUS3](https://www.dectris.com/detectors/x-ray-detectors/pilatus3/) / [EIGER2](https://www.dectris.com/detectors/x-ray-detectors/eiger2/) or [SACLA](https://sacla.xfel.jp/?lang=en) MPCCD Detectors (central hole geometry) but one-module systems like the [Bruker](https://www.bruker.com/en.html) [Photon II](https://www.bruker.com/en/products-and-solutions/diffractometers-and-scattering-systems/single-crystal-x-ray-diffractometers/sc-xrd-components/detectors.html) and [Rayonix](https://www.rayonix.com/) [MX-HS](https://www.rayonix.com/rayonix-mx-hs-series/) are possible as well.
 - It uses [python3](https://www.python.org), [numpy](https://numpy.org), [pyqt6](https://www.riverbankcomputing.com/software/pyqt/), [pyqtgraph](https://pyqtgraph.readthedocs.io/en/latest/), [pyFAI](https://pyfai.readthedocs.io/en/v2023.1/) and [Dans_Diffraction](https://github.com/DanPorter/Dans_Diffraction).

>[!TIP]
>## Short how-to:
> - (python3 -m) pip install xrdPlanner.
> - Type _xrdPlanner_ in a terminal and hit enter.
> - Choose a detector and a model from the _Detector_ menu.
> - Plot diffraction cone contours:
>   - Pick a reference from the _Reference_ menu ([pyFAI](https://pyfai.readthedocs.io/en/v2023.1/)), no hkl tooltip available.
>   - Drop a .cif file onto the window ([Dans_Diffraction](https://github.com/DanPorter/Dans_Diffraction)), click a contour to get a hkl tooltip.
>   - Enter custom unit cell dimensions, click a contour to get a hkl tooltip.
> - Use the units from the _Units_ menu you are the most comfortable with.
> - Hover over the grey line at the top to show the sliders. Click it to make it stay open.
> - Drag the sliders to change energy and geometry.
> - Use [hotkeys](#hotkeys) to streamline your experience.
>
>## Customisation:
>  - Use the [export window](#export-window) from the settings menu.
>    - Build up the available detector and beamstop bank.
>    - Review and change the parameters.
>    - Tooltips should help you along the way.
>  - Use the [detector db editor](#detector-db-editor).
>    - Define the detector you collect your data with.
>    - Add all the different sizes and versions.
>    - Use _Preview_ to play around and _Save_ to update the data bank.
>  - Or edit the _settings.json_ file and the _detector_db.json_ files directly.
>    - Use _Settings_ -> _Edit files_ to edit the _current settings_ or _Detector db_ file.
>    - Check the [settings file documentation](#settings-file-documentation) to get started.
>  - Use pre-set beamline settings files:
>    - Download a settings file from here (e.g. settings/DanMAX.json).
>    - Use the import settings function from the GUI to import.

## Conventions
The geometry is defined with the center of rotation at the sample position, such that the radius of the rotation circle is equal to the sample to detector distance (SDD). That is, the rotation moves the detector along the goniometer circle, keeping the point of normal incidence (PONI) at the same position relative to the detector surface. At 0° the detector is vertical and at 90° the detector is horizontal with the detector normal pointing down.
The tilt angle is defined relative to the detector face and such that the PONI shifts along the detector face, keeping the SDD fixed. The detector face is thus always tangential to the goniometer circle. The $tilt$ can intuitively be described as a rolling motion along the goniometer circle and is considered a convenience function as it is equivalent to a combination of $rotation$ and $y_{offset}$. Consequently, the vertical shift of the PONI position ($\Delta y_{PONI}$) on the detector face is equal to the arclength of a section on the goniometer circle with an angular span equal to the tilt angle.
```math
\Delta y_{PONI} = SDD \cdot tilt
```
The vertical shift of the point of beam incidence (POBI) in the detector reference frame ($\Delta y_{POBI}$) can be described by the side length of the right-angle triangle spanned by the goniometer origin, the PONI, and the POBI, subtracted from the vertical shift of the PONI.
```math
\Delta y_{POBI} = SDD \cdot \left( tilt - tan \left( rotation + tilt \right) \right)
```

[**pyFAI**](https://pyfai.readthedocs.io/en/v2023.1/geometry.html#geometry) uses three rotations, $rot1$, $rot2$, and $rot3$, and three translations, $dist$, $poni1$, and $poni2$, to define the detector position. xrdPlanner uses the same translations ($SDD$, $y_{offset}$, and $x_{offset}$), but only one of the rotations ($rot2$). Apart from the change in sign, the pyFAI $rot2$ and xrdPlanner $rotation$ are equivalent, however, the $tilt$ in pyFAI convention is described by a combination of $rot2$ and shift of $poni1$:
```math
rot2 = -\left( rotation + tilt\right)
```
combined with
```math
\Delta poni1 = SDD \cdot tilt + y_{offset}
```
Additionally, pyFAI places the origin at the lower left corner of the detector, whereas xrdPlanner uses the center of the detector as origin, so one must account for translational shifts corresponding to half the width and height of the detector.

## Polarisation
The polarisation value indicated in xrdPlanner is the apparent scattering intensity, i.e. for P=0.2 only 20 % of the nominal intensity is observed, such that $I_{ij}^{obs}=P_{ij}\cdot I_0$. (Omitting solid angle).
The polarisation value $P$ for the $ij^{th}$ pixel is defined as [[1]](https://doi.org/10.1002/9781119998365.ch1):

```math
P_{ij} = 1 - [p (\vec{v}_{ij} \cdot  \vec{e}_x)^2+(1-p) (\vec{v}_{ij} \cdot  \vec{e}_y)^2] \qquad
```

 - $p$ = $1.0$ *horizontal polarisation*

 - $p$ = $0.5$ *random polarisation*

 - $p$ = $0.0$ *vertical polarisation*

where $p$ is the polarisation factor, $\vec{v}_{ij}$ is the normalised vector coordinate of the $ij^{th}$ pixel, and $\vec{e}_x$ and $\vec{e}_y$ are the horizontal and vertical basis vectors.
**NB:** The polarisation factor differs from the pyFAI convention, which goes from -1 to 1 rather than from 0 to 1, such that $p=(p _{pyFAI}+1)/2$  

## Solid-angle
The solid-angle correction factor $S_{ij}$ is defined as the normalized reciprocal of the cubed length of the vector $v_{ij}$ pointing from the sample position to the detector pixel $p_{ij}$.
```math
S_{ij} = \left(\frac{|\vec{v}_{ij}|}{SDD}\right)^{-3}
```
where $SDD$ is the sample to detector distance.

>[!NOTE]
>## Update to Version 2.0.0 and later (released 23.09.2024)
> - The definition of detector modules was changed from mm to px to be more accurate and consistent.
> - Consequently, custom detectors added to the detector_db file pre 2.0.0 are no longer working.
> - A new [detector db editor](#detector-db-editor) was added (_Settings_ > _Detector db editor_) to make the addition of custom detectors more feasible.
> - Define and add whatever detector you collect your data with, now more easily.
> - If something goes wrong use _Settings_ > _Reset detector db_.
> - There is a backup of your _custom_ detector_db (detector_db.json.bak) in the xrdPlanner folder (_Help_ > _xrdPlanner_).
> - A new parameter was added to allow for detector screen padding (plo.plot_padding), default is 0.

## Latest updates:
  - 2025-03-04 Update: Added pxrd ghost lines to play with and compare the scattering using different parameters.
  - 2025-03-04 Update: Made some improvements to the highlighting.
  - 2025-03-04 Update: Trained most sub-windows to remember their place.
  - 2025-03-04 Update: Started to declutter the code.
  - 2025-03-04 Bugfix: Resynchronized the menu button for the angular grid (Thanks Frederik).
  - 2025-02-12 Bugfix: Fixed f-string SyntaxError crash on startup (python version < 3.12).
  - 2025-02-10 Update: Added the option to estimate Thomas-Cox-Hastings parameters (gaussian) from the estimated instrumental broadening.
  - 2025-02-10 Update: Reworked the hkl display and cone highlighting to make identification and location of hkl more straightforward.
  - 2025-02-10 Update: Highlighting is now linked between the main window and the pxrd viewer.
  - 2024-12-05 Update: Added the option to export the FWHM grid as a compressed numpy .npz file (View > Functions > Export FWHM).
  - 2024-11-08 Update: Added a visual indicators for azimuthal bins (View > Overlay > Grid or hotkey G).
  - 2024-11-08 Bugfix: Fixed invisible menus on windows os.
  - 2024-10-29 Update: Added a menu action to invert the coloring of the resolution and reference conics (_View_ > _Invert cone colors_) to increase visibility, when needed.
  - 2024-10-28 Update: Added a [PXRD plot window](#pxrd-pattern-plot) for cif files to show the diffraction pattern (gaussian peak shape) including [intrumental peak broadening](#instrumental-broadening) (setup FWHM).
  - 2024-10-28 Update: cif files are stored now (name and link-to-file).
  - 2024-10-28 Update: Added somewhat useful doc-strings to all functions (Thanks Github Copilot).
  - 2024-10-01 Bugfix: Fixed an error crashing the detector db editor on python < 3.12 (Thanks Lise).
  - 2024-09-23 Update: Added a [detector db editor](#detector-db-editor) to define and add whatever detector you collect your data with, now more easily.
  - 2024-09-23 Update: Changed the detector module definition from mm to px (pre 2.0.0 custom detectors won't work!).
  - 2024-09-23 Update: Added a preview plot to the FWHM setup window.

<details>
<summary>Older updates</summary>
  
  - 2024-06-20 Update: Added highlighting of a clicked contour.
  - 2024-06-20 Bugfix: Incorrect removal of beamstop contour.
  - 2024-05-30 Update: Added support for custom unit cells to calculate and plot conics.
  - 2024-05-22 Bugfix: Fixed a bug removing conics (Thanks Mads).
  - 2024-05-02 Update: Added the option to calculate the [resolution function](https://doi.org/10.1107/S2053273321007506) (as FWHM, H) for powder diffraction with area detectors (Thanks Stefano).
  - 2024-05-02 Update: Automatic label placement (conic_label_auto, default=True).
  - 2024-05-02 Update: Dynamic cones try to keep the number of visible conics constant (conic_tth_auto, default=True).
  - 2024-05-02 Update: Current colors are displayed in the export window and a colorpicker is used.
  - 2024-02-17 Bugfix: Fixed crash on importing settings files (Thanks Andy).
  - 2024-02-16 Update: Added a unit value at mouse position (hover) to the overlay.
  - 2024-02-09 Update: Updated the About window with links to github and the publication.
  - 2024-02-08 Bugfix: Fixed crash on cif drag/drop (Thanks Mads).
  - 2023-11-28 Update: Added hotkeys to toggle between units/colormaps/overlays.
  - 2023-11-28 Update: Added polarisation and solid angle correction factor overlays.
  - 2023-11-12 Update: Added a new window to export settings to a file.
  - 2023-11-12 Update: Added the option to limit the available detectors for a settings file.
  - 2023-11-12 Update: Upon crash the program will start using the default settings.
  - 2023-09-26 Update: Added a PONI marker.
  - 2023-09-26 Update: Added the option to add custom labels to the sliders.
  - 2023-09-26 Update: Added the option to automatically find a reasonable window size (set plot_size to 0).
  - 2023-09-26 Bugfix: Fixed a bug that prevented the beamstop menu from updating upon changing the available beamstop list.
  - 2023-09-26 Bugfix: Fixed a bug in the calculation of the beamcenter for the combination of rotation and tilt.
  - 2023-08-29 Update: Added a feature to import and switch between settings files.
  - 2023-08-22 Bugfix: Fixed missing symbols and the slider bar on Linux.
  - 2023-08-22 Update: Added a beamstop, define distance to sample with a slider and pick a size from the menu.
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
  - 2023-03-23 Update: Settings are saved to (if doesn't exist) or loaded from (if exists) a 'settings.json' file.
  - 2023-03-23 Update: Added horizontal offset support and slider.
  - 2022-06-07 Update: Added functionality to plot Standard (LaB6, CeO2, ...) contours (needs [pyFAI](https://pyfai.readthedocs.io/en/master/)).
  - 2022-04-28 Update: Changed contour line generation to accept a list of 2-theta values as input.
  - 2022-04-27 Update: Added support for [SACLA](https://sacla.xfel.jp/?lang=en) MPCCD Detectors (central hole geometry).
  - 2022-04-25 Bugfix: Calculation of the beamcenter (rotation and tilt).
</details>

<details>
<summary>After the update</summary>
 
  - Sometimes I might change the name of a parameter and you will get a warning message upon startup looking something like this: _WARNING: "conic_ref_min_int" is not a valid key_! Either that key is no longer in use or its name got changed and is now reset to the default value. The settings file is updated and the warning should no longer appear after restart. Apart from this, your edited settings file will not be altered after updating.
</details>

>[!IMPORTANT]
>## Known bugs and limitations:
>  - On Windows: Switching Dark/Light mode requires restart to change the window frame color.
>  - As of now, there is no consistency check performed on the imported .json file.
>  - The projections fail at a combined angle (rotation + tilt) of 90 degrees and beyond.
>  - If the program crashes (libc++abi) when opening the export window (on a Mac) please update PyQt6 to the latest version.
>  - A crash on startup due to an f-string SyntaxError can be circumvented by upgrading python to version 3.12 or later.

## Examples
#### A PILATUS3 300K detector and a Rubrene sample.
<img src="https://github.com/LennardKrause/xrdPlanner/blob/main/examples/Figure_2_example_light.png" width="50%">

#### A rotated EIGER2 9M detector and a Rubrene sample (darkmode).
<img src="https://github.com/LennardKrause/xrdPlanner/blob/main/examples/Figure_1_example_dark.png" width="50%">

#### Export window
##### Use the export window to select and build up the available detector / beamstop bank and review/change parameters.
<img src="https://github.com/LennardKrause/xrdPlanner/blob/main/examples/Figure_3_export_window.png" width="50%">

#### Detector db editor
##### Define and add the detector you collect your data with.
<img src="https://github.com/LennardKrause/xrdPlanner/blob/main/examples/Figure_4_detector_db_editor.png" width="50%">

#### PXRD pattern plot
##### Look at a 1D PXRD pattern plot of the calculated intensities (cif-file), displayed using gaussian peak shapes and the estimated peak broadening (FWHM).
<img src="https://github.com/LennardKrause/xrdPlanner/blob/main/examples/Figure_5_PXRD.png" width="50%">

#### Instrumental broadening
##### Estimate FWHM (H, in degrees) using the following parameters:
 - Estimated effective thickness of the detector absorption layer depending on sensor material and incident wavelength
 - Detector sensor layer material
 - Estimated beam divergence
 - Expected bandwidth of the incident X-ray beam
 - Intersection of the sample size and the beam size as the diameter

> **The interested user is referred to the article:<br>**
>  <i>On the resolution function for powder diffraction with area detectors</i><br>
>  [D. Chernyshov <i> et al., Acta Cryst.</i> (2021). <b>A77</b>, 497-505](https://doi.org/10.1107/S2053273321007506)
<img src="https://github.com/LennardKrause/xrdPlanner/blob/main/examples/Figure_6_FWHM.png" width="50%">

## Hotkeys

| Key        | Action                    |
| ---------- | ------------------------- |
| F1         | Show About window         |
| F2         | Show Geometry conventions |
| F3         | Show Hotkeys              |
| x          | Show PXRD pattern         |
| f          | Show FWHM window          |
| #          | *Display units*           |
| t          | 2-Theta                   |
| d          | d-spacing                 |
| q          | q-space                   |
| s          | $sin(\theta)/\lambda$     |
| #          | *Label Position*          |
| up arrow   | Top left                  |
| down arrow | Bottom left               |
| #          | *Toggle Overlay*          |
| p          | Show polarisation         |
| a          | Show solid angle          |
| h          | Highlight / Transparency  |
| u          | Toggle unit hover display |
| #          | *Cycle colormaps*         |
| c          | Next                      |
| Shift + c  | Previous                  |

## Settings file documentation
<details>
<summary>geo - startup defaults</summary>
 
    det_type = 'EIGER2'  # [str]  Pilatus3 / Eiger2 / etc.
                         #        -> Detector menu entry
    det_size = '4M'      # [str]  300K 1M 2M 6M / 1M 4M 9M 16M
                         #        -> Detector submenu entry
    ener = 21            # [keV]  Beam energy
    dist = 75            # [mm]   Detector distance
    yoff = 0             # [mm]   Detector offset (vertical)
    xoff = 0             # [mm]   Detector offset (horizontal)
    rota = 25            # [deg]  Detector rotation
    tilt = 0             # [deg]  Detector tilt
    bssz = 'None'        # [mm]   Current beamstop size (or 'None')
    bsdx = 40            # [mm]   Beamstop distance
    unit = 1             # [0-3]  Contour legend
                         #         0: 2-Theta
                         #         1: d-spacing
                         #         2: q-space
                         #         3: sin(theta)/lambda
    reference = 'None'   # [str]  Plot reference contours
                         #          pick from pyFAI
    darkmode = False     # [bool] Darkmode
    colormap = 'viridis' # [cmap] Contour colormap
    bs_list = [1.5,      # [list] Available beamstop sizes
               2.0,
               2.5,
               3.0,
               5.0]
</details>

<details>
<summary>plo - plot settings</summary>
 
    # - geometry contour section - 
    conic_tth_min = 5               # [int]    Minimum 2-theta contour line
    conic_tth_max = 100             # [int]    Maximum 2-theta contour line
    conic_tth_num = 15              # [int]    Number of contour lines
    beamcenter_marker = 'o'         # [marker] Beamcenter marker
    beamcenter_size = 6             # [int]    Beamcenter size
    poni_marker = 'x'               # [marker] Poni marker
    poni_size = 8                   # [int]    Poni size
    conic_linewidth = 2.0           # [float]  Contour linewidth
    conic_label_size = 14           # [int]    Contour label size
    
    # - reference contour section - 
    conic_ref_linewidth = 2.0       # [float]  Reference contour linewidth
    conic_ref_num = 100             # [int]    Number of reference contours
    conic_ref_cif_int = 0.01        # [float]  Minimum display intensity (cif)
    conic_ref_cif_kev = 10.0        # [float]  Energy [keV] for intensity calculation
    conic_ref_cif_irel = True       # [bool]   Linewidth relative to intensity
    conic_ref_cif_lw_min = 0.1      # [float]  Minimum linewidth when using irel
    conic_ref_cif_lw_mult = 3.0     # [float]  Linewidth multiplier when using irel
    conic_hkl_show_int = False      # [bool]   Show intensity in hkl tooltip
    conic_hkl_label_size = 14       # [int]    Font size of hkl tooltip
    
    # - module section - 
    det_module_alpha = 0.20         # [float]  Detector module alpha
    det_module_width = 1            # [int]    Detector module border width
    
    # - general section - 
    conic_steps = 100               # [int]    Conic resolution
    plot_size = 0                   # [int]    Plot size, px (0 for auto)
    plot_size_fixed = True          # [bool]   Fix window size
    unit_label_size = 16            # [int]    Label size, px
    polarisation_fac = 0.99         # [float]  Horizontal polarisation factor
    show_polarisation = False       # [bool]   Show polarisation overlay
    show_solidangle = False         # [bool]   Show solid angle overlay
    overlay_resolution = 300        # [int]    Overlay resolution
    overlay_toggle_warn = True      # [bool]   Overlay warn color threshold
    
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
    enable_slider_bsdx = True       # [bool]   Show beamstop distance slider
    
    # - slider labels - 
    slider_label_ener = 'Energy\n[keV]'            # [str] Label for energy slider
    slider_label_dist = 'Distance\n[mm]'           # [str] Label for distance slider
    slider_label_rota = 'Rotation\n[\u02da]'       # [str] Label for rotation slider
    slider_label_voff = 'Vertical\noffset\n[mm]'   # [str] Label for vertical offset slider
    slider_label_hoff = 'Horizontal\noffset\n[mm]' # [str] Label for horizontal offset slider
    slider_label_tilt = 'Tilt\n[\u02da]'           # [str] Label for tilt slider
    slider_label_bsdx = 'Beamstop\ndistance\n[mm]' # [str] Label for beamstop distance slider
            
    # - update/reset - 
    update_settings = True          # [bool]   Update settings file after load
    update_det_bank = True          # [bool]   Update detector bank after load
    reset_settings = False          # [bool]   Reset settings file
    reset_det_bank = False          # [bool]   Reset detector bank
    
    # - debug/testing -
    set_debug = False               # [bool]   Debug mode

</details>

<details>
<summary>thm - theme</summary>
 
    color_dark = '#404040'                # [color]  Global dark color
    color_light = '#EEEEEE'               # [color]  Global light color
    
    # light mode
    light_conic_label_fill = '#FFFFFF'    # [color]  Contour label fill color
    light_conic_ref_color = '#DCDCDC'     # [color]  Reference contour color
    light_beamstop_color = '#FF000080'    # [color]  Beamstop color
    light_beamstop_edge_color = '#FF0000' # [color]  Beamstop edge color
    light_det_module_color = '#404040'    # [color]  Detector module border color
    light_det_module_fill = '#404040'     # [color]  Detector module background color
    light_plot_bg_color = '#FFFFFF'       # [color]  Plot background color
    light_unit_label_color = '#808080'    # [color]  Label color
    light_unit_label_fill = '#FFFFFF'     # [color]  Label fill color
    light_slider_border_color = '#808080' # [color]  Slider frame border color
    light_slider_bg_color = '#AAC0C0C0'   # [color]  Slider frame background color
    light_slider_bg_hover = '#C0C0C0'     # [color]  Slider frame hover color
    light_slider_label_color = '#000000'  # [color]  Slider frame label color
    
    # dark mode
    dark_conic_label_fill = '#000000'     # [color]  Contour label fill color
    dark_conic_ref_color = '#505050'      # [color]  Reference contour color
    dark_beamstop_color = '#FF000080'     # [color]  Beamstop color
    dark_beamstop_edge_color = '#FF0000'  # [color]  Beamstop edge color
    dark_det_module_color = '#EEEEEE'     # [color]  Detector module border color
    dark_det_module_fill = '#EEEEEE'      # [color]  Detector module background color
    dark_plot_bg_color = '#000000'        # [color]  Plot background color
    dark_unit_label_color = '#C0C0C0'     # [color]  Label color
    dark_unit_label_fill = '#000000'      # [color]  Label fill color
    dark_slider_border_color = '#202020'  # [color]  Slider frame border color
    dark_slider_bg_color = '#AA303030'    # [color]  Slider frame background color
    dark_slider_bg_hover = '#303030'      # [color]  Slider frame hover color
    dark_slider_label_color = '#C0C0C0'   # [color]  Slider frame label color

</details>

<details>
<summary>lmt - limits</summary>

    # energy [keV]
    ener_min =  5    # [int] Energy minimum [keV]
    ener_max =  100  # [int] Energy maximum [keV]
    ener_stp =  1    # [int] Energy step size [keV]
    
    # distance [mm]
    dist_min =  40   # [int] Distance minimum [mm]
    dist_max =  1000 # [int] Distance maximum [mm]
    dist_stp =  1    # [int] Distance step size [mm]
    
    # X-offset [mm]
    xoff_min = -150  # [int] Horizontal offset minimum [mm]
    xoff_max =  150  # [int] Horizontal offset maximum [mm]
    xoff_stp =  1    # [int] Horizontal offset step size [mm]
    
    # Y-offset [mm]
    yoff_min = -250  # [int] Vertical offset minimum [mm]
    yoff_max =  250  # [int] Vertical offset maximum [mm]
    yoff_stp =  1    # [int] Vertical offset step size [mm]
    
    # rotation [deg]
    rota_min = -45   # [int] Rotation minimum [deg]
    rota_max =  45   # [int] Rotation maximum [deg]
    rota_stp =  1    # [int] Rotation step size [deg]
    
    # tilt [deg]
    tilt_min = -40   # [int] Tilt minimum [deg]
    tilt_max =  40   # [int] Tilt maximum [deg]
    tilt_stp =  1    # [int] Tilt step size [deg]

    bsdx_min = 5     # [int] Beamstop distance minimum [mm]
    bsdx_max = 1000  # [int] Beamstop distance maximum [mm]
    bsdx_stp = 1     # [int] Beamstop distance step size [mm]

</details>

<details>
<summary>Detector db entries</summary>

The detector_db.json is stored as a dictionary
- key: detector name e.g. PILATUS3
- value: dictionary {Entry:Value,}, see table below
 
 | Entry |  Value  | Hint |
 |-------|---------|------|
 |  pxs  | 172e-3  | [mm] Pixel size
 |  hmp  | 487     | [px] Module size (horizontal)
 |  vmp  | 195     | [px] Module size (vertical)
 |  hgp  | 7       | [px] Gap between modules (horizontal)
 |  vgp  | 17      | [px] Gap between modules (vertical)
 |  cbh  | 0       | [px] Central beam hole

The size Entry is a dictionary {key:value,}
- key: detector size / type, e.g. 300K
- value: list [hmn, vmn]
  - hmn: [int]  Number of modules (horizontal)
  - vmn: [int]  Number of modules (vertical)
 
#### Example of the PILATUS3 entry
    "PILATUS3": {
        "pxs": 0.172,
        "hmp": 487,
        "vmp": 195,
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

</details>

## Example code
<details>
<summary>Example code for adding xrdPlanner as a widget into an existing GUI</summary>
 
#### xrdPlanner uses its own menu bar, setting the GUI as the parent for xrdPlanner makes it add its menus to the parents menu bar, and likely more in the future.

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
</details>

## I hope this turns out to be useful for someone!

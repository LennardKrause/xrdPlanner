from setuptools import setup
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / 'README.md').read_text()

setup(
    name='xrdPlanner',
    version='1.0.5',
    description='A tool to project X-ray diffraction cones on a detector screen at different geometries (tilt, rotation, offset) and X-ray energies.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/LennardKrause/xrdPlanner/',
    author='Lennard Krause',
    author_email='lkrause@chem.au.dk',
    license='GNU GENERAL PUBLIC LICENSE',
    entry_points = {'console_scripts': ['xrdPlanner=xrdPlanner.run_xrdPlanner:main'],},
    install_requires=['numpy>=1.24',
                      'pyqtgraph>=0.13',
                      'Dans_Diffraction>=2.2.3',
                      'PyQt6>=6.5',
                      'pyFAI>=2023.5',
                      ],

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Topic :: Education',
        'Topic :: Scientific/Engineering',
        'Operating System :: OS Independent',
    ],
)

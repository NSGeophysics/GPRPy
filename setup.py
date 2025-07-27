import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

    
setuptools.setup(
    name="gprpy",
    version="1.0.14",
    author="Alain Plattner",
    author_email="plattner@alumni.ethz.ch",
    description="GPRPy - open source ground penetrating radar processing and visualization",
    entry_points={'console_scripts': ['gprpy = gprpy.__main__:main']},
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NSGeophysics/GPRPy",
    packages=['gprpy'],
    package_data={'gprpy': ['exampledata/GSSI/*.DZT',
                            'exampledata/GSSI/*.txt',
                            'exampledata/SnS/ComOffs/*.xyz',
                            'exampledata/SnS/ComOffs/*.DT1',
                            'exampledata/SnS/ComOffs/*.HD',
                            'exampledata/SnS/WARR/*.DT1',
                            'exampledata/SnS/WARR/*.HD',
                            'exampledata/pickedSurfaceData/*.txt',
                            'examplescripts/*.py',
                            'toolbox/splashdat/*.png',
                            'toolbox/*.py',
                            'irlib/*.py',
                            'irlib/external/*.py']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=['tqdm','numpy','scipy','matplotlib','Pmw','pyevtk']
)

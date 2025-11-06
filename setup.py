import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gprpy-lib",
    version="2.0.0",
    author="Alain Plattner",
    author_email="plattner@alumni.ethz.ch",
    description="GPRPy - a Python library for ground penetrating radar processing and visualization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NSGeophysics/GPRPy",
    packages=['gprpy_lib'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=['numpy','scipy','matplotlib','pyevtk']
)

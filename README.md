# GPRPy
Ground Penetrating Radar processing and visualization software for python 3.

![alt text](https://github.com/NSGeophysics/GPRPy/blob/master/gprPy.png)

## Simplemost installation

1) Download the GPRPy software for example by clicking on the green "clone or download" 
button and then "Download ZIP". Save the file somewhere on your computer and extract the 
zip folder (typically by double-clicking on it).

2) Install anaconda python (the python 3.6 version) from https://www.anaconda.com/download/

3) Once the installation finished, open Anaconda Prompt (on Windows: click on Start, 
then enter "Anaconda Prompt", without the quotation marks into the 
"Search programs and files" field).

4) In the Anaconda Prompt, type the following and press enter afterward:
`pip install numpy scipy matplotlib tqdm Pmw pyevtk`
This will install the dependencies, if they are not already installed

5) In the Anaconda Prompt, change into the directory where you downloaded the GPRPy files.
This is usually through a command like for example
`cd Desktop\GPRPy`
if you downloaded GPRPy directly onto your desktop

6) In Anaconda Prompt, type and press Enter afterward:
`python gprpyGUI.py`
This should open the graphical user interface


## General comments

#### Requires:

numpy, scipy, matplotlib, pickle, tkinter, struct, re, tqdm, Pmw, pyevtk

#### To run the GUI:

`python gprpyGUI.py`

or

`python3 gprpyGUI.py`





# GPRPy
Ground Penetrating Radar processing and visualization software for python 3.

![alt text](https://github.com/NSGeophysics/GPRPy/blob/master/profileGUI.png)

![alt text](https://github.com/NSGeophysics/GPRPy/blob/master/CWGUI.png)


## Simplemost installation

1) Download the GPRPy software for example by clicking on the green "clone or download" 
button and then "Download ZIP". Save the file somewhere on your computer and extract the 
zip folder (typically by double-clicking on it).

2) Install Python 3.7 for example from https://conda.io/miniconda.html

3) Once the installation finished, open a command prompt that can run Python
(on Windows: click on Start, then enter "Anaconda Prompt", without the 
quotation marks into the "Search programs and files" field).

4) In the command prompt, change to the directory  where you downloaded the GPRPy files.
This is usually through a command like for example\
`cd Desktop\GPRPy`\
if you downloaded GPRPy directly onto your desktop. Then type the following and press enter
afterward:\
`python installMigration.py`\
Then type the following and press enter afterward:\
`pip install .`\
(don't forget the period at the end).

## Running the software
To run the **profile graphical user interface**, you can either run the included shell script 
by, in the command prompt in the GPRPy folder, typing and pressing enter\
`profile.sh`

Alternatively, you can open the profile graphical user interface from any folder by running in the command prompt:\
`python -c "import gprpy.__main__" p`
 
To run the **CMP / WARR graphical user interface**, you can either run the included shell script 
by, in the command prompt in the GPRPy folder, typing and pressing enter\
`cmpwarr.sh`

Alternatively, you can open the CMP / WARR graphical user interface from any folder by running in the command prompt:\
`python -c "import gprpy.__main__" c`

If you just run\
`python -c "import gprpy.__main__"`\
then GPRPy will ask you if you want to run it in the profile or CMP/WARR mode.


## In case of trouble

If you have several versions of python installed, for example on a Mac or Linux system, 
replace in the previous commands\
`python` with `python3`\
and\
`pip` with `pip3`

If you have any troubles getting the software running, please send me an email or open an issue
on GitHub and I will help you getting it running.



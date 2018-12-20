# GPRPy
Open-source Ground Penetrating Radar processing and visualization software.

![Profile GUI](profileGUI.png)

![CMP/WARR GUI](https://github.com/NSGeophysics/GPRPy/blob/master/CWGUI.png)


## Simplemost installation
1) Download the GPRPy software for example by clicking on the green "clone or download" 
button and then "Download ZIP". Save the file somewhere on your computer and extract the 
zip folder.

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
After installation, you can run the script from the Anaconda Prompt (or your Python-enabled prompt) by running
either\
`gprpy`\
or\
`python -m gprpy`

The first time you run GPRPy it could take a while to initialize.
GPRPy will ask you if you want to run the profile [p] or WARR / CMP [c] user interface.
Type\
`p`\
and then enter for profile, or\
`c`\
and then enter for CMP / WARR.

You can also directly select one by running either\
`gprpy p`\
or\
`gprpy c`\
or\
`python -m gprpy p`\
or\
`python -m gprpy c`


## Running automatically generated scripts
To run automatically generated scripts, open the command prompt that can run python (for example Anaconda Prompt), switch to the folder with the automatically generated script and run\
`python myscriptname.py`\
where myscriptname.py is the name of your automatically generated script.  


## In case of trouble
If you have several versions of python installed, for example on a Mac or Linux system, 
replace, in the commands shown earlier,
`python` with `python3`\
and\
`pip` with `pip3`

If you have any troubles getting the software running, please send me an email or open an issue
on GitHub and I will help you getting it running.


## Uninstalling GPRPy
To uninstall GPRPy, simply run, in the (Anaconda) command prompt\
`pip uninstall gprpy`


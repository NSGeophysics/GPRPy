https://github.com/NSGeophysics/GPRPy/archive/master.zip
# GPRPy
Open-source Ground Penetrating Radar processing and visualization software.

![Profile GUI](profileGUI.png)

![CMP/WARR GUI](CWGUI.png)


## Simplemost installation
1) Download the GPRPy software for example from [https://github.com/NSGeophysics/GPRPy/archive/master.zip](https://github.com/NSGeophysics/GPRPy/archive/master.zip). Save the file somewhere on your computer and extract the 
zip folder. <br/>
As an **alternative**, you can install git from [https://git-scm.com/](https://git-scm.com/), then run in a command prompt:<br/>
`git clone https://github.com/NSGeophysics/GPRPy.git`<br/>
The advantage of the latter is that you can easily update your software by running from the GPRPy folder in a command prompt:<br/>
`git pull origin master`

2) Install Python 3.7 for example from [https://conda.io/miniconda.html](https://conda.io/miniconda.html)

3) Once the installation finished, open a command prompt that can run Python
(on Windows: click on Start, then enter "Anaconda Prompt", without the 
quotation marks into the "Search programs and files" field).

4) In the command prompt, change to the directory  where you downloaded the GPRPy files.
This is usually through a command like for example<br/>
`cd Desktop\GPRPy`<br/>
if you downloaded GPRPy directly onto your desktop. Then type the following and press enter
afterward:<br/>
`python installMigration.py`<br/>
Then type the following and press enter afterward:<br/>
`pip install .`<br/>
(don't forget the period at the end).


## Running the software
After installation, you can run the script from the Anaconda Prompt (or your Python-enabled prompt) by running
either<br/>
`gprpy`<br/>
or<br/>
`python -m gprpy`

The first time you run GPRPy it could take a while to initialize.
GPRPy will ask you if you want to run the profile [p] or WARR / CMP [c] user interface.
Type<br/>
`p`<br/>
and then enter for profile, or<br/>
`c`<br/>
and then enter for CMP / WARR.

You can also directly select one by running either<br/>
`gprpy p`<br/>
or<br/>
`gprpy c`<br/>
or<br/>
`python -m gprpy p`<br/>
or<br/>
`python -m gprpy c`


## Running automatically generated scripts
To run automatically generated scripts, open the command prompt that can run python (for example Anaconda Prompt), switch to the folder with the automatically generated script and run<br/>
`python myscriptname.py`<br/>
where myscriptname.py is the name of your automatically generated script.  


## In case of trouble
If you have several versions of python installed, for example on a Mac or Linux system, 
replace, in the commands shown earlier,
`python` with `python3`<br/>
and<br/>
`pip` with `pip3`

If you have any troubles getting the software running, please send me an email or open an issue
on GitHub and I will help you getting it running.


## Uninstalling GPRPy
To uninstall GPRPy, simply run, in the (Anaconda) command prompt<br/>
`pip uninstall gprpy`


# GPRPy
Open-source Ground Penetrating Radar processing and visualization software.

![Profile GUI](profileGUI.png)

![CMP/WARR GUI](CWGUI.png)


## Simplemost installation

**In the following instructions, if you use Windows, use the comands `python` and `pip`. If you use Mac or Linux, use the commands `python3` and `pip3` instead.**

1) Download the GPRPy software from 
   [https://github.com/NSGeophysics/GPRPy/archive/master.zip](https://github.com/NSGeophysics/GPRPy/archive/master.zip). <br/>
   Save the file somewhere on your computer and extract the zip folder. <br/>
   As an **alternative**, you can install git from [https://git-scm.com/](https://git-scm.com/), then run in a command prompt:<br/>
   `git clone https://github.com/NSGeophysics/GPRPy.git`<br/>
   The advantage of the latter is that you can easily update your software by running from the GPRPy folder in a command prompt:<br/>
   `git pull origin master`

2) Install Python 3.7 for example from [https://conda.io/miniconda.html](https://conda.io/miniconda.html)

3) Once the installation finished, open a command prompt that can run Python <br/>
   On Windows: click on Start, then enter "Anaconda Prompt", without the quotation marks into the "Search programs and files" field. On Mac or Linux, open the regular terminal.

4) In the command prompt, change to the directory  where you downloaded the GPRPy files.
   This is usually through a command like for example<br/>
   `cd Desktop\GPRPy`<br/>
   if you downloaded GPRPy directly onto your desktop. Then type the following and press enter afterward:<br/>
   `python installMigration.py`<br/>
   Then type the following and press enter afterward:<br/>
   `pip install .`<br/>
   **don't forget the period "." at the end of the `pip install` command**


## Running the software
After installation, you can run the script from the Anaconda Prompt (or your Python-enabled prompt) by running either<br/>
`gprpy`<br/>
or<br/>
`python -m gprpy`

The first time you run GPRPy it could take a while to initialize. GPRPy will ask you if you want to run the profile [p] or WARR / CMP [c] user interface. Type<br/>
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
If you have several versions of python installed, for example on a Mac or Linux system, replace, in the commands shown earlier,
`python` with `python3`<br/>
and<br/>
`pip` with `pip3`

If you have any troubles getting the software running, please send me an email or open an issue on GitHub and I will help you getting it running.


## Uninstalling GPRPy
To uninstall GPRPy, simply run, in the (Anaconda) command prompt<br/>
`pip uninstall gprpy`

## News
Follow [@GPRPySoftware](https://twitter.com/GPRPySoftware) on twitter to hear about news and updates.
Recent tweets:

<blockquote class="twitter-tweet" data-lang="en"><p lang="en" dir="ltr">Fixed small issue that led to multiples when picking points in profile. Thanks Marcus Pacheco for pointing it out! If you use picking in profile mode, please update to version 1.0.3 (uninstall the old version before).</p>&mdash; GPRPy (@GPRPySoftware) <a href="https://twitter.com/GPRPySoftware/status/1139243564469313536?ref_src=twsrc%5Etfw">June 13, 2019</a></blockquote>

<blockquote class="twitter-tweet" data-lang="en"><p lang="en" dir="ltr">I will post updates, changes, and GPRPy news here.</p>&mdash; GPRPy (@GPRPySoftware) <a href="https://twitter.com/GPRPySoftware/status/1089246592786485251?ref_src=twsrc%5Etfw">January 26, 2019</a></blockquote>

<blockquote class="twitter-tweet" data-lang="en"><p lang="en" dir="ltr">GPRPy is a free ground penetrating radar processing and visualization software developed at the University of Alabama. You can download it and install it following the instructions here: <a href="https://nsgeophysics.github.io/GPRPy/">nsgeophysics.github.io/GPRPy/</a></p>&mdash; GPRPy (@GPRPySoftware) <a href="https://twitter.com/GPRPySoftware/status/1088806792191197188?ref_src=twsrc%5Etfw">January 25, 2019</a></blockquote>



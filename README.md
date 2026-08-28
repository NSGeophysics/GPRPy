# GPRPy
Open-source Ground Penetrating Radar processing and visualization software. Supported by the National Science Foundation under grants [EAR-1550732](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1550732) and [EAR-2022671](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2022671)

Please cite https://library.seg.org/doi/10.1190/tle39050332.1

![Profile GUI](profileGUI.png)

![CMP/WARR GUI](CWGUI.png)


## Simplemost installation using anaconda python

1) Download the GPRPy software from 
   [https://github.com/NSGeophysics/GPRPy/archive/master.zip](https://github.com/NSGeophysics/GPRPy/archive/master.zip). <br/>
   Save the file somewhere on your computer and extract the zip folder. <br/>
   As an **alternative**, you can install git from [https://git-scm.com/](https://git-scm.com/), then run in a command prompt:<br/>
   `git clone https://github.com/NSGeophysics/GPRPy.git`<br/>
   The advantage of the latter is that you can easily update your software by running from the GPRPy folder in a command prompt:<br/>
   `git pull origin master`

2) Install anaconda or miniconda [https://www.anaconda.com/download](https://www.anaconda.com/download)

3) Once the installation finished, open a command prompt:
   On Windows: click on Start, then enter "Anaconda Prompt", without the quotation marks into the "Search programs and files" field.
   On Mac or Linux, open the regular terminal. In the command prompt/terminal, run<br/>
   `conda activate`<br/>

5) Change to the directory  where you downloaded the GPRPy files.
   This is usually through a command like for example<br/>
   `cd Desktop\GPRPy`<br/>
   if you downloaded GPRPy directly onto your desktop. Then type the following and press enter afterward:<br/>
   `python installMigration.py`<br/>
   Then type the following and press enter afterward:<br/>
   `pip install .`<br/>
   **don't forget the period "." at the end of the `pip install .` command**


## Running the software
After installation, you can run the script from the Anaconda Prompt (or terminal after running `conda activate`) by running<br/>
`gprpy`<br/>

The first time you run GPRPy it could take a while to initialize. GPRPy will ask you if you want to run the profile [p] or WARR / CMP [c] user interface. Type<br/>
`p`<br/>
and then enter for profile, or<br/>
`c`<br/>
and then enter for CMP / WARR.

You can also directly select one by running either<br/>
`gprpy p`<br/>
or<br/>
`gprpy c`<br/>


## Running automatically generated scripts
To run automatically generated scripts, open Anaconda Prompt (or the terminal and run `conda activate`), switch to the folder with the automatically generated script and run<br/>
`python myscriptname.py`<br/>
where myscriptname.py is the name of your automatically generated script. Note that the data need to be in the directories referenced by `myscriptname.py`. 


## Uninstalling GPRPy
To uninstall GPRPy, simply run, in the Anaconda Prompt (or terminal after `conda activate`)<br/>
`pip uninstall gprpy`

## News
Follow [@GPRPySoftware](https://twitter.com/GPRPySoftware) on twitter to hear about news and updates.
Recent tweets:

<blockquote class="twitter-tweet" data-lang="en"><p lang="en" dir="ltr">Fixed small issue that led to multiples when picking points in profile. Thanks Marcus Pacheco for pointing it out! If you use picking in profile mode, please update to version 1.0.3 (uninstall the old version before).</p>&mdash; GPRPy (@GPRPySoftware) <a href="https://twitter.com/GPRPySoftware/status/1139243564469313536?ref_src=twsrc%5Etfw">June 13, 2019</a></blockquote>

<blockquote class="twitter-tweet" data-lang="en"><p lang="en" dir="ltr">I will post updates, changes, and GPRPy news here.</p>&mdash; GPRPy (@GPRPySoftware) <a href="https://twitter.com/GPRPySoftware/status/1089246592786485251?ref_src=twsrc%5Etfw">January 26, 2019</a></blockquote>

<blockquote class="twitter-tweet" data-lang="en"><p lang="en" dir="ltr">GPRPy is a free ground penetrating radar processing and visualization software developed at the University of Alabama. You can download it and install it following the instructions here: <a href="https://nsgeophysics.github.io/GPRPy/">nsgeophysics.github.io/GPRPy/</a></p>&mdash; GPRPy (@GPRPySoftware) <a href="https://twitter.com/GPRPySoftware/status/1088806792191197188?ref_src=twsrc%5Etfw">January 25, 2019</a></blockquote>



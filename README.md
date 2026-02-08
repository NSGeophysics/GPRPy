# GPRPy-lib

GPRPy-lib is a Python library for processing and visualizing Ground Penetrating Radar (GPR) data. It is designed to be used in a Jupyter Notebook environment, providing a flexible and interactive way to work with GPR data.

It's a fork of the original [GPRPy software](https://github.com/NSGeophysics/GPRPy), with some modifications and improvements. 

## Installation

**In the following instructions, if you use Windows, use the comands `python` and `pip`. If you use Mac or Linux, use the commands `python3` and `pip3` instead.**

1) Download the GPRPy software from 
   [https://github.com/allanspadini/GPRPy-lib](https://github.com/allanspadini/GPRPy-lib/archive/refs/heads/master.zip). <br/>
   Save the file somewhere on your computer and extract the zip folder. <br/>
   `git clone https://github.com/allanspadini/GPRPy-lib.git`<br/>
   The advantage of the latter is that you can easily update your software by running from the GPRPy folder in a command prompt:<br/>
   `git pull origin master`

2) Install Python 3.7 or higher. You can obtain it for example from [https://conda.io/miniconda.html](https://conda.io/miniconda.html)

3) Once the installation finished, open a command prompt that can run Python <br/>
   On Windows: click on Start, then enter "Anaconda Prompt", without the quotation marks into the "Search programs and files" field. On Mac or Linux, open the regular terminal.

4) In the command prompt, change to the directory  where you downloaded the GPRPy files.
   This is usually through a command like for example<br/>
   `cd Desktop\GPRPy-lib`<br/>
   if you downloaded GPRPy-lib directly onto your desktop. Then type the following and press enter afterward:<br/>
   `pip install .`<br/>
   **don't forget the period "." at the end of the `pip install` command**


## Usage Example

Here's a quick example of how to use the GPRPy library in a Jupyter Notebook:

```python
import gprpy_lib as gpr
import numpy as np

# 1. Loading Data
filepath = 'gprpy/exampledata/GSSI/FILE____032.DZT'
data, info = gpr.gprIO_DZT.readdzt(filepath)

# 2. Initial Data Visualization
profilePos = info['rhf_position'] + np.linspace(0.0, data.shape[1] / info['rhf_spm'], data.shape[1])
twtt = np.linspace(0, info['rhf_range'], info['rh_nsamp'])
gpr.plot_profile(data, profilePos, twtt)

# 3. Data Processing
# Apply dewow filter
data_dewowed = gpr.dewow(data, window=50)

# Apply smoothing
data_smoothed = gpr.smooth(data_dewowed, window=5)

# Remove mean trace
data_processed = gpr.remMeanTrace(data_smoothed, ntraces=10)

# 4. Processed Data Visualization
gpr.plot_profile(data_processed, profilePos, twtt)
```

For a more detailed example, please see the Jupyter Notebook in the `examples` directory.

## Uninstalling GPRPy
To uninstall GPRPy, simply run, in the (Anaconda) command prompt<br/>
`pip uninstall gprpy`



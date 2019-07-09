# ePhysDataPlotta

is a collection of mostly high-level Matlab routines for plotting electrophysiological data. The code makes use of in-house formats for time information ('time stamp lists', tsl).

Please note that the code in this repository is not self-sufficient, you'll additionally need the following repositories:
* fileIO
* etslfunc
* graphics
* sampledSeries
* utilities
 
Matlab toolboxes required (for some code files):
* Signal Processing
* Statistics and Machine Learning
* Curve Fitting

Some of the routines:

### gfplotta
'Gap-free' data plotter - produces plots of multi-channel, continuously recorded electrophysiological data in a configurable and ressource-saving way. 

### swplotta
Plotting function which creates plots of 'sweeps' of raw time series plus optionally derived data, currently from current clamp recordings with current injection 

![snapshot](/doc/IVCurve.png)

### PSCplotta
Plotting function for post-synaptic currents (PSCs) which is designed to create i) raw data plots and ii) overlay plots of cutouts thereof 

### periburc 
Performs analysis of and plots current clamp data with hyperpolarizing current injection, computing cell parameters (R, C). It is able to do so for raw data partitioned into bursts and silent periods by threshdetgui. 

### peribuexcit
Analyzes and plots responses of cells in current clamp mode to sinusodal ramp current injection

## General note on repositories in the ExpAnesth organization
The code in these repositories provides basic tools for the analysis of electrophysiological time series to members of the Section of Experimental Anesthesiology, Department of Anesthesiology, University Hospital of Tuebingen. Except where noted, code was written by Harald Hentschke. It has been designed primarily for in-house use by individuals who were instructed on its scope and limitations. Also, a substantial proportion of the code has been developed and extended over a time span of >10 years. In detail,

* the implementation of algorithms reflects the evolution of Matlab itself, that is, code that had been initially developed on older versions of Matlab does not necessarily feature newer techniques such as the new automatic array expansion as introduced in Matlab Release 2016b
* nonetheless, all code has been tested to run on Matlab R2018b
* while most m-files contain ample comments, documentation exists only for a few repositories
* checks of user input are implemented to varying degrees
* the code will be improved, updated and documented when and where the need arises

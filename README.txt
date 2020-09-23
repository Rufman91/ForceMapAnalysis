QUICKSTART GUIDE Force Map Analysis
%%% by Manuel Rufin, program version 23.09.2020%%%

This data analysis program is based around the MATLAB Command Window.

To set up the program, add the force map analysis folder to the MATLAB search path.
The following Toolboxes might be needed:
- Statistics and Machine Learning Toolbox
- Deep Learning Toolbox
- Signal Processing Toolbox
- Parallel Computing Toolbox
- Image Processing Toolbox
- Curve Fitting Toolbox

To create an Experiment-file, which will be the handle class object in your MATLAB workspace you actually work with, call the constructor method

>> ExperimentNameInWorkspace = Experiment();

choose a very short name for ExperimentNameInWorkspace, as every function/data will be a method/property of this object (e.g. 'E','Ex').
The programm will lead you through the process of creating your experiment: 
what kind of measurements have been done? (Force Maps, KPFM) --> Number of specimen/ multiple measurements per specimen/ name of experiment --
--> choose your *.jpk-force-map files (currently only works for the newest jpk force map files (23.09.2020))

The data will then be extracted from the jpk files, which could take a few minutes, and your experiment will be saved to the location you chose at the beginning.

You can now already call several informations, called properties, from you experiment class. Assuming the experiment is called 'E' in the workspace , examples would be:

>> E.ForceMapNames		% List of all force map names
>> E.NumFiles			% Overall number of files
>> E.ForceMapFolders		% List of force maps .mat-file folders
>> E.ExperimentName 		% Name of the experiment

ect. 
As a general rule of thump: functions, in the context of object classes called 'methods', are written in lower case letter, multiple words connected by underlines

>> E.save_experiment

while object properties, such as lists, arrays, structs, or even neural networks are called with capitalized names, without additional spacing

>> E.ExperimentName

Force maps and surface potential maps are themselves object classes in this programm and are saved under the property name E.FM and E.SPM, respectivley.
Single instances of these objects have to be called with curved brackets:

>> E.FM{7}
>> E.SPM{2}

Both of these subclasses come with their own methods and properties. For example:

>> E.FM{3}.FibDiam			% gives the calculated fibril diameter (if it has been calculated already)

or

>> E.FM{1}.quality_control_oliver_pharr		% creates a few plots to validate Oliver-Pharrian analysis of indentation modulus (if force map analysis has been run already)

Most properties will be empty and most methods will produce errors, though, until you run data analysis on the data of interest.
This is done for, respectively, force maps and surface potential maps by

>> E.force_map_analysis_fibril(CPOption, EModOption)
>> E.surface_potential_analysis_fibril()

the options for the force map analysis are explained in the functions header comments. At the moment, before you can run Oliver-Pharrian analysis, you need to manually assign the reconstructed tip data
to the right experimet property. The tip deconvolution has to be done in the old script by Orestis Andriotis.

Load the struct named 'Tip', that is saved to the *.mat file from the old tip deconvolution script, into the workplace
and assign it to the Experiment classes 'CantileverTip' property

>> E.CantilerTip = Tip

You should now be able to run the force analysis. Note that the E-Mod of every force curve will be calculated.

The methods for statistical analysis have been carelessly slapped together for a specific use case and are subject to generalization in the near future.
How E-Moduli from just the fibril apex can be extracted from the FM class can be nearly identically copied from those methods (replacing 'obj.' with 'ExperimentNameInWorkspace.' e.g. 'E.')




QUICKSTART GUIDE Force Map Analysis
%%% by Manuel Rufin, program version 23.11.2020%%%



1.
How to clone and pull from GitHub:

If you wanna clone the program to a ilsb workstation, you first have to set the proxy for git:

$ git config --global http.proxy http://proxy.ilsb.tuwien.ac.at:8080

you can check if you have set the right proxy to the config with

$ cat .gitconfig

then, create a folder where the git repository can be cloned into.

$ mkdir MyFolder
$ cd MyFolder

you can now clone the git repository using the link 'https://github.com/Rufman91/ForceMapAnalysis.git':

$ git clone https://github.com/Rufman91/ForceMapAnalysis.git

the programm is now downloaded to your directory. If there is a new version available, you can pull it from the repository.
'cd ' into your folder and use:

$ git pull

if you have made changes to the programm, you need to 'stash' them and then pull.
Beware that the changes are overwritten, if you stash and pull!

$ git stash
$ git pull

For further information or if you want to contribute to the programm contact rufman@ilsb.tuwien.ac.at



2. Set up programm

To set up the program, add the force map analysis folder to the MATLAB search path under Home > Set Path.
You need the 7z.exe software (https://www.7-zip.org/download.html) and also check in which filedirectory 
the exe file is located (if you picked the default location, everything should work without intervention from your side).

The following MATLAB Toolboxes might be needed:
- Statistics and Machine Learning Toolbox
- Deep Learning Toolbox
- Signal Processing Toolbox
- Parallel Computing Toolbox
- Image Processing Toolbox
- Curve Fitting Toolbox

Note that the contact point determination using one of the Neural Networks will run significantly faster on a device with a good GPU and the Parallel Computing Toolbox!



3. Create and work with Experiment

To create an Experiment-file, which will be the handle class object in your MATLAB workspace you actually work with, call the constructor method

>> ExperimentNameInWorkspace = Experiment();

choose a very short name for ExperimentNameInWorkspace, as every function/data will be a method/property of this object (e.g. 'E','Ex').
The programm will lead you through the process of creating your experiment: 
what kind of measurements have been done? (Force Maps, KPFM) --> Number of files / name of experiment --
--> choose your *.jpk-force-map files (currently only works for the newest two generations of jpk force map files (23.09.2020))

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

>> E.force_map_analysis_fibril(CPOption, EModOption)      % for collagen fibril force maps
>> E.force_map_analysis_general(CPOption, EModOption)     % for more general cases of force maps (e.g. no fibril masking, no FibDiam calculation etc.)
>> E.surface_potential_analysis_fibril()

the options for the force map analysis are explained in the functions header comments.

The tip deconvolution was just adapted from the old script by Orestis Andriotis. It is run in the context of the force_map_analysis_fibril method right at the beginning.

You should now be able to run the force analysis. Note that the E-Mod of every force curve will be calculated.

The methods for statistical analysis have been carelessly slapped together for a specific use case and are subject to generalization in the near future.

The fibrils overall moduli, as well as the apex moduli and standard deviations are saved into the property struct E.EMod


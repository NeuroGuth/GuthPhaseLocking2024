# GuthPhaseLocking2024
Analysis code for Guth et al., bioRxiv, 2024: "Theta-phase locking of single neurons during human spatial memory".

# TreasureHunt
This folder contains the scripts for all analyses in the paper. Different folders contain the scripts for different analyses.
For example, the subfolder "PhaseAnalysis_20230921" contains the script for the analysis of neuronal phase locking.

# Demo
The live script "TG_Demo_20240729" provides an overview of how we assessed neuronal theta-phase locking. This script, along with example data, is available in the subfolder "Demo_20240729".

To run the script, follow these steps:

   1. Download the folder.
   2. Open the live script in MATLAB.
   3. Press "Run" in the MATLAB Live Editor.

The expected runtime is less than one minute. You can further reduce this time by lowering the value of the parameter "params.nSur".


# System requirements and installation guide
The code is written in MATLAB R2021a and R2023a (MathWorks, Natick, MA, USA) and is compatible with these versions. The scripts have been tested on Windows 10 and Windows 11 (Microsoft Corporation, Redmond, WA, USA).

# Dependencies
The scripts utilize the following external functions and toolboxes. Ensure these are downloaded and added to your MATLAB environment for the scripts to execute properly:
[Fieldtrip](https://github.com/fieldtrip/fieldtrip) (Tested version: 20210912), 
[WaveClus](https://github.com/csn-le/wave_clus) (Tested version: 3.0.3), 
[Generalized Phase](https://github.com/mullerlab/generalized-phase) (Tested version: 20201007), 
[CircStat](https://github.com/circstat/circstat-matlab) (Tested version: 1.21.0.0), 
[SPRiNT](https://github.com/lucwilson/SPRiNT) (Tested version: 1.0.3), 
[myBinomTest](https://de.mathworks.com/matlabcentral/fileexchange/24813-mybinomtest-s-n-p-sided) (Tested version: 2.0.0.0), 
[cline](https://de.mathworks.com/matlabcentral/fileexchange/14677-cline) (Tested version: 1.0.0.0), 
[rude](https://de.mathworks.com/matlabcentral/fileexchange/6436-rude-a-pedestrian-run-length-decoder-encoder) (Tested version: 1.0.0.0), and 
[chi2test](https://de.mathworks.com/matlabcentral/fileexchange/16177-chi2test?requestedDomain=) (Tested version: 1.0.0.0).

# Contact
tim.guth@ukbonn.de
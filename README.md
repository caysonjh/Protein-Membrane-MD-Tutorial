# Protein-Membrane-MD-Tutorial

## Introduction

### The Task

These instructions/tutorial will explain the process for running an atomistic simulation of a protein system embedded in a custom lipid bilayer membrane. We will first use the `memb_builder.py` python script to insert the protein into the membrane, specifying which lipids will compose the membrane and in what ratios. Then, we will prepare the protein-membrane system for molecular dynamics (MD). Each of the preparation steps and the simulation themselves will be run using the GROMACS simulation software package. After preparing the system for MD, we will first run an energy minimization step so as to ensure that the atoms are in the lowest-possible energy configuration, so that no infinite or unrealistic forces arise from overlapping atoms in the production MD run. Following the energy minimization, we will run an equilibration step that will stabilize the system parameters (temperature, pressure, volume) and ensure that the system is in a physiologically reasonable configuration. Finally, we can run the production MD run. It is from this run that we will get a trajectory file, which we can use to visualize what is happening in our simulation. The visualization of our simulation will be done using ChimeraX, the molecular visualization and analysis software.

### Tutorial Purpose

The purpose of this tutorial is to provide an example workflow for running a protein-membrane MD simulation, and not to be comprehensive for everything that can be accomplished using these software. Each specific use-case may require slight adaptations to the general pipeline. Likewise, each use-case may also encourage certain additions to the basic simulation, such as adding external forces or introducing applied constraints. These instructions will seek to explain the most stripped-down base case, though can be used as scaffolding for more complex cases. 

### Prerequisites/Disclaimers

This tutorial assumes you are using a unix-based operating system (Linux or MacOS). The instructions may work on Windows, but have not been tested. The instructions also assume that you are familiar with the command line interface (CLI) and basic bash commands. If you are not familiar with the command line, it is recommended to familiarize yourself with the basics of bash before proceeding.

The tutorial will embed the protein DotG from *Legionella pneumophila*, but the process can be completed with any protein in pdb ([protein databank](https://www.rcsb.org/)) file format. The lipids to be used to generate the membrane are POPC, POPE, and POPG. As with the protein, any combination of lipid files can be used, the lipid pdb files can be obtained from the [Charmm GUI Lipid Library](https://www.charmm-gui.org/?doc=archive&lib=lipid)
### Software Requirements

Installation of GROMACS and ChimeraX should precede beginning the instructions. Instructions for installation of these software can be found at the following links: [GROMACS Installation](https://manual.gromacs.org/2024.4/install-guide/index.html), [ChimeraX Installation (click "Other releases" for Linux/Windows)](https://www.cgl.ucsf.edu/chimerax/download.html). Another software that will be required is the python package manager anaconda, so [install anaconda](https://docs.anaconda.com/anaconda/install/) for your specific operating system if you don't already have it installed. 

To ensure that GROMACS is installed correctly, run the following command in the terminal:
```
gmx --version
```
If GROMACS is installed correctly, you should see a message that includes the version and specifics of your GROMACS build. 

Similarly, if conda is installed correctly you should be able to run the following: 
```
conda --version
```
If you see a message including your version of anaconda. 

ChimeraX should be installed as an application on your computer, and you should be able to open it from your applications folder.

---

## Instructions

### 1. Prepare the Environment

To begin, clone this repository to your local machine. You can do this by running the following command in your terminal: 
```
git clone https://github.com/caysonjh/Protein-Membrane-MD-Tutorial.git
```
Navigate to the directory where you cloned the repository. The `environment.yml` file is included in the repository, and includes the specs to create the anaconda environment, specifically defining the package [biopython](https://biopython.org/), which includes utility for manipulating pdb files. You can create the environment and activate it using the following set of commands. 
```
conda env create -f environment.yml
conda activate memb_builder
```
---
### 2. Insert the Protein into the Membrane with `memb_builder.py`

The first step will be to insert the protein into the membrane. We will do this using the provided python script `memb_builder.py`. This script is best run through the command line and includes multiple different command line options in order to customize the embedding. When complete, the entire command will look something like the following, and is the command we will use in this tutorial. 
```
python memb_builder.py --protein DotG.pdb --lipids POPC POPE POPG --lipid_ratios 3 5 2 --output output.pdb --box_size 300 --z -90 --buffer 1 --z-buffer 4
```
For specific details on the command line options, see [Command Line Options + Customization](#command-line-options--customization-md_builderpy) in the [Appendix](#appendix) at the bottom of the tutorial.

For the purpose of this tutorial, the parameters are defined as follows: 

- `--protein`: The name of our protein pdb file. In this case, we are using `DotG.pdb`, which is a 18 member polymer that functions as the opening to the Type IV Secretion System (T4SS) in *Legionella pneumophila*. 
- `--lipids`: The names of the lipids we want to use to compose the membrane. In this case, we are using a combination of `POPC`, `POPE`, and `POPG`, which are the lipids that compose the outer membrane of *Legionella pneumophila*.
- `--lipid_ratios`: The ratios of the lipids we want to use to compose the membrane. In this case, we are using a ratio of 3:5:2 for POPC:POPE:POPG.
- `--output`: The name of the output file. In this case, we are using `output.pdb`, but this can be any name you choose.
- `--box_size`: The size of the simulation box to be created. In this case, we are using a box size of 300 &Aring;ngstr&ouml;ms, because this is about twice the width of DotG and gives ample space on either side of the protein.
- `--z`: The z-coordinate where to begin inserting the membrane. In this case, we are using a z-coordinate of -90 &Aring;ngstr&ouml;ms, which aligns the membrane with the polar regions of DotG. This number was obtained through trial and error with various numbers until the membrane was inserted in the right place. 
- `--buffer`: The amount of buffer space to leave between each lipid during insertion. In this case, we are using a buffer of 1 &Aring;ngstr&ouml;ms, which avoids major atom clashes but gives us a realistic lipid density.
- `--z-buffer`: The amount of buffer space to leave between the upper and lower leaflet in the membrane. In this case, we are using a z-buffer of 4 &Aring;ngstr&ouml;ms, which is a good starting point for separating the two leaflets. This number was also obtained through trial and error with various numbers until the two leaflets were separated enough to avoid clashes but close enough so the two leaflets don't separate during simulation.

Once you've run the command and generated the output file, you should be able to view it in ChimeraX. You can either open ChimeraX and use the built in GUI to select the file, or just run `open output.pdb` in the terminal. Ensure that the system looks as desired before proceeding to the next step.

### 3. Prepare the System for MD








## Appendix

#### Command Line Options + Customization (md_builder.py)

 `--protein`: The name of the protein pdb file. 
  - This file should be in the same directory as the script. 
  - To use a different protein, simply specify the name of a valid pdb file in the same directory as the script. 
  - This option is also configured to accept multiple proteins, and can thus embed multiple proteins into a membrane together. 
  - The protein pdb files *must* be in your desired location (i.e. centered vs located at other coordinates) prior to running the script. If you are unsure of the configuration, run `gmx editconf -f your_protein.pdb -center 0`. The output will include a line that indicates the center of your protein. It will automatically center the protein and output a file called `out.pdb`, which you can delete if centering the protein is not needed for what you are doing. 
  - If you are including multiple proteins, or a single protein with a large hollow space, you will need to split the proteins into multiple files and provide them separately to the `--protein` option. The script is not sophisticated enough to recognize empty space within a protein's area and will not place lipids there. 

`--lipids`: The names of the lipids the membrane will be composed of.
  - The lipid files corresponding to the names of provided here should be in the same directory as the script (i.e. if `--lipid POPC` is specified then `POPC.pdb` should be in the directory)
  - Any number of lipids can be specified, and can be inserted in any ratio, specified by `--lipid_ratios`.
  - The script will automatically center and prepare the lipid file, so no additional preparation is needed here.

`--lipid_ratios`: The ratios of the lipids specified in `--lipids`.
  - The ratios should be specified in the same order as the lipids in `--lipids`.
  - The ratios should be integers, and there should be the same number specified as there are lipids in `--lipids`.
  - This will not affect the total number of lipids inserted into the membrane, but rather allows for membranes composed of multiple lipids.

`--output`: The name of the output file.

`--box_size`: The size of the simulation box to be created (in &Aring;ngstr&ouml;ms).
  - This primarily determines the size of the membrane, and should be large enough to encompass the entirety of the protein. It will also directly determine how many lipids are inserted.

`--z`: The z-coordinate where to begin inserting the membrane. 
  - This is the z-coordinate where lipid insertion will begin, specifically the lower leaflet.
  - This option can be adjusted to change where in relation to your protein the membrane is placed. The default will insert the membrane at the center, but often membranes should be inserted at the polar region of the protein. 
  - This is often best figured out by trial and error, as the script will not automatically determine the best z-coordinate for the membrane.

`--buffer`: The amount of buffer space to leave between each lipid during insertion (in &Aring;ngstr&ouml;ms). 
  - This will directly affect the density of the membrane. 
  - Adding a smaller buffer will increase membrane density, but will also increase the risk of atoms overlapping and causing infinite forces and other issues during simulation. 

`--z-buffer`: The amount of buffer space to leave between the upper and lower leaflet in the membrane (in &Aring;ngstr&ouml;ms). 
  - This will directly impact how close the two membrane leaflets are to each other.
  - The script will automatically separate the two leaflets by the average height of each of the lipids included, but this is often not enough to separate the two leaflets. 
  - Because we want there to be as little gap between the leaflets as possible, this is often best figured out by trial and error. 



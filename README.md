# GEOS-Chem-adjoint-v35
The GEOS-Chem adjoint model with the capability of multiple species chemical data assimilation

This code is based on the official GEOS-chem adjoint model V35l code (www.geos-chem.org). The major modifications comparing to the
standard code are:
  1. The capability of optimizing initial conditions and emissions simultaneously. 
  2. Addition of observation operator including IASI CO, IASI O3 (profile/column), MLS O3, OSIRIS O3, OSIRIS NO2 and OMI HCHO.
  3. Implementation of weak constraint 4D-var (already included in standard code after V35m)
  4. Updates on emission inventories (NEI11 and MEGAN 4.2)
  
The following is the instructions on how to compile the code on fujita-machine environment in Department of Physics, University of Toronto and on SCINET environment.

This document will provide instructions on downloading GEOS-Chem adjoint code, modifying configuration files (e.g. Makefile, objects.*, Object.* etc.) as well as setting up run files (e.g. input.gcadj and input.geos in directory ./gcadj\
_std.git/runs/v8-02-01/geos5). Please note that this is not a stand-alone document. For complete information of GEOS-Chem adjoint model, please read GEOS-Chem Online User's Guide available at: http://acmg.seas.harvard.edu/geos/doc/man/i\
ndex.html and GEOS-Chem Adjoint User's Guide available at: http://adjoint.colorado.edu/~yanko/gcadj_std/GC_adj_man.pdf.

What is GEOS-Chem and GEOS-Chem adjoint?
In one sentence, GEOS-Chem model is a global 3D Chemical Transport Model (CTM) driven by the assimilated meteorological observations from the Goddard Earth Observing System (GEOS). GEOS-Chem adjoint model refers to an adjoint model deri\
ved from the GEOS-Chem CTM, which focuses on inverse modelling studies (e.g. sensitivity tests, data assimilation). In terms of product properties, the GEOS-Chem code package includes the GEOS-Chem model (which is a forward model). The \
GEOS-Chem adjoint code package includes both the GEOS-Chem forward model and its derivative adjoint model. If you would like to further understand what is GEOS-Chem and GEOS-Chem adjoint, here are some sources that are helpful:

Mathematical Modeling of Atmospheric Chemistry by Brasseur G.P. and Jacob D.J. (draft in September, 2015)
This book provides theoretical backgrounds on CTMs in general. If you are interested in physical and chemical processes, the numerical methods applied, model evaluation as well as inverse modelling studies of a CTM, then this book would\
 be strongly recommended as your first reading. The electronic version of the draft of this book is available at: http://acmg.seas.harvard.edu/education/brasseur_jacob/index.html. Currently, this book is not officially published. The pa\
ssword to access each chapter is "ctm". In the future when the book is available, it is recommended that you obtain a hard copy.

Narrative description of the GEOS-Chem Model by Atmospheric Chemistry Modelling Group at Harvard University (last updated in May, 2015)
This page provides a quick reference for GEOS-Chem features and capabilities. Some details on meteorological fields, model resolution, emission inventories, chemistry scheme, parameterization scheme etc. are introduced. If you would lik\
e to understand the theorectical backgrounds of these details, some published references are also available at the bottom of the page. This description is online at: http://acmg.seas.harvard.edu/geos/geos_chem_narrative.html

GEOS-Chem Harvard website Wiki Main Page by Atmospheric Chemistry Modelling Group at Harvard University
This page provides a more detailed list on GEOS-Chem features and capabiities. You can also find some coding backgrounds (e.g. how to read each code, how to modify and debug the codes etc.) related to GEOS-Chem and analysis software (e.\
g. IDL, python). In fact, if you encounter any trouble when create/modify/compile codes within GEOS-Chem, you can always google "keywords" (e.g. geos_chem_adj_mod) + "GEOS-Chem", and the corresponding wiki page explaining the code shoul\
d appear on top of the google list. This wiki page is available at: http://wiki.seas.harvard.edu/geos-chem/index.php/Main_Page

GEOS-Chem Online User's Guide by Atmospheric Chemistry Modelling Group at Harvard University (last updated in June, 2015)
This guide provides an overview on installing, compiling and how to use the forward model code. Make sure you read the Online User's Guide before reading the following instructions.
This guide is available at: http://acmg.seas.harvard.edu/geos/doc/man/index.html

GEOS-Chem Adjoint User's Guide by Henze Group at Colorado Unversity, Boulder (last updated in June, 2015)
This guide provides an overview on installing, compiling and how to use the adjoint model code. Make sure you read the Adjoint User's Guide before reading the following instructions.
This guide is available at: http://adjoint.colorado.edu/~yanko/gcadj_std/GC_adj_man.pdf

Gitlab Tutorial by Henze Group at Colorado University, Boulder
This tutorial shows some features of Gitlab, which are helpful for version control of GEOS-Chem adjoint model. This is essential if you need to turn on/off particular features within the code, but still prefer the code to be stored in o\
ne directory (instead of creating multiple versions of GEOS-Chem adjoint codes in different directories). For example, let's suppose you have a piece of code in general, naming it version A. If you modify the code to enable function \al\
pha, naming it version B. After that, you modify the code to enable function \beta, naming it version C. Gitlab helps to track all the changes you made, and you can recall the code with the designated version as long as you save the ver\
sion using git. Otherwise, you have to store different versions of the code into different directories, or need to make those changes once again if you need to change from one version to another.
Full tutorial on how to use gitlab for GEOS-Chem is available at: http://adjoint.colorado.edu/~yanko/gcadj_std/GitLab_Tutorial.pdf

Guide on how to get started on scinet or fujitaXX platform
For both scinet and fujitaXX users, the first thing to do is to contact your supervisor to create a scinet/fujita account. After that,

For scinet users: Here are the online guides of SCINet user tutorial: https://wiki.scinet.utoronto.ca/wiki/images/5/54/SciNet_Tutorial.pdf
                  In particular, SCINet help wiki on GPC Quickstart : https://wiki.scinet.utoronto.ca/wiki/index.php/GPC_Quickstart

For fujitaXX users: You need to go to Physics Computer Services (PCS) to create your accounts under fujin server user list (including animus and fujitaXX). After that, you could ask/email Gregory Wu on whichever question you are interes\
ted in on how to use the machine, and he will give you some instructions. (No printing guides for fujitaXX so far.)

Download GEOS-Chem adjoint model (V35)

For SciNet users:


For local linux machine (e.g. fujitaXX) users:
To download the standard GEOS-Chem adjoint model, it is recommanded you use animus or any other linux machine that has git installed. SSH to your account, and then type the command:
git clone ssh://git@adjoint.colorado.edu:2222/yanko.davila/gcadj_std.git

This will download a project directory "gcadj_std" containing code and run subdirectories. To explain the file name, "gc" stands for GEOS-Chem. "adj" stands for adjoint code. "std" stands for standard code. "git" means that this code is\
 actively tracked by Gitlab. So all together, you end up getting a GEOS-Chem adjoint model standard code. Recall that for adjoint code you have both the forward GEOS-Chem model and the adjoint model. For its subdirectories, gcadj_std.gi\
t/code section includes all the programs that contribute to GEOS-Chem forward and adjoint model. gcadj_std.git/run section includes all the option/configuration files that allow you to choose how to run your model.

 you have already been provided with a GEOS-Chem adjoint model V35 code, you are strongly recommended to contact the original code owner (the person who wrote/modified the code or the person who passed the code to you). If the code ow\
ner uses Gitlab to track their modifications/updates on the code, please visit http://adjoint.colorado.edu:8080/ to check their branch and see what has been modified. If the code owner does not use Gitlab, then you must contact with the\
 code owner personally, and let the owner tell you what has been changed.

For the worst case scenario, if you cannot contact the code owner, and the owner does not use Gitlab, then the "primitive" way to track what the coder has done is to go to each directory (e.g. ./gcadj_std_modified/code/), and see which \
file has been modified. You can either look for files with "*.f/f90~", "#*.f/f90#" patterns or type emacs on a folder (NOT ON A FILE) (e.g. emacs ./gcadj_std_modified/code/obs_operators) to observe each file's properties (e.g. creating/\
modifying time etc.). After that, for all the modified codes, type command diff for each of them in order to see what has been changed. For example, you found there was mopitt_co_obs_mod.f~ under ./gcadj_std_modified/code/obs_operators.\
 Then, you should download the standard code, compare the standard code with the modified code by typing:
diff ./gcadj_std.git/code/obs_operators/mopitt_co_obs_mod.f ./gcadj_std_modified/code/obs_operators/mopitt_co_obs_mod.f

The output will give you what has been modified. Then, move to the next modified file and repeat the same...

As you can see, such "primitive" method to track a code takes times and effort. It is also casual, and easy to miss out key changes. If you are going to work with GEOS-Chem for some time, you should use Gitlab ASAP.

Configure GEOS-Chem adjoint model (V35)

Whether you have a standard code or a modified code, these gcadj_std codes do not require any installation. However, some other programs are required to install before running these codes. These programs are called library. (For instanc\
e, in your local computer, if you have all the pdfs downloaded, but without proper pdf reader to open it, then nothing can be read or modified. A proper pdf reader in this case is the library for reading pdfs.)

For SciNet users:
Since the users in SciNet are unable to install system files (including libraries) themselves, the technicians of Compute Canada has already installed some libraries that should be good enough to support GEOS-Chem adjoint model (V35). H\
owever, once you encounter any issue with the libraries due to version updates etc., you should email your supervisor ASAP and let them contact the techicians from Compute Canada.

For local linux machine (e.g. fujitaXX) users:
Depending on your project, we provide different routines to install different versions of libraries. For further details on installing library, please read GEOS-Chem library installation guide, which is available at XXX.

Now, libraries should be ready. And the GEOS-Chem adjoint code is ready as well. What to do next is to link the installed libraries to the GEOS-Chem adjoint code so that the gcadj_std code can be compiled successfully.

For SciNet users:
To set the environment properly, you should go to your code's run directory (e.g. /scratch/d/djones/USER_NAME/gcadj_std.git/runs/v8-02-01/geos5/) and open the run file. Add the following lines to the run script to load the appropriate l\
ibraries:

****************************************
module load intel/14.0.1
module load hdf5/1811-v18-serial-intel
module load netcdf/4.2.1.1_serial-intel
ROOT_LIBRARY_DIR=/scinet/gpc/Libraries/netcdf-4.2.1.1/serial_intel/
GC_BIN=$ROOT_LIBRARY_DIR/bin
GC_INCLUDE=$ROOT_LIBRARY_DIR/include
GC_LIB=$ROOT_LIBRARY_DIR/lib
export GC_LIB
export GC_BIN
export GC_INCLUDE
ulimit -s unlimited

For local linux machine (e.g. fujitaXX) users:
You first need to modify your .bashrc file, which should be available in diretory: /home/USER_NAME/ .bashrc file includes a series of configurations for the terminal session, including setting up your working directories, enabling colou\
ring, the shell history, library path selection, compiler selection etc. This file excecuted whenver a new terminal session starts in the interactive mode. In other words, let's say if you modify any commands in .bashrc file, the impact\
 will take place after you restart this terminal session (exit, and then ssh to login again) or create another terminal session.

To configure your .bashrc file, here is an example to look at:
To configure your .bashrc file, here is an example to look at:

******************************************************************
[source,perl]
# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

#-------Intel Compiler environment --------

if [ "$HOSTNAME" = "fujita01.atmosp.physics.utoronto.ca" ]; then
  source /opt/intel/Compiler/11.1/080/bin/ifortvars.sh intel64
  export OMP_NUM_THREADS=8
  export F_UFMTENDIAN=big
  export KMP_STACKSIZE=32m
  ulimit -s 32000000
fi

if [ "$HOSTNAME" = "fujita02.atmosp.physics.utoronto.ca" ]; then
  source /opt/intel/composer_xe_2013.0.079/bin/compilervars.sh intel64
  export OMP_NUM_THREADS=32
  export F_UFMTENDIAN=big
  export KMP_STACKSIZE=32m
  ulimit -s 32000000
fi

if [ "$HOSTNAME" = "fujita03.atmosp.physics.utoronto.ca" ]; then
  source /opt/intel/composer_xe_2013.0.079/bin/compilervars.sh intel64
  export OMP_NUM_THREADS=32
  export F_UFMTENDIAN=big
  export KMP_STACKSIZE=32m
  ulimit -s 32000000
fi
if [ "$HOSTNAME" = "fujita04.atmosp.physics.utoronto.ca" ]; then
  source /opt/intel/composer_xe_2013.0.079/bin/compilervars.sh intel64
  export OMP_NUM_THREADS=32
  export F_UFMTENDIAN=big
  export KMP_STACKSIZE=32m
  ulimit -s 32000000
fi

if [ "$HOSTNAME" = "fujita05.atmosp.physics.utoronto.ca" ]; then
  source /opt/intel/Compiler/11.1/080/bin/ifortvars.sh intel64
  export OMP_NUM_THREADS=8
  export F_UFMTENDIAN=big
  export KMP_STACKSIZE=32m
fi

if [ "$HOSTNAME" = "fujita06.atmosp.physics.utoronto.ca" ]; then
  source /opt/intel/composer_xe_2013.0.079/bin/compilervars.sh intel64
  export OMP_NUM_THREADS=16
  export F_UFMTENDIAN=big
  export KMP_STACKSIZE=32m
  ulimit -s 32000000
fi

if [ "$HOSTNAME" = "fujita07.atmosp.physics.utoronto.ca" ]; then
  source /opt/intel/composer_xe_2013.0.079/bin/compilervars.sh intel64
  export OMP_NUM_THREADS=16
  export F_UFMTENDIAN=big
  export KMP_STACKSIZE=32m
  ulimit -s 32000000
fi

# Python environment
export PATH=/usr/local/package/science/epd/epd/bin:$PATH

# IDL environment
IDL_STARTUP=/home/xzhang/IDL/idl_startup.pro
export IDL_STARTUP

# NetCDF environment (GEOS-Chem library)
export PATH=/home/xzhang/libraries_new/bin/:$PATH
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/xzhang/libraries_new/lib
export LD_LIBRARY_PATH

# Specify NetCDF environment subdirecotries
export GC_BIN=/home/xzhang/libraries_new/bin
export GC_LIB=/home/xzhang/libraries_new/lib
export GC_INCLUDE=/home/xzhang/libraries_new/include
**********************************************************************

In the .bashrc file, 3 components are included: global definitions of .bashrc file, compiler specification, and library specification.

For the compiler setup, fujita01 is used as system node. You should not run any GEOS-Chem code in this node. fujita02-04 are used for running GEOS-Chem codes for advanced users. fujita05-fujita07 are used for running GEOS-Chem codes for\
 primary users. All these nodes can be compiled using either Intel Fortran Compiler 11.1 or Intel Composer XE 2013. Make sure you contact your supervisor and Physics Computer Service (PCS) to activate one of these accounts. The followin\
g information shows the technical properties of each fujita node:

fujita01:
Intel Xeon E5410 2.33 GHz: 16GB RAM, 8 cores, no hyper-threading
fujita02-04 & 08:
Intel Xeon E5-2650v2, 2.6 GHz: 128GB RAM, 16 cores, 32 threads
fujita05:
Intel Xeon E5410 2.33 GHz: 8GB RAM, 8 cores, no hyper-threading
fujita06-07:
Intel Xeon E5620 2.4 GHz: 32GB RAM, 8 cores, and 16 threads

For the software setup, you need to specify NetCDF and related libraries to support GEOS-Chem, IDL library to support gamap (GEOS-Chem output analysis software), and python library to support python (coding software).

Second part of the configuration requires you to look at Makefile in directory ./gcadj_std.git/code/. Nowadays, no matter whether you use scinet or fujitaXX machines, you DO NOT need to modify this file as long as you follow the standar\
d procedure above. Your Makefile should be the same as the default standard code. If your Makefile is non-standard, please contact the original code owner for more information.

GEOS-Chem adjoint model code subdirectory:
Recall that code subdirectory is the section computing both forward model and adjoint model. In general, here are several files and folders in gcadj_std.git/code/ that might require modification based on your project:
Dependencies.mk        Dependencies listing for each file. If you are creating new file in ./gcadj_std.git/code/, then you will need to create a dependency statement here.
Makefile               Compilation environment setup file. If you are creating a new file (mostly for adjoint calculation), then you will need to create a compilation statement here.
Objects.default        All the default compiled object files (*.o) used in both forward and adjoint model. If you are creating new file for forward modelling calculation, then you will need to create an object statement here.
Objects.mk             All the compiled object files (*.o) used in both forward and adjoint model. Your modification should be similar as how Objects.default is modified.
Objects.mkl            All the compiled object files (*.o) used in both forward and adjoint model. Usually you don't need to modify this file

objects.sh             Compiling option file. You will notice that in this file, different modes (HDF, SAT_NETCDF, etc.) compile different set of code. You will modify this file after you create new adjoint-related functions (e.g. new o\
bservation operators) to be added in your run. For instance, if you create a new OMI HCHO observation operator, which has HDF-related functions in the code, then you need to add omi_hcho_obs_mod.o under HDF compiling options. When you t\
rigger "HDF" compiling options in your run file (under gcadj_std.git/runs/v8-02-01/geos5_test), the OMI HCHO observation operator will be compiled. Here are some brief explanations on several known modes:

DEFAULT: Other than the standard forward and adjoint codes, no additional functions (e.g. observation operator, radiative transfer model (RTM)) are compiled.
LIDORT: Other than the standard forward and adjoint codes, LIDORT (a RTM) is also compiled.
HDF: Other than the standard forward and adjoint codes, the HDF4/5-enabled codes are also compiled. (HDF4/5 refers to a dataset format applied to various satellite retrieval data (e.g. OMI, MLS).
SAT_NETCDF: Other than the standard forward and adjoint codes, the NetCDF-enabled codes are also compiled. (NetCDF refers to a dataset format applied to various satellite retrieval data (e.g. TES).
HDF_NETCDF: Other than the standard forward and adjoint codes, both HDF4/5-enabled and NetCDF-enabled codes are also compiled. It is a mandatory setup when you assimilate various datasets at the same time.

Note: If you use "off-diagonal" or "Compute BFGS inverse Hessian" modes, please contact your supervisor for further information.
define.h              Model resolution and compiler type selection file. For example, if you would like to compile a GEOS-Chem model with GEOS5 global meteorological field using linux fortran compiler (which is the compiler for both sci\
net and fujitaXX), then you need to uncomment
#define GRID4x5  'GRID4x5'
#define GRIDREDUCED 'GRIDREDUCED'
#define LINUX_IFORT 'LINUX_IFORT'
#define IN_CLOUD_OD  'IN_CLOUD_OD' (Note: this comment is for geos5 only)

If you need to shift from a global model to a regional model, please follow section X for instructions.

modified/             The folder including all forward model calculation. In particular, geos_chem_mod.f in /gcadj_std.git/code/modified/ is the main code that supports forward calculation. If you would like to create code that modify o\
r diagnose the forward model calculation, then you need to create file under this folder.
obs_operators/        The folder including all the observation operator files for satellite retrieval inverse modelling/data assimilation calculation. If you would like to create new observation operator, then you need to create file un\
der this folder.
adjoint/              The folder including all the adjoint model calculation, In particular, geos_chem_adj_mode.f and inverse_driver.f in /gcadj_std.git/code/modified/ are the main codes that support adjoint calculation. If you have cre\
ated a new observation operator, then you need to modify these two files. Please learn from a similar observation operator (and see how it is written in those files) as your example.
adjoint/define_adj.h  Observation operator selection file. Uncomment particular observation operator(s) if you would like to use it(them). Here is an example of using MOPITT CO observation operator:
(...)
!----------CO observations------
! pick any combination
! => MOPITT CO
!    => MOPITT v3
!    => MOPITT v4
!    => MOPITT v5
! => AIRS CO
! => SCIA Bremen CO
!#define  MOPITT_V3_CO_OBS 'MOPITT_V3_CO_OBS'
!#define  MOPITT_V4_CO_OBS 'MOPITT_V4_CO_OBS'
#define  MOPITT_V5_CO_OBS 'MOPITT_V5_CO_OBS'
!#define  AIRS_CO_OBS     'AIRS_CO_OBS'
!#define  SCIA_BRE_CO_OBS   'SCIA_BRE_CO_OBS'
(...)

So far, here are the most common codes that require modification before you execute your code. If you need to modify other codes, please use GEOS-Chem model user guide (mentioned in Section X) for further details.

Now, after you setup all the dependencies, the compiling environment, observation operators as well as the model resolution, you are ready to move on to the run subdirectories.

GEOS-Chem adjoint model run subdirectory
Run subdirectory is the setup section that determines how to run the code section. In this section, we specifically focus on GEOS-Chem adjoint model run directory used for GEOS5 meteorogical field. Here are several files and folders in \
gcadj_std.git/run/v8-02-01/geos5 that might require modification based on your project:

restart.geos5.4x5.2005070100             Restart file that pre-defines the initial conditions of various tracer concentrations. You must have this file with the correct resolution (in this case: 4x5) as well as the correct date (in this\
 case, your starting period of your model should be 20050701. Otherwise, let's say you need to run the model from 20100101 to 20110101, whereas you don't have the restart file on that day, then please ask your supervisor to provide you \
the restart file. (You should only generate your own restart file by running the forward model from 20050701 to 20100101, and save the restart file after you have contacted your supervisor.)

input.geos                               GEOS-Chem forward model input file. You can select the period of your runs, whether to save the restart file, emission inventory (including specify the inventory directory as well as diagnostic o\
ptions. You can also save the boundary conditions here in order to be used for regional studies. Several menu that you need to pay attention:

First of all, under most conditions, the following menu does not require any modification:
TRACER MENU
GAMAP MENU
BENCHMARK MENU
UNIX CMDS MENU

For following menu, use T(F) to turn on(off) particular option(s).
TRANSPORT MENU
CONVECTION MENU
EMISSIONS MENU
AEROSOL MENU
DEPOSITION MENU
CHEMISTRY MENU
For following menu, you need to pay a bit more attention...

SIMULATION MENU: When you set your start time, make sure you already have the restart file on that date. Somehow the model also prefers the start hour to be only 000000.

OUTPUT MENU: By entering "3", you set a tick to start or stop your time window. "3" also means that an output will be saved by the end of the time window. Hence, make sure you put "3" at the time period you set in Simulation menu. By en\
tering "0", the model would not save any output until a "3" shows up. In fact, any "3" or "0" outside of the time window is not going to be read by the model. Here is an example, if you would like to start at July 1st and end on July 10\
th, by setting your output menu as
Schedule output for JUL : 3000000003000000000000000000000
You will get outputs that average the data (e.g. O3 concentrations) from July 1st to July 10th. Instead, if you would like outputs for everyday from the same time period, then your output menu should like like this:
Schedule output for JUL : 3333333333000000000000000000000
If you would like two averages (one from 1st to 3rd, one from 3rd to 10th), then your output menu will be:
Schedule output for JUL : 3030000003000000000000000000000

DIAGNOSTIC MENU: all the basic dignostics (concentrations, mass fluxes, sources, depositions) except for adjoint diagnostics (e.g. scaling factor, adjoint forcing gradient, cost function, satellite diagnostics) will be saved in ctm.bpch
Under /gcadj_std.git/runs/v8-02-01/runs, there are also ctm.bpch.X being saved while the code is running. These are diagnostic files for each iteration. ctm.bpch, in particular, is the updated diagnostics based on the latest iteration. \
In general, it your run needs 30 iterations, then ctm.bpch.30 or ctm.bpch (only if your run has finished up to 30 iterations) will be your posteriori outputs. To get the a priori, you need to run forward model only in order to retrieve \
ctm.bpch.0. Go to run instructions on how to run forward model only. For diagnostic menu, L means how much vertical level you would like to save. "0" means no level (nothing) is saved. "1" means surface level is saved. "47" means all ve\
rtical level is saved. "For Tracers to print out for each diagnostics", "all" means all tracer involved (based on your tracer and chemistry menu) will be saved for the corresponding diagnostic(s). "N" refers to the tracer index that can\
 be referenced in the tracer menu. For example, if "4" appears, then CO will be saved in this diagnostic.

(These menu only need to be turned on if necessary. Check the adjoint model guide for further information.
PLANEFLIGHT MENU
ND48/49/50/51 MENU
PROD & LOSS MENU

ARCHIVED OH/ O3 P/L MENU: Usually you don't need to change anything here. For OH, you may notice that the model seems to read verion 5 OH despite the forward model has a version 9/10. This is because version 5 produces the most realisti\
c mean OH concentration so far. If you would like to do experiments involving varying mean OH profiles, these menu should be modified.

NESTED GRID MENU: You only turn them on when you do regional studies. A detailed instruction on how to run nested-grid GEOS-Chem is available at:
http://wiki.seas.harvard.edu/geos-chem/index.php/Setting_up_GEOS-Chem_nested_grid_simulations


input.gcadj                              GEOS-Chem adjoint model input file. You can select the inverse modelling/data assimilation algorithms, turn on/off chemistry or aerosols, the a priori, specify the output directory, whether to op\
timize emission/initial condition/reaction rate of a certain tracer, what to observe in adjoint forcing arrays, as well as several diagnostic options. You can also specify whether do a finite difference test here for sensitvity studies.

For following menu, use T(F) to turn on(off) particular functions,

ADJOINT SIMULATION MENU
FORWARD MODEL OPTIONS
ADJOINT MODEL OPTIONS

N.B. No matter what studies you do (even if you are doing a forward model run only), you have to turn LADJ on. Besides that, choose proper work you prefer for your inverse modelling/data assmilation studies.
     In some versions of this file, you may notice there are typos such as "Selecet" (Line 5) or "Invese" (Line 6). In fact, they will not have an impact on your model. In other words, if you change them to "Select" or "Inverse" respect\
ively, they should produce identical results as the typo version.

DIRECTORIES:
OptData/                                Folder saves the cost function, adjoint forcing gradient outputs, scaling factors outputs as well as satellite diagnostics.
adjtmp/                                 Folder saves the temporary file after the forward model run in order to be read by the adjoint model.
diagadj/                                Folder saves the adjoint diagnostics.

CONTROL VARIABLE MENU:
At first, you need to turn on/off initial conditions/emissions/strat prod/loss to select which tracer you would like to optimize. If you turn off any option, then the corresponding detailed options displayed underneath would not matter.\
 For example, if LADJ_EMS is F, then any settings under "FOR LADJ_EMS" will all be turned off.

For LICS, when you speficy the optimized tracer, the tracer index is according to the tracer menu in input.geos For example, if you would like to optimize O3 initial conditions with the default scaling factor being 1 (if you want half o\
f the O3 concentration, set SF to be 0.5. Usually SF not being "1" is assigned for pseudo-obs or OSSE studies), you should end up with
>------------------------------<
          FOR LICS             :
NSOPT: number of tracers opt   : 1
  => opt these tracers------>  : TRC# trc_name SF_DEFAULT REG_PARAM ERROR
Tracer  #1                     : 2    Ox       1          1         1

N.B. For now, we don't quite understand how REG_PARAM or ERROR is going to impact the assimilation results.

For LADJ_EMS, only the emission type shown underneath is available. You cannot create your own emission type since the emission will be unavailable in the model. Once you turn on an existing emission type, follow the similar instruction\
s as FOR LICS. For example, if you would like to turn on anthropogenic CO emissions, you should have:
%%% CONTROL VARIABLE MENU %%%
Initial conditions LICS        : F
... OR emissions   LADJ_EMS    : T
(...)
>------------------------------<
          FOR LADJ_EMS         :
NNEMS: ems groups implemented  : 33
(...)
Emission #28                   : 28   IDADJ_ECO_an   T    1          1         0.2       100     100
(...)

N.B. This section is directly related to "NOBS" as well as "NOBS_CSPEC". As you can see from the files in /gcadj_std.git/code/obs_operators/, some observation operators require functions STT_ADJ or CSPEC_AFTER_CHEM. The principle is tha\
t as long as STT_ADJ exists in the observation operator you need to use, you have to add the tracer to NOBS. Similarly, if your observation operator has CSPEC_AFTER_CHEM, you have to add the species to NOBS_CSPEC. For example, if I turn\
 out MOPITT CO as well as TES O3 observation operator, then your NOBS and NOBS_CSPEC will look like this:
>------------------------------<
NOBS: number of tracers to obs : 1
  => obs these tracers------>  : TRC# tracer_name
Tracer  #1                     : 2    Ox
>------------------------------<
NOBS_CSPEC: # of species to obs: 2
  => obs these species------>  : species_name
Species #1                     : CO
Species #2                     : O3

N.B. Another understanding of the principle is that STT_ADJ is mostly used for stratospheric model observation comparison, whereas CSPEC_AFTER_CHEM is mostly used for tropospheric model observation comparison. Hence, you will assgin tra\
cer under NOBS(NOBS_CSPEC) if the observation applies for stratosphere(troposphere). But keep in mind that this is not always true.

WEAK CONSTRAINT MENU
This is an optional menu (in standard code download from gitlab, you might not see this menu in this file) that applies weak constraint 4D-var algorithm unstead of the regular 4D-var. Please contact your supervisor if you need to turn o\
n this menu.

DIAGNOSTICS MENU
General mode is going to print all the debugging info on the screen, and give you diagnostics (adjoint forcing gradient etc.) for each time steps throughout the time window.
CO satellite diagnostics are only for users who turn on MOPITT CO observation operators.

For following menu, please read the GEOS-Chem adjount model user guide for more information.
OBSERVATION MENU
FINITE DIFFERENCE MENU





run                                      GEOS-Chem adjoint model executables (similar to the "*.exe" you double click to access a software in windows/mac). You can select the starting and finishing iterations, specify run directories, a\
s well as choose which mode you would like to run the codes (the mode refers to the settings in object.sh in code subdirectory).

Iteration number setup:
# Set the start (or current ) iteration number
X=1

# Set the stopping iteration number
XSTOP=4
# Give every run a unique name (default is $PBS_JOBID)
RNAME=gcadj_std.git

# Specify Type of Run "DEFAULT, HDF, SAT_NETCDF, LIDORT"
TYPE=SAT_NETCDF

Run forward model only, specify X=XSTOP=0.
Run only interation, specify X=XSTOP=1.
Run 15 full interations, specify X=1, XSTOP=15.
Run 15 interations in two steps, specify X=1, XSTOP=N (1<N<15, N is integer). Then specify X=N, XSTOP=15 in the same run subdirectories.
Output the diagnostic for interation N (make sure you already have interation N calculated before), specify X=N, XSTOP=N.
Make sure which sets of library you would like to use so that you know which type of runs you do. (More reference on objects.sh)

Run directories setup
(...)
# Directory in the package where the executable runs
DRUNDIR=runs/v8-02-01/geos5_test

# Directory in the package with the source code
DCODE=code

# Package directory name
DPACK=$RNAME

# Directory where run packages are unpacked and run
#DRUN=/scratch/d/djones/USER_NAME/  #for scinet
DRUN=/users/jk/09/USER_NAME/        #for fujitaXX

# Directory where run packages are stored and saved (if not saved locally)
DSAVE=

# Directory where run packages are backed up (if not saved locally)
DARCHIVE=
(...)
#  Set number of threads
export OMP_NUM_THREADS=8 #16 for scinet, 16 or 24 for fujita02-04, 8 for fujita05, 8 or 16 for fujita06-07.

Starting from here, you remain rest of the code as default. The GEOS-Chem adjoint model should be ready to start. If you would like to do multiple species data assimilation, and your data format includes more than one type (e.g. if I wo\
uld like to use MOPITT CO and TES O3 at the same time, and they are HDF5 and NetCDF format respectively), then you need to use HDF_NETCDF mode. This mode is not included in the standard code. To create this mode, you need to change obje\
ct.sh under /gcadj_std.git/code/ and the run file mentioned here. The detailed changes are shown as follows:

objects.sh

(...)

if [ $1 = "DEFAULT" ]; then
mv Objects.mkl Objects.mk
fi

#starting from here, you create mode HDF_NETCDF
if [ $1 = "HDF_NETCDF" ]; then

find="rpmares_mod.o"
replace="rpmares_mod.o                 \
gosat_co2_mod.o               \
tes_nh3_mod.o                 \
tes_o3_mod.o                  \
tes_o3_irk_mod.o"

sed -e "s/$find/$replace/g" Objects.mkl > output1
find="tes_ch4_mod.o"
replace="tes_ch4_mod.o                 \
scia_ch4_mod.o"

sed -e "s/$find/$replace/g" output1 > output2
#rm output1
#mv output Objects.mk

find="getifsun.o"
replace="getifsun.o                    \
gvchsq.o "
sed -e "s/$find/$replace/g" output2 > output3

find="input_mod.o"
replace="input_mod.o            \
He4IncludeModule.o            \
He4ErrorModule.o              \
He4GridModule.o               \
He4SwathModule.o              \
findinv.o                     \
airsv5_mod.o                  \
airs_co_obs_mod.o             \
HdfIncludeModule.o            \
HdfSdModule.o                 \
HdfVdModule.o                 \
omi_no2_obs_mod.o             \
interp.o                      \
gaussj.o                      \
mopitt_obs_mod.o"
sed -e "s/$find/$replace/g" output3 > output
rm output1
rm output2
rm output3
mv output Objects.mk
fi

run

(...)
#starting from here, you create HDF_NETCDF mode.
        ./objects.sh $TYPE

            rm INPUT_FOLDER

        if   [ $TYPE = 'DEFAULT' ]; then
            IFORT_OPT="$IFORT_OPT "
        elif [ $TYPE = 'HDF' ]; then
            IFORT_OPT="$IFORT_OPT HDF=yes"
        elif [ $TYPE = 'SAT_NETCDF' ]; then
            IFORT_OPT="$IFORT_OPT SAT_NETCDF=yes"
        elif [ $TYPE = 'HDF_NETCDF' ]; then
            IFORT_OPT="$IFORT_OPT SAT_NETCDF=yes HDF=yes"
        elif [ $TYPE = 'LIDORT' ]; then
            IFORT_OPT="$IFORT_OPT LIDORT=yes"
        fi



So far we introduce the most common codes that require modification before you execute your version of the GEOS-Chem adjoint model. Please use GEOS-Chem adjoint model user guide (mentioned in Section X) for further details.
Execute the model:

After everything has been put into the right place, now you can go ahead to run your model. For both scinet and fujitaXX, to run the code in interactive mode, simply type ./run to execute the model. Please note that in scinet, interacti\
ve mode can only run for an hour or so. This is not enough for most of your runs (unless you are doing a test run). It isn't a good idea to run the interactive mode in fujitaXX, either. Since if your model needs more than 24 hrs to run,\
 you have to turn on your local computer for this much time as well.

In fact, there are alternative ways to run the code so that your model will still be running even after you logout your ssh session or shut down your local computer. We name this scenario as the prevent disturbing mode. But before we in\
troduce this mode, it is still recommended when you test the code, choose a shorter time period with 3 iterations, and run the model in the interactive mode. If the compiling is successfull without any error messages shown throughout th\
e first 3 iterations, then you can run the model with the full time period using prevent disturbing mode.

To run the prevent disturbing mode,

For scinet:      type qsub -l nodes=1:ppn=8,walltime=1:00:00 run

This allocates a total of 1 hour for the run, using 1 node and 8 CPUs. You can only run GEOS-Chem on 1 node since we don’t use MPI. However, in the run script change the number of threads to:

export OMP_NUM_THREADS=16

so that you are running with 16 “processors” in the queue.

For fujitaXX:   You could either use nohup or GNU screen function to run the model in prevent disturbing mode.

For nohup       type nohup ./run to run the model. Now, you can logout your ssh session or shut down your local computer. If you would like to stop running the model at some point, ssh to your account, type top to launch the task manage\
r, then find the ID number (PID) for your job. After that, type Ctrl + Z to exit the tast manager. Finally, type kill -9 PID to stop.

N.B. please use top command everytime you login to fujitaXX so that you can check if other user is running any model (usually GEOS-Chem related runs will be shown with the task name "geos"). The idea is that for each fujita nodes, the c\
omputing environment will be optimal for 2 GEOS-Chem runs. More than 2 GEOS-Chem runs (whether you triggers them all, or both you and your collegues trigger them together) in the same node will slow down the computing environment for ev\
ery user. In this case, logout and try a different node. If all the nodes are busy, and you need to run your model right away, please contact your supervisor for assistance.

For screen GNU  type screen Now the system is in screen GNU interface. Simply type ./run to run the code, then type Ctrl + A + D to detach (shut down) the screen. Now, you can logout your ssh session or shut down your local computer. If\
 you would like to retach (revisit) your task at some point, ssh to your account, type screen -r  You will be able to see the screen back again.

As you can see, nohup is the simpliest way to guarantee the model will run in the remote computer regardless of the disturbance of the local computer. Screen GNU is more powerful since it saves the screen of the remote computer so that \
you can access it on a later time, no matter whether you logout your ssh or not. Feel free to use either of the two commands. But keep in mind, make sure you kill the job or end the screen session properly after you finish running your \
model. Otherwise, the model will undergo a "runaway" problem, which the model keeps running and occupying CPU usage as well as storage even after your run is complete.

Overall, this guide provides you a simplified approach to set up your GEOS-Chem adjoint model, and run the code using your local machine. However, you may need to make necessary changes to get the code to work. Some other linux programm\
ing skills are essential to play with the giant toy as well. Whenever you encounter a new problem, make sure you look them up on the internet or ask your collegues. These are the quickest way to solve your problem. Good luck.

Appendix: Common GEOS-Chem error messages

When you run the model for the first time, most likely your screen will stop and say "THE COMPILE DID NOT FINISH". This is a sign that the model has something wrong that you need to fix. Everytime when GEOS-Chem model broke down while r\
unning, it will also give you some error messages. 50% of the chance you can figure those out simply by reading the error messages. The other 50% require you to debug yourself. Here, we introduce some quick references/checklists that wi\
ll help you identify quickly what is wrong.

1. GEOS-Chem community has collected some common GEOS-Chem error messages so that users can understand what is wrong. The link is:
http://wiki.seas.harvard.edu/geos-chem/index.php/Common_GEOS-Chem_error_messages

2. Compilation error
When you run the code for the first time, make sure you go over every single file the manual mentioned above. There is hardly any quicker way to setup the model properly than following this guide (unless you start with some existing cod\
e that has worked well in the past). Moreover, by default, there should be no "*.o" or "*.mod" in your code directory before you your initial run (our after you switch/update your library). Make sure delete these module files if you hav\
e them. Once you make sure all above are all set, you execute the run file. The code will start to compile the code with the screen printing looking like this:
ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
       R U N   F O R   N _ C A L C _ S T O P  =  1

 run: Removing old files
  - checking for old core files
ls: core.*: No such file or directory
removed `ctm.bpch'
  - checking for *.chk.* file
ls: adjtmp/*.chk.*: No such file or directory
ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
ITER file updated
./objects.sh: line 8: [: too many arguments
./objects.sh: line 19: [: too many arguments
ifort -cpp -w -O3 -auto -noalign -convert big_endian -fp-model precise -vec-report0  -traceback -openmp -Dmultitask -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include  -c fjx_acet_mo\
d.f
ifort -cpp -w -O3 -auto -noalign -convert big_endian -fp-model precise -vec-report0  -traceback -openmp -Dmultitask -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include  -c charpak_mod\
.f
ifort -cpp -w -O3 -auto -noalign -convert big_endian -fp-model precise -vec-report0  -traceback -openmp -Dmultitask -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include  -c error_mod.f
ifort -cpp -w -O3 -auto -noalign -convert big_endian -fp-model precise -vec-report0  -traceback -openmp -Dmultitask -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include  -c -r8 ./new/n\
etcdf_util_mod.f
ifort -cpp -w -O3 -auto -noalign -convert big_endian -fp-model precise -vec-report0  -traceback -openmp -Dmultitask -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include  -c logical_mod\
.f
ifort -cpp -w -O3 -auto -noalign -convert big_endian -fp-model precise -vec-report0  -traceback -openmp -Dmultitask -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include  -c directory_m\
od.f
ifort -cpp -w -O3 -auto -noalign -convert big_endian -fp-model precise -vec-report0  -traceback -openmp -Dmultitask -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include  -c unix_cmds_m\
od.f
ifort -cpp -w -O3 -auto -noalign -convert big_endian -fp-model precise -vec-report0  -traceback -openmp -Dmultitask -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include  -c -r8 modifie\
d/tracer_mod.f
ifort -cpp -w -O3 -auto -noalign -convert big_endian -fp-model precise -vec-report0  -traceback -openmp -Dmultitask -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include  -c -r8 modifie\
d/julday_mod.f
ifort -cpp -w -O3 -auto -noalign -convert big_endian -fp-model precise -vec-report0  -traceback -openmp -Dmultitask -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include  -c file_mod.f
ifort -cpp -w -O3 -auto -noalign -convert big_endian -fp-model precise -vec-report0  -traceback -openmp -Dmultitask -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include  -c -r8 modifie\
d/grid_mod.f
ifort -cpp -w -O3 -auto -noalign -convert big_endian -fp-model precise -vec-report0  -traceback -openmp -Dmultitask -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include  -c -r8 modifie\
d/time_mod.f
...
ifort -cpp -w -O3 -auto -noalign -convert big_endian -fp-model precise -vec-report0  -traceback -openmp -Dmultitask -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include -I/home/ydavila/GC_Library/include   fjx_acet_mod.\
o charpak_mod.o error_mod.o netcdf_util_mod.o logical_mod.o directory_mod.o unix_cmds_mod.o tracer_mod.o julday_mod.o file_mod.o grid_mod.o time_mod.o logical_adj_mod.o directory_adj_mod.o bpch2_mod.o regrid_1x1_mod.o pressure_mod.o tra\
nsfer_mod.o future_emissions_mod.o lai_mod.o tracerid_mod.o benchmark_mod.o comode_mod.o diag_mod.o dao_mod.o tropopause_mod.o gckpp_adj_Precision.o gckpp_adj_Parameters.o adj_arrays_mod.o gckpp_adj_Function.o gckpp_adj_Global.o gckpp_a\
dj_Monitor.o gckpp_adj_Util.o gckpp_adj_HessianSP.o gckpp_adj_Hessian.o gckpp_adj_Initialize.o gckpp_adj_JacobianSP.o gckpp_adj_Jacobian.o gckpp_adj_LinearAlgebra.o gckpp_adj_Rates.o gckpp_adj_StoichiomSP.o gckpp_adj_Stoichiom.o gckpp_a\
dj_Integrator.o gckpp_adj_Model.o checkpoint_mod.o pbl_mix_mod.o pbl_mix_adj_mod.o diag03_mod.o diag04_mod.o diag41_mod.o diag42_mod.o diag48_mod.o diag49_mod.o diag50_mod.o diag51_mod.o diag56_mod.o diag59_mod.o diag_oh_mod.o diag_pl_m\
od.o ocean_mercury_mod.o drydep_mod.o scale_anthro_mod.o edgar_mod.o bravo_mod.o emep_mod.o nei2005_anthro_mod.o epa_nei_mod.o streets_anthro_mod.o icoads_ship_mod.o arctas_ship_emiss_mod.o cac_anthro_mod.o vistas_anthro_mod.o geia_mod.\
o global_oh_mod.o global_hno3_mod.o global_no3_mod.o global_nox_mod.o global_o1d_mod.o global_o3_mod.o hippo_mod.o uvalbedo_mod.o RnPbBe_mod.o Kr85_mod.o acetone_mod.o aerosol_mod.o aircraft_nox_mod.o retro_mod.o biofuel_mod.o gc_biomas\
s_mod.o gfed2_biomass_mod.o gfed3_biomass_mod.o biomass_mod.o global_ch4_mod.o global_ch4_adj_mod.o c2h6_mod.o ch3i_mod.o a3_read_mod.o a6_read_mod.o i6_read_mod.o gcap_read_mod.o gwet_read_mod.o xtra_read_mod.o megan_mod.o carbon_mod.o\
 carbon_adj_mod.o optdepth_mod.o planeflight_mod.o restart_mod.o checkpt_mod.o population_mod.o lightning_nox_mod.o rpmares_mod.o                 gosat_co2_mod.o               tes_nh3_mod.o                 tes_o3_mod.o                  \
tes_o3_irk_mod.o rpmares_adj_mod.o isoropiaIIcode_adj.o isoropiaII_adj_mod.o wetscav_mod.o wetscav_adj_mod.o seasalt_mod.o sulfate_mod.o sulfate_adj_mod.o hcn_ch3cn_mod.o tagged_co_mod.o tagged_co_adj_mod.o tagged_ox_mod.o tagged_ox_adj\
_mod.o h2_hd_mod.o gcap_convect_mod.o fvdas_convect_mod.o convection_mod.o fvdas_convect_adj_mod.o convection_adj_mod.o pjc_pfix_mod.o pjc_pfix_geos5_window_mod.o dust_dead_mod.o dust_mod.o dust_adj_mod.o co2_mod.o co2_adj_mod.o mercury\
_mod.o toms_mod.o tpcore_bc_mod.o tpcore_fvdas_mod.o tpcore_mod.o tpcore_window_mod.o tpcore_geos5_window_mod.o transport_mod.o linoz_mod.o linoz_adj_mod.o upbdflx_adj_mod.o upbdflx_mod.o strat_chem_mod.o strat_chem_adj_mod.o chemistry_\
mod.o chemistry_adj_mod.o emissions_mod.o emissions_adj_mod.o weak_constraint_mod.o gamap_mod.o input_mod.o            He4IncludeModule.o            He4ErrorModule.o              He4GridModule.o               He4SwathModule.o           \
   findinv.o                     airsv5_mod.o                  airs_co_obs_mod.o             HdfIncludeModule.o            HdfSdModule.o                 HdfVdModule.o                 mls_o3_obs_mod.o              mls_hno3_obs_mod.o     \
       omi_no2_obs_mod.o             omi_ch2o_obs_mod.o            interp.o                      gaussj.o                      mopitt_obs_mod.o improve_bc_mod.o geos_chem_mod.o ErrorModule.o sciabr_co_obs_mod.o tes_ch4_mod.o            \
     scia_ch4_mod.o mem_ch4_mod.o leo_ch4_mod.o geocape_ch4_mod.o osiris_obs_mod.o geos_chem_adj_mod.o inv_hessian_mod.o input_adj_mod.o inverse_mod.o inverse_driver.o CO_strat_pl.o CO_strat_pl_adj.o airmas.o anthroems.o arsl1k.o adBuff\
er.o adStack.o backsub.o biofit.o boxvl.o calcrate.o calcrate_adj.o chemdr.o chemdr_adj.o cleanup.o cleanup_adj.o decomp.o diag1.o diag3.o diag_2pm.o diagoh.o emf_scale.o emfossil.o emisop.o emisop_grass.o emisop_mb.o emissdr.o emmonot.\
o fertadd.o findmon.o fyrno3.o fyhoro.o gasconc.o get_global_ch4.o getifsun.o                    gvchsq.o initialize.o jsparse.o ksparse.o lump.o lump_adj.o ndxx_setup.o ohsave.o partition.o partition_adj.o pderiv.o physproc.o precipfra\
c.o pulsing.o rdisopt.o rdlai.o rdland.o rdlight.o rdmonot.o rdsoil.o readchem.o reader.o readlai.o routines.o ruralbox.o schem.o schem_adj.o setbase.o setemdep.o setemis.o setemis_adj.o setmodel.o sfcwindsqr.o smvgear.o soilbase.o soil\
crf.o soilnoxems.o soiltemp.o soiltype.o subfun.o sunparam.o tcorr.o tropopause.o update.o xltmmp.o                        ifort_errmsg.o BLKSLV.o CLDSRF.o EFOLD.o FLINT.o GAUSSP.o GEN.o JRATET.o JVALUE.o LEGND0.o MATIN4.o MIESCT.o NOAB\
S.o OPMIE.o RD_TJPL.o SPHERE.o XSEC1D.o XSECO2.o XSECO3.o fast_j.o fjfunc.o inphot.o jv_index.o mmran_16.o photoj.o rd_js.o rd_prof.o set_aer.o set_prof.o                     \
            -L/home/ydavila/GC_Library/lib -L/home/ydavila/GC_Library/lib  -L/home/ydavila/GC_Library/lib -lhdfeos -lGctp -lmfhdf -ldf -ldl -lrt -ljpeg -lsz -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lnetcdff -lnetcdf -lz -lm -L \
/libmkl_solver_lp64.a -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -openmp -lpthread -o geos
/opt/intel/composer_xe_2013.0.079/compiler/lib/intel64/libimf.so: warning: warning: feupdateenv is not implemented and will always fail
===============================================================================
G E O S - C H E M   U S E R   I N P U T

READ_INPUT_FILE: Reading input.geos

SIMULATION MENU
---------------
Start time of run           : 20091102 000000
End time of run             : 20091115 000000
...

If the model fails during this process, then you must encounter a compilation failure. A compilation failure usually indicates either error in the coding format or error in compilation. For error in the code format, the line of the code\
 which has problem will be displayed. This is usually easy to correct. For error in compilation, either you have module files in your code directories not been cleaned up or your library is not compatible with the code. In particular, i\
f the compilation stops when dealing with NetCDF-enabled code/command (e.g./new/netcdf_util_mod.f) or HDF-enabled code/command (, then you should consider change your library. Ask your supervisor for further information.

3. Segmentation fault

There are different error causing the segmentation fault. The most common ones are the data retrieving error or index out of shape error, which indicates either the dataset is corrupted/unrecognized/missing or the loop index is running \
out of the dimenisions of your arrays. There are also segmentation fault caused by compilation error, you need to check your code format and library in this case.

4. Error occured during the running period

Most of the error should occur within the first iteration. However, there are error occured when the model was running after several iterations. These error are mostly due to computer setups or model setups. Keep in mind that idealy, yo\
ur computing condition and your code structure should not be changed when the model (the code section) is actively running. In other words, if your computer server is down/reset or your code is modified while your model is running, then\
 most likely you will have this error. Error may also occur during the running period if your computing server is running out of memory or storage. Moreover, if the cost function is converging to a certain threshold, the model will also\
 stop immediately.

5. The model is running very slow

As you run the model several times, you would have a sense on how fast the model runs. If the model stucks at some point or runs very slowly (under the condition that you have never modified anything before this run), then this problem \
should be relaetd to the computing server being busy. Keep in mind that for fujitaXX machines in particular, usually the server will nomrally (run in normal speed) support two GEOS-Chem runs. If there are multiple users submit 3/3+ runs\
, the server will run very slowly. The server will also fail to react if you have ran the code several times a day, which is a common situation when you test/debug the code in the interactive mode. In all the conditions mentioned above,\
 stop submitting more jobs into the server, and run them at a later time. If you have to submit a large number of runs within a week, make sure you tell your group members/supervisor about this.

6. The model runs without error, but the result seems to be way off.

This type of error can be caused by bugs in the code (e.g. calculation of your observation operator is incorrect),improper run setups (e.g. CO optimization disabled when you need to optimize CO), or bugs in your analyzing tool (e.g. bug\
s in your python codes when plotting/importing the model output). Make sure you go over what you have created/modified, and check if there is any error inside. You can also write some diagnostics to help you identify the issues.

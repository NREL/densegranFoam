# densegranFoam

All of these files corresponds to the latest version of the densgranFoam, created at the National Renewable Energy Laboratory (NREL).  The files are based on the OpenFOAM software. Please see the included OpenFOAM readme file ("[README.OpenFOAM](https://github.com/NREL/densegranFoam/blob/main/README.OpenFOAM)") and the GPL licence information ("[COPYING](https://github.com/NREL/densegranFoam/blob/main/COPYING)"). Access to and use of densgranFoam imposes obligations on the user, as set forth in the DATA USE DISCLAIMER AGREEMENT ("[AGREEMENT](https://github.com/NREL/densegranFoam/blob/main/AGREEMENT)").

Required OpenFOAM Version:  [5.x](https://github.com/OpenFOAM/OpenFOAM-5.x)

Solvers and Codes Included:

A.  Solvers
    densGranEulerFoam - This solver is used for simulating dense granular flow based on multiPhaseEulerFoam. The solver at current version uses mu(I)-rheology model for transport parameters. In later versions other models may be incorporated. 
    
    granularEulerFoam - This solver is used for simulating dense granular flow based on twoPhaseEulerFoam. The solver at current version uses mu(I)-rheology model for transport parameters. In later versions other models may be incorporated. 

B.  Tutorials
    densGranEulerFoam - conicalHopper, damBreak4phase, inclinedPlane
    granularEulerFoam - biomassInclinedFlow, dataProcess

Installation/Compiling:
OpenFOAM:
   The included codes work only with the OpenFOAM CFD Toolbox. The densegranFoam repository does not include OpenFOAM. Please visit [OF-5.x](https://github.com/OpenFOAM/OpenFOAM-5.x) to download and install OpenFOAM.  


densegranFoam compilation:
Once OpenFOAM is installed and densegranFoam downloaded, please follow these steps:
   1. Copy the densegranFoam solver codes from applications to WM_PROJECT_USER_DIR.
   2. Run ./Allwclean
   3. Run ./Allwmake
   4. Optionally run ./Allwclean again to clean up any extraneous files.
   5. Make sure that no error messages appeared.  

Running Tutorials:
   1. Copy the tutorial case to the anticipated run location.
   2. Run ./Allwclean to clean up case and copy necessay initial/boundary conditions for volume fraction fields.
   3. Run ./Allrun to run mesh utility, set initial volume fraction fields, necessary decomposition for paraller run, and densegranFoam solver.

The totorial is created to familiarize the user with the general case file structure, parameters setting, and numerical schemes for densegranFoam solver.

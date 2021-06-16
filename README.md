# densegranFoam

All of these files corresponds to the latest version of the densgranFoam, created at the National Renewable Energy Laboratory (NREL).  The files are based on the OpenFOAM software. Please see the included OpenFOAM readme file ("README.OpenFOAM") and the GPL licence information ("COPYING"). Access to and use of densgranFoam imposes obligations on the user, as set forth in the DATA USE DISCLAIMER AGREEMENT ("AGREEMENT").

Required OpenFOAM Version:  [5.x](https://github.com/OpenFOAM/OpenFOAM-5.x)

Solvers and Codes Included:
A.  Solvers
    densegranFoam - A  solver, primarily for performing large-eddy simulation    (LES), for computing atmospheric boundary layer turbulent flow with
        the ability to specify surface roughness, stability, and wind speed        and direction. It must be used with hexahedral meshes (like those
        created by OpenFOAM's blockMesh utility) in order for it to perform       the planar averaging it does.


B.  Tutorials
    setFieldsABL - A utility to initialize the flow field for performing        atmospheric boundary layer LES.  With the utility, you can specify
        an initial mean profile and perturbations that accelerate the        development of turbulence will be superimposed.  You may also 
        specify the initial temperature profile and location and strength        of the capping inversion.



Installation/Compiling:
OpenFOAM:
   The included codes work only with the OpenFOAM CFD Toolbox.  The SOWFA
   repository does not include OpenFOAM.  Please visit 
   github.com/OpenFOAM/OpenFOAM-2.4.x and github.com/OpenFOAM/ThirdParty-2.4.x
   to download and install OpenFOAM.  This release of SOWFA is meant to work
   with up to OpenFOAM-2.4.x.  Making SOWFA work with versions 3.0 and higher
   may come in the near future.  


SOWFA:
   Once OpenFOAM and OpenFAST (optional) are installed and SOWFA downloaded, 
   please follow these steps:
   1.  Put the SOWFA directory where you want it using git clone.
   2   The SOWFA repository includes a .bash_profile example file with a
       function to load the environment variables for your version of OpenFOAM
       and SOWFA.  Please adapt this and use it as necessary and get the 
       environment variables properly loaded.  One point is that the variable
       WM_PROJECT_USER_DIR will be where SOWFA gets compiled.  For example,
       WM_PROJECT_USER_DIR might get set to SOWFA-2.4.x.  Then you will have
       your original SOWFA clean repository and a compilation repository
       specific to the version of OpenFOAM you are using. (This may be helpful
       when we make SOWFA compatible with 5.0.x or v1712 within the near
       future).  Or, you can just set WM_PROJECT_USER_DIR to SOWFA so that
       compilation happens in the SOWFA repository that you downloaded from
       github.
   3.  Run ./Allwclean
   4.  Run ./Allwmake
   5.  Optionally run ./Allwclean again to clean up any extraneous files.
   6.  Make sure that no error messages appeared and that all libraries and
       applications are listed as "up to date."  



Running Tutorials:
Tutorial example cases are provide for each solver. The tutorials are
as follows
1.  example.ABL.flatTerrain.neutral
    example.ABL.flatTerrain.unstable
    example.ABL.flatTerrain.stable
    -Example cases for computing flat-terrain laterally periodic precursor
     atmospheric large-eddy simulations under the full range of stability.
    -Uses ABLSolver

2.  example.ALM
    example.ALMAdvanced
    example.ADM
    -Example cases that use the actuator turbine models using the 
     windPlantSolver.<X> solvers.  These cases assume that a precursor
     case has been run to generate inflow velocity and temperature 
     data.  Because the inflow data is large, we did not include it here.
    -Uses windPlantSolver.ALM
          windPlantSolver.ADM
          windPlantSolver.ALMAdvanced

3.  example.NREL5MW.ALMAdvancedFAST8
    -Example case that uses the FAST 8 coupled actuator line turbine model.
     The example is a single rotor in uniform, non-turbulent inflow.
    -Uses pisoFoamTurbine.ALMAdvancedFAST8 
 
4.  example.UAE_PhaseVI.ALMAdvanced/
    -A case with the setup of the NREL Unsteady Aerodynamics Experiment
     Phase VI, in the NASA Ames 80'x120' wind tunnel.
    -Uses pisoFoamTurbine.ALMAdvanced

5.  example.mesoscaleInfluence.SWiFTsiteLubbock.11Nov2013Diurnal
    -An example of using mesoscale influence to drive the atmospheric
     boundary layer LES.  This case is driven by WRF model output
     for the DOE/Sandia Scaled Wind Farm Technology (SWiFT) site in 
     Lubbock, Texas.  The terrain is flat with laterally periodic
     boundaries in this example.  The case has enough WRF output
     contained in the "forcing" directory to simulate two diurnal
     cycles starting November 11, 2013 at 00:00:00 UTC.   

To run a tutorial, change to that tutorial directory and run the
"runscript.preprocess" script to set up the mesh, etc.  Then run the
"runscript.solve.1" script to run the solver.

These are basic tutorials meant to familiarize the user with the
general file structure of a case and the various input files.  They
can be run on a small amount of processors with coarse meshes, but if 
that is done, they will generate poor results. 

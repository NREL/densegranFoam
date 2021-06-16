"""
script to create blockMeshDict and setFieldsDict for a hopper (2D)
"""

## Define dimensions of the geometry and then write out what is needed to
## blockMeshDict, setFieldsDict, and hRef.


#### TODO ####

## it would be straightforward to make a module of functions (possibly
## object-oriented) to write out the various parts (vertices, blocks, patches,
## etc.)

## Also write new or adjust initial files (in ../0/) ?


import numpy as np

# set dimensions of reservoir (rectangular section above wedge section)
h1 = 0.6 # m
w1 = 0.494 # m

# set dimensions of wedge
semiAngle = 24.4
thetaDeg = 90-semiAngle # degrees; angle of ramp
theta = thetaDeg*np.pi/180
wo = 0.05
print(wo)
#wo = 0.2 # m; width of orifice

# number of gridpoints
Nx = 30 # in x
Ny = 60 # in y

# set the zone where there is material
matlow = 0.0 # m; lower level -- 0 if level with orifice, positive values if
             # above
#matTotLength = 0.68
#matup = matTotLength-wedgeDepth
#matup = 0.68 # m; upper level

# calculate other dimensions
wb = 0.5*(w1-wo) # base of ramps in the wedge
h2 = wb*np.tan(theta) # height of wedge
matup = h2+h1/2
#matup = matTotLength-h2
# arbitrary z dimension
z0 = 0.1*h1

# number of y gridpoints per zone
htot = h1 + h2
NY = np.empty(2, dtype=int)
NY[0] = int(Ny*h1/htot)
NY[1] = int(Ny*h2/htot)

# define the points of the geometry
npoints = 12
points = np.zeros((npoints,3))
points[0,:2] = [wb, 0] # z value already 0
points[1,:2] = [0, h2]
points[2,:2] = [w1, h2]
points[3,:2] = [wb+wo, 0]
points[4:8] = points[:4]
points[4:8,2] = z0
points[8,:2] = [0, htot]
points[9,:2] = [w1, htot]
points[10:] = points[8:10]
points[10:,2] = z0

# define the blocks
# (using a list of arrays in anticipation of different-length arrays for
# various geometries)
blocks = []
blocks.append( np.array([0, 3, 2, 1, 4, 7, 6, 5], dtype=int) )
blocks.append( np.array([1, 2, 9, 8, 5, 6, 11, 10], dtype=int) )

# define the patches
top = []
top.append( np.array([8, 9, 11, 10]) )
walls = []
walls.append( np.array([8, 1, 5, 10]) )
walls.append( np.array([9, 2, 6, 11]) )
walls.append( np.array([1, 0, 4, 5]) )
walls.append( np.array([2, 3, 7, 6]) )
orifice = []
orifice.append( np.array([0, 3, 7, 4]) )


# function to form a string of the numbers in a 1-d array (or list) separated
# only by space
def strarr(arr):
    return " ".join(str(i) for i in arr)


# now construct blockMeshDict
import os, sys
path = os.path.dirname(os.path.abspath(sys.argv[0])) + "/"
bmd = open(path + "blockMeshDict", 'w')
bmd.write("""\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5                                     |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
""")
for i in range(npoints):
    bmd.write(" (%g %g %g)\n"%tuple(points[i]))
bmd.write("""\
);

blocks
(
""")
for i in range(len(blocks)):
    bmd.write(" hex (" + strarr(blocks[i])
              +") (%i %i 1) simpleGrading (1 1 1)\n"%(Nx, NY[i]))
bmd.write("""\
);

edges
(
);

patches
(
 wall walls
 (
""")
for i in range(len(walls)):
    bmd.write("  (" + strarr(walls[i]) + ")\n")
bmd.write("""\
 )
 atmosphere top
 (
""")
for i in range(len(top)):
    bmd.write("  (" + strarr(top[i]) + ")\n")
bmd.write("""\
 )
 atmosphere orifice
 (
""")
for i in range(len(top)):
    bmd.write("  (" + strarr(orifice[i]) + ")\n")
bmd.write("""\
 )
);

mergePatchPairs
(
);

// ************************************************************************* //
""")
bmd.close()


# write setFieldsDict
sfd = open(path + "/setFieldsDict", 'w')
sfd.write("""\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5                                     |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      setFieldsDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defaultFieldValues
(
 volScalarFieldValue alpha.air 1
 volScalarFieldValue alpha.particles 0
);

regions
(
 boxToCell
 {
""")
sfd.write("   box (%g %g %g)"%(0, matlow, 0))
sfd.write(" (%g %g %g);\n"%(w1, matup, z0))
sfd.write("""\
   fieldValues
     (
      volScalarFieldValue alpha.air 0.55
      volScalarFieldValue alpha.particles 0.45
     );
 }
);


// ************************************************************************* //
""")
sfd.close()


# write href
href = open(path + "/../constant/hRef", 'w')
href.write("""\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5                                     |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       uniformDimensionedScalarField;
    location    "constant";
    object      hRef;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 0 0 0 0 0];
""")
href.write("value           %g;\n"%matup)
href.write("""\
// ************************************************************************* //
""")

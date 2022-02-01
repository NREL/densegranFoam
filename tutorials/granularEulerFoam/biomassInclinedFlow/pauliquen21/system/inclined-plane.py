"""
script to create blockMeshDict and setFieldsDict for an inclined plane (2D)
"""

## Define dimensions of the geometry and then write out what is needed to
## blockMeshDict and setFieldsDict files. Also write new or adjust initial files
## (in ../0/) ?

## it would be straightforward to make a module of functions (and possible
## object-oriented) to write out the various parts (vertices, blocks, patches,
## etc.)

import numpy as np

## set dimensions of ramp
theta = 23 # degrees; angle of ramp from horizontal
L = 1 # m; length of ramp
depth = 0.05 # m; trivial value for 2D
ha = 0.5 # m; height above left edge of ramp 

# number of gridpoints
Nx = 1000 # in x
Ny = 500 # in y

# set the zone where there is material
#matwidth = 0.3 # m, from left edge of ramp
matheight = 1.0 # m, above left edge of ramp; will fill to top of domain if set
                # large enough

# calculate other dimensions
b = L*np.cos(theta*(np.pi/180)) # m; base of ramp
h = L*np.sin(theta*(np.pi/180)) # m; height of ramp
h2 = h + ha # m; height of domain

# make material extend to lower edge of ramp
matwidth = b

bmat = -b + matwidth # m, right edge of material 
hmat = h + matheight # m, top edge of material (maybe outside domain?)

# number of x gridpoints per block
Nxb = Nx//2

# define the points of the geometry
npoints = 12
points = np.zeros((npoints,3))
points[1] = [-b, h, 0]
points[2] = [-b, h2, 0]
points[3] = [0, h2, 0]
points[4] = [b, h2, 0]
points[5] = [b, 0, 0]
points[6:] = points[:6]
points[6:,2] = depth

# define the blocks
# (using a list of arrays in anticipation of different-length arrays for
# different geometries)
blocks = []
blocks.append( np.array([1, 0, 3, 2, 7, 6, 9, 8], dtype=int) )
blocks.append( np.array([0, 5, 4, 3, 6, 11, 10, 9], dtype=int) )

# define the patches
left = np.array([2, 1, 7, 8])
ramp = np.array([1, 0, 6, 7])
top = []
top.append( np.array([2, 3, 9, 8]) )
top.append( np.array([3, 4, 10, 9]) )
right = np.array([4, 5, 11, 10])
bottom = np.array([0, 5, 11, 6])

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
              +") (%i %i 1) simpleGrading (1 1 1)\n"%(Nxb, Ny))
bmd.write("""\
);

edges
(
);

patches
(
 wall left
  (
""")
bmd.write("  (" + strarr(left) + ")\n")
bmd.write("""\
 )
 wall ramp
  (
""")
bmd.write("  (" + strarr(ramp) + ")\n")
bmd.write("""\
 )
 atmosphere top
 (
""")
for i in range(len(top)):
    bmd.write("  (" + strarr(top[i]) + ")\n")
bmd.write("""\
 )
 atmosphere right
 (
""")
bmd.write("  (" + strarr(right) + ")\n")
bmd.write("""\
 )
 atmosphere bottom
 (
""")
bmd.write("  (" + strarr(bottom) + ")\n")
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
 volScalarFieldValue alpha.air 0
 volScalarFieldValue alpha.particles 0
);

regions
(
 boxToCell
 {
""")
sfd.write("   box (%g %g %g)"%(-b, 0, 0))
sfd.write(" (%g %g %g);\n"%(bmat, hmat, depth))
sfd.write("""\
   fieldValues
     (
      volScalarFieldValue alpha.air 0
      volScalarFieldValue alpha.particles 0.55
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
href.write("value           %g;\n"%h2)
href.write("""\
// ************************************************************************* //
""")

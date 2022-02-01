/*--------------------------------*- C++ -*----------------------------------*\
| =========                |                                                 |
| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \    /   O peration     | Version:  5                                     |
|   \  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])


   convertToMeters 0.001;

   define(PI, 3.14159265)
   define(angD, calc(PI/180)) 
   define(mu, calc(10)) //semi angle 24.4
   define(Lw, 600)
   define(D1, 50)//openning 50 //5 mm column diameter
   define(D2, calc(2*(D1/2+Lw*sin(angD*mu))))
   define(L1, calc(Lw*cos(angD*mu))) //10 cm length
   define(L2, calc(L1+600))
   
   define(R1, calc(D1/2))
   define(CW1, calc(D1/4)) //Width of middle square section
   
   define(CX1, calc(R1*cos(angD*45)))
   define(CZ1, calc(R1*sin(angD*45)))

   define(R2, calc(D2/2))
   define(CW2, calc(D2/4)) //Width of middle square section

   define(CX2, calc(R2*cos(angD*45)))
   define(CZ2, calc(R2*sin(angD*45)))
   
   define(NPS, 4) //how many cells in the square section
   define(NPD, 4) //how many cells from square section to perimeter
   define(NPY, 20) // how many cells from top to bottom

   vertices
   (
    ( CW1 0.0  CW1) vlabel(fiveoclocksqb)
    (-CW1 0.0  CW1) vlabel(sevenoclocksqb)
    (-CW1 0.0 -CW1) vlabel(elevenoclocksqb)
    ( CW1 0.0 -CW1) vlabel(oneoclocksqb)
   
    ( CX1 0.0  CZ1) vlabel(fiveoclockcb)
    (-CX1 0.0  CZ1) vlabel(sevenoclockcb)
    (-CX1 0.0 -CZ1) vlabel(elevenoclockcb)
    ( CX1 0.0 -CZ1) vlabel(oneoclockcb)

    ( CW2 L1  CW2) vlabel(fiveoclocksqm)
    (-CW2 L1  CW2) vlabel(sevenoclocksqm)
    (-CW2 L1 -CW2) vlabel(elevenoclocksqm)
    ( CW2 L1 -CW2) vlabel(oneoclocksqm)
   
    ( CX2 L1  CZ2) vlabel(fiveoclockcm)
    (-CX2 L1  CZ2) vlabel(sevenoclockcm)
    (-CX2 L1 -CZ2) vlabel(elevenoclockcm)
    ( CX2 L1 -CZ2) vlabel(oneoclockcm)

    ( CW2 L2  CW2) vlabel(fiveoclocksqt)
    (-CW2 L2  CW2) vlabel(sevenoclocksqt)
    (-CW2 L2 -CW2) vlabel(elevenoclocksqt)
    ( CW2 L2 -CW2) vlabel(oneoclocksqt)

    ( CX2 L2  CZ2) vlabel(fiveoclockct)
    (-CX2 L2  CZ2) vlabel(sevenoclockct)
    (-CX2 L2 -CZ2) vlabel(elevenoclockct)
    ( CX2 L2 -CZ2) vlabel(oneoclockct)
   );				

   blocks
   (
    //square block
    hex (
       sevenoclocksqb fiveoclocksqb oneoclocksqb elevenoclocksqb
       sevenoclocksqm fiveoclocksqm oneoclocksqm elevenoclocksqm
       )
    (NPS NPS NPY)
    simpleGrading (1 1 1)

    //slice1
    hex (
       sevenoclockcb fiveoclockcb fiveoclocksqb sevenoclocksqb
       sevenoclockcm fiveoclockcm fiveoclocksqm sevenoclocksqm
       )
    (NPS NPD NPY)
    simpleGrading (1 1 1)

    //slice2
    hex (
       sevenoclocksqb elevenoclocksqb elevenoclockcb sevenoclockcb 
       sevenoclocksqm elevenoclocksqm elevenoclockcm sevenoclockcm 
       )
   (NPS NPD NPY)
simpleGrading (1 1 1)

   //slice3
   hex (
         elevenoclocksqb oneoclocksqb oneoclockcb elevenoclockcb
         elevenoclocksqm oneoclocksqm oneoclockcm elevenoclockcm
       )
   (NPS NPD NPY)
simpleGrading (1 1 1)

   //slice4
   hex (
         oneoclocksqb fiveoclocksqb fiveoclockcb oneoclockcb
         oneoclocksqm fiveoclocksqm fiveoclockcm oneoclockcm
       )
   (NPS NPD  NPY)
simpleGrading (1 1 1)

    //******************************//


        //square block
    hex (
       sevenoclocksqm fiveoclocksqm oneoclocksqm elevenoclocksqm
       sevenoclocksqt fiveoclocksqt oneoclocksqt elevenoclocksqt
       )
    (NPS NPS NPY)
    simpleGrading (1 1 1)

    //slice1
    hex (
       sevenoclockcm fiveoclockcm fiveoclocksqm sevenoclocksqm
       sevenoclockct fiveoclockct fiveoclocksqt sevenoclocksqt
       )
    (NPS NPD NPY)
    simpleGrading (1 1 1)

    //slice2
    hex (
       sevenoclocksqm elevenoclocksqm elevenoclockcm sevenoclockcm
       sevenoclocksqt elevenoclocksqt elevenoclockct sevenoclockct
       )
   (NPS NPD NPY)
simpleGrading (1 1 1)

       //slice3
   hex (
         elevenoclocksqm oneoclocksqm oneoclockcm elevenoclockcm
         elevenoclocksqt oneoclocksqt oneoclockct elevenoclockct
       )
   (NPS NPD NPY)
simpleGrading (1 1 1)

   //slice4
   hex (
         oneoclocksqm fiveoclocksqm fiveoclockcm oneoclockcm
         oneoclocksqt fiveoclocksqt fiveoclockct oneoclockct
       )
   (NPS NPD  NPY)
simpleGrading (1 1 1)

    

   );


   //create the quarter circles
   edges
   (
    arc fiveoclockcb sevenoclockcb (0.0 0.0 R1)
    arc sevenoclockcb elevenoclockcb (-R1 0.0 0.0)
    arc elevenoclockcb oneoclockcb (0.0 0.0 -R1)
    arc oneoclockcb fiveoclockcb (R1 0.0 0.0)

    arc fiveoclockcm sevenoclockcm (0.0 L1 R2)
    arc sevenoclockcm elevenoclockcm (-R2 L1 0.0)
    arc elevenoclockcm oneoclockcm (0.0 L1 -R2)
    arc oneoclockcm fiveoclockcm (R2 L1 0.0)

    arc fiveoclockct sevenoclockct (0.0 L2 R2)
    arc sevenoclockct elevenoclockct (-R2 L2 0.0)
    arc elevenoclockct oneoclockct (0.0 L2 -R2)
    arc oneoclockct fiveoclockct (R2 L2 0.0)

   );

   patches
   (
    patch outlet
    (
     (fiveoclocksqb oneoclocksqb elevenoclocksqb sevenoclocksqb)
     (fiveoclocksqb fiveoclockcb oneoclockcb oneoclocksqb)
     (fiveoclockcb fiveoclocksqb sevenoclocksqb sevenoclockcb)
     (sevenoclocksqb elevenoclocksqb elevenoclockcb sevenoclockcb)
     (oneoclocksqb oneoclockcb elevenoclockcb elevenoclocksqb)
    )

    patch inlet
    (
     (fiveoclocksqt oneoclocksqt elevenoclocksqt sevenoclocksqt)
     (fiveoclocksqt fiveoclockct oneoclockct oneoclocksqt)
     (fiveoclockct fiveoclocksqt sevenoclocksqt sevenoclockct)
     (sevenoclocksqt elevenoclocksqt elevenoclockct sevenoclockct)
     (oneoclocksqt oneoclockct elevenoclockct elevenoclocksqt)
    )

    wall walls
    (
     (sevenoclockcb fiveoclockcb fiveoclockcm sevenoclockcm)
     (sevenoclockcb sevenoclockcm elevenoclockcm elevenoclockcb)
     (elevenoclockcb elevenoclockcm oneoclockcm oneoclockcb)
     (oneoclockcb oneoclockcm fiveoclockcm fiveoclockcb)
     (sevenoclockcm fiveoclockcm fiveoclockct sevenoclockct)
     (sevenoclockcm sevenoclockct elevenoclockct elevenoclockcm)
     (elevenoclockcm elevenoclockct oneoclockct oneoclockcm)
     (oneoclockcm oneoclockct fiveoclockct fiveoclockcm)
    )

);

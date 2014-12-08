FoamFile
{
    version         2.0;
    format          ascii;
    class           dictionary;
    object          blockMeshDict;
}

//original file provided by Chalmers University

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

//Mathematical constants:
m4_define(pi, 3.1415926536)

//Geometry
m4_define(rmax, 1.0) // pipe radius

//Number of cells:
//Radial direction:
m4_define(rNumberOfCells, 5)  // in the O-grid belt
//Tangential direction:
m4_define(tNumberOfCells, 7)
//Axial direction:
m4_define(zIQnumberOfCells, 8)


//Radial grading:
m4_define(rGrading, 3)

//Axial grading:
m4_define(zGradingIQ, 3)

//Plane I:
m4_define(zI, 0.00)
m4_define(rI, rmax)
m4_define(rRelI, 0.724)

//Plane Q:
m4_define(zQ, 10.00)
m4_define(rQ, rmax)
m4_define(rRelQ, 0.724)


convertToMeters 1.0;

vertices
(
//Plane I:
(calc(rRelI*rI*cos(pi/4)) -calc(rRelI*rI*sin(pi/4)) zI) vlabel(I0)
(calc(rRelI*rI*cos(pi/4)) calc(rRelI*rI*sin(pi/4)) zI) vlabel(I1)
(calc(-rRelI*rI*cos(pi/4)) calc(rRelI*rI*sin(pi/4)) zI) vlabel(I2)
(calc(-rRelI*rI*cos(pi/4)) -calc(rRelI*rI*sin(pi/4)) zI) vlabel(I3)
(calc(rI*cos(pi/4)) -calc(rI*sin(pi/4)) zI) vlabel(I12)
(calc(rI*cos(pi/4)) calc(rI*sin(pi/4)) zI) vlabel(I13)
(calc(-rI*cos(pi/4)) calc(rI*sin(pi/4)) zI) vlabel(I14)
(calc(-rI*cos(pi/4)) -calc(rI*sin(pi/4)) zI) vlabel(I15)


//Plane Q:
(calc(rRelQ*rQ*cos(pi/4)) -calc(rRelQ*rQ*sin(pi/4)) zQ) vlabel(Q0)
(calc(rRelQ*rQ*cos(pi/4)) calc(rRelQ*rQ*sin(pi/4)) zQ) vlabel(Q1)
(calc(-rRelQ*rQ*cos(pi/4)) calc(rRelQ*rQ*sin(pi/4)) zQ) vlabel(Q2)
(calc(-rRelQ*rQ*cos(pi/4)) -calc(rRelQ*rQ*sin(pi/4)) zQ) vlabel(Q3)
(calc(rQ*cos(pi/4)) -calc(rQ*sin(pi/4)) zQ) vlabel(Q12)
(calc(rQ*cos(pi/4)) calc(rQ*sin(pi/4)) zQ) vlabel(Q13)
(calc(-rQ*cos(pi/4)) calc(rQ*sin(pi/4)) zQ) vlabel(Q14)
(calc(-rQ*cos(pi/4)) -calc(rQ*sin(pi/4)) zQ) vlabel(Q15)


);
// Defining blocks:
blocks
(
    //Blocks between plane I and plane Q:
    // block0 - positive x O-grid block belt
    hex (I13 I1 I0 I12 Q13 Q1 Q0 Q12 ) IQ
    (rNumberOfCells tNumberOfCells zIQnumberOfCells)
    simpleGrading (rGrading 1 zGradingIQ)
    // block1 - positive y O-grid block belt
    hex (I14 I2 I1 I13 Q14 Q2 Q1 Q13 ) IQ
    (rNumberOfCells tNumberOfCells zIQnumberOfCells)
    simpleGrading (rGrading 1 zGradingIQ)
    // block2 - negative x O-grid block belt
    hex (I15 I3 I2 I14 Q15 Q3 Q2 Q14 ) IQ
    (rNumberOfCells tNumberOfCells zIQnumberOfCells)
    simpleGrading (rGrading 1 zGradingIQ)
    // block3 - negative y O-grid block belt
    hex (I12 I0 I3 I15 Q12 Q0 Q3 Q15 ) IQ
    (rNumberOfCells tNumberOfCells zIQnumberOfCells)
    simpleGrading (rGrading 1 zGradingIQ)
    // block4 - central O-grid block 
    hex (I0 I1 I2 I3 Q0 Q1 Q2 Q3 ) IQ
    (tNumberOfCells tNumberOfCells zIQnumberOfCells)
    simpleGrading (1 1 zGradingIQ)
);

edges
(
    //Plane I:
    arc I0  I1  (calc(rRelI*rI) 0 zI)
    arc I1  I2  (0 calc(rRelI*rI) zI)
    arc I2  I3  (-calc(rRelI*rI) 0 zI)
    arc I3  I0  (0 -calc(rRelI*rI) zI)
    arc I12 I13 (rI 0 zI)
    arc I13 I14 (0 rI zI)
    arc I14 I15 (-rI 0 zI)
    arc I15 I12 (0 -rI zI)
    
    //Plane Q:
    arc Q0  Q1  (calc(rRelQ*rQ) 0 zQ)
    arc Q1  Q2  (0 calc(rRelQ*rQ) zQ)
    arc Q2  Q3  (-calc(rRelQ*rQ) 0 zQ)
    arc Q3  Q0  (0 -calc(rRelQ*rQ) zQ)
    arc Q12 Q13 (rQ 0 zQ)
    arc Q13 Q14 (0 rQ zQ)
    arc Q14 Q15 (-rQ 0 zQ)
    arc Q15 Q12 (0 -rQ zQ)

);

//patches done 
// Defining patches:
patches
(
    patch inlet
    (
       (I1 I13 I12 I0)
       (I2 I14 I13 I1)
       (I3 I15 I14 I2)
       (I0 I12 I15 I3)
       (I3 I2 I1 I0)
    )
    wall wall
    (
       (I12 I13 Q13 Q12)
       (I13 I14 Q14 Q13)
       (I14 I15 Q15 Q14)
       (I15 I12 Q12 Q15)
    )
    patch outlet
    (
       (Q0 Q12 Q13 Q1)
       (Q1 Q13 Q14 Q2)
       (Q2 Q14 Q15 Q3)
       (Q3 Q15 Q12 Q0)
       (Q0 Q1 Q2 Q3)
    )
);

mergePatchPairs 
(
);

// ************************************************************************* //

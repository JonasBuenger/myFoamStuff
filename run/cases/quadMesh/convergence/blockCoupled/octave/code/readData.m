clear;
clc;


 path2data = '/home/jonas/foam/jonas-3.0/run/cases/quadMesh/convergence/blockCoupled/';

 [cc1 cv1 s1 T1] = getData([path2data '1/']);
 [cc2 cv2 s2 T2] = getData([path2data '2/']);
 [cc3 cv3 s3 T3] = getData([path2data '3/']);
 [cc4 cv4 s4 T4] = getData([path2data '4/']);
 [cc5 cv5 s5 T5] = getData([path2data '5/']);
 [cc6 cv6 s6 T6] = getData([path2data '6/']);

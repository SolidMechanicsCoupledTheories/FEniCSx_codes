// Gmsh project created on Wed Jun 12 11:41:13 2024
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0.0, 0.0, 0.0, 2.5, 0, Pi/2, Pi/2};
//+
Characteristic Length{ PointsOf{ Volume{1}; } } = 0.25;
//+
Physical Surface("yBot", 8) = {3}; //yBot
//+
Physical Surface("zBot", 9) = {2}; //zBot
//+
Physical Surface("xBot", 10) = {4}; //xBot
//+
Physical Surface("Outer", 11) = {1}; //Outer
//+
Physical Volume("totVol", 12) = {1}; //totVol

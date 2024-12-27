// Gmsh project created on Wed Jun 12 10:01:30 2024
SetFactory("OpenCASCADE");
//+

// Sphere radius
r = 2.5;

Sphere(1) = {0, 0, 0, r, 0, Pi/2, Pi/2};
//+
Characteristic Length{ PointsOf{ Volume{1}; } } = 0.75;
//+
Curve Loop(5) = {4, -7, 6};
//+
Plane Surface(5) = {5};
//+
Curve Loop(7) = {5, -3, -6};
//+
Plane Surface(6) = {7};
//+
Curve Loop(8) = {5, 2, -7};
//+
Plane Surface(7) = {8};
//+
Curve Loop(9) = {3, 2, -4};
//+
Plane Surface(8) = {9};
//+
Physical Surface(12, "xBot") = {5};
//+
Physical Surface(13, "yBot") = {6};
//+
Physical Surface(14, "zBot") = {7};
//+
//Physical Surface(15, "outer") = {8};
//+
Physical Volume("totVold", 11) = {1};

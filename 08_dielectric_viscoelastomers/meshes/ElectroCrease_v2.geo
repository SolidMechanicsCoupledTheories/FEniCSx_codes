// Gmsh project created on Thu Jun 13 22:34:38 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {20, 0, 0, 1.0};
//+
Point(3) = {20, 120, 0, 1.0};
//+
Point(4) = {0.02, 120, 0, 1.0};
//+
Point(5) = {0, 119.88, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 1};
//+
Curve Loop(1) = {1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};//+
MeshSize {1, 2} = 5;
//+
MeshSize {3} = 1;
//+
MeshSize {4, 5} = 0.1;
//+
Physical Curve("Bottom", 6) = {1};
//+
Physical Curve("Right", 7) = {2};
//+
Physical Curve("Left", 8) = {5};
//+
Physical Curve("TopImperfection", 9) = {4};
//+
Physical Curve("Top", 10) = {3};
//+
Physical Surface("WholeBody", 11) = {1};

// Gmsh project created on Thu Jun 13 22:34:38 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {20, 0, 0, 1.0};
//+
Point(3) = {20, 20*6, 0, 1.0};
//+
Point(4) = {0, 20*6-1.e-3*20*6, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Physical Curve("Top", 5) = {3};
//+
Physical Curve("Bottom", 6) = {1};
//+
Physical Curve("Left", 7) = {4};
//+
Physical Curve("Right", 8) = {2};
//+
MeshSize {4} = 0.1;
//+
MeshSize {1, 2} = 1;
//+
Curve Loop(1) = {2, 3, 4, 1};
//+
Plane Surface(1) = {1};
//+
Physical Surface("All", 9) = {1};
//+
MeshSize {2, 1} = 2.0;
//+
MeshSize {1, 2} = 5.0;
//+
Split Curve {2} Point {};
//+
Split Curve {2} Point {};

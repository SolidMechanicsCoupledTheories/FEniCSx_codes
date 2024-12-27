// Gmsh project created on Tue Jun 18 22:37:51 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, -0.5, 0, 1.0};
//+
Point(2) = {0, -50, 0, 1.0};
//+
Point(3) = {50, -50, 0, 1.0};
//+
Point(4) = {50, 0, 0, 1.0};
//+
Point(5) = {50, 50, 0, 1.0};
//+
Point(6) = {0, 50, 0, 1.0};
//+
Point(7) = {0, 0.5, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Point(8) = {0, 0., 0, 1.0};
//+
Circle(7) = {7, 8, 1};
//+
Physical Curve("DropRadius", 8) = {7};
//+
Physical Curve("UpperLeft", 9) = {6};
//+
Physical Curve("LowerLeft", 10) = {1};
//+
Physical Curve("Top", 11) = {5};
//+
Physical Curve("Bottom", 12) = {2};
//+
Physical Curve("LowerRight", 13) = {3};
//+
Physical Curve("UpperRight", 14) = {4};
//+
Physical Point("RightGround", 15) = {4};
//+
MeshSize {1, 7} = 0.1;
//+
MeshSize {6, 5, 3, 2} = 5;
//+
Curve Loop(1) = {6, 7, 1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Surface("WholeBody", 16) = {1};
//+
MeshSize {7, 1} = 0.05;

// Gmsh project created on Wed Jun  5 13:27:27 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {50, 0, 0, 1.0};
//+
Point(3) = {50, 2.5, 0, 1.0};
//+
Point(4) = {50, 5, 0, 1.0};
//+
Point(5) = {0, 2.5, 0, 1.0};
//+
Point(6) = {0, 5, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 5};
//+
Line(4) = {5, 1};
//+
Line(5) = {3, 4};
//+
Line(6) = {4, 6};
//+
Line(7) = {6, 5};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, -7, -6, -5};
//+
Plane Surface(2) = {2};
//+
Physical Curve("Left_edge", 8) = {7, 4};
//+
Physical Curve("Bottom_edge", 9) = {1};
//+
Physical Curve("Right_edge", 10) = {2, 5};
//+
Physical Curve("Top_edge", 11) = {6};
//+
Physical Curve("Middle_edge", 12) = {3};
//+
Physical Surface("Bottom_layer", 13) = {1};
//+
Physical Surface("Top_layer", 14) = {2};

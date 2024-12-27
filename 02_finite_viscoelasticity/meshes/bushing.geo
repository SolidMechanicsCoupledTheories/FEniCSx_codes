// Gmsh project created on Fri Feb 23 15:13:20 2024
SetFactory("OpenCASCADE");
//+
l = 4;
Point(1) = {0, 0, 0, l};
//+
Point(2) = {0, 19.3, 0, l};
//+
Point(3) = {22.5, 19.3, 0, l};
//+
Point(4) = {22.5, 0, 0, l};
//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Point(5) = {29.3, 9.65, 0, l};
//+
Circle(4) = {3, 5, 4};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Extrude {{0, 1, 0}, {0, 0, 0}, 2*Pi} {
  Surface{1}; 
}
//+
Physical Surface("yTop", 7) = {2};
//+
Physical Surface("yBot", 8) = {4};
//+
Physical Volume(9) = {1};

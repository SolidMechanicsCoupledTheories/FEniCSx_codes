// Gmsh project created on Mon Jan 15 11:31:36 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0};
//+
Point(2) = {10, 0, 0};
//+
Point(3) = {10, 10, 0};
//+
Point(4) = {0, 10, 0};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 10} {
  Surface{1}; 
}
//+
Sphere(2) = {0, 0, 0, 5, 0, Pi/2, Pi/2};
//+
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; }
//+
Field[1] = Ball;
//+
Field[1].Radius = 5;
//+
Field[1].VIn = 1;
//+
Field[1].VOut = 5;
//+
Field[1].Thickness = 1;
//+
Background Field = 1;
//+
Physical Volume("inclusion", 32) = {2};
//+
Physical Volume("matrix", 33) = {1};

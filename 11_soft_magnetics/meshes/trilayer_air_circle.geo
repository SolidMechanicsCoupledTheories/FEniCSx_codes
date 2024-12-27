// Gmsh project created on Wed May 22 12:33:36 2024
SetFactory("OpenCASCADE"); 
//+
Point(1) = {-10, 0.25, 0, 1.0};
//+
Point(2) = {-10, -0.25, 0, 1.0};
//+
Point(3) = {-10, 0.50, 0, 1.0};
//+
Point(4) = {-10, -0.50, 0, 1.0};
//+
Point(5) = {10, -0.50, 0, 1.0};
//+
Point(6) = {10, -0.25, 0, 1.0};
//+
Point(7) = {10, 0.25, 0, 1.0};
//+
Point(8) = {10, 0.50, 0, 1.0};
//+
Point(9) = {10, 0, 0, 1.0};
//+
Point(10) = {-10, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 3};
//+
Line(8) = {3, 1};
//+
Line(9) = {1, 7};
//+
Line(10) = {2, 6};
//+
Circle(11) = {0, 0, 0, 50, 0, 2*Pi};
//+
Curve Loop(1) = {7, 8, 9, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, -5, -10, -1};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {3, 4, -10, 2};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {11};
//+
Curve Loop(5) = {7, 8, 1, 2, 3, 4, 5, 6};
//+
Plane Surface(4) = {4, 5};
//+
Physical Curve(12) = {11};
//+
Physical Curve(13) = {8, 2, 1};
//+
Physical Curve(14) = {9, 10, 5};
//+
Physical Surface(15) = {4};
//+
Physical Surface(16) = {1, 3};
//+
Physical Surface(17) = {2};
//+
MeshSize {1, 2, 3, 4, 5, 6, 7, 8} = .3;
Mesh.MeshSizeFromCurvature = 50;

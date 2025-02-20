// Gmsh project created on Wed Jul 10 16:41:32 2024
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 50, 50, 50};
//+
MeshSize {8, 7, 3, 4, 6, 2, 1, 5} = 10;
//+
Field[1] = Box;
//+
Field[1].Thickness = 50;
//+
Field[1].VIn  = 1;
Field[1].VOut = 10;
//+
Field[1].XMax = 10;
Field[1].XMin = 0;
//+
Field[1].YMax = 10;
Field[1].YMin = 0;
//+
Field[1].ZMax = 50;
Field[1].ZMin = 40;
//+
Background Field = 1;
//+
Physical Surface(13) = {1};
//+
Physical Surface(14) = {2};
//+
Physical Surface(15) = {3};
//+
Physical Surface(16) = {4};
//+
Physical Surface(17) = {5};
//+
Physical Surface(18) = {6};
//+
Physical Volume(19) = {1};

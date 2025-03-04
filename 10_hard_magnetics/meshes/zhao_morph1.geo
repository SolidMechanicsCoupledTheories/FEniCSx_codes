// Gmsh project created on Fri Oct 27 09:55:10 2023
//SetFactory("OpenCASCADE");
res=1.0;
//+
Point(1) = {21, 0, 0, res};
//+
Point(2) = {27, 0, 0, res};
//+
Point(3) = {27, 11, 0, res};
//+
Point(4) = {11, 11, 0, res};
//+
Point(5) = {11, 27, 0, res};
//+
Point(6) = {0, 27, 0, res};
//+
Point(7) = {0, 21, 0, res};
//+
Point(8) = {5, 21, 0, res};
//+
Point(9) = {5, 5, 0, res};
//+
Point(10) = {21, 5, 0, res};
//+
Point(11) = {16, 5, 0, res};
//+
Point(12) = {16, 11, 0, res};
//+
Point(13) = {5, 16, 0, res};
//+
Point(14) = {11, 16, 0, res};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 12};
//+
Line(4) = {12, 4};
//+
Line(5) = {4, 14};
//+
Line(6) = {14, 5};
//+
Line(7) = {5, 6};
//+
Line(8) = {6, 7};
//+
Line(9) = {7, 8};
//+
Line(10) = {8, 13};
//+
Line(11) = {13, 9};
//+
Line(12) = {9, 11};
//+
Line(13) = {11, 10};
//+
Line(14) = {10, 1};
//+
Line(15) = {10, 3};
//+
Line(16) = {11, 12};
//+
Line(17) = {9, 4};
//+
Line(18) = {13, 14};
//+
Line(19) = {8, 5};

//+
Curve Loop(1) = {1, 2, -15, 14};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {13, 15, 3, -16};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {12, 16, 4, -17};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {11, 17, 5, -18};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {10, 18, 6, -19};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {9, 19, 7, 8};
//+
Plane Surface(6) = {6};
//+
Extrude {0, 0, 0.84} {
  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Layers {2}; 
}
//+
Physical Surface("left_top", 152) = {150};
//+
Physical Volume("volume1", 153) = {1};
//+
Physical Volume("volume2", 154) = {2};
//+
Physical Volume("volume3", 155) = {3};
//+
Physical Volume("volume4", 156) = {4};
//+
Physical Volume("volume5", 157) = {5};
//+
Physical Volume("volume6", 158) = {6};
//+
Physical Surface("right_bot", 159) = {28};

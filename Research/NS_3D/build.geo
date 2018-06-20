// Gmsh project created on Tue Jun 19 14:55:08 2018
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, -0.125, 0.8, 0, 2*Pi};
//+
Circle(2) = {0, 0, -0.125, 0.3, 0, 2*Pi};
//+
Circle(3) = {0, 0, -0.125, 0.1, 0, 2*Pi};
//+
Extrude {{0, 1, 0}, {2, 0, 0}, Pi/4} {
  Curve{1,2,3}; Layers{24};
}//+
Curve Loop(4) = {1};
//+
Curve Loop(5) = {2};
//+
Curve Loop(6) = {2};
//+
Curve Loop(7) = {3};
//+
Curve Loop(8) = {3};
//+
Curve Loop(9) = {5};
//+
Curve Loop(10) = {7};
//+
Curve Loop(11) = {7};
//+
Curve Loop(12) = {9};
//+
Curve Loop(13) = {9};
//+
Plane Surface(4) = {4, 5};
//+
Plane Surface(5) = {6, 7};
//+
Plane Surface(6) = {8};
//+
Plane Surface(7) = {9, 10};
//+
Plane Surface(8) = {11, 12};


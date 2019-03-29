// Gmsh project created on Wed Mar 27 08:02:23 2019
//+
SetFactory("OpenCASCADE");
Cone(1) = {0, 0, 0, 0, -1, 0, 0.5, 0.1, 2*Pi};
//+
Characteristic Length {1} = 1;
//+
Characteristic Length {2} = 2;
//+
Transfinite Surface {1};
//+
Cylinder(2) = {0, 0, 0, 0, 1, 0, 0.5, 2*Pi};
//+
Characteristic Length {3, 4} = 2;

//+
Box(3) = {-0.5, -1.4, -0.5, 1, 0.3, 1};
//+
Characteristic Length {9, 10, 6, 5, 11, 7, 8, 12} = 1;

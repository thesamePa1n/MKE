Mesh.MshFileVersion = 2.0;

p = 0.1;

Point(1) = {0, 0, 0, p}; 
Point(2) = {2, 0, 0, p}; 
Point(3) = {2, 1, 0, p}; 
Point(4) = {0, 1, 0, p}; 

Point(5) = {0, 0, 1, p};
Point(6) = {2, 0, 1, p};
Point(7) = {2, 1, 1, p};
Point(8) = {0, 1, 1, p};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};//+
Line(5) = {5, 1};
//+
Line(6) = {8, 4};
//+
Line(7) = {7, 3};
//+
Line(8) = {6, 2};
//+
Line(9) = {6, 7};
//+
Line(10) = {7, 8};
//+
Line(11) = {8, 5};
//+
Line(12) = {5, 6};
//+
Curve Loop(1) = {6, 4, -5, -11};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {10, 11, 12, 9};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7, -2, -8, 9};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {2, 3, 4, 1};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {10, 6, -3, -7};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {12, 8, -1, -5};
//+
Plane Surface(6) = {6};
//+
Physical Surface(1) = {1};
//+
Physical Surface(2) = {3};
//+
Physical Surface(3) = {2};
//+
Physical Surface(4) = {4};
//+
Physical Surface(5) = {6};
//+
Physical Surface(6) = {5};
//+
Surface Loop(1) = {5, 2, 1, 4, 3, 6};
//+
Volume(1) = {1};
//+
Physical Volume(1) = {1};

SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 2.0;

LX = 2;
LY = 2;
r = 0.2;
p1 = LY / 60;
p2 = r / 10;
Point(1) = {0, 0, 0, p1};
Point(2) = {0, LY, 0, p1};
Point(3) = {LX, LY, 0, p1};
Point(4) = {LX, 0, 0, p1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Point(5) = {1, 1, 0, p2};
Point(6) = {1, 1.5, 0, p2};
Point(7) = {1.5, 1, 0, p2};
Point(8) = {1, 0.5, 0, p2};
Point(9) = {0.5, 1, 0, p2};

Circle(20) = {6, 5, 7};
Circle(21) = {7, 5, 8};
Circle(22) = {8, 5, 9};
Circle(23) = {9, 5, 6};
//+
Curve Loop(1) = {2, 3, 4, 1};
//+
Curve Loop(2) = {23, 20, 21, 22};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve(1) = {1, 2, 3, 4};
//+
Physical Curve(2) = {23, 20, 21, 22};
//+
Physical Surface(1) = {1};

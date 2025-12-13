SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 2.0;

LX = 2;
LY = 1;
r = 0.2;
x0 = LX / 2;
y0 = LY / 2;
p1 = LY / 60;
p2 = r / 20;
Point(1) = {0, 0, 0, p1};
Point(2) = {0, LY, 0, p1};
Point(3) = {LX, LY, 0, p1};
Point(4) = {LX, 0, 0, p1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//+
Point(5) = {0.5, 0.5, 0, p2};
//+
Point(6) = {0.5, 0.7, 0, p2};
//+
Point(7) = {0.7, 0.5, 0, p2};
//+
Point(8) = {0.5, 0.3, 0, p2};
//+
Point(9) = {0.3, 0.5, 0, p2};
//+
Point(10) = {1, 0.5, 0, p2};
//+
Point(11) = {0.8, 0.5, 0, p2};
//+
Point(12) = {1, 0.3, 0, p2};
//+
Point(13) = {1.2, 0.5, 0, p2};
//+
Point(14) = {1, 0.7, 0, p2};
//+
Point(15) = {1.5, 0.5, 0, p2};
//+
Point(16) = {1.5, 0.3, 0, p2};
//+
Point(17) = {1.7, 0.5, 0, p2};
//+
Point(18) = {1.5, 0.7, 0, p2};
//+
Point(19) = {1.3, 0.5, 0, p2};

Circle(20) = {6, 5, 7};
Circle(21) = {7, 5, 8};
Circle(22) = {8, 5, 9};
Circle(23) = {9, 5, 6};

Circle(24) = {11, 10, 12};
Circle(25) = {12, 10, 13};
Circle(26) = {11, 10, 14};
Circle(27) = {14, 10, 13};

Circle(28) = {18, 15, 19};
Circle(29) = {18, 15, 17};
Circle(30) = {17, 15, 16};
Circle(31) = {16, 15, 19};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {20, 21, 22, 23};
Line Loop(3) = {24, 25, 26, 27};
Line Loop(4) = {28, 29, 30, 31};

Plane Surface(1) = {1, 2, 3, 4};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

//+
Physical Curve(1) = {1, 2, 4, 3};
//+
Physical Curve(2) = {20, 23, 22, 21, 26, 24, 27, 25, 28, 31, 29, 30};
//+
Physical Surface(1) = {1, 2, 4};

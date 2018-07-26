// Gmsh project created on Sun Apr 29 18:25:09 2018
Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {1.0, 0.0, 0.0, 1.0};
Point(3) = {0.0, 1.0, 0.0, 1.0};
Point(4) = {1.0, 1.0, 0.0, 1.0};
Point(5) = {0.45, 0.25, 0.0, 1.0};
Point(6) = {0.45, 0.45, 0.0, 1.0};
Point(7) = {0.45, 0.65, 0.0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};

Circle(5) = {5, 6, 7};
Circle(6) = {7, 6, 5};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5};
Physical Line(6) = {6};

Plane Surface(1) = {1, 2};
Plane Surface(2) = {2};

Physical Surface(1) = {1};
Physical Surface(2) = {2};

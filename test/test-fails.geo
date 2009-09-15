Function SET3D

  Point(1) = {-L/2,   0,-D/2};
  Point(2) = {-L/2+l1,0,-D/2};
  Point(3) = { L/2-l1,0,-D/2};
  Point(4) = { L/2,   0,-D/2};

  Point(8) = {-L/2,   0, D/2};
  Point(7) = {-L/2+l1,0, D/2};
  Point(6) = { L/2-l1,0, D/2};
  Point(5) = { L/2,   0, D/2};

  Line(1) = {1,2};
  Line(2) = {2,7};
  Line(3) = {7,8};
  Line(4) = {8,1};

  Line Loop(1) = {1,2,3,4};

  Line(5) = {2,3};
  Line(6) = {3,6};
  Line(7) = {6,7};
  Line(8) = {7,2};

  Line Loop(2) = {5,6,7,8};

  Line(9)  = {3,4};
  Line(10) = {4,5};
  Line(11) = {5,6};
  Line(12) = {6,3};

  Line Loop(3) = {9,10,11,12};

  // Little dense box for molecule
  Point(9)  = {-L/2+l1,0,-boxd/2};
  Point(10) = { L/2-l1,0,-boxd/2};
  Point(11) = { L/2-l1,0, boxd/2};
  Point(12) = {-L/2+l1,0, boxd/2};

  Line(13) = {9,10};
  Line(14) = {10,11};
  Line(15) = {11,12};
  Line(16) = {12,9};

  Line Loop(4) = {13,14,15,16};

  Plane Surface(1) = {1};
  Plane Surface(2) = {2};
  Plane Surface(3) = {3};
  Plane Surface(4) = {4};

  Transfinite Line{1:12}=divisions Using Progression 1;
  Transfinite Line{13:16}=boxdivisions Using Progression 1;

  Transfinite Surface{1,2,3,4};
  Recombine Surface{1,2,3,4};

  Extrude {0,-h1,0}{
  	  Surface{1,2,3}; Layers{4}; Recombine;
  }

  Extrude {0,h2,0}{
  	  Surface{2}; Layers{4}; Recombine;
  }
  Extrude {0,boxh,0}{
  	  Surface{4}; Layers{4}; Recombine;
  }

  Return

divisions=3;
D  = 100;
h1 =  95;
h2 = 100;
L  = 150;
l2 =  68;
l1 = (L-l2)/2.0;

boxd = 13.3;
boxh = 3.89;
boxdivisions=8;

Call SET3D;
// Neumann boundary conditions
Physical Surface(127) = {37, 25, 33, 55, 73, 77, 69, 47, 99, 91, 104};
// Ground voltage; Dirichlet V0
Physical Surface(128) = {60, 38, 82};
// Left gate; Dirichlet V1
Physical Surface(129) = {1, 103};
// Right gate; Dirichlet V2
Physical Surface(130) = {95, 3};

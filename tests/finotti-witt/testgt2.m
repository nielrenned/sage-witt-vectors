load 'gt.m';

p:=7;

F:=GF(p);

// f:=[ [* [ F!1, F!0, F!2 ], 1, 2 *], [* [ F!1, F!1, F!1 ], 3, 2 *] ];

f:=[ [* [ F!1, F!0, F!2, F!1 ], 1, 0 *], [* [ F!0, F!2, F!1, F!0 ], 2 , 2 *], 
 [* [ F!1, F!1, F!1, F!2 ], 3, 2 *] ];


// pols:=etapols(p,2);
// bintab:=BinTab(p,2);

time pols:=etapols(p,3);
time bintab:=BinTab(p,3);


ResetMaximumMemoryUsage();
time g2 := GT2(f : bintab := bintab);
print "Memory Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);

ResetMaximumMemoryUsage();
time g1 := GT(f : pols := pols);
print "Memory Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);

PR:=Parent(g1[1]);

g1 eq [ PR!term : term in g2];

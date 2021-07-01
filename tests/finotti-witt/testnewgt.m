P<x,y>:=PolynomialRing(Integers(),2);

f:=x^10-5*x^2+x-1;
//f:=x^3*y + y^2 -5*x + 2*y -1;

p:=7;
k:=3;
n:=7;
ff:=Pol_GT_Form(f,p,n);
F:=GF(p^k);
v:=[ Random(F) : i in [0..n] ];
w:=[ Random(F) : i in [0..n] ];

//time pols:=etapols(p,n);

//time gt1:=GT(ff : vvars := v cat w, pols:=pols );
//time gt2:=GT_1(ff : vvars := v cat w, pols:=pols );

/* time gt1:=GT(ff : pols:=pols ); */
/* time gt2:=GT_1(ff : pols:=pols ); */



time tab:=bin_tab(10,p^(n+1));

/* time gt2:=GT2_1(ff : tab:=tab ); */
/* time gt1:=GT2(ff : tab:=tab ); */

time gt1:=GT2(ff : tab:=tab, vvars:=v cat w );
time gt2:=GT3(ff : tab:=tab, vvars:=v cat w );


PP:=Parent(gt1[1]);
gt1 eq [ PP!xx : xx in gt2 ];

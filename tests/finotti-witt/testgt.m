load 'gt.m';

p:=5;

F:=GF(p);

f:=[ [* [ F!1, F!0, F!2 ], 1, 2 *], [* [ F!1, F!1, F!1 ], 3, 2 *] ];


pols:=etapols(p,2);

PR<x0,x1,x2,y0,y1,y2>:=PolynomialRing(F,6);
x:=[ x0, x1, x2 ];
y:=[ y0, y1, y2 ];
ff:=newWittProd([ F!1, F!0, F!2 ], x : pols:=pols);                
y2:=newWittProd(y,y : pols:=pols);                                 
ff:=newWittProd( ff, y2 : pols:=pols);                             
x3:=newWittProd(x,x : pols:=pols);                                 
x3:=newWittProd(x3,x : pols:=pols);
tmp:=newWittProd( [ F!1, F!1, F!1 ], x3 : pols:=pols);
tmp:=newWittProd( tmp, y2 : pols:=pols);
ff:=newWittSum(ff, tmp : pols:=pols);
gg:=GT(f);
ff eq [ PR!term : term in gg ];



f:=[ [* [ F!1, F!0, F!2, F!1 ], 1, 0 *], [* [ F!0, F!2, F!1, F!0 ], 2 , 2 *], 
[* [ F!1, F!1, F!1, F!2 ], 3, 2 *] ];

gg:=GT(f);

pols:=etapols(p,3);
PR<x0,x1,x2,x3,y0,y1,y2,y3>:=PolynomialRing(F,8);
x:=[ x0, x1, x2, x3 ];
y:=[ y0, y1, y2, y3 ];
ff:=newWittProd([ F!1, F!0, F!2, F!1 ], x : pols:=pols);
x2:=newWittProd(x,x : pols:=pols);
y2:=newWittProd(y,y : pols:=pols);
tmp:=newWittProd( x2, y2 : pols:=pols);
tmp:=newWittProd([ F!0, F!2, F!1, F!0 ], tmp : pols:=pols);
ff:=newWittSum(ff, tmp : pols:=pols);
tmp:=newWittProd(x, x2 : pols:=pols);
tmp:=newWittProd(tmp, y2 : pols:=pols);
tmp:=newWittProd(tmp, [ F!1, F!1, F!1, F!2 ] : pols:=pols);
ff:=newWittSum(ff, tmp : pols:=pols);
ff eq [ PR!term : term in gg ];

/*


////////////////////////////

p:=5;
n:=3;

load 'witt_gen.m';

vecS:=WittSumPol(n,p);
vecP:=WittProdPol(n,p);
primr:=GenPrimRoots(p,n+1);

F:=GF(p);
P:=PolynomialRing(F,2*n+2);
AssignNames(~P,
            [ "x" cat IntegerToString(i) : i in [0..n] ] cat 
            [ "y" cat IntegerToString(i) : i in [0..n] ] );


// load 'load_js.m';

PP<X,Y>:=PolynomialRing(Integers(),2);
PPP:=PolynomialRing(ResidueClassRing(p^(n+2)),2);

mp:=PP!(PPP!ClassicalModularPolynomial(p));


vx:=[ P.i : i in [1..(n+1)] ];
vy:=[ P.i : i in [(n+2)..(2*n+2)] ];

vxs:=<[P!1] cat [P!0 : i in [1..n] ], vx>;


maxdeg:=Max(Degree(mp,1), Degree(mp,2));
// compute powers of vx
for i in [2..maxdeg] do
    Append(~vxs,WittProd(vxs[i],vx : vecP:=vecP));
end for;

evalv:=vy cat [ F!0 : i in [1..(n+1)] ];
// compute powers of vy
vys:=[ [ Evaluate(t,evalv) : t in vxs[i] ] : i in [1..(#vxs)] ];

res:=[ P!0 : i in [0..n] ];

// compute GT
for t in Terms(mp) do
    coef:=LeadingCoefficient(t);
    degX:=Degree(t,1);
    degY:=Degree(t,2);
    tmp:=IntToWitt(coef,p,n : primr:=primr);
    tmp:=WittProd(tmp,vxs[degX+1] : vecP:=vecP);
    tmp:=WittProd(tmp,vys[degY+1] : vecP:=vecP);
    
    res:=WittSum(res,tmp : vecS:=vecS);
    
end for;

// NOW, NEW WAY

load 'etas.m';
load 'mp_gt.m';


pols:=etapols(p,n);


mmp:=Pol_GT_Form(mp,p,n : primr:=primr );

res2:=GT(mmp : pols:=pols, primr:=primr);

res eq [ P!x : x in res2 ];

// Now, using the new form mmp

res3:=[ P!0 : i in [0..n] ];

// compute GT
for t in mmp do
    coef:=t[1];
    degX:=t[2];
    degY:=t[3];
    tmp:=newWittProd(coef,vxs[degX+1] : pols:=pols);
    tmp:=newWittProd(tmp,vys[degY+1] : pols:=pols);
    
    res3:=newWittSum(res3,tmp : pols:=pols);
    
end for;



*/

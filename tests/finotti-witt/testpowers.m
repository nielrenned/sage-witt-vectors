p:=3;
m:=5;
F:=GF(p^m);
n:=5;
v:=[ Random(F) : i in [0..n] ];
k:=4;

pols:=etapols(p,n);
bintab:=BinTab(p,n);

v1:=WittPower(v,k : choice:=3,bintab:=bintab);

v2:=v;
for i in [2..k] do
    v2:=WittProd(v2,v : pols:=pols);
end for;

v1 eq v2;

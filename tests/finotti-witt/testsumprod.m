load 'etas.m';

p:=7;PR<X0,X1,X2,X3,X4,Y0,Y1,Y2,Y3,Y4>:=PolynomialRing(GF(p),10);

ResetMaximumMemoryUsage();
time w1:=newWittProd([X0,X1,X2,X3,X4],[Y0,Y1,Y2,Y3,Y4]);
print "Memory Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);


ResetMaximumMemoryUsage();
time w2:=newWittProd2([X0,X1,X2,X3,X4],[Y0,Y1,Y2,Y3,Y4]);
print "Memory Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);


ResetMaximumMemoryUsage();
time w3:=newWittProd3([X0,X1,X2,X3,X4],[Y0,Y1,Y2,Y3,Y4]);
print "Memory Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);

w1 eq w2;
w2 eq w3;

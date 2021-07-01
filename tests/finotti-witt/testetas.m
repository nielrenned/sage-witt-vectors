load 'etas.m';

p:=7;F:=GF(p);P<X0,X1,X2,X3,X4,Y0,Y1,Y2,Y3,Y4>:=PolynomialRing(F,10);

ResetMaximumMemoryUsage();
time v1:=vetav(p,3,[X0,X1,X2,X3]);
print "Memory Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);

ResetMaximumMemoryUsage();
time v2:=vetav2(p,3,[X0,X1,X2,X3]);
print "Memory Usage (in MB): ", GetMaximumMemoryUsage()/((1024.0)^2);

v1 eq v2;


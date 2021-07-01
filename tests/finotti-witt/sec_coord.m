Coef1 := function(n,p,i)
    return Binomial(n*p,i) div p;
end function;


Coef2 := function(n,p,i)
    return (Binomial(n,i)-Binomial(n*p,i*p)) div p;
end function;



fSpcFctn := function(x,y,n)

    p:=Characteristic(Parent(x));
    term:= -&+[ Coef1(n,p,i)*x^i*y^(n*p-i) : i in [1..(n*p-1)] | (i mod p) ne 0 ];

    if (p mod n) eq 0 then
        return term;
    end if;

    term2:=&+[ Coef2(n,p,i)*x^(p*i)*y^(p*(n-i)) : i in [1..(n-1)]];

    return term2 + term;

end function;



SpcFctn := function(f,n)

    if f eq 0 then
        return 0;
    end if;

    P:=Parent(f);
    p:=Characteristic(P);
    term:=Terms(f)[1];
    g:=f-term;

    if g eq P!0 then
        return 0;
    end if;

    res1:=fSpcFctn(term,g,n);
    
    if #Terms(g) eq 1 then
        return res1;
    end if;

    if n gt 1 then 
        res2:= &+[ Binomial(n,i)*term^(p*i)*$$(g,n-i) : i in [1..(n-1)]];
    else
        res2:=0;
    end if;
        
    return $$(g,n) + res1 + res2;
    
end function;


SpcFct := function(f)
    return SpcFctn(f,1);
end function;

fSpcFct := function(x,y)
    return fSpcFctn(x,y,1);
end function;


/*

TESTED!!!

RR<X,Y,Z,W>:=PolynomialRing(Integers(),4);

p:=17;R<x,y,z,w>:=PolynomialRing(GF(p),4);
n:=1;

p1:=R!(((X^p+Y^p+Z^p+W^p)^n-(X+Y+Z+W)^(p*n)) div p);
q1:=SpcFctn(x+y+z+w,n);

p1 eq q1;


*/


SpcFctv := function(v)

    if (#v eq 1) or (#v eq 0) then
        return 0;
    end if;

    term:=v[1];
    P:=Parent(term);
    p:=Characteristic(P);
    w:=Remove(v,1);

    g:=&+w;

    return fSpcFct(term,g) + $$(w);

end function;


NextCoord := function(a,p)
    return (a^p - a) div p;
end function;


Coef3:= function(p,i)
    return NextCoord(Binomial(p,i) div p,p);
end function;


// mu
Spc2Fct := function(f)
    
    P:=Parent(f);
    p:=Characteristic(P);
    term:=Terms(f)[1];
    g:=f-term;
    if g eq P!0 then
        return 0;
    end if;

    res1:= &+[ Coef3(p,i)*term^(p*i)*g^(p*(p-i)) : i in [1..(p-1)] ];

    v:=[ -Coef1(1,p,i)*term^i*g^(p-i) : i in [1..(p-1)] ];
    
    // sf:=fSpcFct(term,g);
    sf := &+v;
    res2:= SpcFctv(v);
    
    if #Terms(g) eq 1 then
        return res1 + res2;
    end if;

    res3:=fSpcFct(SpcFct(g),sf);

    res4:= -&+[ Coef1(1,p,i)*term^(p*i)*SpcFctn(g,p-i) : i in [1..(p-1)] ];
    
    return $$(g) + res1 + res2 + res3 + res4;

end function;



/*

TESTED!!!!!

RR<X,Y,Z,W>:=PolynomialRing(Integers(),4);

p:=11;R<x,y,z,w>:=PolynomialRing(GF(p),4);

time p1:=R!( (((X^(p^2)+Y^(p^2)+Z^(p^2)+W^(p^2) - (X^p+Y^p+Z^p+W^p)^p) div p) - ((X^p+Y^p+Z^p+W^p - (X+Y+Z+W)^p) div p)^p) div p );
time q1:=Spc2Fct(x+y+z+w);

p1 eq q1;


*/

Coef4 := function(p,i)
    return (Binomial(p^2,i) div p^2);
end function;


// psi_2
Spc3Fct := function(f)

    P:=Parent(f);
    p:=Characteristic(P);
    term:=Terms(f)[1];
    g:=f-term;

    if g eq P!0 then
        return 0;
    end if;

    res1:= -&+[ Coef4(p,i)*term^i*g^(p^2-i) : i in [1..(p^2-1)] | (i mod p) ne 0 ]; 

    if #Terms(g) eq 1 then
        return res1;
    end if;

    res2:= &+[ Coef1(1,p,i)*term^(p*i)*SpcFctn(g,p-i) : i in [1..(p-1)] ];

    return $$(g) + res1 + res2;

end function;

/*


TESTED!!!!

RR<X,Y,Z,W>:=PolynomialRing(Integers(),4);

p:=11;R<x,y,z,w>:=PolynomialRing(GF(p),4);

time p1:=R!(((X^p+Y^p+Z^p+W^p)^p-(X+Y+Z+W)^(p^2)) div p^2);
time q1:=Spc3Fct(x+y+z+w);

p1 eq q1;





*/


/*

TESTED!!!!

Now test second coord.

load 'witt_gen.m';

p:=31;R<x0,x1,x2,y0,y1,y2>:=PolynomialRing(GF(p),6);

time p1:=R!WittSumPol(2,p)[3];

time q1:=x2+y2+fSpcFct(x1,y1)+fSpcFct(x1+y1,SpcFct(x0+y0))+Spc3Fct(x0+y0)+Spc2Fct(x0+y0);

p1 eq q1;

*/


/*

Now, test the eq. of the elliptic curve

load 'witt_gen.m';


p:=11;P<a0,a1,a2,b0,b1,b2,x0,x1,x2>:=PolynomialRing(GF(p),9);
PP<A0,A1,A2,B0,B1,B2,X0,X1,X2>:=PolynomialRing(Integers(),9);


x:=[x0,x1,x2];
a:=[a0,a1,a2];
b:=[b0,b1,b2];

f:=x0^3+a0*x0+b0;
sf:=SpcFct(f);

rhs1:=WittPwr(x,3);
rhs1:=WittSum(rhs1,WittProd(a,x));
rhs1:=WittSum(rhs1,b);

rhs1[2] eq (3*x0^(2*p)+a0^p)*x1 + a1*x0^p + b1 + sf;

rhs1:=rhs1[3];


ev:=[ 3*x0^(2*p)*x1 , a0^p*x1, a1*x0^p , b1]; 

term1:= (3*x0^(2*p^2)+a0^(p^2))*x2 + 3*x0^(p^2)*x1^(2*p) +
    (-NextCoord(3,p)*x0^(2*p^2)+a1^p)*x1^p + a2*x0^(p^2) + b2;

term2:= SpcFctv(ev);

term3:= fSpcFct( (3*x0^(2*p)+a0^p)*x1 + a1*x0^p + b1 , sf );

term4:= Spc2Fct(f);

term5:= Spc3Fct(f);

rhs2:=term1+term2+term3+term4+term5;

rhs2 eq rhs1;


Term1:=(3*X0^(2*p^2)+A0^(p^2))*X2 + 3*X0^(p^2)*X1^(2*p) +
    (-NextCoord(3,p)*X0^(2*p^2)+A1^p)*X1^p + A2*X0^(p^2) + B2;

Term2:=(A0^(p^2)*X1^p + A1^p*X0^(p^2) + B1^p + (3*X0^(2*p)*X1)^p -
    (A0^p*X1 + A1*X0^p + B1 + 3*X0^(2*p)*X1 +
        ((X0^(3*p) + A0^p*X0^p + B0^p - (X0^3 + A0*X0 + B0)^p) div p))^p);

Term3:=(X0^(3*p^2) + A0^(p^2)*X0^(p^2) + B0^(p^2) - (X0^3 + A0*X0 + B0)^(p^2));

rhs3:= Term1 + (( p* Term2 + Term3 ) div p^2);

rhs3:= P!rhs3;

rhs3 eq rhs1;

*/


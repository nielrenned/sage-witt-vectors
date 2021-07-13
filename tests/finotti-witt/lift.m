// computes canonical and minimal degree liftings
// using new methods (etas, formulas for GT, etc.)
// It works for n>2 when computing MINIMAL DEGREE LIFTS!!!!
// But it still needs work for canonical liftings (for n>2)...


// we need the functions for the Greenberg transform
load "gt.m"; // which loads etas.m, and witt_gen.m


function my_int(f)
    // formal integral of polynomials
    deg:=Degree(f);
    P:=Parent(f);
    p:=Characteristic(P);
    res:=0;
    for d in [0..deg] do
        if ((d+1) mod p) ne 0 then
            res+:=Coefficient(f,d)/(d+1)*P.1^(d+1);
        end if;
    end for;
    return res;

end function;


function lift(a0,b0,n : pols:=[], minimal:=false)
    // lenght (n+1) -- up to a_n, b_n
    // returns vectors va, vb, vF, vH
    // pols are etapols from etas.m
    // computes canonical liftins by default
    // set minimal:=true to compute minimal degree liftings

    F:=Parent(a0);
    p:=Characteristic(F);

    if #pols eq 0 then
        pols:=etapols(p,n);
    end if;

    // ring for the resulting polynomials in x0
    Pres<x0>:=PolynomialRing(F);


    // results
    resa:=[a0];
    resb:=[b0];
    resF:=[x0];
    resH:=[Pres!1];


    f:=x0^3+a0*x0+b0;
    ff:=f^((p-1) div 2);

    // Hasse Invariant
    HI:=F!(Coefficient(ff,p-1));


    // main loop
    for i in [1..n] do

        //polynomial ring for this coordinate
        M:=(3*p^(i-1)-1) div 2; // c_i's
        // N:=((i+3)*p^i - i*p^(i-1) - 3) div 2; // d_i's
        Pi:=PolynomialRing(F,2 + 2 + (M+1) + 1);
        // x_0, y_0, a_n, b_n, c_i's, Hn

        // a temporary power of f to help compute F_i's:
        // f^((p^i-1)/2)
        if i eq 1 then
           tmppf:=ff;
        else
           tmppf:=tmppf^p*ff;
        end if;

        Fi:= HI^(-(p^i-1) div (p-1))*tmppf - x0^(p^i-1);
        if i gt 1 then
            Fi-:=&+[ resF[j+1]^(p^(i-j)-1)*Derivative(resF[j+1]) :
                     j in [1..(i-1)] ];
        end if;
        Fi:=my_int(Fi);
        Fi:=Evaluate(Fi,Pi.1);

        // will make c_(p^(n-1))=0
        Fi+:=&+[ Pi.(5+j)*(Pi.1)^(p*j) : j in [0..M] | j ne p^(i-1) ];

        // //////////////////////////////////////////////////////////
        // ADD CONDITION IF NOT MINIMAL!!!!!!!

        if (not minimal) and (i eq 2) then
            tmp:= F!(3/4)*resF[2]^2;
            Fi+:= &+[(Coefficient(tmp,Integers()!(i/p+p)))^p*Pi.1^i :
                     i in [Integers()!((3*p^2+p)/2) .. 2*p^2-p by p]];
        end if;

        // Condition for i ge 3
        // if (not minimal) and (i ge 3) then
        //
        // end if;


        // //////////////////////////////////////////////////////////

        va:=[ Pi!x : x in resa ] cat [Pi.3];
        vb:=[ Pi!x : x in resb ] cat [Pi.4];
        vF:=[ Evaluate(x,Pi.1) : x in resF ] cat [Fi];
        vG:=[ Pi.2*Evaluate(x,Pi.1) : x in resH ] cat [Pi.2*Pi.(M+6)];
        vone:=[Pi!1] cat [ Pi!0 : j in [1..i]];

        vvars:=vF cat vG;


        GTx:=GT( [[* vone,3,0 *], [* va,1,0 *], [* vb,0,0 *]] : pols:=pols, vvars:=vvars);

        GTy:=GT( [[* vone,0,2 *]] : pols:=pols, vvars:=vvars);

        RHS:=GTx[i+1];
        delete GTx;
        LHS:=GTy[i+1];
        delete GTy;
        LHS:=Coefficient(LHS,M+6,0); // terms without Hn

        // now, replace y0^2 by f
        deg:=Degree(LHS,2);
        tmppf2:=1;
        tmpLHS:=Coefficient(LHS,2,0);
        for d in [j : j in [2..deg] | (j mod 2) eq 0 ] do
           tmppf2*:=Evaluate(f,Pi.1);
           tmpLHS+:=Coefficient(LHS,2,d)*tmppf2;
        end for;

        RHS-:=tmpLHS;

        delete tmpLHS, LHS, tmppf2;

        tmppf2:=Evaluate(tmppf*f,Pi.1); //f^((p^i+1)/2)
        // long division

        RHS*:=Pi!(1/2);
        deg1:=Degree(RHS,Pi.1);
        deg2:=(3*(p^i+1) div 2);
        quo:=0;

        while deg1 ge deg2 do
            lterm:=Coefficient(RHS,1,deg1)*((Pi.1)^(deg1-deg2));
            RHS-:=lterm*tmppf2;
            quo+:=lterm;
            deg1:=Degree(RHS,1);
        end while;
        // RHS is now the remainder

        delete tmppf2;

        vrem:=Coefficients(RHS,1);

        // matrix of coefficients
        neqts:=#vrem;
        Mat:=Matrix(F,2+M,neqts,
             [[F!(Coefficient(vrem[j],k+2,1)) : j in [1..neqts]]
                 : k in [ ii : ii in [1..(2+(M+1))] | ii ne 3+p^(i-1) ]]);

        vec:=Vector(F,[ -Evaluate(vrem[j],[ 0 : k in [1..(2+2+(M+1)+1)]]) : j in [1..neqts]]);

        vsol:=Solution(Mat,vec);

        // to convert solutoions to Pres/F
        evalvec:=[x0,0] cat [ vsol[j] : j in [1..2+p^(i-1)] ] cat [0] cat
                 [vsol[j] : j in [2+p^(i-1)+1..(2+M)]] cat [0];

        //print evalvec;

        Append(~resa,vsol[1]);
        Append(~resb,vsol[2]);
        Append(~resF,Evaluate(Fi,evalvec));
        Append(~resH,Evaluate(quo,evalvec));


    end for;

    return resa, resb, resF, resH;

end function;

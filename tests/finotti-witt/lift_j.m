// Compute the formulas formulas for the j-invariant
// of the canonical lifting


// we need the Greenberg transform
load "gt.m";

// reduces polynomial modulo p^n
// since the modular polynomial is huge...
function reducepol(f,p,n)
    res := 0;
    for mon in Monomials(f) do
        coef := MonomialCoefficient(f,mon) mod p^n;
        res +:=coef*mon;
    end for;
    return res;
end function;


// we will need to take p-th roots!
function pthroot(x : n:=1)
    // takes the p^n-th root
    F:=Parent(x);
    p:=Characteristic(F);

    // let's check if we are in a function field
    if IsField(F) then
	num:=Numerator(x);
	den:=Denominator(x);
	R:=Parent(num);
    else
	num:=x;
	den:=F!1;
	R:=F;
    end if;

    r:=Rank(R);


    newnum:=0;
    for term in Terms(num) do
        tterm:=term;
	for j in [1..r] do
            if r eq 1 then
                deg:=Degree(tterm);
            else
		deg:=Degree(tterm,j);
            end if;

	    error if (deg mod p^n) ne 0,
		  "Exponent of ", tterm, " not divisible by ", p^n;

	    if r eq 1 then
                tterm:=Evaluate(tterm,1)*(R.1)^(deg div p^n);
            else
                tterm:=Evaluate(tterm,j,1)*(R.j)^(deg div p^n);
            end if;

	end for;
	newnum+:=tterm;
    end for;

    if Degree(den) eq 0 then
        newden:=den;
    else
        newden:=0;
        for term in Terms(den) do
            tterm:=term;
	    for j in [1..r] do
                if r eq 1 then
                    deg:=Degree(tterm);
                else
		    deg:=Degree(tterm,j);
                end if;

		error if (deg mod p^n) ne 0,
		      "Exponent of ", tterm, " not divisible by ", p^n;

		if r eq 1 then
                    tterm:=Evaluate(tterm,1)*(R.1)^(deg div p^n);
                else
                    tterm:=Evaluate(tterm,j,1)*(R.j)^(deg div p^n);
                end if;

	    end for;
	    newden+:=tterm;

        end for;
    end if;

    return F!newnum/F!newden;

end function;



// formulas for J_1, J_2, ... , J_n
function lift_j(p,n : choice:=1, pols:=[], bintab:=[], tab:=[])
    // n is length - 1.
    F<j0> := RationalFunctionField(GF(p));
    P<ji> := PolynomialRing(F);
    mpol := reducepol(ClassicalModularPolynomial(p),p,n+1);
    vj := [ P!j0 ];

    if (choice eq 1) and (#pols eq 0) then
        pols := etapols(p,n);
    end if;

    if (choice ne 1) and (#bintab eq 0) then
        bintab := BinTab(p,n);
    end if;

    // maxdeg = p+1
    if #tab eq 0 then
        tab:=bin_tab(p+1,p^(n+1)); // table of binomials to use
    end if;

    for i in [1..n] do
        mp := Pol_GT_Form(mpol,p,i);
        v1 := [ P!x : x in vj ] cat [ ji ];
        v2 := [ x^p : x in v1 ];
        eqtn := (GT(mp : choice:=choice, pols:=pols, bintab:=bintab, tab:=tab, vvars := v1 cat v2))[i+1];
        vji := -Coefficient(eqtn,0)/Coefficient(eqtn,p);

        delete eqtn;

        Append(~vj,F!(pthroot(vji)));

        delete vji;

    end for;

    return [ F!x : x in vj ];

end function;



// This one lifts a concrete j0, not formulas.
// It probably just works for finite fields, due to the
// use of Root( , p) to take p-th roots.
// j0 cannot be in GF(p^2)!!!
function lift_j0(j0,n : choice:=1, pols:=[], bintab:=[], tab:=[])
    // n is length - 1.

    F := Parent(j0); // hopefully a finite field
    p := Characteristic(F);

    if j0 in GF(p^2) then
        print("j0 cannot be in GF(p^2)");
        return [j0];
    end if;

    P<ji> := PolynomialRing(F);
    mpol := reducepol(ClassicalModularPolynomial(p),p,n+1);
    vj := [ P!j0 ];

    if (choice eq 1) and (#pols eq 0) then
        pols := etapols(p,n);
    end if;

    if (choice ne 1) and (#bintab eq 0) then
        bintab := BinTab(p,n);
    end if;

    // maxdeg = p+1
    if #tab eq 0 then
        tab:=bin_tab(p+1,p^(n+1)); // table of binomials to use
    end if;

    for i in [1..n] do
        mp := Pol_GT_Form(mpol,p,i);
        v1 := [ P!x : x in vj ] cat [ ji ];
        v2 := [ x^p : x in v1 ];
        eqtn := (GT(mp : choice:=choice, pols:=pols, bintab:=bintab, tab:=tab, vvars := v1 cat v2))[i+1];
        vji := -Coefficient(eqtn,0)/Coefficient(eqtn,p);

        delete eqtn;

        Append(~vj,Root(vji,p));

        delete vji;

    end for;

    return [ F!x : x in vj ];

end function;


// deg (not exact) of Gi (denom of j_i)
// needed for interpolation
function degGi(p,i)
    r := (p-1) div 2;
    r1 := Ceiling(r/3);
    r2 := Floor(r/2);

    if (p mod 6) eq 5 then
        si := Max(0,((i-3)*p^i + i*p^(i-1))/3);
    else
        si := 0;
    end if;

    if (p mod 4) eq 3 then
        ri := (i-1)*p^(i-1);
    else
        ri := 0;
    end if;

    return (r2-r1)*(i*p^(i-1)+(i-1)*p^(i-2)) + ri + si;
end function;




// Computes formulas using interpolation.
// Needs p >= 5!
// Much slower, but uses a lot less memory.
function lift_j_int(p,n: choice:=1, pols:=[], bintab:=[], tab:=[])

    if p lt 5 then
        print("Characteristic must be > 3.");
        return [];
    end if;

    // precompute stuff!
    if (choice eq 1) and (#pols eq 0) then
        pols := etapols(p,n);
    end if;

    if (choice ne 1) and (#bintab eq 0) then
        bintab := BinTab(p,n);
    end if;

    // maxdeg = p+1
    if #tab eq 0 then
        tab:=bin_tab(p+1,p^(n+1)); // table of binomials to use
    end if;


    // enough elements to do all interpolations!
    k := Ceiling(Log(p,degGi(p,n)+p^n+1));
    F<a> := GF(p^k);
    P<j0> := RationalFunctionField(F);

    r := (p-1) div 2;
    r1 := Ceiling(r/3);
    r2 := Floor(r/2);

    Spol := F!(-2/9)^r *
            &+[ Binomial(r,i) * Binomial(i,3*i-r) * F!(27/4)^i * j0^(i-r1) * (j0-1728)^(r2-i) : i in [r1..r2] ];

    if p eq 31 then
        iota := 2;
    else
        iota := 1;
    end if;


    I := []; // inputs (for interpolation)
    M := []; // values (for interpolation)

    degFn := degGi(p,n) + p^n - iota;

    // compute enough "numerical" liftings for interpolation
    j := 1; // power of gen of F
    while #I lt degFn+1 do
        x := a^j;
        // make sure not in GF(p^2)!
        if x^(p^2) - x ne 0 then
            Append(~I,x);
            Append(~M,lift_j0(x,n : choice:=choice, pols:=pols, bintab:=bintab, tab:=tab));
        end if;
        j +:= 1;
    end while;

    vj := [ j0 ];

    // interpolate each coordinate
    for i in [1..n] do
        degFi := degGi(p,i) + p^i - iota;

        if i eq 1 then
            Gi := Spol;
        else
            Gi := Spol^(i*p^(i-1)+(i-1)*p^(i-2));
        end if;

        if (p mod 6) eq 5 then
            si := Max(0,((i-3)*p^i + i*p^(i-1)) div 3);
            Gi *:= j0^si;
        end if;

        if (p mod 4) eq 3 then
            ri := (i-1)*p^(i-1);
            Gi *:= (j0-1728)^ri;
        end if;

        vI := [ I[j] : j in [1..degFi + 1] ];
        vV := [ M[j][i+1]*Evaluate(Gi,I[j]) : j in [1..degFi + 1] ];

        Append(~vj,P!Interpolation(vI,vV)/Gi);

    end for;

    PP<j0>:=RationalFunctionField(GF(p));

    return [ PP!x : x in vj];

end function;


// Computes the NON-universal Weierstrass coefficients using
// the lift of the j-invariant via:
//    a = lamba^4 * 27 * j/(4*(1728 - j))
//    b = lamba^6 * 27 * j/(4*(1728 - j))
// where
//    lambda^2 = (b0/a0, 0, 0, ...)
function jweier(p,n : choice:=1, pols:=[], bintab:=[], tab:=[])

    if (choice eq 1) and (#pols eq 0) then
        pols := etapols(p,n);
    end if;

    if (choice ne 1) and (#bintab eq 0) then
        bintab := BinTab(p,n);
    end if;

    // maxdeg = p+1
    if #tab eq 0 then
        tab:=bin_tab(p+1,p^(n+1)); // table of binomials to use
    end if;

    vj := lift_j(p,n : choice:=choice, pols:=pols, bintab:=bintab, tab:=tab);
    P := Universe(vj);

    v1728 := IntToWitt(1728,p,n);
    v1728 := [ P!x : x in v1728 ];

    v27 := IntToWitt(27,p,n);
    v27 := [ P!x : x in v27 ];

    v4 := IntToWitt(4,p,n);
    v4 := [ P!x : x in v4 ];


    vv := WittDiff(v1728,vj : choice:=choice, pols:=pols, bintab:=bintab);
    vv := WittInv(vv : choice:=choice, pols:=pols, bintab:=bintab);
    vv := WittProd(vj,vv : choice:=choice, pols:=pols, bintab:=bintab);
    vv := WittProd(v27,vv: choice:=choice, pols:=pols, bintab:=bintab);
    vv := WittDiv(vv,v4 : choice:=choice, pols:=pols, bintab:=bintab);

    PP<a0,b0>:=RationalFunctionField(GF(p),2);

    j:=1728*4*a0^3/(4*a0^3+27*b0^2);

    vv := [ Evaluate(term,j) : term in vv ];

    va := [ vv[i]*(b0^2/a0^2)^(p^(i-1)) : i in [1..(n+1)]  ];
    vb := [ vv[i]*(b0^3/a0^3)^(p^(i-1)) : i in [1..(n+1)]  ];

    return [va, vb];
end function;

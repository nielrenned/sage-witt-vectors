# Witt Vectors and Canonical Liftings #

This contains various routines to perform computations with Witt
 vectors and canonical liftings on
 [MAGMA](http://magma.maths.usyd.edu.au/magma/).

**Note:** These routines do not have many checks for consistency
of the input!  So it might output errors or, even worse, garbage, if
the input is not consistent with what is expected it to be.


## Files ##

* `etas.m`: Contains the auxiliary "eta functions" for computations
  with Witt vectors.

* `witt.m`: Contains functions to perform operations with Witt
  vectors, except for powers.  (Depends on `etas.m`.)

* `gt.m`: Contains routines to compute the Greenberg transform of
  polynomials in two variables, allowing in particular to evaluate
  these polynomials at a pair of Witt vectors.  It also contains
  functions to compute powers of Witt vectors.  (Depends on `witt.m`.)

* `lift.m`: Contains routines to compute canonical liftings and the
  elliptic Teichmüller lift (up to 3 coordinates) and minimal degree
  liftings of ordinary elliptic curves.  (Depends on `gt.m`.)

* `lift_j.m`: Contains routines to find the coordinates of the
  j-invariant of the canonical lifting.  (Depends on `gt.m`.)  It also
  contains routines to compute (non-universal) Weierstrass
  coefficients of the canonical lifting.

* `test*.m`: These are test files and can be ignored.  (I leave them
  for when I want to test changes.)


## Theory ##

The routines are based on the methods developed in [Computations with
Witt Vectors and the Greenberg Transform][comp].  It provides more
efficient ways to perform computations with Witt vectors.  In
particular, one does not need the (*huge*!) polynomials that define
the sum and product of Witt vectors.

The algorithm to compute the canonical lifting is described in
[Degrees of the Elliptic Teichmuller Lift][canlift], and the one to
compute minimal degree liftings is described in [Minimal Degree
Liftings of Hyperelliptic Curves][minlift].

For the liftings of the j-invariant, see [Coordinates of the
j-Invariant of the Canonical Lifting][jinv].

For the Weierstrass coefficients using the lift of the j-invariant,
see [Weierstrass Coefficients of the Canonical Lifting][wcoef].


[comp]: http://www.math.utk.edu/~finotti/papers/witt.pdf

[canlift]: http://www.math.utk.edu/~finotti/papers/degs.pdf

[minlift]: http://www.math.utk.edu/~finotti/papers/mindeg.pdf

[jinv]: http://www.math.utk.edu/~finotti/papers/jn.pdf

[wcoef]: http://www.math.utk.edu/~finotti/papers/wcoef.pdf



## Usage ##

### Eta Polynomials ###

The "eta polynomials" (described in Section 5 of [Computations with Witt
Vectors and the Greenberg Transform][comp]) that are used in all
computations can be computed in three different ways:

1. We can store some simpler version (in two variables) in memory and
   used them in the computations.

2. Store some more basic data that take much less memory and compute
   the eta polynomials "on the fly".  (This is describe in Section 8 of
   [Computations with Witt Vectors and the Greenberg
   Transform][comp]).  This is usually slower, but requires less
   memory.  In some cases, where the eta polynomials are too large, and
   one does not need many computations, it might be actually faster
   than method 1.
   
3. Use method 2, but store some values of resulting computations of
   the eta polynomial, so that we don't need to recompute them, which is
   likely to happen with method 2.  Of course, it uses more memory
   than method 2, but likely less memory than method 1.


Most of the functions below, which involve computations with Witt
vectors, use these eta polynomials, and therefore have options for the
user to choose the method to be used.  The choice of method is done
via an optional argument `choice`.  Setting, for instance, `choice:=2`
tells the function to use the method 2 above.  The default choice (so
no need to be specified) is method 1, which is probably the best
choice "in average".

As mentioned above, method 1 stores some shorter eta polynomials,
while methods 2 and 3 store much smaller data, which are basically
some binomials.  In principle the corresponding data, eta polynomials
or binomials, is computed every time a function is called.  But it can
also be computed separately, saved, and then given to functions via
optional arguments.  An optional argument `pols` gives the eta
polynomials to be used in method 1, while the optional argument
`bintab` gives the "binomial table" that is used in methods 2 and 3.
(This will be clearer with the examples in the sections below.)

To compute the eta polynomials necessary for computations with Witt
vectors of length `n+1` in characteristic `p` we call `etapols(p,n)`.
To compute the table of binomials we call `BinTab(p,n)`.
```
> etapols(2,4);
[
    x*y,
    x^3*y + x*y^3,
    x^7*y + x^5*y^3 + x^3*y^5 + x*y^7,
    x^15*y + x^13*y^3 + x^12*y^4 + x^11*y^5 + x^10*y^6 + x^9*y^7 + x^7*y^9 + x^6*y^10 + x^5*y^11 + x^4*y^12 + x^3*y^13 + x*y^15
]

> BinTab(3,4);
[
    [ 2, 2 ],
    [ 2, 2, 0, 1, 1, 0, 2, 2 ],
    [ 2, 2, 0, 1, 1, 2, 2, 2, 0, 1, 1, 1, 2, 2, 1, 1, 1, 0, 2, 2, 2, 1, 1, 0, 2, 2 ],
    [ 2, 2, 0, 1, 1, 2, 2, 2, 0, 1, 1, 2, 2, 2, 2, 1, 1, 0, 2, 2, 2, 1, 1, 0, 2, 2, 2, 1, 1, 2, 2, 2, 0, 1, 1, 1, 2, 2, 2, 1, 1, 2, 
    2, 2, 1, 1, 1, 0, 2, 2, 2, 1, 1, 2, 2, 2, 0, 1, 1, 2, 2, 2, 0, 1, 1, 2, 2, 2, 2, 1, 1, 0, 2, 2, 2, 1, 1, 0, 2, 2 ]

```

**Important:** As you can see above, the argument that is passed for
length `n+1` is actually `n`.  This is because we often see vectors of
length `n+1` as `[a0,a1,...,an]`, so `n` corresponds to the last
coordinate (when counting from 0.)  *We use this convention through
out!*


### Basic Computations with Witt Vectors ###



The following functions are provided in `witt.m`:

* `WittSum`: adds two Witt vectors.
* `WittProd`: multiplies two Witt vectors.
* `WittNeg`: gives the additive inverse (i.e., negative) of a Witt
  vector.
* `WittInv`: gives the multiplicative inverse of a Witt vector (if
  exists).
* `WittDiff`: gives the difference of two Witt vectors.
* `WittDiv`: gives the quotient of two Witt vectors (if exists).
  
The file `gt.m` also provides `WittPower` to compute powers of Witt
vectors.
  
As explained above, each of these functions have an optional argument
`choice` that allows you to choose which method (1, 2 or 3, as above)
to use in the computations, and it defaults to method 1.

**Note:** If any of the above functions are called to compute with
Witt vectors over finite fields, then the operations are performed by
converting them to elements of unramified extensions of p-adic
integers, performing the operation in that ring, and then converting
back to Witt vectors.  This method is *a lot* more efficient than any
of the methods above!
  
Here is one simple example:
```
> load "gt.m";
Loading "gt.m"
Loading "witt.m"
Loading "etas.m"

> p:=5; d:=3; F<a>:=GF(p^d); P<x,y>:=PolynomialRing(F,2);
> v := [ x, x*y, P!2 ];
> w := [ x^2, x-y, 1+y^2];
> WittSum(v,w);
[
    x^2 + x,
    4*x^9 + 3*x^8 + 3*x^7 + 4*x^6 + x*y + x + 4*y,
    4*x^49 + 3*x^48 + 3*x^47 + 4*x^46 + 3*x^44 + 2*x^43 + x^41 + 4*x^38 + 4*x^37*y + 3*x^37 + 3*x^36*y + 2*x^36 + x^35*y + 3*x^35 + 
        3*x^34*y + 2*x^34 + 4*x^33*y + 3*x^32*y + x^31 + 3*x^30 + 2*x^29*y^2 + x^29*y + x^29 + 3*x^28*y^2 + x^28*y + x^28 + 
        4*x^27*y^2 + 2*x^27 + 2*x^26*y + x^26 + 2*x^25*y^2 + 4*x^25*y + 4*x^24*y^2 + x^24*y + x^24 + 4*x^23*y^2 + 3*x^23*y + 
        2*x^22*y^2 + 2*x^22*y + x^22 + 3*x^21*y^3 + 4*x^21*y^2 + x^21*y + 3*x^20*y^3 + 2*x^20*y^2 + 2*x^20*y + 4*x^20 + 2*x^19*y^3 + 
        2*x^19*y^2 + 2*x^19*y + 4*x^19 + x^18*y^3 + 4*x^18*y^2 + 3*x^18*y + 4*x^17*y^3 + 4*x^17*y^2 + 2*x^17*y + 4*x^17 + x^16*y^3 + 
        2*x^16*y^2 + 4*x^16*y + 2*x^16 + 4*x^15*y^3 + 4*x^15*y^2 + 3*x^15*y + 3*x^15 + 3*x^14*y^3 + 3*x^14*y^2 + x^14*y + x^13*y^4 + 
        x^13*y^3 + 4*x^13*y + x^13 + 3*x^12*y^4 + 3*x^12*y^3 + 4*x^12*y + 2*x^12 + x^11*y^3 + 4*x^11*y^2 + 2*x^11 + x^10*y^4 + 
        4*x^10*y^2 + x^10*y + x^10 + x^9*y^4 + 4*x^9*y^3 + x^9*y + 4*x^8*y^3 + x^8*y^2 + 3*x^7*y^4 + x^7*y^3 + x^6*y^4 + 4*x^5*y^4 + 
        3*x^5*y^3 + 3*x^5*y^2 + 4*x^5*y + x^4*y^5 + 4*x^4*y^4 + x^4*y^3 + 4*x^4*y^2 + 3*x^3*y^5 + 4*x^3*y^4 + 4*x^3*y^3 + 2*x^2*y^5 +
        4*x^2*y^4 + 4*x*y^5 + y^2 + 3
]

> WittProd(v,w);
[
    x^3,
    x^11*y + x^6 + 4*x^5*y,
    4*x^50*y^4 + 2*x^50 + x^49*y^5 + 3*x^45*y^3 + 4*x^44*y^4 + 3*x^43*y^5 + 3*x^40*y^2 + x^39*y^3 + 4*x^38*y^4 + 2*x^37*y^5 + 
        4*x^35*y + 4*x^34*y^2 + 4*x^33*y^3 + 4*x^32*y^4 + 4*x^31*y^5 + x^25*y^2 + x^25 + x^10*y^5 + 4*x^5*y^10
]

> WittNeg(v); 
[
    4*x,
    4*x*y,
    3
]

> WittPower(v,5);
[
    x^5,
    0,
    x^105*y^5
]

```
Note that we cannot invert `v` (or `w`) since they are not units.
```
> WittInv(v);

WittInv(
    v: [ x, x*y, 2 ]
)
WittInv1(
    v: [ x, x*y, 2 ]
)
In file "witt.m", line 325, column 19:
>>         w2 := [ PR!res[j] : j in [1..i] ] cat [x];
                     ^
Runtime error in '!': Illegal coercion
```
But:
```
> WittInv([P!2,x,y]);
[
    3,
    x,
    2*x^10 + y
]
```

As explained above, every time these functions are called, the eta
polynomials are computed, wasting time.  If one needs to make many
computations for longer length or large characteristic, it is better
to save these and reuse them.  For example, continuing the example
above, we could save the eta polynomials (for prime `p`, defined to be
5 in the example, and length 3) with:
```
epols:=etapols(p,2);
```
(Note again that we use `2` as the second argument to give length of
3.)  Then, when we need to compute with Witt vectors of length 3 with
entries in characteristic `p`, we can pass these polynomials as an
optional argument `pols`.  The routines then use those, instead of
computing them each time:
```
> WittSum(v,w : pols:=epols);
[
    x^2 + x,
    4*x^9 + 3*x^8 + 3*x^7 + 4*x^6 + x*y + x + 4*y,
    4*x^49 + 3*x^48 + 3*x^47 + 4*x^46 + 3*x^44 + 2*x^43 + x^41 + 4*x^38 + 
        4*x^37*y + 3*x^37 + 3*x^36*y + 2*x^36 + x^35*y + 3*x^35 + 3*x^34*y + 
        2*x^34 + 4*x^33*y + 3*x^32*y + x^31 + 3*x^30 + 2*x^29*y^2 + x^29*y + 
        x^29 + 3*x^28*y^2 + x^28*y + x^28 + 4*x^27*y^2 + 2*x^27 + 2*x^26*y + 
        x^26 + 2*x^25*y^2 + 4*x^25*y + 4*x^24*y^2 + x^24*y + x^24 + 4*x^23*y^2 +
        3*x^23*y + 2*x^22*y^2 + 2*x^22*y + x^22 + 3*x^21*y^3 + 4*x^21*y^2 + 
        x^21*y + 3*x^20*y^3 + 2*x^20*y^2 + 2*x^20*y + 4*x^20 + 2*x^19*y^3 + 
        2*x^19*y^2 + 2*x^19*y + 4*x^19 + x^18*y^3 + 4*x^18*y^2 + 3*x^18*y + 
        4*x^17*y^3 + 4*x^17*y^2 + 2*x^17*y + 4*x^17 + x^16*y^3 + 2*x^16*y^2 + 
        4*x^16*y + 2*x^16 + 4*x^15*y^3 + 4*x^15*y^2 + 3*x^15*y + 3*x^15 + 
        3*x^14*y^3 + 3*x^14*y^2 + x^14*y + x^13*y^4 + x^13*y^3 + 4*x^13*y + x^13
        + 3*x^12*y^4 + 3*x^12*y^3 + 4*x^12*y + 2*x^12 + x^11*y^3 + 4*x^11*y^2 + 
        2*x^11 + x^10*y^4 + 4*x^10*y^2 + x^10*y + x^10 + x^9*y^4 + 4*x^9*y^3 + 
        x^9*y + 4*x^8*y^3 + x^8*y^2 + 3*x^7*y^4 + x^7*y^3 + x^6*y^4 + 4*x^5*y^4 
        + 3*x^5*y^3 + 3*x^5*y^2 + 4*x^5*y + x^4*y^5 + 4*x^4*y^4 + x^4*y^3 + 
        4*x^4*y^2 + 3*x^3*y^5 + 4*x^3*y^4 + 4*x^3*y^3 + 2*x^2*y^5 + 4*x^2*y^4 + 
        4*x*y^5 + y^2 + 3
]

> WittProd(v,w : pols:=epols);
[
    x^3,
    x^11*y + x^6 + 4*x^5*y,
    4*x^50*y^4 + 2*x^50 + x^49*y^5 + 3*x^45*y^3 + 4*x^44*y^4 + 3*x^43*y^5 + 
        3*x^40*y^2 + x^39*y^3 + 4*x^38*y^4 + 2*x^37*y^5 + 4*x^35*y + 4*x^34*y^2 
        + 4*x^33*y^3 + 4*x^32*y^4 + 4*x^31*y^5 + x^25*y^2 + x^25 + x^10*y^5 + 
        4*x^5*y^10
]
```

Remember that by default, all these functions use method 1.  You can
specify a different method with the optional argument `choice`.  To
use method 2:
```
> WittProd(v,w : choice:=2);
[
    x^3,
    x^11*y + x^6 + 4*x^5*y,
    4*x^50*y^4 + 2*x^50 + x^49*y^5 + 3*x^45*y^3 + 4*x^44*y^4 + 3*x^43*y^5 + 
        3*x^40*y^2 + x^39*y^3 + 4*x^38*y^4 + 2*x^37*y^5 + 4*x^35*y + 4*x^34*y^2 
        + 4*x^33*y^3 + 4*x^32*y^4 + 4*x^31*y^5 + x^25*y^2 + x^25 + x^10*y^5 + 
        4*x^5*y^10
]
```
For method 3:
```
> WittProd(v,w : choice:=3);
[
    x^3,
    x^11*y + x^6 + 4*x^5*y,
    4*x^50*y^4 + 2*x^50 + x^49*y^5 + 3*x^45*y^3 + 4*x^44*y^4 + 3*x^43*y^5 + 
        3*x^40*y^2 + x^39*y^3 + 4*x^38*y^4 + 2*x^37*y^5 + 4*x^35*y + 4*x^34*y^2 
        + 4*x^33*y^3 + 4*x^32*y^4 + 4*x^31*y^5 + x^25*y^2 + x^25 + x^10*y^5 + 
        4*x^5*y^10
]
```

As observed above, methods 2 and 3 do not use the eta polynomials, but some
smaller data, which are in fact just some binomials.  To save time
then, we can precompute these data and reuse it.  For length 3 and
characteristic `p`, we can save the date with:
```
bt:=BinTab(p,2);
```
And then we can use it by using the optional argument `bintab`:
```
> WittProd(v,w : choice:=2, bintab:=bt);
[
    x^3,
    x^11*y + x^6 + 4*x^5*y,
    4*x^50*y^4 + 2*x^50 + x^49*y^5 + 3*x^45*y^3 + 4*x^44*y^4 + 3*x^43*y^5 + 
        3*x^40*y^2 + x^39*y^3 + 4*x^38*y^4 + 2*x^37*y^5 + 4*x^35*y + 4*x^34*y^2 
        + 4*x^33*y^3 + 4*x^32*y^4 + 4*x^31*y^5 + x^25*y^2 + x^25 + x^10*y^5 + 
        4*x^5*y^10
]

> WittProd(v,w : choice:=3, bintab:=bt);
[
    x^3,
    x^11*y + x^6 + 4*x^5*y,
    4*x^50*y^4 + 2*x^50 + x^49*y^5 + 3*x^45*y^3 + 4*x^44*y^4 + 3*x^43*y^5 + 
        3*x^40*y^2 + x^39*y^3 + 4*x^38*y^4 + 2*x^37*y^5 + 4*x^35*y + 4*x^34*y^2 
        + 4*x^33*y^3 + 4*x^32*y^4 + 4*x^31*y^5 + x^25*y^2 + x^25 + x^10*y^5 + 
        4*x^5*y^10
]
```

**Important:** If you are performing many computations with Witt
vectors is highly recommended that you store and reuse either these
eta polynomials (for method 1) or the binomials (for methods 2 and 3)
as above.  Otherwise, they will be computed every time you call
perform a computation, and for longer lengths or higher
characteristics, these alone can take some time.


### Conversions

The file `etas.m` also provide a few basic functions.

* `IntToWitt`: Converts an integer to a Witt vector of a given length
  and entries in a specified prime field.

* `WittVToSeries`: Converts a Witt vector over a finite field to a
  p-adic (truncated) power series in an unramified extension of the
  p-adic integers.

* `SeriesToWittV`: Converts a p-adic (truncated) power series in an
  unramified extension of the p-adic integers to a Witt vector.


Here are some examples.

Let's convert 1728 into a Witt vector of length 3 over the field
with 5 elements.
```
> IntToWitt(1728,5,2);
[ 3, 2, 0 ]
```
(Again, note we use second argument as `2` to give length 3.)

Converting from Witt vector to power series:
```
> F<a> := GF(5^3);
> s := WittVToSeries([a,F!2,a^7,a^100]); s;
-190*aa^2 - 219*aa - 305 + O(5^4)

> Parent(s);
Unramified extension defined by the polynomial (1 + O(5^4))*x^3 + O(5^4)*x^2 + 
    (3 + O(5^4))*x + 3 + O(5^4)
 over 5-adic ring

> F!s;
a
```


Now, converting from series to Witt vector:
```
> Z5 := pAdicRing(5 : Precision:=4);     
> Z125<aa> := ext<Z5 | 3>;
> Z125;
Unramified extension defined by the polynomial (1 + O(5^4))*x^3 + O(5^4)*x^2 + 
    (3 + O(5^4))*x + 3 + O(5^4)
 over 5-adic ring

> F!aa;
a

> x := Random(Z125); x;
-120*aa^2 + 59*aa - 241 + O(5^4)

> SeriesToWittV(x);
[ a^34, 4, a^40, a^65 ]
```

### Greenberg Transform


The file `gt.m` provides the following functions:

* `GT`: Computes the Greenberg transform of polynomials (over Witt
  vectors) in two variables.  In particular, it provides an efficient
  way to evaluate such polynomials.

* `Pol_GT_Form`: Converts a polynomial over the integers in two
  variables to the form that is needed for the input of `GT`.

The input for `GT` (a polynomial in two variables with coefficients in
the ring of Witt vectors) must be given in the form ` [ m1, m2, ... ]`
where each `mi` corresponds to a monomial and has the form `[* vi, ki,
li *]`, where `vi` is the coefficient of the monomial, and `ki` and `li` are the
powers of the first and second variable, respectively, for the
corresponding monomial.

Here are some examples:
```
> p:=3; d:=4; F<a>:=GF(p^d);
> pol := [ [* [a,F!2,a^5], 2,0 *], [* [a,F!0,F!1], 1, 1 *], [* [F!1,a^7,F!0], 0, 0 *] ];
> GT(pol);
[
    a*x0^2 + a*x0*y0 + 1,
    2*x0^6 + a^43*x0^5*y0 + a^43*x0^4*y0^2 + a^42*x0^4 + a^43*x0^3*x1 + a^2*x0^3*y0 + a^3*x0^3*y1 + a^42*x0^2*y0^2 + a^41*x0^2 + 
        a^41*x0*y0 + a^3*x1*y0^3 + a^7,
    a^5*x0^18 + a^38*x0^17*y0 + a^72*x0^16*y0^2 + a^37*x0^16 + a^3*x0^15*x1 + a^46*x0^15*y0^3 + a^31*x0^15*y0 + a^43*x0^15*y1 + 
        a^46*x0^14*x1*y0 + a^37*x0^14*y0^4 + a^2*x0^14*y0^2 + a^6*x0^14*y0*y1 + a^70*x0^14 + a^77*x0^13*x1*y0^2 + a^45*x0^13*x1 + 
        a^49*x0^13*y0^5 + a^48*x0^13*y0^3 + a^37*x0^13*y0^2*y1 + a*x0^13*y0 + a^5*x0^13*y1 + a^6*x0^12*x1^2 + a^57*x0^12*x1*y0^3 + 
        a^36*x0^12*x1*y0 + a^6*x0^12*x1*y1 + a^49*x0^12*y0^4 + a^9*x0^12*y0^3*y1 + a^4*x0^12*y0^2 + a^76*x0^12*y0*y1 + a^6*x0^12*y1^2
        + a^34*x0^12 + a^9*x0^11*x1^2*y0 + a^10*x0^11*x1*y0^4 + a^45*x0^11*x1*y0^2 + a^9*x0^11*x1*y0*y1 + a^75*x0^11*x1 + 
        a^49*x0^11*y0^7 + a^48*x0^11*y0^5 + a^49*x0^11*y0^4*y1 + a^47*x0^11*y0^3 + a^5*x0^11*y0^2*y1 + a^9*x0^11*y0*y1^2 + 
        a^9*x0^11*y0 + a^35*x0^11*y1 + a^9*x0^10*x1^2*y0^2 + a^8*x0^10*x1^2 + a^37*x0^10*x1*y0^5 + a^5*x0^10*x1*y0^3 + 
        a^9*x0^10*x1*y0^2*y1 + a^44*x0^10*x1*y0 + a^8*x0^10*x1*y1 + a^49*x0^10*y0^8 + a^35*x0^10*y0^4 + a^9*x0^10*y0^2*y1^2 + 
        a^52*x0^10*y0^2 + a^4*x0^10*y0*y1 + a^8*x0^10*y1^2 + a^16*x0^10 + a^12*x0^9*x1^3 + a^6*x0^9*x1^2*y0^3 + a^48*x0^9*x1^2*y0 + 
        a^49*x0^9*x1^2*y1 + a^9*x0^9*x1*y0^6 + a^45*x0^9*x1*y0^4 + a^46*x0^9*x1*y0^3*y1 + a^7*x0^9*x1*y0^2 + a^48*x0^9*x1*y0*y1 + 
        a^9*x0^9*x1*y1^2 + a^9*x0^9*x1 + a^49*x0^9*x2 + x0^9*y0^9 + a^8*x0^9*y0^7 + a^8*x0^9*y0^4*y1 + a^73*x0^9*y0^3 + 
        a^47*x0^9*y0^2*y1 + a^48*x0^9*y0*y1^2 + a^35*x0^9*y0 + a^49*x0^9*y1 + a^9*x0^9*y2 + a^9*x0^8*x1^2*y0^4 + a^8*x0^8*x1^2*y0^2 +
        a^7*x0^8*x1^2 + a^49*x0^8*x1*y0^7 + a^5*x0^8*x1*y0^5 + a^49*x0^8*x1*y0^4*y1 + a^35*x0^8*x1*y0^3 + a^8*x0^8*x1*y0^2*y1 + 
        a^13*x0^8*x1*y0 + a^7*x0^8*x1*y1 + a^48*x0^8*y0^8 + a^28*x0^8*y0^4 + a^8*x0^8*y0^2*y1^2 + a^16*x0^8*y0^2 + a^53*x0^8*y0*y1 + 
        a^7*x0^8*y1^2 + a^50*x0^8 + a^9*x0^7*x1^2*y0^5 + a^8*x0^7*x1^2*y0^3 + a^7*x0^7*x1^2*y0 + a^49*x0^7*x1*y0^5*y1 + 
        a^8*x0^7*x1*y0^4 + a^48*x0^7*x1*y0^3*y1 + a^13*x0^7*x1*y0^2 + a^7*x0^7*x1*y0*y1 + a^67*x0^7*x1 + a^47*x0^7*y0^7 + 
        a^46*x0^7*y0^5 + a^47*x0^7*y0^4*y1 + a^45*x0^7*y0^3 + a^53*x0^7*y0^2*y1 + a^7*x0^7*y0*y1^2 + a^7*x0^7*y0 + a^27*x0^7*y1 + 
        a^49*x0^6*x1^3*y0^3 + a^6*x0^6*x1^2*y0^6 + a^48*x0^6*x1^2*y0^4 + a^49*x0^6*x1^2*y0^3*y1 + a^53*x0^6*x1^2 + a^8*x0^6*x1*y0^7 +
        a^47*x0^6*x1*y0^5 + a^8*x0^6*x1*y0^4*y1 + a^49*x0^6*x1*y0^3*y1^2 + a^50*x0^6*x1*y0^3 + a^27*x0^6*x1*y0 + a^53*x0^6*x1*y1 + 
        a^29*x0^6*y0^4 + a^6*x0^6*y0^3*y1 + a^51*x0^6*y0^2 + a^67*x0^6*y0*y1 + a^53*x0^6*y1^2 + a*x0^6 + a^9*x0^5*x1^2*y0^7 + 
        a^8*x0^5*x1^2*y0^5 + a^7*x0^5*x1^2*y0^3 + a^48*x0^5*x1*y0^5*y1 + a^53*x0^5*x1*y0^4 + a^47*x0^5*x1*y0^3*y1 + a^67*x0^5*x1*y0^2
        + a^11*x0^5*x1 + a^45*x0^5*y0^5 + a^44*x0^5*y0^3 + a^27*x0^5*y0^2*y1 + a^17*x0^5*y0 + a^51*x0^5*y1 + a^9*x0^4*x1^2*y0^8 + 
        a^8*x0^4*x1^2*y0^6 + a^7*x0^4*x1^2*y0^4 + a^47*x0^4*x1*y0^7 + a^53*x0^4*x1*y0^5 + a^47*x0^4*x1*y0^4*y1 + a^27*x0^4*x1*y0^3 + 
        a^11*x0^4*x1*y0 + a^26*x0^4*y0^4 + a^17*x0^4*y0^2 + a^51*x0^4*y0*y1 + a^45*x0^4 + a^9*x0^3*x1^3*y0^6 + a^48*x0^3*x1^2*y0^7 + 
        a^49*x0^3*x1^2*y0^6*y1 + a^53*x0^3*x1^2*y0^3 + a^6*x0^3*x1*y0^6 + a^67*x0^3*x1*y0^4 + a^13*x0^3*x1*y0^3*y1 + a^17*x0^3*x1 + 
        a^10*x0^3*y0^3 + a^5*x0^3*y0 + a^57*x0^3*y1 + a^8*x0^2*x1^2*y0^8 + a^7*x0^2*x1^2*y0^6 + a^27*x0^2*x1*y0^5 + a^51*x0^2*x1*y0^3
        + a^45*x0^2*y0^2 + a^7*x0^2 + a^7*x0*x1^2*y0^7 + a^51*x0*x1*y0^4 + a^7*x0*y0 + a^9*x1^6 + a^9*x1^3*y1^3 + a^53*x1^2*y0^6 + 
        a^57*x1*y0^3 + a^9*x2*y0^9
]
```

As with the other functions, you can specify the method for the
computations with the optional argument `choice`.  For method 1
(the default), again as above, you can use the optional argument
`pols` to give the function precomputed eta polynomials.  For methods
2 or 3, you can use the optional argument `bintab` to pass along the
table of the used table of binomials used with those methods.

For example:
```
> epols := etapols(p,2);
> GT(pol : pols:=epols);
[
    a*x0^2 + a*x0*y0 + 1,
    2*x0^6 + a^43*x0^5*y0 + a^43*x0^4*y0^2 + a^42*x0^4 + a^43*x0^3*x1 + a^2*x0^3*y0 + a^3*x0^3*y1 + a^42*x0^2*y0^2 + a^41*x0^2 + 
        a^41*x0*y0 + a^3*x1*y0^3 + a^7,
    a^5*x0^18 + a^38*x0^17*y0 + a^72*x0^16*y0^2 + a^37*x0^16 + a^3*x0^15*x1 + a^46*x0^15*y0^3 + a^31*x0^15*y0 + a^43*x0^15*y1 + 
        a^46*x0^14*x1*y0 + a^37*x0^14*y0^4 + a^2*x0^14*y0^2 + a^6*x0^14*y0*y1 + a^70*x0^14 + a^77*x0^13*x1*y0^2 + a^45*x0^13*x1 + 
        a^49*x0^13*y0^5 + a^48*x0^13*y0^3 + a^37*x0^13*y0^2*y1 + a*x0^13*y0 + a^5*x0^13*y1 + a^6*x0^12*x1^2 + a^57*x0^12*x1*y0^3 + 
        a^36*x0^12*x1*y0 + a^6*x0^12*x1*y1 + a^49*x0^12*y0^4 + a^9*x0^12*y0^3*y1 + a^4*x0^12*y0^2 + a^76*x0^12*y0*y1 + a^6*x0^12*y1^2
        + a^34*x0^12 + a^9*x0^11*x1^2*y0 + a^10*x0^11*x1*y0^4 + a^45*x0^11*x1*y0^2 + a^9*x0^11*x1*y0*y1 + a^75*x0^11*x1 + 
        a^49*x0^11*y0^7 + a^48*x0^11*y0^5 + a^49*x0^11*y0^4*y1 + a^47*x0^11*y0^3 + a^5*x0^11*y0^2*y1 + a^9*x0^11*y0*y1^2 + 
        a^9*x0^11*y0 + a^35*x0^11*y1 + a^9*x0^10*x1^2*y0^2 + a^8*x0^10*x1^2 + a^37*x0^10*x1*y0^5 + a^5*x0^10*x1*y0^3 + 
        a^9*x0^10*x1*y0^2*y1 + a^44*x0^10*x1*y0 + a^8*x0^10*x1*y1 + a^49*x0^10*y0^8 + a^35*x0^10*y0^4 + a^9*x0^10*y0^2*y1^2 + 
        a^52*x0^10*y0^2 + a^4*x0^10*y0*y1 + a^8*x0^10*y1^2 + a^16*x0^10 + a^12*x0^9*x1^3 + a^6*x0^9*x1^2*y0^3 + a^48*x0^9*x1^2*y0 + 
        a^49*x0^9*x1^2*y1 + a^9*x0^9*x1*y0^6 + a^45*x0^9*x1*y0^4 + a^46*x0^9*x1*y0^3*y1 + a^7*x0^9*x1*y0^2 + a^48*x0^9*x1*y0*y1 + 
        a^9*x0^9*x1*y1^2 + a^9*x0^9*x1 + a^49*x0^9*x2 + x0^9*y0^9 + a^8*x0^9*y0^7 + a^8*x0^9*y0^4*y1 + a^73*x0^9*y0^3 + 
        a^47*x0^9*y0^2*y1 + a^48*x0^9*y0*y1^2 + a^35*x0^9*y0 + a^49*x0^9*y1 + a^9*x0^9*y2 + a^9*x0^8*x1^2*y0^4 + a^8*x0^8*x1^2*y0^2 +
        a^7*x0^8*x1^2 + a^49*x0^8*x1*y0^7 + a^5*x0^8*x1*y0^5 + a^49*x0^8*x1*y0^4*y1 + a^35*x0^8*x1*y0^3 + a^8*x0^8*x1*y0^2*y1 + 
        a^13*x0^8*x1*y0 + a^7*x0^8*x1*y1 + a^48*x0^8*y0^8 + a^28*x0^8*y0^4 + a^8*x0^8*y0^2*y1^2 + a^16*x0^8*y0^2 + a^53*x0^8*y0*y1 + 
        a^7*x0^8*y1^2 + a^50*x0^8 + a^9*x0^7*x1^2*y0^5 + a^8*x0^7*x1^2*y0^3 + a^7*x0^7*x1^2*y0 + a^49*x0^7*x1*y0^5*y1 + 
        a^8*x0^7*x1*y0^4 + a^48*x0^7*x1*y0^3*y1 + a^13*x0^7*x1*y0^2 + a^7*x0^7*x1*y0*y1 + a^67*x0^7*x1 + a^47*x0^7*y0^7 + 
        a^46*x0^7*y0^5 + a^47*x0^7*y0^4*y1 + a^45*x0^7*y0^3 + a^53*x0^7*y0^2*y1 + a^7*x0^7*y0*y1^2 + a^7*x0^7*y0 + a^27*x0^7*y1 + 
        a^49*x0^6*x1^3*y0^3 + a^6*x0^6*x1^2*y0^6 + a^48*x0^6*x1^2*y0^4 + a^49*x0^6*x1^2*y0^3*y1 + a^53*x0^6*x1^2 + a^8*x0^6*x1*y0^7 +
        a^47*x0^6*x1*y0^5 + a^8*x0^6*x1*y0^4*y1 + a^49*x0^6*x1*y0^3*y1^2 + a^50*x0^6*x1*y0^3 + a^27*x0^6*x1*y0 + a^53*x0^6*x1*y1 + 
        a^29*x0^6*y0^4 + a^6*x0^6*y0^3*y1 + a^51*x0^6*y0^2 + a^67*x0^6*y0*y1 + a^53*x0^6*y1^2 + a*x0^6 + a^9*x0^5*x1^2*y0^7 + 
        a^8*x0^5*x1^2*y0^5 + a^7*x0^5*x1^2*y0^3 + a^48*x0^5*x1*y0^5*y1 + a^53*x0^5*x1*y0^4 + a^47*x0^5*x1*y0^3*y1 + a^67*x0^5*x1*y0^2
        + a^11*x0^5*x1 + a^45*x0^5*y0^5 + a^44*x0^5*y0^3 + a^27*x0^5*y0^2*y1 + a^17*x0^5*y0 + a^51*x0^5*y1 + a^9*x0^4*x1^2*y0^8 + 
        a^8*x0^4*x1^2*y0^6 + a^7*x0^4*x1^2*y0^4 + a^47*x0^4*x1*y0^7 + a^53*x0^4*x1*y0^5 + a^47*x0^4*x1*y0^4*y1 + a^27*x0^4*x1*y0^3 + 
        a^11*x0^4*x1*y0 + a^26*x0^4*y0^4 + a^17*x0^4*y0^2 + a^51*x0^4*y0*y1 + a^45*x0^4 + a^9*x0^3*x1^3*y0^6 + a^48*x0^3*x1^2*y0^7 + 
        a^49*x0^3*x1^2*y0^6*y1 + a^53*x0^3*x1^2*y0^3 + a^6*x0^3*x1*y0^6 + a^67*x0^3*x1*y0^4 + a^13*x0^3*x1*y0^3*y1 + a^17*x0^3*x1 + 
        a^10*x0^3*y0^3 + a^5*x0^3*y0 + a^57*x0^3*y1 + a^8*x0^2*x1^2*y0^8 + a^7*x0^2*x1^2*y0^6 + a^27*x0^2*x1*y0^5 + a^51*x0^2*x1*y0^3
        + a^45*x0^2*y0^2 + a^7*x0^2 + a^7*x0*x1^2*y0^7 + a^51*x0*x1*y0^4 + a^7*x0*y0 + a^9*x1^6 + a^9*x1^3*y1^3 + a^53*x1^2*y0^6 + 
        a^57*x1*y0^3 + a^9*x2*y0^9
]

> bt := BinTab(p,2);
> GT(pol : choice :=2, bintab:=bt);
[
    a*x0^2 + a*x0*y0 + 1,
    2*x0^6 + a^43*x0^5*y0 + a^43*x0^4*y0^2 + a^42*x0^4 + a^43*x0^3*x1 + a^2*x0^3*y0 + a^3*x0^3*y1 + a^42*x0^2*y0^2 + a^41*x0^2 + 
        a^41*x0*y0 + a^3*x1*y0^3 + a^7,
    a^5*x0^18 + a^38*x0^17*y0 + a^72*x0^16*y0^2 + a^37*x0^16 + a^3*x0^15*x1 + a^46*x0^15*y0^3 + a^31*x0^15*y0 + a^43*x0^15*y1 + 
        a^46*x0^14*x1*y0 + a^37*x0^14*y0^4 + a^2*x0^14*y0^2 + a^6*x0^14*y0*y1 + a^70*x0^14 + a^77*x0^13*x1*y0^2 + a^45*x0^13*x1 + 
        a^49*x0^13*y0^5 + a^48*x0^13*y0^3 + a^37*x0^13*y0^2*y1 + a*x0^13*y0 + a^5*x0^13*y1 + a^6*x0^12*x1^2 + a^57*x0^12*x1*y0^3 + 
        a^36*x0^12*x1*y0 + a^6*x0^12*x1*y1 + a^49*x0^12*y0^4 + a^9*x0^12*y0^3*y1 + a^4*x0^12*y0^2 + a^76*x0^12*y0*y1 + a^6*x0^12*y1^2
        + a^34*x0^12 + a^9*x0^11*x1^2*y0 + a^10*x0^11*x1*y0^4 + a^45*x0^11*x1*y0^2 + a^9*x0^11*x1*y0*y1 + a^75*x0^11*x1 + 
        a^49*x0^11*y0^7 + a^48*x0^11*y0^5 + a^49*x0^11*y0^4*y1 + a^47*x0^11*y0^3 + a^5*x0^11*y0^2*y1 + a^9*x0^11*y0*y1^2 + 
        a^9*x0^11*y0 + a^35*x0^11*y1 + a^9*x0^10*x1^2*y0^2 + a^8*x0^10*x1^2 + a^37*x0^10*x1*y0^5 + a^5*x0^10*x1*y0^3 + 
        a^9*x0^10*x1*y0^2*y1 + a^44*x0^10*x1*y0 + a^8*x0^10*x1*y1 + a^49*x0^10*y0^8 + a^35*x0^10*y0^4 + a^9*x0^10*y0^2*y1^2 + 
        a^52*x0^10*y0^2 + a^4*x0^10*y0*y1 + a^8*x0^10*y1^2 + a^16*x0^10 + a^12*x0^9*x1^3 + a^6*x0^9*x1^2*y0^3 + a^48*x0^9*x1^2*y0 + 
        a^49*x0^9*x1^2*y1 + a^9*x0^9*x1*y0^6 + a^45*x0^9*x1*y0^4 + a^46*x0^9*x1*y0^3*y1 + a^7*x0^9*x1*y0^2 + a^48*x0^9*x1*y0*y1 + 
        a^9*x0^9*x1*y1^2 + a^9*x0^9*x1 + a^49*x0^9*x2 + x0^9*y0^9 + a^8*x0^9*y0^7 + a^8*x0^9*y0^4*y1 + a^73*x0^9*y0^3 + 
        a^47*x0^9*y0^2*y1 + a^48*x0^9*y0*y1^2 + a^35*x0^9*y0 + a^49*x0^9*y1 + a^9*x0^9*y2 + a^9*x0^8*x1^2*y0^4 + a^8*x0^8*x1^2*y0^2 +
        a^7*x0^8*x1^2 + a^49*x0^8*x1*y0^7 + a^5*x0^8*x1*y0^5 + a^49*x0^8*x1*y0^4*y1 + a^35*x0^8*x1*y0^3 + a^8*x0^8*x1*y0^2*y1 + 
        a^13*x0^8*x1*y0 + a^7*x0^8*x1*y1 + a^48*x0^8*y0^8 + a^28*x0^8*y0^4 + a^8*x0^8*y0^2*y1^2 + a^16*x0^8*y0^2 + a^53*x0^8*y0*y1 + 
        a^7*x0^8*y1^2 + a^50*x0^8 + a^9*x0^7*x1^2*y0^5 + a^8*x0^7*x1^2*y0^3 + a^7*x0^7*x1^2*y0 + a^49*x0^7*x1*y0^5*y1 + 
        a^8*x0^7*x1*y0^4 + a^48*x0^7*x1*y0^3*y1 + a^13*x0^7*x1*y0^2 + a^7*x0^7*x1*y0*y1 + a^67*x0^7*x1 + a^47*x0^7*y0^7 + 
        a^46*x0^7*y0^5 + a^47*x0^7*y0^4*y1 + a^45*x0^7*y0^3 + a^53*x0^7*y0^2*y1 + a^7*x0^7*y0*y1^2 + a^7*x0^7*y0 + a^27*x0^7*y1 + 
        a^49*x0^6*x1^3*y0^3 + a^6*x0^6*x1^2*y0^6 + a^48*x0^6*x1^2*y0^4 + a^49*x0^6*x1^2*y0^3*y1 + a^53*x0^6*x1^2 + a^8*x0^6*x1*y0^7 +
        a^47*x0^6*x1*y0^5 + a^8*x0^6*x1*y0^4*y1 + a^49*x0^6*x1*y0^3*y1^2 + a^50*x0^6*x1*y0^3 + a^27*x0^6*x1*y0 + a^53*x0^6*x1*y1 + 
        a^29*x0^6*y0^4 + a^6*x0^6*y0^3*y1 + a^51*x0^6*y0^2 + a^67*x0^6*y0*y1 + a^53*x0^6*y1^2 + a*x0^6 + a^9*x0^5*x1^2*y0^7 + 
        a^8*x0^5*x1^2*y0^5 + a^7*x0^5*x1^2*y0^3 + a^48*x0^5*x1*y0^5*y1 + a^53*x0^5*x1*y0^4 + a^47*x0^5*x1*y0^3*y1 + a^67*x0^5*x1*y0^2
        + a^11*x0^5*x1 + a^45*x0^5*y0^5 + a^44*x0^5*y0^3 + a^27*x0^5*y0^2*y1 + a^17*x0^5*y0 + a^51*x0^5*y1 + a^9*x0^4*x1^2*y0^8 + 
        a^8*x0^4*x1^2*y0^6 + a^7*x0^4*x1^2*y0^4 + a^47*x0^4*x1*y0^7 + a^53*x0^4*x1*y0^5 + a^47*x0^4*x1*y0^4*y1 + a^27*x0^4*x1*y0^3 + 
        a^11*x0^4*x1*y0 + a^26*x0^4*y0^4 + a^17*x0^4*y0^2 + a^51*x0^4*y0*y1 + a^45*x0^4 + a^9*x0^3*x1^3*y0^6 + a^48*x0^3*x1^2*y0^7 + 
        a^49*x0^3*x1^2*y0^6*y1 + a^53*x0^3*x1^2*y0^3 + a^6*x0^3*x1*y0^6 + a^67*x0^3*x1*y0^4 + a^13*x0^3*x1*y0^3*y1 + a^17*x0^3*x1 + 
        a^10*x0^3*y0^3 + a^5*x0^3*y0 + a^57*x0^3*y1 + a^8*x0^2*x1^2*y0^8 + a^7*x0^2*x1^2*y0^6 + a^27*x0^2*x1*y0^5 + a^51*x0^2*x1*y0^3
        + a^45*x0^2*y0^2 + a^7*x0^2 + a^7*x0*x1^2*y0^7 + a^51*x0*x1*y0^4 + a^7*x0*y0 + a^9*x1^6 + a^9*x1^3*y1^3 + a^53*x1^2*y0^6 + 
        a^57*x1*y0^3 + a^9*x2*y0^9
]
```

One can also use the function `GT` to evaluate polynomials in two
variables using the optional argument `vvars`.  For instance, if given
Witt vectors `v` and `w` we want to compute, say, `v^2*w + 2*v*w+
w^3`, one can do:
```
> P<x>:=PolynomialRing(F);
> v:=[ a^29, x + a^11, a^68 ];
> w:=[ a^25, a^48, a^9*x^2 ]; 


> vone := [ F!1, F!0, F!0 ];
> vtwo := IntToWitt(2,p,2); vtwo;
[ 2, 1, 0 ]

> pol := [ [* vone, 2, 1 *], [* vtwo, 1, 1 *], [* vone, 0, 3 *] ];
> GT(pol : pols:=epols, vvars:=v cat w);
[
    a^34,
    a^17*x + a^67,
    a^65*x^6 + a^30*x^3 + a^61*x^2 + a^16*x + 1
]
```

Note that we could also have used `Pol_GT_Form` to produce the `pol`
above:
```
> PP<X,Y>:=PolynomialRing(Integers(),2);
> zpol := X^2*Y + 2*X*Y + Y^3;
> Pol_GT_Form(zpol,p,2);
[ [*
    [ 1, 0, 0 ],
    2, 1
*], [*
    [ 2, 1, 0 ],
    1, 1
*], [*
    [ 1, 0, 0 ],
    0, 3
*] ]
```

### Canonical and Minimal Degree Liftings


The file `lift.m` has routines to compute the canonical and minimal
degree liftings of ordinary elliptic curves, along with their
corresponding lifting of points: the elliptic Teichmüller and minimal
degree respectively.  The curves are given by their Weierstrass
coefficients.

The canonical liftings can be computed up to length 3, while minimal
degree liftings can have any length.

The output is four vectors, say, `va`, `vb`, `vF`, and `vH`.  `va` and
`vb` are the Weierstrass coefficients of the canonical lifting, while
the Teichmüller/minimal degree lift is given by `(vF, y0*vH)`.  (See
[Weierstrass Coefficients of the Canonical Lifting][wcoef] for more
information, including the description of the algorithm.)

[wcoef]: http://www.math.utk.edu/~finotti/papers/wcoef.pdf

For example, to find some generic formulas for the canonical lifting
of `y0^2 = x0^3 + a0*x0 + b0` in characteristic 5 for length 2:
```
> load 'lift.m';
Loading "lift.m"
Loading "gt.m"
Loading "witt.m"
Loading "etas.m"

> p:=5; F<a0,b0>:=RationalFunctionField(GF(p),2);
> lift(a0,b0,1);
[
    a0,
    (a0^3*b0^2 + b0^4)/a0
]
[
    b0,
    4*a0^6*b0 + a0^3*b0^3 + b0^5
]
[
    x0,
    4/a0*x0^7 + 4*b0/a0*x0^4 + a0*x0^3 + 3*b0*x0^2 + 3*b0^2/a0*x0 + a0*b0
]
[
    1,
    1/a0*x0^8 + 2*x0^6 + 3*b0/a0*x0^5 + 2*a0*x0^4 + 3*b0*x0^3 + a0^2*x0^2 + 3*a0*b0*x0 + 3*a0^3 + 3*b0^2
]
>
```

For minimal degree, we set the optional argument `minimal:=true`.
```
> lift(a0,b0,2 : minimal:=true);
[
    a0,
    (a0^3*b0^2 + b0^4)/a0,
    (4*a0^36 + a0^33*b0^2 + a0^30*b0^4 + 3*a0^27*b0^6 + 2*a0^24*b0^8 + a0^18*b0^12 + 4*a0^12*b0^16 + 3*a0^9*b0^18 + 4*a0^6*b0^20 + 
        4*a0^3*b0^22 + 4*b0^24)/a0^11
]
[
    b0,
    4*a0^6*b0 + a0^3*b0^3 + b0^5,
    a0^36*b0 + 4*a0^33*b0^3 + 3*a0^27*b0^7 + 4*a0^21*b0^11 + 4*a0^15*b0^15 + a0^12*b0^17 + 3*a0^6*b0^21 + 4*b0^25
]
[
    x0,
    4/a0*x0^7 + 4*b0/a0*x0^4 + a0*x0^3 + 3*b0*x0^2 + 3*b0^2/a0*x0 + a0*b0,
    2/a0^6*x0^37 + (a0^3*b0^2 + b0^4)/a0^11*x0^35 + 2*b0/a0^6*x0^34 + 3/a0^4*x0^33 + (4*a0^3 + 4*b0^2)/a0^6*x0^31 + 
        3*b0^5/a0^10*x0^30 + 4*b0^2/a0^5*x0^29 + (a0^3 + 3*b0^2)/a0^4*x0^27 + 2*b0/a0^2*x0^26 + (4*a0^3*b0 + 3*b0^3)/a0^4*x0^24 + 
        (4*a0^6 + 3*a0^3*b0^2 + b0^4)/a0^5*x0^23 + (2*a0^3*b0^3 + 4*b0^5)/a0^6*x0^22 + (a0^3*b0^2 + b0^4)/a0^4*x0^21 + (2*a0^12*b0 + 
        4*a0^9*b0^3 + a0^3*b0^7 + b0^9)/a0^11*x0^20 + (a0^9 + 3*a0^6*b0^2 + 4*b0^6)/a0^6*x0^19 + (2*a0^6*b0 + a0^3*b0^3 + 
        4*b0^5)/a0^4*x0^18 + (2*a0^6 + a0^3*b0^2 + b0^4)/a0^2*x0^17 + (3*a0^6*b0^3 + 4*a0^3*b0^5 + 3*b0^7)/a0^6*x0^16 + (2*a0^15 + 
        4*a0^12*b0^2 + 3*a0^6*b0^6 + 4*b0^10)/a0^10*x0^15 + (4*a0^6*b0 + 2*b0^5)/a0^2*x0^14 + (3*a0^9 + 4*a0^6*b0^2 + 3*a0^3*b0^4 + 
        2*b0^6)/a0^3*x0^13 + (3*a0^9*b0 + 2*a0^6*b0^3 + a0^3*b0^5 + 2*b0^7)/a0^4*x0^12 + (2*a0^9*b0^2 + 3*a0^6*b0^4 + 4*a0^3*b0^6 + 
        2*b0^8)/a0^5*x0^11 + (2*a0^6*b0^3 + 3*a0^3*b0^5 + 3*b0^7)/a0^3*x0^10 + (3*a0^9*b0^2 + 2*a0^6*b0^4 + 4*a0^3*b0^6 + 
        4*b0^8)/a0^4*x0^9 + (2*a0^9*b0^3 + 2*a0^6*b0^5 + 4*a0^3*b0^7 + b0^9)/a0^5*x0^8 + (3*a0^9*b0^4 + 3*a0^6*b0^6 + a0^3*b0^8 + 
        2*b0^10)/a0^6*x0^7 + (3*a0^9*b0^3 + 3*a0^3*b0^7 + 2*b0^9)/a0^4*x0^6 + (4*a0^21 + 3*a0^15*b0^4 + a0^12*b0^6 + 2*a0^9*b0^8 + 
        2*a0^3*b0^12 + 2*b0^14)/a0^11*x0^5 + (a0^9*b0^5 + 3*a0^6*b0^7 + 4*a0^3*b0^9 + 2*b0^11)/a0^6*x0^4 + (4*a0^9*b0^4 + 4*a0^6*b0^6
        + a0^3*b0^8 + 3*b0^10)/a0^4*x0^3 + (2*a0^9*b0^5 + 2*a0^6*b0^7 + 4*b0^11)/a0^5*x0^2 + (2*a0^9*b0^6 + 4*b0^12)/a0^6*x0 + 
        (4*a0^21*b0 + a0^18*b0^3 + 4*a0^15*b0^5 + a0^12*b0^7 + a0^9*b0^9 + 3*a0^6*b0^11 + b0^15)/a0^10
]
[
    1,
    1/a0*x0^8 + 2*x0^6 + 3*b0/a0*x0^5 + 2*a0*x0^4 + 3*b0*x0^3 + a0^2*x0^2 + 3*a0*b0*x0 + 3*a0^3 + 3*b0^2,
    1/a0^10*x0^56 + 2/a0^9*x0^54 + 2*b0/a0^10*x0^53 + 1/a0^8*x0^52 + 2*b0/a0^9*x0^51 + b0^2/a0^10*x0^50 + 3/a0^6*x0^48 + (4*a0^3*b0^2
        + 4*b0^4)/a0^11*x0^46 + 4*b0/a0^6*x0^45 + (3*a0^6 + 3*a0^3*b0^2 + 3*b0^4)/a0^10*x0^44 + (3*a0^3*b0^3 + 3*b0^5)/a0^11*x0^43 + 
        (2*a0^6 + 4*a0^3*b0^2 + 4*b0^4)/a0^9*x0^42 + (3*a0^6*b0 + 3*a0^3*b0^3 + 4*b0^5)/a0^10*x0^41 + (a0^9 + 3*a0^6*b0^2 + 
        4*a0^3*b0^4 + 4*b0^6)/a0^11*x0^40 + (2*a0^6*b0 + 2*b0^5)/a0^9*x0^39 + (2*a0^9 + a0^6*b0^2 + 2*b0^6)/a0^10*x0^38 + (2*a0^6*b0 
        + 4*a0^3*b0^3 + b0^5)/a0^8*x0^37 + (4*a0^9 + a0^6*b0^2 + 4*a0^3*b0^4 + 2*b0^6)/a0^9*x0^36 + (4*a0^6*b0^3 + b0^7)/a0^10*x0^35 
        + (3*a0^6 + 3*b0^4)/a0^5*x0^34 + (3*a0^6*b0 + 4*a0^3*b0^3 + 3*b0^5)/a0^6*x0^33 + (a0^6 + 2*b0^4)/a0^4*x0^32 + (4*a0^12*b0 + 
        a0^9*b0^3 + 3*a0^6*b0^5 + 2*a0^3*b0^7 + 2*b0^9)/a0^11*x0^31 + (4*a0^9 + 4*a0^6*b0^2 + 4*b0^6)/a0^6*x0^30 + (2*a0^12*b0 + 
        2*a0^9*b0^3 + 3*a0^6*b0^5 + 4*a0^3*b0^7 + 4*b0^9)/a0^10*x0^29 + (2*a0^12*b0^2 + a0^9*b0^4 + 4*a0^6*b0^6 + 4*a0^3*b0^8 + 
        4*b0^10)/a0^11*x0^28 + (4*a0^12*b0 + a0^9*b0^3 + 4*a0^6*b0^5 + 2*a0^3*b0^7 + 2*b0^9)/a0^9*x0^27 + (2*a0^15 + 2*a0^12*b0^2 + 
        4*a0^3*b0^8 + 4*b0^10)/a0^10*x0^26 + (4*a0^15*b0 + a0^12*b0^3 + 2*a0^9*b0^5 + 2*a0^6*b0^7 + 2*a0^3*b0^9 + 
        2*b0^11)/a0^11*x0^25 + (4*a0^3*b0^2 + 2*b0^4)*x0^24 + (4*a0^9*b0 + 3*a0^6*b0^3 + 2*b0^7)/a0^4*x0^23 + (a0^9 + 3*a0^6*b0^2 + 
        4*a0^3*b0^4 + 2*b0^6)/a0^2*x0^22 + (3*a0^9*b0 + 3*a0^3*b0^5 + b0^7)/a0^3*x0^21 + (a0^9*b0^2 + 4*a0^3*b0^6 + 
        2*b0^8)/a0^4*x0^20 + (4*a0^7*b0 + a0^4*b0^3 + a0*b0^5)*x0^19 + (3*a0^6*b0^2 + 3*a0^3*b0^4 + 3*b0^6)*x0^18 + (3*a0^9*b0 + 
        4*a0^6*b0^3 + a0^3*b0^5 + 3*b0^7)/a0*x0^17 + (4*a0^15 + 3*a0^12*b0^2 + a0^6*b0^6 + 2*b0^10)/a0^5*x0^16 + (a0^12*b0 + 
        2*a0^9*b0^3 + 4*a0^6*b0^5 + 2*a0^3*b0^7 + 4*b0^9)/a0^3*x0^15 + (2*a0^15 + a0^9*b0^4 + 2*a0^3*b0^8 + 4*b0^10)/a0^4*x0^14 + 
        (3*a0^12*b0^3 + 3*a0^9*b0^5 + a0^6*b0^7 + a0^3*b0^9 + 4*b0^11)/a0^5*x0^13 + (3*a0^15 + 4*a0^9*b0^4 + 2*a0^6*b0^6 + 
        2*a0^3*b0^8 + 2*b0^10)/a0^3*x0^12 + (3*a0^21*b0 + 2*a0^18*b0^3 + 4*a0^15*b0^5 + a0^12*b0^7 + a0^9*b0^9 + 4*a0^6*b0^11 + 
        4*b0^15)/a0^10*x0^11 + (a0^18 + 3*a0^15*b0^2 + 3*a0^9*b0^6 + 2*a0^6*b0^8 + 4*a0^3*b0^10 + 2*b0^12)/a0^5*x0^10 + (a0^21*b0 + 
        4*a0^18*b0^3 + a0^15*b0^5 + 4*a0^12*b0^7 + 4*a0^9*b0^9 + 3*b0^15)/a0^9*x0^9 + (4*a0^24 + a0^21*b0^2 + 2*a0^18*b0^4 + 
        2*a0^15*b0^6 + 4*a0^12*b0^8 + a0^9*b0^10 + 3*b0^16)/a0^10*x0^8 + (a0^21*b0 + a0^18*b0^3 + 4*a0^15*b0^5 + 3*a0^12*b0^7 + 
        2*a0^9*b0^9 + 4*b0^15)/a0^8*x0^7 + (3*a0^24 + 2*a0^21*b0^2 + 3*a0^18*b0^4 + 4*a0^15*b0^6 + 4*a0^12*b0^8 + 2*a0^9*b0^10 + 
        3*b0^16)/a0^9*x0^6 + (2*a0^24*b0 + 2*a0^21*b0^3 + a0^18*b0^5 + 2*a0^15*b0^7 + 2*a0^9*b0^11 + 4*b0^17)/a0^10*x0^5 + 
        (a0^13*b0^2 + 4*a0^4*b0^8 + 2*a0*b0^10)*x0^4 + (4*a0^15*b0 + 3*a0^12*b0^3 + 4*a0^6*b0^7 + 2*a0^3*b0^9 + 4*b0^11)*x0^3 + 
        (3*a0^17 + 4*a0^14*b0^2 + 4*a0^11*b0^4 + a0^8*b0^6 + 4*a0^5*b0^8 + a0^2*b0^10)*x0^2 + (4*a0^16*b0 + 4*a0^13*b0^3 + 
        4*a0^10*b0^5 + 4*a0^7*b0^7 + 3*a0^4*b0^9 + 4*a0*b0^11)*x0 + 4*a0^18 + 4*a0^15*b0^2 + a0^12*b0^4 + 3*a0^9*b0^6 + a0^6*b0^8 + 
        2*b0^12
]
```
(Note that for length 2 the minimal and canonical liftings coincide.)

If also allows us to pass the eta polynomials to be used in the
computations using the optional argument `pols` as with the functions
above.

For minimal degree liftings, we can get higher lengths:
```
> epols:=etapols(5,3);

> lift(GF(5)!1,GF(5)!2,3 : minimal:=true, pols:=epols);
[ 1, 0, 4, 0 ]
[ 2, 3, 3, 0 ]
[
    x0,
    4*x0^7 + 3*x0^4 + x0^3 + x0^2 + 2*x0 + 2,
    2*x0^37 + 4*x0^34 + 3*x0^33 + x0^30 + x0^29 + 3*x0^27 + 4*x0^26 + 2*x0^24 + 2*x0^23 + 4*x0^22 + x0^20 + 4*x0^19 + 2*x0^17 + x0^16
        + x0^15 + 2*x0^14 + 4*x0^11 + x0^10 + 4*x0^9 + 4*x0^8 + 4*x0^7 + 2*x0^6 + 3*x0^5 + 3*x0^3 + 2*x0^2 + 2*x0 + 1,
    x0^187 + 3*x0^185 + 2*x0^184 + 4*x0^183 + 3*x0^180 + 4*x0^179 + 2*x0^177 + 3*x0^176 + x0^175 + 2*x0^174 + 4*x0^173 + 3*x0^172 + 
        4*x0^171 + 3*x0^170 + 2*x0^169 + 3*x0^167 + x0^166 + 2*x0^165 + 2*x0^164 + 3*x0^163 + 2*x0^162 + x0^161 + 3*x0^159 + x0^158 +
        2*x0^156 + 3*x0^155 + 2*x0^154 + 2*x0^153 + 2*x0^152 + 3*x0^149 + 4*x0^147 + 3*x0^143 + x0^142 + 2*x0^141 + 4*x0^139 + 
        4*x0^137 + x0^136 + 3*x0^134 + 2*x0^133 + 3*x0^132 + 4*x0^131 + 2*x0^130 + 4*x0^129 + 2*x0^128 + 2*x0^126 + 2*x0^123 + 
        4*x0^122 + 4*x0^121 + x0^119 + x0^118 + 4*x0^117 + 4*x0^116 + 3*x0^115 + 3*x0^113 + 2*x0^112 + 3*x0^111 + 3*x0^110 + 3*x0^109
        + x0^108 + 2*x0^106 + 3*x0^105 + 3*x0^104 + 2*x0^103 + 4*x0^102 + 3*x0^101 + 3*x0^100 + 2*x0^99 + x0^98 + 3*x0^97 + 2*x0^96 +
        3*x0^95 + 4*x0^94 + 3*x0^93 + 3*x0^92 + x0^91 + 3*x0^90 + x0^88 + 3*x0^87 + 4*x0^85 + x0^84 + x0^83 + x0^81 + x0^80 + 2*x0^79
        + 3*x0^78 + 4*x0^77 + x0^76 + 3*x0^74 + x0^73 + 2*x0^72 + 3*x0^71 + x0^69 + 2*x0^68 + x0^67 + x0^65 + 3*x0^64 + 4*x0^63 + 
        x0^62 + 4*x0^61 + 4*x0^60 + 3*x0^59 + x0^58 + 3*x0^57 + 2*x0^56 + x0^55 + 4*x0^54 + x0^53 + 3*x0^52 + 4*x0^51 + 4*x0^50 + 
        3*x0^49 + x0^48 + 2*x0^47 + 4*x0^46 + x0^45 + 3*x0^43 + 3*x0^42 + 4*x0^41 + x0^39 + 4*x0^38 + x0^37 + 2*x0^36 + 4*x0^35 + 
        3*x0^33 + 4*x0^31 + x0^30 + 3*x0^28 + 4*x0^27 + 2*x0^26 + 4*x0^24 + 2*x0^23 + 3*x0^22 + 2*x0^21 + 2*x0^20 + x0^19 + 2*x0^18 +
        4*x0^16 + 4*x0^14 + 3*x0^13 + 2*x0^11 + 3*x0^9 + 4*x0^7 + 2*x0^5 + 4*x0^2 + 4*x0 + 3
]
[
    1,
    x0^8 + 2*x0^6 + x0^5 + 2*x0^4 + x0^3 + x0^2 + x0,
    x0^56 + 2*x0^54 + 4*x0^53 + x0^52 + 4*x0^51 + 4*x0^50 + 3*x0^48 + 3*x0^45 + 3*x0^44 + 2*x0^42 + 3*x0^41 + 3*x0^40 + 3*x0^39 + 
        4*x0^38 + 3*x0^37 + x0^34 + 4*x0^33 + 3*x0^32 + 2*x0^31 + x0^30 + x0^29 + 4*x0^27 + x0^25 + 3*x0^24 + 3*x0^23 + 2*x0^20 + 
        3*x0^19 + 2*x0^18 + 4*x0^17 + 3*x0^16 + x0^14 + 2*x0^13 + 4*x0^11 + 4*x0^8 + 3*x0^7 + 2*x0^5 + x0^4 + 3*x0 + 1,
    x0^336 + 2*x0^334 + 4*x0^333 + x0^332 + 4*x0^331 + 4*x0^330 + 2*x0^326 + 4*x0^324 + 3*x0^323 + 2*x0^322 + 2*x0^321 + 3*x0^320 + 
        3*x0^319 + x0^318 + 4*x0^317 + 2*x0^316 + x0^315 + 2*x0^314 + 4*x0^313 + x0^312 + 3*x0^311 + 4*x0^310 + 3*x0^309 + x0^308 + 
        4*x0^307 + x0^305 + 3*x0^304 + x0^303 + 4*x0^302 + x0^301 + x0^300 + x0^296 + 2*x0^294 + 4*x0^293 + x0^292 + 4*x0^291 + 
        4*x0^290 + 2*x0^286 + 4*x0^284 + 3*x0^283 + 2*x0^282 + 4*x0^281 + 3*x0^280 + 2*x0^279 + x0^278 + x0^277 + x0^276 + x0^275 + 
        x0^274 + 3*x0^273 + x0^271 + 2*x0^270 + 4*x0^269 + 4*x0^267 + 4*x0^265 + 4*x0^263 + 4*x0^262 + 3*x0^261 + x0^260 + 4*x0^259 +
        3*x0^258 + x0^257 + 3*x0^255 + 2*x0^254 + 4*x0^253 + x0^252 + 2*x0^251 + x0^250 + 3*x0^249 + 3*x0^248 + x0^247 + 4*x0^246 + 
        x0^245 + 3*x0^244 + x0^242 + 4*x0^241 + x0^240 + x0^239 + 4*x0^238 + x0^237 + 4*x0^236 + 2*x0^234 + 3*x0^233 + 4*x0^229 + 
        2*x0^228 + 4*x0^227 + x0^226 + 4*x0^224 + 2*x0^223 + 3*x0^222 + 2*x0^221 + x0^220 + x0^219 + x0^217 + x0^216 + 3*x0^215 + 
        x0^214 + x0^213 + 3*x0^212 + 3*x0^211 + 4*x0^210 + 4*x0^208 + 3*x0^207 + 4*x0^206 + 3*x0^205 + 3*x0^204 + 4*x0^203 + 4*x0^202
        + 2*x0^201 + 2*x0^199 + 2*x0^198 + x0^196 + 4*x0^195 + x0^194 + 3*x0^193 + 4*x0^192 + 3*x0^191 + 4*x0^190 + x0^188 + 4*x0^187
        + 4*x0^186 + 2*x0^185 + 3*x0^184 + 2*x0^182 + 4*x0^181 + 2*x0^180 + x0^178 + 2*x0^177 + 4*x0^176 + x0^175 + 3*x0^174 + 
        4*x0^173 + 2*x0^172 + 3*x0^170 + 3*x0^169 + 4*x0^168 + x0^166 + 2*x0^165 + 3*x0^164 + x0^163 + 3*x0^162 + 4*x0^161 + 2*x0^159
        + x0^158 + x0^157 + 2*x0^156 + 2*x0^154 + 4*x0^152 + 2*x0^151 + 4*x0^150 + x0^149 + 4*x0^148 + 3*x0^146 + x0^145 + 3*x0^144 +
        x0^142 + 2*x0^141 + 3*x0^140 + 4*x0^139 + 2*x0^138 + 2*x0^137 + 4*x0^136 + 4*x0^135 + x0^134 + 3*x0^133 + x0^131 + x0^129 + 
        x0^128 + 2*x0^126 + x0^125 + x0^124 + 2*x0^123 + 4*x0^122 + 2*x0^121 + 3*x0^120 + 3*x0^119 + 3*x0^118 + x0^117 + 3*x0^116 + 
        2*x0^114 + 4*x0^113 + 3*x0^111 + 3*x0^110 + 3*x0^109 + 4*x0^108 + x0^107 + 4*x0^106 + 2*x0^104 + 3*x0^103 + x0^101 + 3*x0^100
        + 4*x0^99 + 3*x0^98 + 4*x0^97 + 4*x0^95 + 3*x0^94 + 2*x0^93 + 2*x0^92 + 2*x0^91 + 3*x0^90 + 4*x0^89 + 4*x0^88 + x0^87 + x0^86
        + 3*x0^85 + 2*x0^84 + x0^82 + 4*x0^81 + 2*x0^80 + 3*x0^78 + 4*x0^75 + 3*x0^74 + x0^73 + 4*x0^72 + 2*x0^70 + x0^69 + 4*x0^68 +
        3*x0^66 + 3*x0^64 + 2*x0^63 + 3*x0^62 + x0^61 + x0^59 + 4*x0^58 + 4*x0^57 + x0^56 + 3*x0^55 + x0^54 + x0^53 + 3*x0^52 + 
        3*x0^51 + 2*x0^50 + 3*x0^49 + 4*x0^48 + 2*x0^47 + 4*x0^46 + x0^44 + 4*x0^43 + x0^42 + 3*x0^40 + 2*x0^39 + 3*x0^38 + x0^37 + 
        3*x0^36 + 3*x0^35 + 3*x0^34 + x0^33 + 4*x0^32 + 3*x0^31 + 4*x0^29 + 4*x0^28 + x0^27 + x0^26 + 4*x0^24 + 2*x0^22 + 4*x0^19 + 
        x0^18 + 4*x0^17 + 2*x0^16 + 2*x0^14 + 3*x0^13 + 4*x0^11 + 4*x0^10 + 3*x0^9 + 4*x0^8 + 2*x0^7 + x0^5 + x0^3 + 4*x0^2 + 4*x0
]
```

### Lifting the j-Invariant

The file `lift_j.m` has routines to compute the formulas for the
j-invariant of the canonical lifting as described in [Coordinates of the
j-Invariant of the Canonical Lifting][jinv].

To find the formulas in characteristic 5 for length 4:
```
> lift_j(5,3);
[
    j0,
    j0^4 + 3*j0^3,
    j0^24 + 3*j0^23 + j0^20 + j0^19 + 4*j0^18 + j0^17 + 4*j0^16 + 4*j0^15 + 4*j0^14 + 2*j0^13 + 2*j0^10 + 3*j0^5,
    (j0^149 + 3*j0^148 + j0^145 + j0^144 + 4*j0^143 + j0^142 + 3*j0^141 + 3*j0^140 + 4*j0^139 + j0^138 + 2*j0^137 + 4*j0^136 + j0^135
        + 2*j0^134 + j0^133 + j0^132 + 4*j0^130 + j0^129 + 4*j0^128 + 3*j0^126 + 3*j0^125 + 2*j0^123 + 2*j0^122 + j0^120 + 3*j0^119 +
        3*j0^118 + 4*j0^117 + j0^116 + 3*j0^114 + 2*j0^113 + 2*j0^112 + 3*j0^111 + j0^110 + 4*j0^109 + 3*j0^107 + 2*j0^106 + 4*j0^104
        + 4*j0^103 + 3*j0^102 + 4*j0^100 + 4*j0^99 + 3*j0^98 + j0^97 + j0^95 + 2*j0^94 + 2*j0^93 + 3*j0^92 + 4*j0^91 + 4*j0^90 + 
        3*j0^89 + 2*j0^88 + 2*j0^87 + j0^85 + 3*j0^84 + 2*j0^83 + 3*j0^82 + 3*j0^80 + 2*j0^79 + j0^78 + 3*j0^77 + j0^76 + 3*j0^75 + 
        2*j0^73 + j0^72 + 3*j0^70 + 2*j0^69 + 3*j0^68 + 3*j0^65 + 4*j0^63 + 4*j0^62 + j0^61 + j0^59 + 3*j0^58 + 3*j0^55 + 4*j0^50 + 
        3*j0^45 + 4*j0^40 + j0^25 + 1)/j0^25
]
```

As it uses `GT`, you can again make a choice of method used with the
optional argument `choice` and pass along eta polynomials with `pols`
and binomials with `bintab`.
```
> bt:=BinTab(2,5);
> lift_j(2,5 : choice:=2, bintab:=bt);
[
    j0,
    0,
    0,
    0,
    j0^8,
    j0^24 + j0^16
]
```

One can also compute the j-invariant of the canonical lifting of
curves over finite fields, instead of the general the formulas above,
using `lift_j0`.  In this case, note it is necessary that the
j-invariant given, say `j0`, is not in `GF(p^2)`.  (In particular,
this guarantees that the elliptic curve is not supersingular!)

```
> p:=7; d:=5; F<a>:=GF(p^d);
> j0 := Random(F);
> (j0^(p^2) - j0) eq F!0;
false

> lift_j0(j0,4);
[ a^4843, a^1676, a^6953, a^597, a^9271 ]
```


Using information from [Coordinates of the j-Invariant of the
Canonical Lifting][jinv], we can also compute the formulas by
interpolation (using `lift_j0` sufficiently many times) with
`lift_j_int`.  This is generally much slower than `lift_j`, but it
uses *a lot* less memory.  Also, note that `lift_j_int` is only
implemented for characteristic 5 or more, which `lift_j` works for any
characteristic.
```
> lift_j_int(5,3);
[
    j0,
    j0^4 + 3*j0^3,
    j0^24 + 3*j0^23 + j0^20 + j0^19 + 4*j0^18 + j0^17 + 4*j0^16 + 4*j0^15 + 4*j0^14 + 2*j0^13 + 2*j0^10 + 3*j0^5,
    (j0^149 + 3*j0^148 + j0^145 + j0^144 + 4*j0^143 + j0^142 + 3*j0^141 + 3*j0^140 + 4*j0^139 + j0^138 + 2*j0^137 + 4*j0^136 + j0^135
        + 2*j0^134 + j0^133 + j0^132 + 4*j0^130 + j0^129 + 4*j0^128 + 3*j0^126 + 3*j0^125 + 2*j0^123 + 2*j0^122 + j0^120 + 3*j0^119 +
        3*j0^118 + 4*j0^117 + j0^116 + 3*j0^114 + 2*j0^113 + 2*j0^112 + 3*j0^111 + j0^110 + 4*j0^109 + 3*j0^107 + 2*j0^106 + 4*j0^104
        + 4*j0^103 + 3*j0^102 + 4*j0^100 + 4*j0^99 + 3*j0^98 + j0^97 + j0^95 + 2*j0^94 + 2*j0^93 + 3*j0^92 + 4*j0^91 + 4*j0^90 + 
        3*j0^89 + 2*j0^88 + 2*j0^87 + j0^85 + 3*j0^84 + 2*j0^83 + 3*j0^82 + 3*j0^80 + 2*j0^79 + j0^78 + 3*j0^77 + j0^76 + 3*j0^75 + 
        2*j0^73 + j0^72 + 3*j0^70 + 2*j0^69 + 3*j0^68 + 3*j0^65 + 4*j0^63 + 4*j0^62 + j0^61 + j0^59 + 3*j0^58 + 3*j0^55 + 4*j0^50 + 
        3*j0^45 + 4*j0^40 + j0^25 + 1)/j0^25
]
```

Finally, the file `lift_j.m` also provides the function `jweier` which
gives formulas for the Weierstrass coefficients of the canonical
lifting of the curve `y0^2 = x0^3 + a0*x0 + b0` in terms of `a0` and
`b0`.  These are *not* universal (as defined in [Weierstrass
Coefficients of the Canonical Lifting][wcoef]), as it might contain
extra powers of `a0` and `b0` in the denominator.  (See also
[Denominators of the Weierstrass Coefficients of the Canonical
Lifting][delong] for more details.)

[delong]: http://www.math.utk.edu/~finotti/papers/delong.pdf

To use just give the prime and the index:
```
> jweier(5,2);
[
    [
        a0,
        (2*a0^12 + 3*a0^9*b0^2 + 3*a0^6*b0^4 + 3*a0^3*b0^6 + 3*b0^8)/(a0*b0^4),
        (3*a0^120 + 4*a0^105*b0^10 + 3*a0^96*b0^16 + 2*a0^93*b0^18 + a0^90*b0^20 + 4*a0^87*b0^22 + 4*a0^81*b0^26 + 4*a0^78*b0^28 + a0^75*b0^30 + 
            4*a0^72*b0^32 + 3*a0^66*b0^36 + 4*a0^63*b0^38 + 4*a0^60*b0^40 + 4*a0^57*b0^42 + 2*a0^54*b0^44 + a0^51*b0^46 + 2*a0^45*b0^50 + 4*a0^36*b0^56
            + 4*a0^33*b0^58 + a0^30*b0^60 + 2*a0^27*b0^62 + 2*a0^24*b0^64 + a0^15*b0^70 + 3*b0^80)/(a0^35*b0^40)
    ],
    [
        b0,
        (2*a0^12*b0 + 3*a0^9*b0^3 + 3*a0^6*b0^5 + 3*a0^3*b0^7 + 3*b0^9)/a0^6,
        (3*a0^120 + 4*a0^105*b0^10 + 3*a0^96*b0^16 + 2*a0^93*b0^18 + a0^90*b0^20 + 4*a0^87*b0^22 + 4*a0^81*b0^26 + 4*a0^78*b0^28 + a0^75*b0^30 + 
            4*a0^72*b0^32 + 3*a0^66*b0^36 + 4*a0^63*b0^38 + 4*a0^60*b0^40 + 4*a0^57*b0^42 + 2*a0^54*b0^44 + a0^51*b0^46 + 2*a0^45*b0^50 + 4*a0^36*b0^56
            + 4*a0^33*b0^58 + a0^30*b0^60 + 2*a0^27*b0^62 + 2*a0^24*b0^64 + a0^15*b0^70 + 3*b0^80)/(a0^60*b0^15)
    ]
]
```

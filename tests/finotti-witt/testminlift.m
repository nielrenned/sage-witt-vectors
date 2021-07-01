p:=11;
F<a0,b0>:=RationalFunctionField(GF(p),2);

load "fctWitt-GTnew2.m";

time va, vb, vF, vH := lift(a0,b0,2);


load "minlift.m";

time wa, wb, wF, wH := minlift(a0,b0,2);

va eq wa;
vb eq wb;
vF eq wF;
vH eq wH;

/* va[2] eq a1; */
/* va[3] eq a2; */
/* vb[2] eq b1; */
/* vb[3] eq b2; */

//P:=Parent(F1);

/* Evaluate(vF[2],P.1) eq F1; */
/* Evaluate(vF[3],P.1) eq F2; */
/* Evaluate(vH[2],P.1) eq H1; */
/* Evaluate(vH[3],P.1) eq H2; */

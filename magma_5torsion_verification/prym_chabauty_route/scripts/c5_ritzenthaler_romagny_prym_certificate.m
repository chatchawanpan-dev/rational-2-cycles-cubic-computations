/*
    Ritzenthaler--Romagny Prym certificate for the X_0(5) curve C5.

    This script supplies the missing Section 6 certification point: it
    specializes Theorem 1.1 of Ritzenthaler--Romagny, "On the Prym variety
    of genus 3 covers of genus 1 curves", to the plane quartic C5.

    The theorem applies to a smooth non-hyperelliptic genus-3 quartic written
    as

        y^4 - h(x,z)y^2 + f(x,z)g(x,z) = 0,

    where f, g, h are binary quadratics over the base field and the
    3-by-3 coefficient matrix A with rows f,h,g is invertible.  It gives
    an explicit genus-2 curve X_RR whose Jacobian is the Prym factor, up to
    Q-isogeny:

        Jac(C5) ~ Jac(D) x Jac(X_RR).

    The script verifies:

      1. the normalized C5 equation has the required form over Q;
      2. det(A) is nonzero;
      3. the Ritzenthaler--Romagny formula gives a genus-2 curve X_RR;
      4. X_RR is Q-isomorphic to the proposed P5 used in the rank/torsion
         computations.
*/

SetColumns(0);

Q := Rationals();

procedure Header(s)
    print "";
    print s;
    print &cat["-" : i in [1..#s]];
end procedure;

Header("1. Normalize C5 to the Ritzenthaler--Romagny form");

P2<X,Y,W> := ProjectiveSpace(Q, 2);
Fcan :=
    729*X^4 + 864*X^2*Y^2 + 256*Y^4
    + 20736*X^3*W - 24576*X*Y^2*W
    + 73728*X^2*W^2 + 65536*Y^2*W^2
    - 786432*X*W^3;
C5 := Curve(P2, Fcan);

Fmonic := Fcan/256;

f := X*(3*X + 64*W);
h := -Q!27/8*X^2 + 96*X*W - 256*W^2;
g := Q!3/256*(81*X^2 + 576*X*W - 4096*W^2);

print "Fcan/256 =";
print Fmonic;
print "f =", f;
print "h =", h;
print "g =", g;
print "Check Fcan/256 = Y^4 - h*Y^2 + f*g:", Fmonic eq Y^4 - h*Y^2 + f*g;
assert Fmonic eq Y^4 - h*Y^2 + f*g;
assert IsNonsingular(C5);
assert Genus(C5) eq 3;
assert not IsHyperelliptic(C5);

Header("2. Coefficient matrix and Ritzenthaler--Romagny data");

/*
   Rows are coefficients of f, h, g in the basis X^2, XW, W^2.
*/
A := Matrix(Q, 3, 3, [
    3,       64,       0,
    -Q!27/8, 96,    -256,
    Q!243/256, Q!27/4, -48
]);

Ainv := A^-1;
print "A =";
print A;
print "det(A) =", Determinant(A);
print "A^-1 =";
print Ainv;
assert Determinant(A) ne 0;

R<t> := PolynomialRing(Q);

/*
   If A^-1 = [a_i b_i c_i]_{i=1..3}, then
      a(t)=a1+2*a2*t+a3*t^2,
      b(t)=b1+2*b2*t+b3*t^2,
      c(t)=c1+2*c2*t+c3*t^2.
*/
arr := Eltseq(Ainv);
a1 := arr[1]; b1 := arr[2]; c1 := arr[3];
a2 := arr[4]; b2 := arr[5]; c2 := arr[6];
a3 := arr[7]; b3 := arr[8]; c3 := arr[9];

a := a1 + 2*a2*t + a3*t^2;
b := b1 + 2*b2*t + b3*t^2;
c := c1 + 2*c2*t + c3*t^2;
fRR := b*(b^2 - a*c);

print "a(t) =", a;
print "b(t) =", b;
print "c(t) =", c;
print "Ritzenthaler--Romagny genus-2 polynomial b*(b^2-a*c) =";
print fRR;

XRR := HyperellipticCurve(fRR);
print "Genus(X_RR) =", Genus(XRR);
assert Genus(XRR) eq 2;

Header("3. Compare X_RR with the proposed P5");

fP5 := -5*(6*t^2 - 2*t + 1)*
           (324*t^4 + 144*t^3 + 14*t^2 - 6*t - 1);
P5 := HyperellipticCurve(fP5);

print "P5 polynomial =";
print fP5;
print "Genus(P5) =", Genus(P5);
assert Genus(P5) eq 2;

is_iso, iso := IsIsomorphic(XRR, P5);
print "Is X_RR Q-isomorphic to P5? =", is_iso;
print "One Q-isomorphism X_RR -> P5 is:";
print iso;
assert is_iso;

isos := Isomorphisms(XRR, P5);
print "All Q-isomorphisms returned by Magma:";
print isos;

Header("4. Consequence");

print "The Ritzenthaler--Romagny formula applied over Q gives X_RR.";
print "Magma proves X_RR is Q-isomorphic to the displayed P5.";
print "Therefore the Prym factor attached to C5 -> C5/<Y -> -Y>";
print "is Q-isogenous to J(P5).  This is the missing correspondence";
print "certificate needed by the Section 6 Prym-torsion argument.";

quit;

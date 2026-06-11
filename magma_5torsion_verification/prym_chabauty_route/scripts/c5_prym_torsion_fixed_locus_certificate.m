/*
    Prym-torsion and fixed-locus certificate for the X_0(5) curve C5.

    Feedback 6 observes that a direct Chabauty computation on the
    non-hyperelliptic genus-3 quartic is unnecessary if the Prym P attached
    to C5 -> C5/<Y -> -Y> has rank zero and rational torsion of exponent
    at most 2.

    This script verifies the computational pieces needed for that argument:

      1. C5 is smooth, genus 3, and non-hyperelliptic.
      2. The quotient by Y -> -Y has elliptic model
             E: y^2 = x^3 - x^2 + 12*x + 72.
      3. The proposed genus-2 Prym candidate P5 has rank-zero Jacobian.
      4. #J(P5)(F_7), #J(P5)(F_11), #J(P5)(F_13) are 44, 108, 134,
         with gcd 2.
      5. The fixed locus Y=0 on C5 has exactly the two rational points
             (0:0:1), (-64/3:0:1).

    The companion script c5_ritzenthaler_romagny_prym_certificate.m verifies
    the explicit Ritzenthaler--Romagny decomposition and the Q-isomorphism
    between the resulting genus-2 curve and P5.
*/

SetColumns(0);

Q := Rationals();
Zint := Integers();

procedure Header(s)
    print "";
    print s;
    print &cat["-" : i in [1..#s]];
end procedure;

function SquarePointCount(rhs)
    if rhs eq 0 then
        return 1;
    elif IsSquare(rhs) then
        return 2;
    else
        return 0;
    end if;
end function;

function CountP5(Fq)
    N := 0;
    for s in Fq do
        rhs := -5*(6*s^2 - 2*s + 1)*
                  (324*s^4 + 144*s^3 + 14*s^2 - 6*s - 1);
        N +:= SquarePointCount(rhs);
    end for;

    lc := Fq!(-9720);
    if lc ne 0 and IsSquare(lc) then
        N +:= 2;
    end if;

    return N;
end function;

function JacOrderGenus2FromCounts(p)
    Fp := GF(p);
    Fp2 := GF(p^2);

    N1 := CountP5(Fp);
    N2 := CountP5(Fp2);

    a1 := p + 1 - N1;
    a2 := (N2 - (p^2 + 1) + a1^2) div 2;

    return p^2 + 1 + a2 - (p + 1)*a1, N1, N2, a1, a2;
end function;

Header("1. The genus-3 curve C5");

P2<X,Y,W> := ProjectiveSpace(Q, 2);
Fcan :=
    729*X^4 + 864*X^2*Y^2 + 256*Y^4
    + 20736*X^3*W - 24576*X*Y^2*W
    + 73728*X^2*W^2 + 65536*Y^2*W^2
    - 786432*X*W^3;
C5 := Curve(P2, Fcan);

print "Defining polynomial of C5:";
print Fcan;
print "Genus(C5) =", Genus(C5);
print "IsNonsingular(C5) =", IsNonsingular(C5);
print "IsHyperelliptic(C5) =", IsHyperelliptic(C5);
assert Genus(C5) eq 3;
assert IsNonsingular(C5);
assert not IsHyperelliptic(C5);

Header("2. Quotient elliptic curve");

/*
   As in the earlier probe, write the quartic as a quadratic in Y^2.
   Its discriminant gives the quotient by Y |-> -Y.
*/
R<x> := PolynomialRing(Q);
coef_w2 := Q!256/729;
coef_w1 := Q!32/27*x^2 - Q!8192/243*x + Q!65536/729;
coef_w0 := x^4 + Q!256/9*x^3 + Q!8192/81*x^2 - Q!262144/243*x;
disc_w := coef_w1^2 - 4*coef_w2*coef_w0;

H := HyperellipticCurve(disc_w);
Pbase := H![0, Q!65536/729, 1];
Eraw, HtoEraw := EllipticCurve(H, Pbase);
E, ErawToE := MinimalModel(Eraw);

Eexpected := EllipticCurve([Q!0, Q!-1, Q!0, Q!12, Q!72]);
print "Quotient elliptic curve minimal model =", E;
print "Expected model =", Eexpected;
assert E eq Eexpected;

rElo, rEhi := RankBounds(E);
print "RankBounds(E) =", rElo, rEhi;
print "TorsionSubgroup(E) =", TorsionSubgroup(E);
assert rElo eq 1 and rEhi eq 1;

Header("3. Prym candidate rank and finite-field torsion bound");

S<s> := PolynomialRing(Q);
fP := -5*(6*s^2 - 2*s + 1)*
          (324*s^4 + 144*s^3 + 14*s^2 - 6*s - 1);
P5 := HyperellipticCurve(fP);
JP5 := Jacobian(P5);

print "P5: v^2 =", fP;
print "Genus(P5) =", Genus(P5);
print "Discriminant(fP) factorization =", Factorization(Zint!Discriminant(fP));
assert Genus(P5) eq 2;

rPlo, rPhi := RankBounds(JP5);
print "RankBounds(Jacobian(P5)) =", rPlo, rPhi;
assert rPlo eq 0 and rPhi eq 0;

TP, TPmap := TorsionSubgroup(JP5);
print "TorsionSubgroup(Jacobian(P5)) =", TP;

orders := [];
for p in [7, 11, 13] do
    order, N1, N2, a1, a2 := JacOrderGenus2FromCounts(p);
    Append(~orders, order);
    printf "p=%o: #P5(F_p)=%o, #P5(F_{p^2})=%o, a1=%o, a2=%o, #J(P5)(F_p)=%o\n",
           p, N1, N2, a1, a2, order;

    /*
       Direct finite-field Jacobian order, if available, as an independent
       consistency check.
    */
    Fp := GF(p);
    JP5p := Jacobian(ChangeRing(P5, Fp));
    direct_order := #JP5p;
    print "  direct #Jacobian(ChangeRing(P5, GF(p))) =", direct_order;
    assert direct_order eq order;
end for;

print "Orders at p=7,11,13 =", orders;
print "GCD of these orders =", GCD(orders);
assert orders eq [44, 108, 134];
assert GCD(orders) eq 2;

Header("4. Fixed locus of the involution Y -> -Y");

print "The only projective fixed point with Y != 0 would be (0:1:0).";
print "Fcan(0,1,0) =", Evaluate(Fcan, [Q!0, Q!1, Q!0]);
assert Evaluate(Fcan, [Q!0, Q!1, Q!0]) ne 0;

fixed_poly := Evaluate(Fcan, [X, 0, W]);
expected_fixed :=
    3*X*(3*X + 64*W)*(81*X^2 + 576*X*W - 4096*W^2);

print "Fcan(X,0,W) =", fixed_poly;
print "Factorization of Fcan(X,0,W) =", Factorization(fixed_poly);
assert fixed_poly eq expected_fixed;

Rx<xaff> := PolynomialRing(Q);
quad := 81*xaff^2 + 576*xaff - 4096;
disc_quad := Discriminant(quad);
print "Quadratic factor after W=1:", quad;
print "Discriminant =", disc_quad;
print "576^2 * 5 =", 576^2 * 5;
print "Is discriminant a rational square? =", IsSquare(disc_quad);
assert disc_quad eq 576^2 * 5;
assert not IsSquare(disc_quad);

fixed_pts_Q := {@ C5![0, 0, 1], C5![-Q!64/3, 0, 1] @};
print "Rational fixed points certified by the factorization =", fixed_pts_Q;
for Pt in fixed_pts_Q do
    assert Pt[2] eq 0;
end for;

Header("5. Link to the singular X_0(5) fiber product");

print "The two fixed points correspond to the two visible canonical points";
print "from the X_0(5) fiber-product normalization certificate:";
print "  (0:0:1) and (-64/3:0:1).";
print "The companion normalization certificate avoids singular-model";
print "fiber artifacts and proves directly on Z=1 that";
print "the second canonical coordinate cannot vanish for rational L*T != 0.";

Header("6. Certificate status");

print "Verified computational inputs for the Prym/fixed-locus proof route:";
print "  - C5 is smooth, genus 3, and non-hyperelliptic.";
print "  - The elliptic quotient is y^2 = x^3 - x^2 + 12*x + 72.";
print "  - The Ritzenthaler--Romagny Prym candidate has rank-zero Jacobian.";
print "  - The finite-field Jacobian orders at 7, 11, 13 are 44, 108, 134.";
print "  - Their gcd is 2.";
print "  - The fixed locus Y=0 has exactly two rational points.";
print "Together with the RR certificate and the X_0(5) normalization";
print "certificate, this supplies the computational input cited in Section 6.";

quit;

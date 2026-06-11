/*
    Normalization and degeneracy certificate for the X_0(5) fiber product.

    This script records only the proof-grade data used in Theorem 6.2:
      1. the projective X_0(5) fiber product;
      2. its canonical smooth plane-quartic model;
      3. the canonical map on the affine nondegenerate chart;
      4. a direct contradiction for any rational point with L*T*Z != 0.

    The rational-point determination on the canonical model is proved by the
    Prym certificate in the companion scripts.  Given that determination, the
    factorization below shows that no rational point of the fiber product has
    L*T*Z != 0.  This avoids using scheme-theoretic inverse images of a
    rational map on the singular plane model.
*/

SetColumns(0);

Q := Rationals();

procedure Header(s)
    print "";
    print s;
    print &cat["-" : i in [1..#s]];
end procedure;

Header("1. X_0(5) fiber product");

P2<L,T,Z> := ProjectiveSpace(Q, 2);
A := L^2 + 10*L*Z + 5*Z^2;
F := T^2*(16*Z^2 + 27*T^2)*A^3 + 4096*L*Z^9;
Csing := Curve(P2, F);

singular_pts := SingularPoints(Csing);
expected_singular := {@ Csing![0,1,0], Csing![1,0,0] @};
assert singular_pts eq expected_singular;

print "Projective equation:";
print F;
print "Genus(Csing) =", Genus(Csing);
print "Singular points =", singular_pts;
print "Every singular point has L*T*Z = 0.";

Header("2. Canonical plane-quartic normalization");

phi := CanonicalMap(Csing);
Ccan, is_hyp := CanonicalImage(Csing, phi);
P2c<X,Y,W> := Ambient(Ccan);
Fcan := DefiningPolynomial(Ccan);
Fexpected :=
    729*X^4 + 864*X^2*Y^2 + 256*Y^4
    + 20736*X^3*W - 24576*X*Y^2*W
    + 73728*X^2*W^2 + 65536*Y^2*W^2
    - 786432*X*W^3;

assert Fcan eq Fexpected;
assert Genus(Ccan) eq 3;
assert IsNonsingular(Ccan);
assert not is_hyp;

print "Canonical image is non-hyperelliptic =", not is_hyp;
print "Canonical quartic equation:";
print Fcan;
print "Genus(Ccan) =", Genus(Ccan);
print "IsNonsingular(Ccan) =", IsNonsingular(Ccan);

Header("3. Canonical map on the affine nondegenerate chart");

polys := DefiningPolynomials(phi);
R<l,t> := PolynomialRing(Q, 2, "grevlex");
affine_subs := [l, t, 1];
phi1_aff := Evaluate(polys[1], affine_subs);
phi2_aff := Evaluate(polys[2], affine_subs);
phi3_aff := Evaluate(polys[3], affine_subs);
A_aff := l^2 + 10*l + 5;
phi2_factor := (Q!729/64)*l*(t^2 + Q!8/27)*A_aff;
assert phi2_aff eq phi2_factor;

print "On the affine chart Z=1, lambda=L/Z=l and t=T/Z=t.";
print "The second canonical coordinate is";
print phi2_aff;
print "Factorization:";
print Factorization(phi2_aff);
print "Equivalently, phi_2 = (729/64)*l*(t^2 + 8/27)*(l^2 + 10*l + 5).";
print "If l*t != 0 and l,t are rational, then phi_2 cannot vanish:";
print "  t^2 + 8/27 has no rational root.";
print "  l^2 + 10*l + 5 has discriminant 80, not a rational square.";
assert not IsSquare(Q!80);

Header("4. Direct affine fiber check over the two rational canonical points");

F_aff := Evaluate(F, affine_subs);
I_origin := ideal<R | F_aff, phi1_aff, phi2_aff>;
I_second := ideal<R | F_aff, phi2_aff, 3*phi1_aff + 64*phi3_aff>;
J_origin := Saturation(I_origin, ideal<R | l*t*phi3_aff>);
J_second := Saturation(I_second, ideal<R | l*t*phi3_aff>);
assert 1 in J_origin;
assert 1 in J_second;

print "For canonical point (0:0:1), saturating by l*t*phi_3 gives";
print GroebnerBasis(J_origin);
print "For canonical point (-64/3:0:1), saturating by l*t*phi_3 gives";
print GroebnerBasis(J_second);
print "Thus there is no affine point with l*t != 0 and phi_3 != 0 mapping to either canonical rational point.";

Header("5. Conclusion");

print "Assuming Proposition 6.1, Ccan(Q) has only the two points above.";
print "Both points have second canonical coordinate Y=0.";
print "On the rational affine nondegenerate chart this would force";
print "l*(t^2 + 8/27)*(l^2 + 10*l + 5) = 0, impossible when l*t != 0.";
print "Hence the fiber product has no rational point with lambda*t != 0.";

quit;

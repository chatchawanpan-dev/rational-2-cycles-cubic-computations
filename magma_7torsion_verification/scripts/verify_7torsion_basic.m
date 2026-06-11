Q := Rationals();
Z := Integers();

procedure Section(title)
    print "";
    print "======================================================================";
    print title;
    print "======================================================================";
end procedure;

function RatSquareRoot(q)
    if q lt 0 then
        return false, Q!0;
    end if;

    num := Z!Numerator(q);
    den := Z!Denominator(q);
    bool_num, root_num := IsSquare(num);
    bool_den, root_den := IsSquare(den);
    if bool_num and bool_den then
        return true, (Q!root_num)/(Q!root_den);
    end if;

    return false, Q!0;
end function;

function StringSet(seq)
    return SequenceToSet([ Sprint(x) : x in seq ]);
end function;

Section("Basic Magma checks for the 7-torsion argument");
print "Manuscript target: Section \"Exclusion of rational 7-torsion\"";

P<mu> := PolynomialRing(Q);

A7 := mu^2 - mu + 1;
B7 := mu^6 - 11*mu^5 + 30*mu^4 - 15*mu^3 - 10*mu^2 + 5*mu + 1;
D7 := mu^3 - 8*mu^2 + 5*mu + 1;
N7 := A7^3 * B7^3;
T7 := mu^7 * (mu - 1)^7 * D7;
R7 := mu^12 - 18*mu^11 + 117*mu^10 - 354*mu^9 + 570*mu^8 - 486*mu^7
    + 273*mu^6 - 222*mu^5 + 174*mu^4 - 46*mu^3 - 15*mu^2 + 6*mu + 1;

Section("1. Polynomial identity");
lhs := N7 * (N7 - 1728*T7);
rhs := A7^3 * B7^3 * R7^2;
assert lhs eq rhs;
print "Verified:";
print "  N_7(mu) * (N_7(mu) - 1728*T_7(mu)) = A_7(mu)^3 * B_7(mu)^3 * R_7(mu)^2";
print "Degree(lhs) =", Degree(lhs);

Section("2. Hyperelliptic curve C_7");
f := A7 * B7;
C := HyperellipticCurve(f);
J := Jacobian(C);
disc_f := Discriminant(f);
disc_factors := Factorization(Z!disc_f);
bad_primes := [t[1] : t in disc_factors];

assert IsSquarefree(f);
assert Degree(f) eq 8;
assert Genus(C) eq 3;
assert bad_primes eq [2, 3, 7];
assert #Roots(A7) eq 0;
assert #Roots(B7) eq 0;
assert #Roots(f) eq 0;

print "f(mu) = A_7(mu) * B_7(mu)";
print "Degree(f) =", Degree(f);
print "Genus(C_7) =", Genus(C);
print "IsSquarefree(f) =", IsSquarefree(f);
print "Roots(A_7) over Q =", Roots(A7);
print "Roots(B_7) over Q =", Roots(B7);
print "Roots(A_7*B_7) over Q =", Roots(f);
print "Discriminant(f) =", disc_f;
print "Factorization(Discriminant(f)) =", disc_factors;
print "Bad primes =", bad_primes;
print "RankBound(J(C_7)) =", RankBound(J);

Section("3. Known rational points");
mu_values := [Q!0, Q!1, -(Q!1)/2, (Q!2)/3, Q!3];
known_finite_points := [];
positive_y_values := [];

for muv in mu_values do
    f_val := Evaluate(f, muv);
    bool, yv := RatSquareRoot(f_val);
    assert bool;

    Append(~positive_y_values, yv);
    Append(~known_finite_points, C![muv, yv, 1]);
    Append(~known_finite_points, C![muv, -yv, 1]);

    printf "mu = %o, f(mu) = %o, y = +/- %o\n", muv, f_val, yv;
    printf "  points: %o and %o\n", C![muv, yv, 1], C![muv, -yv, 1];
end for;

known_infinity := SetToSequence(PointsAtInfinity(C));
known_points := known_finite_points cat known_infinity;

print "Points at infinity =", known_infinity;
print "Total number of known points =", #known_points;
print "Known points as a sequence:";
print known_points;

search_points := SetToSequence(Points(C : Bound := 100));
print "Points(C_7 : Bound := 100) =", search_points;
print "Known points match Points(C_7 : Bound := 100):", StringSet(search_points) eq StringSet(known_points);

Section("4. Reduction at p = 5");
p := 5;
Fp := GF(p);
P5<m5> := PolynomialRing(Fp);
C5 := HyperellipticCurve(P5!f);
points_C5 := SetToSequence(Points(C5));

assert #points_C5 eq 12;
print "#C_7(F_5) =", #points_C5;
print "Points(C_7(F_5)) =", points_C5;

reduced_known_points := [];
print "Reduction of the known affine rational points modulo 5:";
for i in [1..#mu_values] do
    muv := mu_values[i];
    yv := positive_y_values[i];
    mu_red := Fp!muv;
    y_red := Fp!yv;

    Pplus := C5![mu_red, y_red, 1];
    Pminus := C5![mu_red, -y_red, 1];
    Append(~reduced_known_points, Pplus);
    Append(~reduced_known_points, Pminus);

    printf "mu = %o, y = +/- %o  -->  %o and %o\n", muv, yv, Pplus, Pminus;
end for;

known_infinity_mod5 := SetToSequence(PointsAtInfinity(C5));
reduced_known_points cat:= known_infinity_mod5;
print "Points at infinity modulo 5 =", known_infinity_mod5;

assert StringSet(reduced_known_points) eq StringSet(points_C5);
print "Known rational points reduce bijectively onto C_7(F_5): true";

all_affine_y_nonzero := true;
for Pt in points_C5 do
    coords := Eltseq(Pt);
    if coords[3] ne 0 and coords[2] eq 0 then
        all_affine_y_nonzero := false;
    end if;
end for;
assert all_affine_y_nonzero;
print "All affine F_5-points have y != 0:", all_affine_y_nonzero;

Section("5. Square obstruction at the noncuspidal mu-values");
Pu<u> := PolynomialRing(Q);
noncuspidal_mu_values := [-(Q!1)/2, (Q!2)/3, Q!3];
expected_roots := SequenceToSet([34992/343, -950272/9261]);

for muv in noncuspidal_mu_values do
    Nv := Evaluate(N7, muv);
    Tv := Evaluate(T7, muv);
    quad := 27*Nv*u^2 + 16*Nv*u + 4096*Tv;
    roots := Roots(quad);

    assert SequenceToSet([r[1] : r in roots]) eq expected_roots;

    print "mu =", muv;
    print "  quadratic in u = t^2:", quad;
    print "  factorization:", Factorization(quad);
    print "  roots:", roots;

    for pair in roots do
        uv := pair[1];
        is_square, sqrt_uv := RatSquareRoot(uv);
        if uv lt 0 then
            printf "    u = %o is negative, so it cannot equal t^2.\n", uv;
        else
            printf "    u = %o is nonnegative.\n", uv;
            print "      numerator factorization =", Factorization(Numerator(uv));
            print "      denominator factorization =", Factorization(Denominator(uv));
            print "      rational square? =", is_square;
            if is_square then
                print "      square root =", sqrt_uv;
            end if;
            assert not is_square;
        end if;
    end for;
end for;

Section("Completed");
print "All basic 7-torsion Magma checks completed successfully.";

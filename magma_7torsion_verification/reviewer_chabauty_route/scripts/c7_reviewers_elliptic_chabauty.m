/*
    Reviewer-recommended elliptic Chabauty verification for C_7.

    This script replaces the invalid proof pattern

        nonempty fake 2-Selmer set covered by known points
        => all rational points are known

    by the standard finite-covering route described in the Magma handbook:

      1. compute the fake 2-Selmer set of C_7;
      2. give explicit representatives delta_i in Q[x]/(f);
      3. for each delta_i, construct the associated elliptic curve
             E_i: Y^2 = gamma_i B_+(X)
         over K = Q(sqrt(21));
      4. compute a finite odd-index Mordell-Weil subgroup of E_i(K);
      5. apply Magma's elliptic Chabauty routine to determine the points
         on E_i(K) whose X-image lies in P^1(Q);
      6. pull the resulting rational X-values back to C_7.

    The crucial certificate is:

      - PseudoMordellWeilGroup returns success = true, so the supplied subgroup
        has finite odd index in E_i(K);
      - Chabauty(... : IndexBound := 2) returns an R_i which is a power of 2;
      - hence gcd(index, R_i) = 1 for every i.
*/

SetVerbose("Selmer", 1);
SetVerbose("EllChab", 1);

Q := Rationals();
Z := Integers();

procedure Section(title)
    print "";
    print "======================================================================";
    print title;
    print "======================================================================";
end procedure;

function IsPowerOfTwo(n)
    n := Z!n;
    if n le 0 then
        return false;
    end if;
    while n mod 2 eq 0 do
        n div:= 2;
    end while;
    return n eq 1;
end function;

function StringSet(seq)
    return SequenceToSet([ Sprint(x) : x in seq ]);
end function;

function RationalSquareRoot(q)
    if q lt 0 then
        return false, Q!0;
    end if;
    ok_num, root_num := IsSquare(Z!Numerator(q));
    ok_den, root_den := IsSquare(Z!Denominator(q));
    if ok_num and ok_den then
        return true, (Q!root_num)/(Q!root_den);
    end if;
    return false, Q!0;
end function;

function ElementToRational(a)
    if Parent(a) cmpeq Q then
        return Q!a;
    end if;

    coeffs := Eltseq(a);
    assert #coeffs ge 1;
    for j in [2..#coeffs] do
        assert coeffs[j] eq 0;
    end for;
    return Q!coeffs[1];
end function;

function P1PointToRationalP1(Pt, P1Target)
    coords := Eltseq(Pt);
    assert #coords eq 2;
    c0 := coords[1];
    c1 := coords[2];

    if c1 eq 0 then
        assert c0 ne 0;
        return P1Target![1, 0];
    end if;

    return P1Target![ElementToRational(c0/c1), 1];
end function;

Section("0. Curve, known points, and notation");

P<x> := PolynomialRing(Q);
A7 := x^2 - x + 1;
B7 := x^6 - 11*x^5 + 30*x^4 - 15*x^3 - 10*x^2 + 5*x + 1;
f := A7 * B7;

C := HyperellipticCurve(f);
assert Genus(C) eq 3;
assert IsSquarefree(f);

P1Q := ProjectiveSpace(Q, 1);
CtoP1 := map< C -> P1Q | [C.1, C.3] >;

known_mu_values := [ Q!0, Q!1, -Q!1/2, Q!2/3, Q!3 ];
known_x_values := { P1Q | P1Q![m, 1] : m in known_mu_values };
known_x_values join:= { P1Q![1, 0] };

known_points := SetToSequence(PointsAtInfinity(C));
for mu in known_mu_values do
    ok, y := RationalSquareRoot(Evaluate(f, mu));
    assert ok;
    Append(~known_points, C![mu, y, 1]);
    Append(~known_points, C![mu, -y, 1]);
end for;

print "Curve C_7:";
print C;
print "Known rational X-values on P^1:";
print known_x_values;
print "Known rational points on C_7:";
print known_points;

Section("1. Fake 2-Selmer set and explicit representatives");

S, AtoS := TwoCoverDescent(C);
A<theta> := Domain(AtoS);

delta_names := [
    "infinity",
    "mu = 0",
    "mu = 1",
    "mu = -1/2",
    "mu = 2/3",
    "mu = 3"
];
delta_representatives := [
    A!1,
    A!(-theta),
    A!(1 - theta),
    A!(-1/2 - theta),
    A!(2/3 - theta),
    A!(3 - theta)
];

assert #S eq 6;
assert { AtoS(d) : d in delta_representatives } eq S;

print "Fake 2-Selmer set size =", #S;
print "Explicit representatives delta_i in A = Q[theta]/(f(theta)):";
for i in [1..#delta_representatives] do
    print "  i =", i, " label =", delta_names[i], " delta =", delta_representatives[i],
          " abstract class =", AtoS(delta_representatives[i]);
end for;
print "Verified: the six displayed delta_i represent the full fake 2-Selmer set.";

Section("2. Splitting field and elliptic-cover construction");

R<z> := PolynomialRing(Q);
K<s> := NumberField(z^2 - 21);
KX<X> := PolynomialRing(K);

Bplus := X^3 - ((11 + s)/2)*X^2 + ((5 + s)/2)*X + 1;
Bminus := X^3 - ((11 - s)/2)*X^2 + ((5 - s)/2)*X + 1;
assert Bplus * Bminus eq KX!B7;

LTHETA<THETA> := quo<KX | Bplus>;
j := hom< A -> LTHETA | THETA >;

print "K = Q(sqrt(21)) with s^2 = 21.";
print "B_+(X) =", Bplus;
print "B_-(X) =", Bminus;
print "Verified: B_+(X) B_-(X) = B_7(X).";

expected_gamma_values := [
    K!1,
    K!1,
    K!-1,
    K!(1/8*(-3*s - 14)),
    K!(1/27*(3*s + 14)),
    K!(-3*s - 14)
];

base_x_values := [
    K!0,
    K!0,
    K!1,
    K!(-1/2),
    K!(2/3),
    K!3
];

P1K := ProjectiveSpace(Q, 1);
all_chabauty_x_values := { P1Q | };

for i in [1..#delta_representatives] do
    Section(Sprintf("3.%o. Elliptic Chabauty for %o", i, delta_names[i]));

    delta := delta_representatives[i];
    gamma := Norm(j(delta));
    assert gamma eq expected_gamma_values[i];

    Ecover := HyperellipticCurve(gamma * Bplus);
    Pbase := Ecover![base_x_values[i], gamma, K!1];
    Eraw, EcoverToEraw := EllipticCurve(Ecover, Pbase);
    Emin, ErawToEmin := MinimalModel(Eraw);

    EcoverToP1 := map< Ecover -> P1K | [Ecover.1, Ecover.3] >;
    EminToP1 := Expand(Inverse(ErawToEmin) * Inverse(EcoverToEraw) * EcoverToP1);

    success, MW, MWtoE := PseudoMordellWeilGroup(Emin);
    assert success;

    V, Rchab := Chabauty(MWtoE, EminToP1 : IndexBound := 2);
    assert IsPowerOfTwo(Rchab);

    pi := Extend(EminToP1);
    x_images := { P1Q | P1PointToRationalP1(pi(MWtoE(v)), P1Q) : v in V };
    assert x_images subset known_x_values;
    all_chabauty_x_values join:= x_images;

    print "delta =", delta;
    print "gamma = Norm(j(delta)) =", gamma;
    print "Elliptic curve E_i/K after minimalization:";
    print Emin;
    print "Mordell-Weil subgroup returned by PseudoMordellWeilGroup:";
    print MW;
    print "Images of Mordell-Weil subgroup generators on E_i(K):";
    for j in [1..Ngens(MW)] do
        print j, "->", MWtoE(MW.j);
    end for;
    print "success =", success, "(finite odd index in E_i(K))";
    print "Chabauty returned R_i =", Rchab;
    print "R_i is a power of 2:", IsPowerOfTwo(Rchab);
    print "Therefore gcd([E_i(K):G_i], R_i) = 1.";
    print "Rational P^1 images found by elliptic Chabauty:";
    print x_images;
end for;

Section("4. Pullback to C_7 and final assertion");

assert all_chabauty_x_values eq known_x_values;

print "Union of rational X-values from all six elliptic Chabauty computations:";
print all_chabauty_x_values;
print "Verified: this union is exactly {infinity, 0, 1, -1/2, 2/3, 3}.";

pulled_back_points := { C | };
for xval in all_chabauty_x_values do
    fiber_points := RationalPoints(xval @@ CtoP1);
    print "Fiber over", xval, "=", fiber_points;
    pulled_back_points join:= { C | P : P in fiber_points };
end for;

assert StringSet(SetToSequence(pulled_back_points)) eq StringSet(known_points);

print "";
print "Final certified rational points:";
print pulled_back_points;
print "";
print "CERTIFIED RESULT:";
print "C_7(Q) consists exactly of the two points at infinity and the affine";
print "points above mu = 0, 1, -1/2, 2/3, 3.";

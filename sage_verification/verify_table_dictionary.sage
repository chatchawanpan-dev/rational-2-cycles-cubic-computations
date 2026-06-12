from sage.all import QQ, EllipticCurve


def row(a, n):
    a = QQ(a)
    E = EllipticCurve(QQ, [0, 4, 0, 0, 16*a^2])
    P = E(0, -4*a)
    Q = n*P
    X, Y = Q.xy()

    b = -X/4 - 2 - 16*a^2/X^2
    x_plus = (-4*a - Y)/(2*X)
    x_minus = (-4*a + Y)/(2*X)

    def f(x):
        return x^3 + b*x + a

    assert X != 0
    assert Y != 0
    assert Y^2 == X^3 + 4*X^2 + 16*a^2
    assert f(x_plus) == x_minus
    assert f(x_minus) == x_plus
    assert x_plus != x_minus

    return a, n, Q, b, x_plus, x_minus


rows = [
    (1, 2),
    (1, 3),
    (1, 4),
    (1, 5),
    (2, 2),
    (2, 3),
    (3, 2),
    (3, 3),
    (QQ(1)/2, 2),
    (QQ(1)/2, 3),
]

print("Verified table rows for Theorem 3.1 / Table 1")
print("a | n | [n]P_a | b | (x_plus, x_minus)")
for a, n in rows:
    print(row(a, n))

print("All displayed rows verified in exact rational arithmetic.")

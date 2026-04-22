/*
 * main.c  —  Driver for polynomials.c
 *
 * Exercises every function declared in polynomials.h.
 * Compile with:   make
 * Run with:       ./poly_demo
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "polynomials.h"

/* ---------------------------------------------------------------
 * Helpers
 * --------------------------------------------------------------- */
static void section(const char *title) {
    printf("\n╔══════════════════════════════════════════════════╗\n");
    printf("  %s\n", title);
    printf("╚══════════════════════════════════════════════════╝\n");
}

static void subsection(const char *title) {
    printf("\n── %s ──\n", title);
}

static void print_complex(const char *label, Complex c) {
    if (fabs(c.imaginary) < 1e-12)
        printf("  %-30s  %.6f\n", label, c.real);
    else
        printf("  %-30s  %.6f %s %.6fi\n", label, c.real,
               c.imaginary >= 0 ? "+" : "-", fabs(c.imaginary));
}

/* ---------------------------------------------------------------
 * 1. Utility
 * --------------------------------------------------------------- */
static void test_factorial(void) {
    section("1. Utility — factorial");
    int cases[] = {0, 1, 5, 10, 20, 21, -1};
    for (int i = 0; i < 7; i++) {
        long long r = factorial(cases[i]);
        printf("  factorial(%3d) = %lld\n", cases[i], r);
    }
}

/* ---------------------------------------------------------------
 * 2. Complex arithmetic
 * --------------------------------------------------------------- */
static void test_complex(void) {
    section("2. Complex arithmetic");

    Complex a = {3.0,  4.0};
    Complex b = {1.0, -2.0};

    printf("  a = %.1f + %.1fi\n", a.real, a.imaginary);
    printf("  b = %.1f + %.1fi\n", b.real, b.imaginary);

    print_complex("a + b",          add_complex(a, b));
    print_complex("a - b",          subtract_complex(a, b));
    print_complex("a * b",          multiply_complex(a, b));
    print_complex("a / b",          divide_complex(a, b));

    Complex zero = {0.0, 0.0};
    Complex q    = divide_complex(a, zero);
    printf("  %-30s  real=%.1f  imag=%.1f  (both NaN expected)\n",
           "a / 0", q.real, q.imaginary);
}

/* ---------------------------------------------------------------
 * 3. Polynomial basics
 * --------------------------------------------------------------- */
static void test_basics(void) {
    section("3. Polynomial basics");

    /* p(x) = 2x^3 - 3x^2 + 0x + 5 */
    double p[] = {2.0, -3.0, 0.0, 5.0};
    int    ps  = 4;

    subsection("print_polynomial");
    printf("  p(x) = ");
    print_polynomial(p, ps);
    printf("\n");

    subsection("degree");
    printf("  degree(p) = %d  (expected 3)\n", degree(p, ps));

    double leading_zero[] = {0.0, 0.0, 4.0, -1.0};
    printf("  degree([0,0,4,-1]) = %d  (expected 2)\n",
           degree(leading_zero, 4));

    subsection("evaluate_polynomial (Horner)");
    double xs[] = {0.0, 1.0, -1.0, 2.0};
    for (int i = 0; i < 4; i++)
        printf("  p(%.0f) = %.4f\n", xs[i], evaluate_polynomial(p, xs[i], ps));
}

/* ---------------------------------------------------------------
 * 4. Calculus
 * --------------------------------------------------------------- */
static void test_calculus(void) {
    section("4. Calculus");

    /* p(x) = x^3 - 6x^2 + 11x - 6  →  roots at 1, 2, 3 */
    double p[] = {1.0, -6.0, 11.0, -6.0};
    int    ps  = 4;

    subsection("derivative");
    double *d = derivative(p, ps);
    printf("  p'(x) = ");
    print_polynomial(d, ps - 1);
    printf("\n");

    subsection("integral");
    double *ig = integral(p, ps);
    printf("  ∫p dx = ");
    print_polynomial(ig, ps + 1);
    printf("\n");

    subsection("evaluate_derivative / evaluate_integral at x=2");
    printf("  p'(2) = %.6f  (exact = %.6f)\n",
           evaluate_derivative(p, ps, 2.0),
           evaluate_polynomial(d, 2.0, ps - 1));
    printf("  ∫p(2) = %.6f\n", evaluate_integral(p, ps, 2.0));

    subsection("numerical derivatives (h = 1e-5)");
    double h = 1e-5;
    printf("  central_difference  p'(2) = %.8f\n",
           central_difference(p, ps, 2.0, h));
    printf("  second_derivative   p''(2)= %.8f\n",
           second_derivative(p, ps, 2.0, h));
    printf("  richardson_extrap   p'(2) = %.8f\n",
           richardson_extrapolation(p, ps, 2.0, h));

    free(d);
    free(ig);
}

/* ---------------------------------------------------------------
 * 5. Polynomial arithmetic
 * --------------------------------------------------------------- */
static void test_arithmetic(void) {
    section("5. Polynomial arithmetic");

    /* p = x^2 + 2x + 1,  q = x + 1 */
    double p[] = {1.0, 2.0, 1.0};
    double q[] = {1.0, 1.0};
    int    ps  = 3, qs = 2;

    subsection("add / subtract");
    double *s = add_polynomials(p, q, ps, qs);
    printf("  p + q = ");  print_polynomial(s, ps); printf("\n");
    free(s);

    double *dif = subtract_polynomials(p, q, ps, qs);
    printf("  p - q = ");  print_polynomial(dif, ps); printf("\n");
    free(dif);

    subsection("multiply (direct convolution)");
    double *m = multiply_polynomials(p, q, ps, qs);
    printf("  p * q = ");  print_polynomial(m, ps + qs - 1); printf("\n");
    free(m);

    subsection("multiply (FFT)");
    double *mf = multiply_fft(p, ps, q, qs);
    printf("  p * q (FFT) = ");  print_polynomial(mf, ps + qs - 1); printf("\n");
    free(mf);

    subsection("divide");
    /* (x^2 + 2x + 1) / (x + 1) = (x + 1) */
    double *qt = divide_polynomials(p, ps, q, qs);
    if (qt) {
        printf("  p / q = ");  print_polynomial(qt, ps - qs + 1); printf("\n");
        free(qt);
    }

    subsection("remainder");
    int rem_sz = 0;
    double *rem = remainder(p, ps, q, qs, &rem_sz);
    if (rem) {
        printf("  p %% q = ");  print_polynomial(rem, rem_sz); printf("\n");
        free(rem);
    }

    subsection("point-evaluation shortcuts at x = 3");
    double x = 3.0;
    printf("  (p+q)(3) = %.4f\n", add(p, q, ps, qs, x));
    printf("  (p-q)(3) = %.4f\n", sub(p, q, ps, qs, x));
    printf("  (p*q)(3) = %.4f\n", multiply(p, q, ps, qs, x));
    printf("  (p/q)(3) = %.4f\n", divide(p, q, ps, qs, x));
    printf("  (p/q)(x=-1)= %.4f  (NaN expected — root of denominator)\n",
           divide(p, q, ps, qs, -1.0));
}

/* ---------------------------------------------------------------
 * 6. Root finding & bounds
 * --------------------------------------------------------------- */
static void test_roots(void) {
    section("6. Root finding & bounds");

    /* p(x) = x^3 - 6x^2 + 11x - 6  →  roots 1, 2, 3 */
    double p[] = {1.0, -6.0, 11.0, -6.0};
    int    ps  = 4;

    subsection("isRoot");
    printf("  isRoot(p, 1.0) = %s\n", isRoot(p, 1.0, ps) ? "true" : "false");
    printf("  isRoot(p, 1.5) = %s\n", isRoot(p, 1.5, ps) ? "true" : "false");

    subsection("Newton–Raphson");
    double r_newton = approximate_root_Newton(p, ps, 0.5);
    printf("  root near x0=0.5 → %.8f  (expected ≈ 1)\n", r_newton);
    r_newton = approximate_root_Newton(p, ps, 2.5);
    printf("  root near x0=2.5 → %.8f  (expected ≈ 3)\n", r_newton);

    subsection("Bisection");
    double r_bisect = bisection_root_method(p, ps, 0.0, 1.5, 1e-9);
    printf("  bisection [0, 1.5] → %.8f  (expected ≈ 1)\n", r_bisect);
    double r_bad = bisection_root_method(p, ps, 4.0, 5.0, 1e-9);
    printf("  bisection [4, 5]   → %.4f  (NaN expected — no sign change)\n", r_bad);

    subsection("Cauchy root bound");
    printf("  bound = %.6f  (all roots lie in [-bound, bound])\n",
           Cauchy_root_bound(p, ps));
}

/* ---------------------------------------------------------------
 * 7. Quadratic helpers
 * --------------------------------------------------------------- */
static void test_quadratic(void) {
    section("7. Quadratic helpers");

    subsection("real roots: x^2 - 5x + 6  (roots 2, 3)");
    double a=1, b=-5, c=6;
    printf("  evaluate at x=2: %.4f\n", evaluate_quadratic(a, b, c, 2.0));
    printf("  real_roots_quadratic: %s\n",
           real_roots_quadratic(a, b, c) ? "true" : "false");
    quadratic_roots rq = find_roots_quadratic(a, b, c);
    print_complex("root1", rq.root1);
    print_complex("root2", rq.root2);

    subsection("complex roots: x^2 + x + 1  (discriminant < 0)");
    quadratic_roots rc = find_roots_quadratic(1, 1, 1);
    print_complex("root1", rc.root1);
    print_complex("root2", rc.root2);

    subsection("minimum / maximum");
    printf("  min of x^2 - 4x + 5  = %.4f  (expected 1)\n",
           minimum_quadratic(1, -4, 5));
    printf("  max of -x^2 + 4x - 3 = %.4f  (expected 1)\n",
           maximum_quadratic(-1, 4, -3));
    printf("  min of -x^2 ...      = %.4f  (NaN expected — a<=0)\n",
           minimum_quadratic(-1, 4, -3));
}

/* ---------------------------------------------------------------
 * 8. Cubic helpers
 * --------------------------------------------------------------- */
static void test_cubic(void) {
    section("8. Cubic helpers");

    printf("  evaluate_cubic(1,-6,11,-6, x=2) = %.4f  (expected 0)\n",
           evaluate_cubic(1, -6, 11, -6, 2.0));

    printf("  is_monotonic_cubic(1,0,1,0): %s  (disc<0 → no extrema → monotone)\n",
           is_monotonic_cubic(1, 0, 1, 0) ? "true" : "false");

    printf("  real_roots_count_cubic(1,-6,11,-6) = %d  (expected 3)\n",
           real_roots_count_cubic(1, -6, 11, -6));
    printf("  real_roots_count_cubic(1,0,1,0)   = %d  (expected 1 — monotone)\n",
           real_roots_count_cubic(1, 0, 1, 0));
}

/* ---------------------------------------------------------------
 * 9. Numerical integration
 * --------------------------------------------------------------- */
static void test_integration(void) {
    section("9. Numerical integration  ∫₀¹ x² dx  (exact = 1/3)");

    /* p(x) = x^2 */
    double p[] = {1.0, 0.0, 0.0};
    int    ps  = 3;

    double simp = Simpson_method(p, ps, 0.0, 1.0);
    double trap  = Trapezoid_rule(p, ps, 0.0, 1.0);
    printf("  Simpson  = %.10f  (error = %.2e)\n", simp, fabs(simp - 1.0/3.0));
    printf("  Trapezoid= %.10f  (error = %.2e)\n", trap, fabs(trap - 1.0/3.0));
}

/* ---------------------------------------------------------------
 * 10. Lagrange interpolation
 * --------------------------------------------------------------- */
static void test_interpolation(void) {
    section("10. Lagrange interpolation");

    /* Interpolate sin at 0, π/6, π/2 */
    double xs[] = {0.0, PI/6.0, PI/2.0};
    double ys[] = {0.0, 0.5,    1.0};

    double val = Lagrange_interpolation(xs, ys, 2, PI/4.0);
    printf("  Lagrange sin(π/4) ≈ %.6f  (sin(π/4) = %.6f)\n",
           val, sin(PI/4.0));

    /* Duplicate x — should return NaN */
    double xdup[] = {1.0, 1.0, 2.0};
    double ydup[] = {1.0, 2.0, 3.0};
    printf("  Duplicate x:       %.4f  (NaN expected)\n",
           Lagrange_interpolation(xdup, ydup, 2, 1.5));
}

/* ---------------------------------------------------------------
 * 11. Taylor series
 * --------------------------------------------------------------- */
static void test_taylor(void) {
    section("11. Taylor series approximations  (k = 10 terms)");
    int k = 10;
    double x = 1.0;

    printf("  sin(1) approx = %.10f  exact = %.10f  err = %.2e\n",
           sine_approximation(x, k),  sin(x),  fabs(sine_approximation(x,k)-sin(x)));
    printf("  cos(1) approx = %.10f  exact = %.10f  err = %.2e\n",
           cosine_approximation(x, k), cos(x), fabs(cosine_approximation(x,k)-cos(x)));
}

/* ---------------------------------------------------------------
 * 12. Special polynomials
 * --------------------------------------------------------------- */
static void test_special(void) {
    section("12. Special polynomials");

    int    n    = 8;
    int   *found_c = calloc(n + 1, sizeof(int));
    double *memo_c  = calloc(n + 1, sizeof(double));
    int   *found_l = calloc(n + 1, sizeof(int));
    double *memo_l  = calloc(n + 1, sizeof(double));

    subsection("Chebyshev T_n(cos θ) = cos(nθ)  at x = 0.5");
    for (int i = 0; i <= n; i++)
        printf("  T_%d(0.5) = %9.6f  (cos(%d·π/3) = %9.6f)\n",
               i,
               implement_chebyshev(0.5, i, found_c, memo_c),
               i, cos(i * PI / 3.0));

    subsection("Legendre P_n at x = 0.5");
    for (int i = 0; i <= n; i++)
        printf("  P_%d(0.5) = %9.6f\n",
               i, implement_Legendre(0.5, i, found_l, memo_l));

    free(found_c); free(memo_c);
    free(found_l); free(memo_l);
}

/* ---------------------------------------------------------------
 * 13. Sturm sequences
 * --------------------------------------------------------------- */
static void test_sturm(void) {
    section("13. Sturm sequences  — x^3 - 6x^2 + 11x - 6  (roots 1, 2, 3)");

    double p[] = {1.0, -6.0, 11.0, -6.0};
    int    ps  = 4;

    int seq_count = 0;
    int *sizes    = NULL;
    double **seq  = sturm_sequences(p, ps, &seq_count, &sizes);

    printf("  Sturm chain has %d polynomials\n", seq_count);
    printf("  real roots in (-∞, +∞) : %d  (expected 3)\n",
           count_real_roots_sturm(seq, sizes, seq_count, -1e9, 1e9));
    printf("  real roots in (0, 1.5) : %d  (expected 1)\n",
           count_real_roots_sturm(seq, sizes, seq_count, 0.0,  1.5));
    printf("  real roots in (1.5, 4) : %d  (expected 2)\n",
           count_real_roots_sturm(seq, sizes, seq_count, 1.5,  4.0));

    sturm_free(seq, sizes, seq_count);
}

/* ---------------------------------------------------------------
 * 14. Companion matrix
 * --------------------------------------------------------------- */
static void test_companion(void) {
    section("14. Companion matrix  — x^3 - 6x^2 + 11x - 6");

    double p[] = {1.0, -6.0, 11.0, -6.0};
    int    ps  = 4;
    int    deg = degree(p, ps);

    double **C = companion_matrix(p, ps);
    if (!C) { printf("  NULL returned (leading coeff == 0)\n"); return; }

    printf("  %dx%d Frobenius companion matrix:\n", deg, deg);
    for (int i = 0; i < deg; i++) {
        printf("  [");
        for (int j = 0; j < deg; j++)
            printf(" %7.3f", C[i][j]);
        printf(" ]\n");
        free(C[i]);
    }
    free(C);
}

/* ---------------------------------------------------------------
 * 15. FFT (exposed for testing)
 * --------------------------------------------------------------- */
static void test_fft(void) {
    section("15. FFT — round-trip test (size 4)");

    /* x = [1, 2, 3, 4] */
    int     sz = 4;
    Complex a[] = {{1,0},{2,0},{3,0},{4,0}};

    printf("  Input:        ");
    for (int i = 0; i < sz; i++) printf("%.1f ", a[i].real);
    printf("\n");

    fft(a, sz, false);
    printf("  After FFT:    ");
    for (int i = 0; i < sz; i++) printf("(%.2f,%.2f) ", a[i].real, a[i].imaginary);
    printf("\n");

    fft(a, sz, true);
    printf("  After IFFT:   ");
    for (int i = 0; i < sz; i++) printf("%.3f ", a[i].real);
    printf("  (should match input)\n");
}

/* ---------------------------------------------------------------
 * main
 * --------------------------------------------------------------- */
int main(void) {
    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║          polynomials.c  —  full demo driver          ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n");

    test_factorial();
    test_complex();
    test_basics();
    test_calculus();
    test_arithmetic();
    test_roots();
    test_quadratic();
    test_cubic();
    test_integration();
    test_interpolation();
    test_taylor();
    test_special();
    test_sturm();
    test_companion();
    test_fft();

    printf("\n✓ All demos complete.\n\n");
    return 0;
}

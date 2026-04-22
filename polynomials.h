#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H

#include <math.h>
#include <stdbool.h>

#define PI 3.14159265358979323846

/* ============================================================
 *  Structs
 * ============================================================ */

typedef struct {
    double real;
    double imaginary;
} Complex;

typedef struct {
    Complex root1;
    Complex root2;
} quadratic_roots;

/* ============================================================
 *  Utility
 * ============================================================ */

/* Returns n!, or -1 if n < 0 or n > 20 (overflow guard) */
long long factorial(int n);

/* ============================================================
 *  Complex arithmetic
 * ============================================================ */

Complex add_complex      (Complex a, Complex b);
Complex subtract_complex (Complex a, Complex b);
Complex multiply_complex (Complex a, Complex b);

/* Returns {NAN, NAN} if b == 0 */
Complex divide_complex   (Complex a, Complex b);

/* ============================================================
 *  Polynomial basics
 *  Polynomials are stored as coefficient arrays in descending
 *  order of degree: coeff[0] is the leading coefficient.
 *  'size' is the number of elements in the array.
 * ============================================================ */

/* Returns the degree of the polynomial (skips leading zeros) */
int    degree             (double *coeff, int size);

/* Evaluates the polynomial at x using Horner's method */
double evaluate_polynomial(double *coeff, double x, int size);

/* Prints the polynomial in human-readable form, e.g. "3.00x^2 + x - 1.00" */
void   print_polynomial   (double *coeff, int size);

/* ============================================================
 *  Calculus
 * ============================================================ */

/* Returns a heap-allocated array of size (size-1) holding the derivative.
 * Returns a single-element zero array for a constant polynomial.
 * Caller must free(). */
double *derivative        (double *coeff, int size);

/* Returns a heap-allocated array of size (size+1) holding the integral.
 * The constant of integration is set to 0. Caller must free(). */
double *integral          (double *coeff, int size);

/* Evaluates the derivative of the polynomial at x */
double  evaluate_derivative(double *coeff, int size, double x);

/* Evaluates the (indefinite) integral of the polynomial at x */
double  evaluate_integral  (double *coeff, int size, double x);

/* Numerical first derivative using the central-difference formula */
double  central_difference (double *coeff, int size, double x, double h);

/* Numerical second derivative using the central-difference formula */
double  second_derivative  (double *coeff, int size, double x, double h);

/* Richardson extrapolation for a more accurate numerical first derivative */
double  richardson_extrapolation(double *coeff, int size, double x, double h);

/* ============================================================
 *  Polynomial arithmetic  (return heap-allocated arrays)
 * ============================================================ */

/* Returns coeff1 + coeff2 as a new array of size max(size1, size2).
 * Caller must free(). */
double *add_polynomials     (double *coeff1, double *coeff2, int size1, int size2);

/* Returns coeff1 - coeff2 as a new array of size max(size1, size2).
 * Caller must free(). */
double *subtract_polynomials(double *coeff1, double *coeff2, int size1, int size2);

/* Returns coeff1 * coeff2 as a new array of size (size1 + size2 - 1).
 * Uses direct convolution. Caller must free(). */
double *multiply_polynomials(double *coeff1, double *coeff2, int size1, int size2);

/* Returns coeff1 * coeff2 using FFT, size (size1 + size2 - 1).
 * Caller must free(). */
double *multiply_fft        (double *coeff1, int size1, double *coeff2, int size2);

/* Returns the quotient of dividend / divisor via polynomial long division.
 * Returns NULL if degree(divisor) > degree(dividend). Caller must free(). */
double *divide_polynomials  (double *dividend, int size1, double *divisor, int size2);

/* Returns the remainder of dividend / divisor.
 * Sets *rem_size to the length of the returned array. Caller must free(). */
double *remainder           (double *coeff1, int size1, double *coeff2, int size2, int *rem_size);

/* ============================================================
 *  Polynomial evaluation shortcuts (evaluate at a point x)
 * ============================================================ */

double add     (double *coeff1, double *coeff2, int size1, int size2, double x);
double sub     (double *coeff1, double *coeff2, int size1, int size2, double x);
double multiply(double *coeff1, double *coeff2, int size1, int size2, double x);

/* Returns NAN if coeff2(x) == 0 */
double divide  (double *coeff1, double *coeff2, int size1, int size2, double x);

/* ============================================================
 *  Root finding & bounds
 * ============================================================ */

/* Returns true if |p(x)| < 1e-9 */
bool   isRoot                (double *coeff, double x, int size);

/* Newton–Raphson root finder starting from x0.
 * Returns NAN if derivative is too small (flat region). */
double approximate_root_Newton(double *coeff, int size, double x0);

/* Bisection method on [a, b]. Returns NAN if f(a)*f(b) > 0. */
double bisection_root_method  (double *coeff, int size, double a, double b, double error);

/* Cauchy bound: all real roots lie in [-bound, bound].
 * Returns NAN if leading coefficient is 0. */
double Cauchy_root_bound      (double *coeff, int size);

/* ============================================================
 *  Quadratic helpers
 * ============================================================ */

double          evaluate_quadratic   (double a, double b, double c, double x);

/* Returns true iff the discriminant is strictly positive */
bool            real_roots_quadratic (double a, double b, double c);

/* Returns both roots (real or complex) via the quadratic formula */
quadratic_roots find_roots_quadratic (double a, double b, double c);

/* Returns the minimum value of ax²+bx+c. Returns NAN if a <= 0. */
double          minimum_quadratic    (double a, double b, double c);

/* Returns the maximum value of ax²+bx+c. Returns NAN if a >= 0. */
double          maximum_quadratic    (double a, double b, double c);

/* ============================================================
 *  Cubic helpers
 * ============================================================ */

double evaluate_cubic         (double a, double b, double c, double d, double x);

/* Returns true iff the cubic has no local extrema (monotone everywhere) */
bool   is_monotonic_cubic     (double a, double b, double c, double d);

/* Returns 1, 2, or 3: the number of distinct real roots of the cubic */
int    real_roots_count_cubic (double a, double b, double c, double d);

/* ============================================================
 *  Numerical integration
 * ============================================================ */

/* Composite Simpson's rule with N = 1,000,000 sub-intervals */
double Simpson_method  (double *coeff, int size, double lower_limit, double upper_limit);

/* Composite trapezoid rule with N = 1,000,000 sub-intervals */
double Trapezoid_rule  (double *coeff, int size, double lower_limit, double upper_limit);

/* ============================================================
 *  Interpolation
 * ============================================================ */

/* Lagrange interpolation through n+1 points (x[0..n], y[0..n]).
 * Returns NAN if any two x-values are equal. */
double Lagrange_interpolation(double *x, double *y, int n, double value);

/* ============================================================
 *  Taylor series approximations
 * ============================================================ */

/* Approximates sin(x) using k terms of the Taylor series */
double sine_approximation        (double x, int k);

/* Approximates cos(x) using k terms of the Taylor series */
double cosine_approximation      (double x, int k);

/* Approximates e^x using k terms of the Taylor series */
double exponential_approximation (double x, int k);

/* ============================================================
 *  Special polynomials  (memoised top-down recurrence)
 *  Both functions require caller-allocated arrays:
 *    found[0..n] zeroed (int)
 *    memo[0..n]  zeroed (double)
 * ============================================================ */

/* Evaluates the nth Chebyshev polynomial T_n(x) */
double implement_chebyshev(double x, int n, int *found, double *memo);

/* Evaluates the nth Legendre polynomial P_n(x) */
double implement_Legendre (double x, int n, int *found, double *memo);

/* ============================================================
 *  Sturm sequences  (for counting real roots in an interval)
 * ============================================================ */

/* Builds the Sturm sequence for the polynomial.
 * Returns a heap-allocated array of polynomial arrays.
 * *sequence_count is set to the number of polynomials.
 * *sizes_out is set to a heap-allocated array of their sizes.
 * Free with sturm_free(). */
double **sturm_sequences     (double *coeff, int size, int *sequence_count, int **sizes_out);

/* Frees memory allocated by sturm_sequences() */
void     sturm_free          (double **seq, int *sizes, int sequence_count);

/* Counts sign changes in the Sturm sequence at a single point x */
int      count_signs_sturm   (double **seq, int *sizes, int sequence_count, double x);

/* Returns the number of distinct real roots in the open interval (a, b) */
int      count_real_roots_sturm(double **seq, int *sizes, int sequence_count, double a, double b);

/* ============================================================
 *  Companion matrix
 * ============================================================ */

/* Returns the deg×deg Frobenius companion matrix as a heap-allocated
 * array of row pointers. Returns NULL if the leading coefficient is 0.
 * Free each row then the array itself. */
double **companion_matrix(double *coeff, int size);

/* ============================================================
 *  FFT (internal — exposed for testing)
 * ============================================================ */

/* In-place Cooley–Tukey FFT. size must be a power of 2.
 * Set invert=true for the inverse transform. */
void fft(Complex *a, int size, bool invert);

#endif /* POLYNOMIALS_H */

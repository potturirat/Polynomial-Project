# polynomials.c

A comprehensive single-file C library for polynomial mathematics. Covers everything from basic arithmetic and calculus to numerical root-finding, orthogonal polynomials, FFT-based multiplication, and Sturm sequences.

---

## Table of Contents

- [Requirements](#requirements)
- [Building](#building)
- [Data Representation](#data-representation)
- [API Reference](#api-reference)
  - [Complex Arithmetic](#complex-arithmetic)
  - [Basic Polynomial Operations](#basic-polynomial-operations)
  - [Calculus](#calculus)
  - [Polynomial Arithmetic](#polynomial-arithmetic)
  - [Root Finding](#root-finding)
  - [Quadratic & Cubic Utilities](#quadratic--cubic-utilities)
  - [Numerical Integration](#numerical-integration)
  - [Numerical Differentiation](#numerical-differentiation)
  - [Interpolation & Approximation](#interpolation--approximation)
  - [Orthogonal Polynomials](#orthogonal-polynomials)
  - [Sturm Sequences](#sturm-sequences)
  - [Companion Matrix](#companion-matrix)
  - [Series Approximations](#series-approximations)
  - [Utility](#utility)
- [Usage Examples](#usage-examples)
- [Notes & Limitations](#notes--limitations)

---

## Requirements

- C99 or later (`-std=c99`)
- Standard math library (`-lm`)
- No external dependencies

## Building

Include `polynomials.c` directly in your project, or compile it alongside your driver:

```bash
gcc -std=c99 -O2 main.c polynomials.c -lm -o poly_program
```

If compiling as a standalone object:

```bash
gcc -std=c99 -c polynomials.c -o polynomials.o -lm
```

---

## Data Representation

Polynomials are represented as **arrays of `double`** in **descending order of degree**:

```
p(x) = a_n * x^n + a_{n-1} * x^{n-1} + ... + a_1 * x + a_0
coeff[] = { a_n, a_{n-1}, ..., a_1, a_0 }
size    = n + 1
```

**Example:** `3x³ - 2x + 5` → `coeff = {3.0, 0.0, -2.0, 5.0}`, `size = 4`

The caller is responsible for freeing any heap-allocated arrays returned by the library.

---

## API Reference

### Complex Arithmetic

The library defines a `Complex` struct and four arithmetic operations, used internally for FFT-based polynomial multiplication.

```c
typedef struct { double real; double imaginary; } Complex;

Complex add_complex(Complex a, Complex b);
Complex subtract_complex(Complex a, Complex b);
Complex multiply_complex(Complex a, Complex b);
Complex divide_complex(Complex a, Complex b);   // returns {NAN, NAN} on divide-by-zero
```

---

### Basic Polynomial Operations

```c
int degree(double *coeff, int size);
```
Returns the degree of the polynomial, ignoring leading zero coefficients.

```c
double evaluate_polynomial(double *coeff, double x, int size);
```
Evaluates p(x) using Horner's method. O(n).

```c
void print_polynomial(double *coeff, int size);
```
Prints the polynomial in human-readable form, e.g. `3.00x^3 - 2.00x + 5.00`.

```c
bool isRoot(double *coeff, double x, int size);
```
Returns `true` if `|p(x)| < 1e-9`.

---

### Calculus

```c
double *derivative(double *coeff, int size);
```
Returns a newly allocated array representing p′(x). Caller must `free()` the result.

```c
double *integral(double *coeff, int size);
```
Returns a newly allocated array representing ∫p(x)dx (constant of integration = 0). Caller must `free()` the result.

```c
double evaluate_derivative(double *coeff, int size, double x);
```
Evaluates p′(x) at a point directly, without allocating.

```c
double evaluate_integral(double *coeff, int size, double x);
```
Evaluates ∫p(x)dx at a point directly, without allocating.

---

### Polynomial Arithmetic

```c
double *add_polynomials(double *coeff1, double *coeff2, int size1, int size2);
double *subtract_polynomials(double *coeff1, double *coeff2, int size1, int size2);
```
Return newly allocated arrays for p₁ ± p₂. Caller must `free()` the result.

```c
double add(double *coeff1, double *coeff2, int size1, int size2, double x);
double sub(double *coeff1, double *coeff2, int size1, int size2, double x);
double multiply(double *coeff1, double *coeff2, int size1, int size2, double x);
double divide(double *coeff1, double *coeff2, int size1, int size2, double x);
```
Evaluate the sum, difference, product, or quotient of two polynomials **at a single point x**. No memory allocation.

```c
double *multiply_polynomials(double *coeff1, double *coeff2, int size1, int size2);
```
Returns a newly allocated coefficient array for p₁ × p₂ via direct convolution. O(n²). Caller must `free()`.

```c
double *multiply_fft(double *coeff1, int size1, double *coeff2, int size2);
```
Returns a newly allocated coefficient array for p₁ × p₂ using the **Fast Fourier Transform**. O(n log n). Preferred for large-degree polynomials. Caller must `free()`.

```c
double *divide_polynomials(double *dividend, int size1, double *divisor, int size2);
```
Returns the **quotient** polynomial array from polynomial long division. Returns `NULL` if the divisor has higher degree. Caller must `free()`.

```c
double *remainder(double *coeff1, int size1, double *coeff2, int size2, int *rem_size);
```
Returns the **remainder** polynomial from dividing coeff1 by coeff2. `rem_size` is set to the length of the returned array. Caller must `free()`.

---

### Root Finding

```c
double approximate_root_Newton(double *coeff, int size, double x0);
```
Newton-Raphson method starting from initial guess `x0`. Runs up to 1000 iterations with tolerance `1e-7`. Returns `NAN` if the derivative vanishes.

```c
double bisection_root_method(double *coeff, int size, double a, double b, double error);
```
Bisection method on the interval `[a, b]`. Requires `p(a)` and `p(b)` to have opposite signs (returns `NAN` otherwise). Terminates when the interval width or function value is within `error`.

```c
double Cauchy_root_bound(double *coeff, int size);
```
Returns the **Cauchy bound** — an upper bound on the absolute value of all real roots:
`bound = 1 + max(|aᵢ|) / |aₙ|`.

---

### Quadratic & Cubic Utilities

```c
typedef struct { Complex root1; Complex root2; } quadratic_roots;

double evaluate_quadratic(double a, double b, double c, double x);
bool real_roots_quadratic(double a, double b, double c);
quadratic_roots find_roots_quadratic(double a, double b, double c);
double minimum_quadratic(double a, double b, double c);   // valid only when a > 0
double maximum_quadratic(double a, double b, double c);   // valid only when a < 0
```

```c
double evaluate_cubic(double a, double b, double c, double d, double x);
bool is_monotonic_cubic(double a, double b, double c, double d);
int real_roots_count_cubic(double a, double b, double c, double d);
```
Counts the number of distinct real roots of a cubic using its discriminant.

---

### Numerical Integration

```c
double Simpson_method(double *coeff, int size, double lower_limit, double upper_limit);
```
Composite Simpson's rule with N = 1,000,000 subintervals.

```c
double Trapezoid_rule(double *coeff, int size, double lower_limit, double upper_limit);
```
Composite Trapezoid rule with N = 1,000,000 subintervals.

Both return the definite integral ∫[lower, upper] p(x) dx.

---

### Numerical Differentiation

```c
double central_difference(double *coeff, int size, double x, double h);
```
Approximates p′(x) using the central difference formula: `[p(x+h) - p(x-h)] / (2h)`.

```c
double second_derivative(double *coeff, int size, double x, double h);
```
Approximates p″(x) using: `[p(x+h) - 2p(x) + p(x-h)] / h²`.

```c
double richardson_extrapolation(double *coeff, int size, double x, double h);
```
Higher-accuracy first derivative using Richardson extrapolation:
`[4·D(h/2) - D(h)] / 3`, where D(h) is the central difference estimate.

---

### Interpolation & Approximation

```c
double Lagrange_interpolation(double *x, double *y, int n, double value);
```
Evaluates the degree-n Lagrange interpolating polynomial through the points `(x[0],y[0]), ..., (x[n],y[n])` at `value`. Returns `NAN` if any two x-values coincide. Note: `n` here is the **number of intervals** (i.e., n+1 points).

---

### Orthogonal Polynomials

Both functions use **memoised top-down recursion** (O(n) time) via caller-supplied arrays.

```c
double implement_chebyshev(double x, int n, int *found, double *memo);
```
Evaluates the Chebyshev polynomial Tₙ(x) using the recurrence Tₙ(x) = 2x·Tₙ₋₁(x) − Tₙ₋₂(x).

```c
double implement_Legendre(double x, int n, int *found, double *memo);
```
Evaluates the Legendre polynomial Pₙ(x) using Bonnet's recurrence.

**Usage pattern:** Allocate and zero-initialise both arrays before the first call:
```c
int    found[N+1]; memset(found, 0, sizeof(found));
double memo[N+1];  memset(memo,  0, sizeof(memo));
double val = implement_chebyshev(x, N, found, memo);
```

---

### Sturm Sequences

Used to count the exact number of distinct real roots in an interval, without needing any initial guesses.

```c
double **sturm_sequences(double *coeff, int size, int *sequence_count, int **sizes_out);
```
Builds the Sturm chain for the polynomial. Returns a 2D array of polynomials. Sets `*sequence_count` to the number of polynomials in the chain, and `*sizes_out` to their respective sizes.

```c
int count_real_roots_sturm(double **seq, int *sizes, int sequence_count, double a, double b);
```
Returns the number of distinct real roots in the open interval (a, b) using Sturm's theorem.

```c
int count_signs_sturm(double **seq, int *sizes, int sequence_count, double x);
```
Counts sign changes in the Sturm sequence evaluated at x. (Used internally by `count_real_roots_sturm`.)

```c
void sturm_free(double **seq, int *sizes, int sequence_count);
```
Frees all memory allocated by `sturm_sequences`. **Always call this after use.**

**Typical usage:**
```c
int seq_count, *sizes;
double **seq = sturm_sequences(coeff, size, &seq_count, &sizes);
int roots = count_real_roots_sturm(seq, sizes, seq_count, -100.0, 100.0);
sturm_free(seq, sizes, seq_count);
```

---

### Companion Matrix

```c
double **companion_matrix(double *coeff, int size);
```
Builds the (size-1) × (size-1) companion matrix of the polynomial. The eigenvalues of this matrix are the roots of the polynomial. Returns `NULL` if the leading coefficient is zero. Each row is a separately allocated `double*`; the caller must free all rows and then the outer array.

---

### Series Approximations

These use Taylor series truncated to k terms.

```c
double sine_approximation(double x, int k);
double cosine_approximation(double x, int k);
double exponential_approximation(double x, int k);
```

Uses the `factorial` helper, which supports integers 0–20 (returns -1 for n > 20).

---

### Utility

```c
long long factorial(int n);
```
Returns n! for 0 ≤ n ≤ 20. Returns -1 for negative input or n > 20 (overflow guard).

---

## Usage Examples

### Evaluate and differentiate a polynomial

```c
// p(x) = 2x^3 - x + 4
double coeff[] = {2.0, 0.0, -1.0, 4.0};
int size = 4;

double val = evaluate_polynomial(coeff, 3.0, size);  // p(3)
printf("p(3) = %.4f\n", val);  // 55.0

double *der = derivative(coeff, size);
double dval = evaluate_polynomial(der, 3.0, size - 1); // p'(3)
printf("p'(3) = %.4f\n", dval);  // 53.0
free(der);
```

### Find a root with Newton-Raphson

```c
// p(x) = x^2 - 2  (root near 1.414...)
double coeff[] = {1.0, 0.0, -2.0};
double root = approximate_root_Newton(coeff, 3, 1.0);
printf("root ≈ %.10f\n", root);  // 1.4142135624
```

### Count real roots using Sturm's theorem

```c
// p(x) = x^3 - 3x + 1  (3 real roots)
double coeff[] = {1.0, 0.0, -3.0, 1.0};
int seq_count, *sizes;
double **seq = sturm_sequences(coeff, 4, &seq_count, &sizes);
int n_roots = count_real_roots_sturm(seq, sizes, seq_count, -10.0, 10.0);
printf("Number of real roots: %d\n", n_roots);  // 3
sturm_free(seq, sizes, seq_count);
```

### Multiply polynomials with FFT

```c
// (x + 1)(x - 1) = x^2 - 1
double p1[] = {1.0,  1.0};
double p2[] = {1.0, -1.0};
double *product = multiply_fft(p1, 2, p2, 2);
print_polynomial(product, 3);  // 1.00x^2 + 0.00x - 1.00
free(product);
```

### Evaluate Legendre polynomial

```c
int found[6]; memset(found, 0, sizeof(found));
double memo[6]; memset(memo, 0, sizeof(memo));
// P5(0.5)
double val = implement_Legendre(0.5, 5, found, memo);
printf("P5(0.5) = %.6f\n", val);  // 0.089844
```

---

## Notes & Limitations

- **Memory management:** All functions returning `double*` or `double**` allocate on the heap. The caller must `free()` these, or use `sturm_free()` for Sturm chains and manually free companion matrix rows.
- **Coefficient convention:** Arrays are always in **descending** degree order. Passing arrays in ascending order will silently produce wrong results.
- **Bisection precondition:** `bisection_root_method` requires `p(a)·p(b) < 0`. It returns `NAN` if this is not satisfied — always check the return value.
- **Factorial overflow:** `factorial` is capped at n = 20. Series approximation functions (`sine_approximation`, etc.) may silently degrade for large k due to this.
- **FFT precision:** `multiply_fft` rounds small coefficients in the result to zero (`< 1e-9`), which prevents noise from floating-point arithmetic but may affect polynomials with legitimately tiny coefficients.
- **Sturm sequences and repeated roots:** Sturm's theorem counts **distinct** real roots. A root of multiplicity m is counted once.
- **Companion matrix:** The caller must free each row `result[i]` individually, then free the outer `result` pointer. The matrix is not freed automatically.
- **Thread safety:** All functions are stateless except for memoised orthogonal polynomial functions, which require caller-managed state arrays (`found`, `memo`). These are safe provided each thread uses its own arrays.

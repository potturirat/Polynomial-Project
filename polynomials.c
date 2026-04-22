#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <stdbool.h>

#define PI 3.14159265358979323846 

long long factorial(int n){
    if(n<0) return -1; 
    if(n==0) return 1;
    if(n==1) return 1;
    if(n > 20) return -1;
    return n*factorial(n-1);
}

typedef struct{
    double real;
    double imaginary;
}Complex;

Complex add_complex(Complex a,Complex b){
    Complex sum;
    sum.real=a.real+b.real;
    sum.imaginary=a.imaginary+b.imaginary;
    return sum;
}

Complex subtract_complex(Complex a, Complex b){
    Complex diff;
    diff.real=a.real-b.real;
    diff.imaginary=a.imaginary-b.imaginary;
    return diff;
}

Complex multiply_complex(Complex a, Complex b){
    Complex prod;
    prod.real=a.real*b.real - a.imaginary*b.imaginary;
    prod.imaginary=a.imaginary*b.real+b.imaginary*a.real;
    return prod;
}

Complex divide_complex(Complex a, Complex b){
    Complex quotient;
    double sq=b.real*b.real+b.imaginary*b.imaginary;
    if(sq == 0.0){ 
        Complex nan = {NAN, NAN}; 
        return nan; 
    }
    else{
    quotient.real=(a.real*b.real - a.imaginary*b.imaginary)/sq;
    quotient.imaginary=(a.imaginary*b.real - b.imaginary*a.real)/sq;
    return quotient;
    }
}

typedef struct{
    Complex root1;
    Complex root2;
}quadratic_roots;

int degree(double *coeff, int size){
    for(int i = 0; i < size; i++){
        if(coeff[i] != 0) {
        return size - 1 - i;
        }
    }
    return 0;
}

double evaluate_polynomial(double *coeff, double x, int size){ //evaluate value of polynomial at a given point 
    double polynomial = 0.0;
    for(int i=0;i<size;i++){
        polynomial = coeff[i] + x*polynomial;
   }
   return polynomial;
}

double* derivative(double *coeff,int size){ // returns the polynomial array of the derivative of the input
    if(size <= 1){
    double *zero = calloc(1, sizeof(double));
    return zero;  // derivative of constant is 0
}
    double *der=malloc((size-1) * sizeof(double));
    for(int i=0;i<size-1;i++){
        der[i]=(size-1-i)*coeff[i];
    }
    return der;
}

double* integral(double *coeff, int size){ // gives the polynomial array of integral of the input
    double *integ=malloc((size+1) * sizeof(double));
    for(int i=0;i<size;i++){
        integ[i]=coeff[i]/(size-i);
    }
    integ[size]=0.0; //a random integral constant
    return integ;
}

double evaluate_derivative(double *coeff, int size, double x){//evaluate derivative at a given point
    double *der=derivative(coeff, size);
    int der_size = (size > 1) ? size - 1 : 1;
    double result = evaluate_polynomial(der, x, der_size); 
    free(der); 
    return result;  
}

double evaluate_integral(double *coeff, int size, double x){ // evaluate integral at a given point
    double *integ=integral(coeff, size);
    double result =evaluate_polynomial(integ,x,size+1);
    free(integ);
    return result;
    
}

double* add_polynomials(double *coeff1, double *coeff2, int size1, int size2){
    if(size1 >= size2){
        double *add = malloc(size1 * sizeof(double));
        for(int i = 0; i < size2; i++){
            add[(size1 - size2) + i] = coeff1[(size1 - size2) + i] + coeff2[i];
        }
        for(int i = 0; i < size1 - size2; i++){
            add[i] = coeff1[i];
        }
        return add;
    }
    else{
    double *add = calloc(size2, sizeof(double));
    int diff = size2 - size1;
    for(int i = 0; i < size2; i++) add[i] = coeff2[i];
    for(int i = 0; i < size1; i++) add[diff + i] += coeff1[i];
    return add;
    }
}

double* subtract_polynomials(double *coeff1, double *coeff2, int size1, int size2){ //gives the difference array
    if(size1>=size2){
        double *sub=malloc(size1 * sizeof(double));
        for(int i=size1-1;i>=size1-size2+1;i--){
            sub[i]=coeff1[i];
        }
        for(int i = 0; i < size2; i++){
        sub[(size1 - size2) + i] = coeff1[(size1 - size2) + i] - coeff2[i];
        }
        for(int i = 0; i < size1 - size2; i++){
        sub[i] = coeff1[i];
        }
        return sub;
    }
    else{
    double *sub = calloc(size2, sizeof(double));
    int diff = size2 - size1;
    for(int i = 0; i < size2; i++) sub[i] = -coeff2[i];
    for(int i = 0; i < size1; i++) sub[diff + i] += coeff1[i];
    return sub;
    }
}

double add(double *coeff1, double *coeff2, int size1, int size2, double x){//evaluates sum of two polynomials
    double *arr = add_polynomials(coeff1, coeff2, size1, size2);
    int sz = (size1 >= size2) ? size1 : size2;
    double result = evaluate_polynomial(arr, x, sz);
    free(arr);
    return result;
}

double sub(double *coeff1, double *coeff2, int size1, int size2, double x){//evaluates difference of two polynomials
    double *arr = subtract_polynomials(coeff1, coeff2, size1, size2);
    int sz = (size1 >= size2) ? size1 : size2;
    double result = evaluate_polynomial(arr, x, sz);
    free(arr);
    return result;
}

bool isRoot(double *coeff,double x,int size){ //check whether x is a root or not
    if(fabs(evaluate_polynomial(coeff,x,size))<1e-9){
        return true;
    }
    else{
        return false;
    }
}

double multiply(double *coeff1, double *coeff2, int size1, int size2, double x){ //multiplies two polynomials at a given value
    return (evaluate_polynomial(coeff1,x,size1))*(evaluate_polynomial(coeff2,x,size2));
}

double divide(double *coeff1, double *coeff2, int size1, int size2, double x){//divides two polynomials at a given value
    if(isRoot(coeff2,x,size2)){
        return NAN;
    }
    else{
        return (evaluate_polynomial(coeff1,x,size1))/(evaluate_polynomial(coeff2,x,size2));
    }
}

double* multiply_polynomials(double *coeff1, double *coeff2,int size1, int size2){//convolution logic
    double *res=calloc(size1+size2-1, sizeof(double));
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            res[i+j]+=coeff1[i]*coeff2[j];
        }
    }
    return res;
} 

void fft(Complex *a,int size, bool invert){ // Fast Fourier Transformation implementation
    if(size==1) return;
    Complex *a0=(Complex*)malloc(size/2*sizeof(Complex));
    Complex *a1=(Complex*)malloc(size/2*sizeof(Complex));
    for(int i=0;i<size/2;i++){
            a0[i]=a[2*i];
            a1[i]=a[2*i+1];
    }
    fft(a0,size/2,invert);
    fft(a1,size/2,invert);
    double angle;
    if(invert==true){
        angle= 2* PI/(size*(-1));
    }
    else{
        angle=2*PI/(size);
    }
    Complex w = {1.0,0.0};
    Complex wn={cos(angle), sin(angle)};
    for(int i=0;i<size/2;i++){
        Complex prod = multiply_complex(w, a1[i]);
        a[i] = add_complex(a0[i], prod);
        a[i+size/2] = subtract_complex(a0[i], prod);
        w = multiply_complex(wn, w);
    }
    
    if(invert){
        for(int i = 0; i < size; i++){
            a[i].real/= size;
            a[i].imaginary/= size;
        }
    }
    free(a0);
    free(a1);
}

double* multiply_fft(double *coeff1, int size1, double *coeff2, int size2){ //multiplying polynomials using FFT
    long long size = 1;
    while(size<size1+size2){
        size=size<<1;
    }
    Complex *fft1=(Complex*)malloc(size*sizeof(Complex));
    Complex *fft2=(Complex*)malloc(size*sizeof(Complex));
    for(int i=0;i<size;i++){
        fft1[i]=(Complex){0.0,0.0};
        fft2[i]=(Complex){0.0,0.0};
    }
    for(int i=0;i<size1;i++){
        fft1[i].real=coeff1[i];
    }
    for(int i=0;i<size2;i++){
        fft2[i].real=coeff2[i];
    }
    fft(fft1,size,false);
    fft(fft2,size,false);
    for(int i=0;i<size;i++){
        fft1[i]= multiply_complex(fft1[i], fft2[i]);
    }
    fft(fft1, size, true);
    free(fft2);
    int result_size = size1 + size2 - 1;
    double *result = malloc(result_size * sizeof(double));
    for(int i = 0; i < result_size; i++){
    result[i] = fft1[i].real;
    }
    free(fft1);
    return result;
}

double* divide_polynomials(double *dividend, int size1, double *divisor, int size2){
    if(size2>size1){
        return NULL;
    }
    int deg_quotient = size1-size2;
    double *quotient = (double*)malloc((deg_quotient+1)*sizeof(double));
    double *temp = (double*)malloc((size1)*sizeof(double));
    for(int i=0;i<size1;i++){
        temp[i]=dividend[i];
    }
    for(int i=0;i<=deg_quotient;i++){
        double num = temp[i]/divisor[0];
        quotient[i] = num;
        for(int j=0;j<size2;j++){
            temp[i+j] = temp[i+j] - num*divisor[j];
        }
    }
    free(temp);
    return quotient;
}

void print_polynomial(double *coeff,int size){ //print the polynomial
    if(fabs(coeff[0]) == 1.0){
    printf("%sx^%d", coeff[0] < 0 ? "-" : "", size-1);
    } 
    else{
    printf("%0.2fx^%d", coeff[0], size-1);
    }
    if(coeff[0]<0){    
    printf("- %0.2fx^%d",fabs(coeff[0]),size-1);
    }
    else if(coeff[0]>0){
        printf("%0.2fx^%d",coeff[0],size-1);
    }
    for(int i=1;i<size-1;i++){
        if(coeff[i]>0){
            printf(" + %0.2fx^%d",coeff[i],size-i-1);
        }
        else if(coeff[i]<0){
            printf(" - %0.2fx^%d",fabs(coeff[i]),size-i-1);
        }
    }
    if(coeff[size-1]>0){
        printf(" + %0.2f",coeff[size-1]);
    }
    else if(coeff[size-1]<0){
        printf(" - %0.2f",fabs(coeff[size-1]));
    }
}

double evaluate_quadratic(double a,double b,double c, double x){
    return a*x*x+b*x+c;
}

bool real_roots_quadratic(double a, double b, double c){ //checking nature of roots of quadratic
    double discriminant=b*b-4*a*c;
    if(discriminant>0){
        return true;
    }
    else{ 
        return false;
    }
}
quadratic_roots find_roots_quadratic(double a, double b, double c){ //evalauting the roots of a quadratic
    quadratic_roots solution;
    double discriminant=b*b-4*a*c;
    if(discriminant>0){
        solution.root1.imaginary=0;
        solution.root2.imaginary=0;
        solution.root1.real=(-b+sqrt(discriminant))/(2*a);
        solution.root2.real=(-b-sqrt(discriminant))/(2*a);
    }
    else if(discriminant<0){
        solution.root1.imaginary=(sqrt(-discriminant))/(2*(a));
        solution.root2.imaginary=(-sqrt(-discriminant))/(2*(a));
        solution.root1.real=-b/(2*a);
        solution.root2.real=-b/(2*a);
    }
    else if(discriminant==0){
        solution.root1.imaginary=solution.root2.imaginary=0;
        solution.root1.real=-b/(2*a);
        solution.root2.real=-b/(2*a);
    }
    return solution;
}

double minimum_quadratic(double a, double b, double c){ //calculate minimum of quadratic
    if(a<=0){
        return NAN;
    }
    else{
        double discriminant=b*b-4*a*c;
        return -discriminant/(4*a);
    }
}

double maximum_quadratic(double a, double b, double c){ //calculate maximum of quadratic
   if(a>=0){
        return NAN;
    }
    else{
        double discriminant=b*b-4*a*c;
        return -discriminant/(4*a);
    }
}

double evaluate_cubic(double a,double b,double c,double d,double x){
    return a*x*x*x+b*x*x+c*x+d;
}

bool is_monotonic_cubic(double a, double b, double c, double d){
    (void)d;
    double disc = (2*b)*(2*b) - 4*(3*a)*c;
    return disc <= 0; 
}

int real_roots_count_cubic(double a,double b,double c,double d){ // counting number of real roots of a cubic
    if(is_monotonic_cubic(a, b, c, d)){
        return 1;
    }
    else{
        quadratic_roots extremum=find_roots_quadratic(3*a,2*b,c);
        double val1=evaluate_cubic(a, b, c, d,extremum.root1.real);
        double val2=evaluate_cubic(a, b, c, d,extremum.root2.real);
        if(val1*val2>0){
            return 1;
        }
        else if(val1*val2<0){
            return 3;
        }
        else {
            return 2;
        }
    }
}

double sine_approximation(double x, int k){ // approxiamate the value of sine upto k terms 
    double sine=0.0;
    for(int i=0;i<k;i++){
        long long fact = factorial(2*i+1);
        if(fact < 0) break; 
        sine = sine + pow(-1,i) * pow(x,2*i+1) / (double)fact;
    }
    return sine;
}

double cosine_approximation(double x, int k){ // approxiamate the value of cosine upto k terms 
    double cosine=0.0;
    for(int i=0;i<k;i++){
        long long fact = factorial(2*i);
        if(fact < 0) break; 
        cosine = cosine + pow(-1,i) * pow(x,2*i) / (double)fact;
    }
    return cosine;
}

double exponential_approximation(double x, int k){ // approxiamate the value of exp function using polynomials to k terms
    double result = 0.0;
    for(int i=0;i<k;i++){
        result = result + pow(x,2*i) / (double)fact;
        
    }
    return result;
}

double approximate_root_Newton(double *coeff, int size, double x0){
    double error = 1e-7; //x0 is the initial assumption of a root
    for(int i=0;i<1000;i++){
        double f_x = evaluate_polynomial(coeff,x0,size);
        double f_der = evaluate_derivative(coeff, size, x0);
        if(fabs(f_der)< 1e-9){
            return NAN;
        }
        double next_value = x0 - f_x/f_der;
        if(fabs(next_value - x0)<error){
            return next_value;
        }
        x0= next_value;
    }
    return x0;
}

double implement_chebyshev(double x, int n,int *found, double *memo){ //evaluate Tn. Initialize memo[] to zeroes. Instead of a time complexity of O(2^n) we are able to do it in O(n).
    if(n < 0) return NAN;
    if(n==0){
        return 1;
    } 
    if(n==1){
        return x;
    }
    if(found[n]){
        return memo[n];
    }
    memo[n] = 2*x*implement_chebyshev(x,n-1,found,memo) - implement_chebyshev(x,n-2,found,memo);
    found[n]=1;
    return memo[n];
}

double implement_Legendre(double x, int n, int *found, double *memo){ // evaluate nth value using top down approach-O(n) approach
    if(n<0) return NAN;
    if(n==0){
        return 1;
    }
    if(n==1){
        return x;
    }
    if(found[n]){
        return memo[n];
    }
    memo[n] = ((2*n-1)*x*implement_Legendre(x,n-1,found,memo) - (n-1)*implement_Legendre(x,n-2,found,memo))/n;
    found[n]=1;
    return memo[n];
}

double Simpson_method(double *coeff, int size, double lower_limit, double upper_limit){//calculating the integral using Simpson's method
    long long N = 1000000;
    double h = (upper_limit - lower_limit)/(2*N);
    double integral =evaluate_polynomial(coeff,lower_limit,size);
    for(long long i=1;i<2*N;i++){
        double x_value = lower_limit+i*h;
        if(i%2!=0){
            integral = integral + 4* evaluate_polynomial(coeff, x_value,size);
        }
        else{
            integral =  integral + 2* evaluate_polynomial(coeff, x_value,size);
    }
}
    integral = integral + evaluate_polynomial(coeff, upper_limit, size);
    integral = integral*(h/3);
    return integral;
}

double Lagrange_interpolation(double *x, double *y,int n, double value){ // Lagrange interpolation
    double result=0.0;
    for(int i=0;i<=n;i++){
        double temp = y[i];
       for(int j=0;j<=n;j++){
        if(j==i) continue;
        if(x[i] == x[j]) return NAN; 
        temp*= (value - x[j])/(x[i] - x[j]);
        }
        result = result + temp;
    }
    return result;
}

double bisection_root_method(double *coeff, int size, double a, double b, double error){ 
    if(evaluate_polynomial(coeff,a,size) * evaluate_polynomial(coeff,b,size) > 0){
    return NAN;
    }
    else{
    double m =(a+b)/2.0;
    if(isRoot(coeff,m,size) || fabs(evaluate_polynomial(coeff,m,size))<=error || (b-a)/2.0<error ){
        return m;
    }
    else{
        if(evaluate_polynomial(coeff,m,size)<0){
            return bisection_root_method(coeff,size,m,b,error);
        }
        else{
            return bisection_root_method(coeff,size,a,m,error);
        }
    }
    }
}

double Cauchy_root_bound(double *coeff, int size){
    if(coeff[0] == 0.0) return NAN;
    double max = 0.0;
    for(int i = 1; i < size; i++){
        if(fabs(coeff[i]) > max){
            max = fabs(coeff[i]);
        }
    }
    return 1.0 + max / fabs(coeff[0]);
}

double Trapezoid_rule(double *coeff,int size, double lower_limit , double upper_limit){
    long long N = 1000000;
    double fa = evaluate_polynomial(coeff,lower_limit,size);
    double fb = evaluate_polynomial(coeff, upper_limit,size);
    double h =(upper_limit-lower_limit)/N;
    double integral = fa+fb;
    for(long long i=1;i<N;i++){
        double x_value = lower_limit+i*h;
        integral = integral + 2*evaluate_polynomial(coeff,x_value,size);
    }
    integral = integral*h/2;
    return integral;
}

double *remainder(double *coeff1, int size1,double *coeff2, int size2,int *rem_size){
    int deg1 = size1 - 1;
    int deg2 = size2 - 1;

    if (deg2 > deg1) {
        *rem_size = size1;
        double *copy = (double*)malloc(size1 * sizeof(double));
        if (!copy) { *rem_size = 0; return NULL; }
        for (int i = 0; i < size1; i++) copy[i] = coeff1[i];
        return copy;
    }

    int quotient_size = deg1 - deg2 + 1;  

    double *quot = divide_polynomials(coeff1, size1, coeff2, size2);
    if (!quot) { *rem_size = 0; return NULL; }
    double *prod = multiply_polynomials(coeff2, size2, quot, quotient_size);
    free(quot);
    if (!prod) { 
        *rem_size = 0; return NULL; 
    }
    double *rem = (double*)calloc(size2, sizeof(double));
    if (!rem) { 
        free(prod); 
        *rem_size = 0; 
        return NULL; 
    }
    for(int i = 0; i < size2; i++){
    rem[i] = coeff1[quotient_size + i] - prod[quotient_size + i];
}
    free(prod);
    int actual_size = size2;
    while (actual_size > 1 && fabs(rem[0]) < 1e-12) {
        for (int i = 0; i < actual_size - 1; i++) rem[i] = rem[i+1];
        actual_size--;
    }

    *rem_size = actual_size;
    return rem;
}

double** sturm_sequences(double *coeff, int size,int *sequence_count,int **sizes_out){
    double **seq  = (double**)malloc((size + 1) * sizeof(double*));
    int    *sizes = (int*)   malloc((size + 1) * sizeof(int));
    if (!seq || !sizes) {
        free(seq); free(sizes);
        *sequence_count = 0;
        *sizes_out = NULL;
        return NULL;
    }
    seq[0] = (double*)malloc(size * sizeof(double));
    if (!seq[0]) { free(seq); free(sizes); return NULL; }
    for (int i = 0; i < size; i++) seq[0][i] = coeff[i];
    sizes[0] = size;
    seq[1]  = derivative(coeff, size);
    sizes[1] = (size > 1) ? size - 1 : 1;
    int k = 2;
    while (k <= size && sizes[k-1] > 1) {
        int rem_size = 0;
        double *rem = remainder(seq[k-2], sizes[k-2],seq[k-1], sizes[k-1],&rem_size);
        if (!rem || rem_size == 0) {
            free(rem);
            break;   
        }
        for (int i = 0; i < rem_size; i++) {
            rem[i] = -rem[i];
        }
        seq[k]   = rem;
        sizes[k] = rem_size;
        k++;
    }
    *sequence_count = k;
    *sizes_out      = sizes;
    return seq;
}

void sturm_free(double **seq, int *sizes, int sequence_count){ 
    if (!seq) return;
    for (int i = 0; i < sequence_count; i++) free(seq[i]);
    free(seq);
    free(sizes);
}

int count_signs_sturm(double **seq, int *sizes, int sequence_count, double x){
    int count = 0;
    double prev_value = evaluate_polynomial(seq[0], x, sizes[0]);
    for (int i = 1; i < sequence_count; i++) {
        double cur = evaluate_polynomial(seq[i], x, sizes[i]);
        if (fabs(cur) > 1e-12) {
            if ((prev_value > 0.0 && cur < 0.0) || (prev_value < 0.0 && cur > 0.0)) {
                count++;
            }
            prev_value = cur;   
        }
    }
    return count;
}

int count_real_roots_sturm(double **seq, int *sizes,int sequence_count, double a, double b){
    int va = count_signs_sturm(seq, sizes, sequence_count, a);
    int vb = count_signs_sturm(seq, sizes, sequence_count, b);
    return abs(va - vb);
}

double central_difference(double *coeff, int size, double x, double h){ //numerical diffrentiation
    return (evaluate_polynomial(coeff,x+h,size)-evaluate_polynomial(coeff,x-h,size))/(2.0*h);
}

double second_derivative(double *coeff, int size, double x, double h){ //numerical differentiation
    return (evaluate_polynomial(coeff,x+h, size) - 2*evaluate_polynomial(coeff,x,size) + evaluate_polynomial(coeff,x-h, size))/(h*h);
}

double richardson_extrapolation(double *coeff, int size, double x, double h){
    double d1 = central_difference(coeff,size,x,h);
    double d2 = central_difference(coeff,size,x,h/2.0);
    double result = (4*d2 - d1)/3.0;
    return result;
}

double** companion_matrix(double *coeff, int size){//building the companion matrix
    int deg = size -1;
    double lead = 0.0;
    for(int i=0;i<size;i++){
    if(coeff[i]){ lead = coeff[i]; break; }
         }
    if(lead == 0.0) return NULL; 

    double **result=(double**)malloc(deg*sizeof(double*));
    for(int i=0;i<deg;i++){
        result[i] = (double*)calloc(deg,sizeof(double));
    }
    for(int i=0;i<deg;i++){
      result[i][deg-1] = -coeff[deg-i]/lead;
    }
    for(int i=1;i<deg;i++){
        result[i][i-1]=1.0;
    }
    return result;
}
























































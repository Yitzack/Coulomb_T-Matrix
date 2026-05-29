#include<iostream>
#include<complex>
#include"Math_func.h"
#include"Numerical_Tables.h"
#include"Around.h"
#include"Elements.h"
#include"Interpolation.h"
using namespace std;

Around<complex<long double>> Int_phi(int, int, int, long double, long double);
Around<complex<long double>> Int_phi_recurs(int, int, int, long double, long double, int, long double, long double);
Around<complex<long double>> Int_theta(int, int, int, long double, long double, long double);
Around<complex<long double>> Int_theta_recurs(int, int, int, long double, long double, long double, int, long double, long double);
Around<complex<long double>> Int_r(int, int, int, long double, long double, long double, long double);
Around<complex<long double>> Int_r_recurs(int, int, int, long double, long double, long double, long double, int, long double, long double);

int main()
{
	return(0);
}

Around<complex<long double>> Int_phi(int l, int lp, int m, long double q, long double qp)
{
	return(Around<complex<long double>>(complex<long double>(0,0),complex<long double>(0,0)));
}

Around<complex<long double>> Int_phi_recurs(int l, int lp, int m, long double q, long double qp, int level, long double a, long double b)
{
	return(Around<complex<long double>>(complex<long double>(0,0),complex<long double>(0,0)));
}

Around<complex<long double>> Int_theta(int l, int lp, int m, long double q, long double qp, long double phi)
{
	return(Around<complex<long double>>(complex<long double>(0,0),complex<long double>(0,0)));
}

Around<complex<long double>> Int_theta_recurs(int l, int lp, int m, long double q, long double qp, long double phi, int level, long double a, long double b)
{
	return(Around<complex<long double>>(complex<long double>(0,0),complex<long double>(0,0)));
}

Around<complex<long double>> Int_r(int l, int lp, int m, long double q, long double qp, long double phi, long double theta)
{
	return(Around<complex<long double>>(complex<long double>(0,0),complex<long double>(0,0)));
}

Around<complex<long double>> Int_r_recurs(int l, int lp, int m, long double q, long double qp, long double phi, long double theta, int level, long double a, long double b)
{
	return(Around<complex<long double>>(complex<long double>(0,0),complex<long double>(0,0)));
}


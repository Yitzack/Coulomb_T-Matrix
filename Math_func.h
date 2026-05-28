//Special math functions
#include<complex>
#include<numbers>
#include<iostream>
#include"Elements.h"
using namespace std;

//Spherical Bessel function wrapper (totally not needed, but I'll probably forget what its called when I want it.
inline long double j(unsigned int l, long double r) {return(sph_bessel(l ,r));}

//Spherical harmonics
complex<long double> Y(unsigned int l, int m, long double theta, long double phi)
{
	using namespace numbers;
	if(m < 0)
		return(pow(-1.l,m)*conj(Y(l, abs(m), theta, phi)));
	return(pow(-1.l,m)*sqrt((2*l+1)*tgammal(l-m+1)/(tgammal(l+m+1)*4.l*pi_v<long double>))*assoc_legendre(l,m,cos(theta))*exp(complex<long double>(0,m*phi)));
}

//Lienard-Wiechert potential in sphere coordinates V^\mu=P^\mu/sqrt(M^2r^2+(r/cdot P)^2)
Elements<long double> Vsphere(Elements<long double> P, Elements<long double> r, long double M)
{
	return(P/sqrt(pow(M,2)*(r*r)+pow(r*P,2)));
}

//Lienard-Wiechert potential in boosted sphereical coordinates
Elements<long double> Vboosted(Elements<long double> P, Elements<long double> rho, long double M)
{
	long double v = P[3]/M;
	return(Elements<long double>(1.l/rho[1],
				     v/rho[1]*cos(rho[2])*sqrt((1.l-pow(v,2))/(1.l-pow(v*cos(rho[2]),2))),
				     v/rho[1]*cos(rho[2])*sqrt((1.l-pow(v,2))/(1.l-pow(v*cos(rho[2]),2))),
				     0));
}

//Convert spherical coordinates to spherical boosted coordinates
Elements<long double> rho_boosted(Elements<long double> P, Elements<long double> r, long double M)
{
	long double v = P[3]/M;
	return(Elements<long double>(r[0],
				     r[1]*sqrt(1.l-pow(v*sin(r[2]),2)),
				     acos(cos(r[2])/sqrt(1.l-pow(v*sin(r[2]),2))),
				     r[3]));
}

//Convert spherical boosted coordinates to spherical coordinates
Elements<long double> rho_boosted(Elements<long double> P, Elements<long double> rho, long double M)
{
	long double v = P[3]/M;
	return(Elements<long double>(rho[0],
				     rho[1]*sqrt((1.l-pow(v*cos(rho[2]),2))/(1.l-pow(v,2))),
				     acos(cos(rho[2])*sqrt((1.l-pow(v,2))/(1.l-pow(v*cos(rho[2]),2))),
				     rho[3]));
}


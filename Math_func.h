//Special math functions
#include<complex>
#include<cmath>
#include<numbers>
#include<iostream>
#include<utility>
#include<limits>
#include"Elements.h"
using namespace std;

long double jj_Yukawa(int, int, long double, long double, long double);
long double jj_Yukawa(int, int, long double, long double, long double, long double);
long double Extrema_jj_on_r(int, int, long double, long double, long double);
long double jj_on_r(int, int, long double, long double, long double);
long double jj_on_r(int, int, long double, long double, long double, long double);
long double j(unsigned int, long double);
complex<long double> Y(unsigned int, int, long double, long double);
Elements<long double> Vboosted(Elements<long double>, Elements<long double>, long double);
Elements<long double> Vsphere(Elements<long double>, Elements<long double>, long double);
Elements<long double> r_spherical(Elements<long double>, Elements<long double>, long double);
Elements<long double> rho_boosted(Elements<long double>, Elements<long double>, long double);
long double Hypergeometric2F1(long double, long double, long double, long double);
long double Si(long double);
complex<long double> Ci(long double);
complex<long double> Ei(complex<long double>);

long double jj_Yukawa(int l, int lp, long double q, long double qp, long double a)
{
	using namespace numbers;
	long double mu = .05;

	if(qp == 0 && lp == 0 && q != 0)
	{
		swap(q,qp);
		swap(l,lp);
	}

	if(q == 0 && l == 0 && qp != 0)
	{
		return(exp(-a*mu)*(qp*cos(lp*pi_v<long double>/2.l-a*qp)-mu*sin(lp*pi_v<long double>/2.l-a*qp))/(qp*(pow(qp,2)+pow(mu,2))));
	}
	if(q == 0 && l == 0 && qp == 0 && lp == 0)
	{
		return(exp(-a*mu)*(1+a*mu)/pow(mu,2));
	}
	if((q == 0 && l != 0) || (qp == 0 && lp != 0))
	{
		return(0);
	}

	long double cos_diff = 0, cos_sum = 0;	//cos((l±lp)*pi/2)
	long double sin_diff = 0, sin_sum = 0;	//sin((l±lp)*pi/2)

	switch((l-lp)%4)
	{
		case 0:
			cos_diff = 1;
			break;
		case -3:
		case 1:
			sin_diff = 1;
			break;
		case -2:
		case 2:
			cos_diff = -1;
			break;
		case -1:
		case 3:
			sin_diff = -1;
			break;
	}
	switch((l+lp)%4)
	{
		case 0:
			cos_sum = 1;
			break;
		case 1:
			sin_sum = 1;
			break;
		case 2:
			cos_sum = -1;
			break;
		case 3:
			sin_sum = -1;
			break;
	}

	complex<long double> phase[4] = {complex<long double>(cos_sum,sin_sum),complex<long double>(-cos_diff,sin_diff),complex<long double>(-cos_diff,-sin_diff),complex<long double>(cos_sum,-sin_sum)};
	complex<long double> ExpIntgrals[4] = {complex<long double>(-a*mu,-a*(q+qp)),complex<long double>(-a*mu,a*(q-qp)),complex<long double>(-a*mu,a*(qp-q)),complex<long double>(-a*mu,a*(q+qp))};

	long double answer = 0;
	for(int i = 0; i < 4; i++)
	{
		answer += (phase[i]*Ei(ExpIntgrals[i])).real();
	}

	//b->inf limit of jj_Yukawa(l,lp,q,qp,a,b)
	if(q-qp < 0)
		answer -= 2.l*pi_v<long double>*(sin_sum+sin_diff);
	else if(q-qp > 0)
		answer -= 2.l*pi_v<long double>*(sin_sum-sin_diff);
	else
		answer -= 2.l*pi_v<long double>*sin_sum;

	answer /= 4.l*q*qp;

	return(answer);
}

long double jj_Yukawa(int l, int lp, long double q, long double qp, long double a, long double b)
{
	using namespace numbers;
	long double mu = .05;

	if(qp == 0 && lp == 0 && q != 0)
	{
		swap(q,qp);
		swap(l,lp);
	}

	if(q == 0 && l == 0 && qp != 0)
	{
		return((exp(-a*mu)*(qp*cos(lp*pi_v<long double>/2.l-a*qp)-mu*sin(lp*pi_v<long double>/2.l-a*qp))+exp(-b*mu)*(-qp*cos(lp*pi_v<long double>/2.l-b*qp)+mu*sin(lp*pi_v<long double>/2.l-b*qp)))/(qp*(pow(qp,2)+pow(mu,2))));
	}
	if(q == 0 && l == 0 && qp == 0 && lp == 0)
	{
		return((exp(-a*mu)*(1+a*mu)-exp(-b*mu)*(1+b*mu))/pow(mu,2));
	}
	if((q == 0 && l != 0) || (qp == 0 && lp != 0))
	{
		return(0);
	}

	long double cos_diff = 0, cos_sum = 0;	//cos((l±lp)*pi/2)
	long double sin_diff = 0, sin_sum = 0;	//sin((l±lp)*pi/2)

	switch((l-lp)%4)
	{
		case 0:
			cos_diff = 1;
			break;
		case -3:
		case 1:
			sin_diff = 1;
			break;
		case -2:
		case 2:
			cos_diff = -1;
			break;
		case -1:
		case 3:
			sin_diff = -1;
			break;
	}
	switch((l+lp)%4)
	{
		case 0:
			cos_sum = 1;
			break;
		case 1:
			sin_sum = 1;
			break;
		case 2:
			cos_sum = -1;
			break;
		case 3:
			sin_sum = -1;
			break;
	}

	complex<long double> phase[8] = {complex<long double>(cos_sum,sin_sum),complex<long double>(-cos_diff,sin_diff),complex<long double>(-cos_diff,-sin_diff),complex<long double>(cos_sum,-sin_sum),complex<long double>(-cos_sum,-sin_sum),complex<long double>(cos_diff,-sin_diff),complex<long double>(cos_diff,sin_diff),complex<long double>(-cos_sum,sin_sum)};
	complex<long double> ExpIntgrals[8] = {complex<long double>(-a*mu,-a*(q+qp)),complex<long double>(-a*mu,a*(q-qp)),complex<long double>(-a*mu,a*(qp-q)),complex<long double>(-a*mu,a*(q+qp)),complex<long double>(-b*mu,-b*(q+qp)),complex<long double>(-b*mu,b*(q-qp)),complex<long double>(-b*mu,b*(qp-q)),complex<long double>(-b*mu,b*(q+qp))};

	long double answer = 0;
	for(int i = 0; i < 8; i++)
		answer += (phase[i]*Ei(ExpIntgrals[i])).real();
	answer /= 4.l*q*qp;

	return(answer);
}

//returns any extrema of jj_on_r for the tail near r.
long double Extrema_jj_on_r(int l, int lp, long double q, long double qp, long double r)
{
	long double xi;
	long double xi1;
	long double h = 1e-4;
	long double extrema = 0;

	for(long double test = r-1; test <= r+1; test += .01)	//Find the most extrem value within .5 of r and use that as the feed for Newton's method.
	{
		if(abs(jj_on_r(l, lp, q, qp, test)) > extrema)
		{
			extrema = jj_on_r(l, lp, q, qp, test);
			xi = test;
		}
	}

	xi1 = xi - h*(jj_on_r(l, lp, q, qp, xi-h)-jj_on_r(l, lp, q, qp, xi+h))/(-2*jj_on_r(l, lp, q, qp, xi-h)+4*jj_on_r(l, lp, q, qp, xi)-2*jj_on_r(l, lp, q, qp, xi+h));
	while(abs(xi1-xi)>1e-4)
	{
		xi = xi1;
		xi1 = xi - h*(jj_on_r(l, lp, q, qp, xi-h)-jj_on_r(l, lp, q, qp, xi+h))/(-2*jj_on_r(l, lp, q, qp, xi-h)+4*jj_on_r(l, lp, q, qp, xi)-2*jj_on_r(l, lp, q, qp, xi+h));
	}

	return(xi1);
}

//int_a^inf j_l(r q)j_l'(r q')/r*r^2dr
long double jj_on_r(int l, int lp, long double q, long double qp, long double a)
{
	//function requires that q>=qp
	if(qp > q)
		swap(q,qp);

	if(qp == 0 && lp == 0 && q != 0)
	{
		return(cos(lp*numbers::pi_v<long double>/2.l-a*q)/pow(q,2));
	}
	if(q == 0 && l == 0 && qp == 0 && lp == 0)
	{
		return(numeric_limits<double>::infinity());
	}
	if((q == 0 && l != 0) || (qp == 0 && lp != 0))
	{
		return(0);
	}

	long double cos_diff = 0, cos_sum = 0;	//cos((l±lp)*pi/2)
	long double sin_diff = 0, sin_sum = 0;	//sin((l±lp)*pi/2)

	switch((l-lp)%4)
	{
		case 0:
			cos_diff = 1;
			break;
		case -3:
		case 1:
			sin_diff = 1;
			break;
		case -2:
		case 2:
			cos_diff = -1;
			break;
		case -1:
		case 3:
			sin_diff = -1;
			break;
	}
	switch((l+lp)%4)
	{
		case 0:
			cos_sum = 1;
			break;
		case 1:
			sin_sum = 1;
			break;
		case 2:
			cos_sum = -1;
			break;
		case 3:
			sin_sum = -1;
			break;
	}

	return((-cos_diff*Ci(a*(q-qp))+cos_sum*Ci(-a*(q+qp))+sin_diff*Si(a*(q-qp))+sin_sum*Si(a*(q+qp))+numbers::pi_v<long double>*(sin_diff-sin_sum)/2.l).real()/(2.l*q*qp));
}

//int_a^b j_l(r q)j_l'(r q')/r*r^2dr
long double jj_on_r(int l, int lp, long double q, long double qp, long double a, long double b)
{
	using namespace numbers;

	//function requires that q>=qp
	if(qp > q)
	{
		swap(q,qp);
		swap(l,lp);
	}

	if(qp == 0 && lp == 0 && q != 0)
	{
		return((cos(lp*pi_v<long double>/2.l-a*q)-cos(lp*pi_v<long double>/2.l-b*q))/pow(q,2));
	}
	if(q == 0 && l == 0 && qp == 0 && lp == 0)
	{
		return((pow(b,2)-pow(a,2))/2.l);
	}
	if((q == 0 && l != 0) || (qp == 0 && lp != 0))
	{
		return(0);
	}

	long double cos_diff = 0, cos_sum = 0;	//cos((l±lp)*pi/2)
	long double sin_diff = 0, sin_sum = 0;	//sin((l±lp)*pi/2)

	switch((l-lp)%4)
	{
		case 0:
			cos_diff = 1;
			break;
		case -3:
		case 1:
			sin_diff = 1;
			break;
		case -2:
		case 2:
			cos_diff = -1;
			break;
		case -1:
		case 3:
			sin_diff = -1;
			break;
	}
	switch((l+lp)%4)
	{
		case 0:
			cos_sum = 1;
			break;
		case 1:
			sin_sum = 1;
			break;
		case 2:
			cos_sum = -1;
			break;
		case 3:
			sin_sum = -1;
			break;
	}

	return((-cos_diff*Ci(a*(q-qp))+cos_diff*Ci(b*(q-qp))+cos_sum*Ci(-a*(q+qp))-cos_sum*Ci(-b*(q+qp))+sin_diff*Si(a*(q-qp))-sin_diff*Si(b*(q-qp))+sin_sum*Si(a*(q+qp))-sin_sum*Si(b*(q+qp))).real()/(2.l*q*qp));
}

//Spherical Bessel (this recursion algorithm runs in 2.75 ns to the standard algorithm's 6.69 us on -O3 and 17.25 us instead of 52.24 ms without optimization)
long double j(unsigned int l, long double r)
{
	//check for r==0
	if(r == 0) return(l==0?1.l:0);

	if(r < .0433371*pow(l,1.73982))	//asymptotic form
		return(pow(r,l)*pow(2,l)*tgammal(l+1)/tgammal(2*l+2)*(1.l-pow(r,2)/(long double)(6+4*l)+pow(r,4)/(long double)(120+128*l+32*pow(l,2))));

	//recursion
	long double jn1 = sin(r)/r;
	long double jn = sin(r)/pow(r,2)-cos(r)/r;
	long double temp;

	if(l == 0)
		return(jn1);
	else if(l == 1)
		return(jn);

	for(int i = 1; i < l; i++)
	{
		temp = (2*i+1)*jn/r-jn1;
		jn1 = jn;
		jn = temp;
	}

	return(jn);
}

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
	return(P/sqrt(pow(M,2)*(r[1]*r[1])+pow(P[1]*r[1]*cos(r[2]),2)));
}

//Lienard-Wiechert potential in boosted sphereical coordinates
Elements<long double> Vboosted(Elements<long double> P, Elements<long double> rho, long double M)
{
	long double v = P[1]/M;
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
Elements<long double> r_spherical(Elements<long double> P, Elements<long double> rho, long double M)
{
	long double v = P[3]/M;
	return(Elements<long double>(rho[0],
				 rho[1]*sqrt((1.l-pow(v*cos(rho[2]),2))/(1.l-pow(v,2))),
				 acos(cos(rho[2])*sqrt((1.l-pow(v,2))/(1.l-pow(v*cos(rho[2]),2)))),
				 rho[3]));
}

long double Hypergeometric2F1(long double a, long double b, long double c, long double z)
{
	long double term = 1;
	long double sum = 1;
	int k = 0;

	for(k; k <= 10000 && abs(term) > 1e-15*abs(sum); k++)
	{
		term *= (a+k)*(b+k)*z/((c+k)*(k+1));
		sum += term;
	}

	return(sum);
}

//Cosine Integral for large positive definite x
complex<long double> Ci(long double x)
{
	if(x > 0)
		return((x*(18+pow(x,2))*cos(x)+(12+22*pow(x,2)+pow(x,4))*sin(x))/(x*(36+24*pow(x,2)+pow(x,4))));
	return(complex<long double>((x*(18+pow(x,2))*cos(x)+(12+22*pow(x,2)+pow(x,4))*sin(x))/(x*(36+24*pow(x,2)+pow(x,4))),numbers::pi_v<long double>));
}

//Sine Integral for large positive definite x
long double Si(long double x)
{
	if(x > 0)
		return(-((12+22*pow(x,2)+pow(x,4))*cos(x)+x*(18+pow(x,2))*sin(x))/(x*(36+24*pow(x,2)+pow(x,4)))+numbers::pi_v<long double>/2.l);
	return(-((12+22*pow(x,2)+pow(x,4))*cos(x)+x*(18+pow(x,2))*sin(x))/(x*(36+24*pow(x,2)+pow(x,4)))-numbers::pi_v<long double>/2.l);
}

//Exponential Integral Ei for large abs(x)
complex<long double> Ei(complex<long double> z)
{
	using namespace numbers;
	complex<long double> Answer;
	complex<long double> Terms(1);
	complex<long double> Temp;
	int n = 1;

	if(abs(z) < 20)
	{
		Answer = 0;
		do
		{
			Temp = pow(z,n)/((long double)(n)*tgammal(n+1));
			Answer += Temp;
			n++;
		}while(abs(Temp)/abs(Answer) > 1e-15 && n <= 100);
		return(Answer+egamma_v<long double>+log(z));
	}

	if(z.imag() > 0)
		Answer = complex<long double>(0,pi_v<long double>);
	else if(z.imag() < 0)
		Answer = complex<long double>(0,-pi_v<long double>);
	else
		Answer = 0;

	do
	{
		Temp = tgammal(n+1)/pow(z,n);
		if(abs(Terms) < abs(Terms+Temp))	//quit because it started diverging
			break;
		Terms += Temp;
		n++;
	}while(n <= 100);

	return(exp(z)*Terms/z+Answer);
}


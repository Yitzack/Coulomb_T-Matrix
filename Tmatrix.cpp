#include<iostream>
#include<complex>
#include<numeric>
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
complex<long double> Integrand(int, int, int, long double, long double, long double, long double, long double);

int main()
{
	return(0);
}

Around<complex<long double>> Int_phi(int l, int lp, int m, long double q, long double qp)	//Phi integral is analytic.
{
	using namespace numbers;
	return(Int_theta(l, lp, m, q, qp, 0.l)*2.l*pi_v<long double>);
}

Around<complex<long double>> Int_theta(int l, int lp, int m, long double q, long double qp, long double phi)
{
	using namespace numbers;
	Around<complex<long double>> Answer(complex<long double>(0,0),complex<long double>(0,0));
	int n = lcm(l, lp);

	for(int i = 0; i < n; i++)
		Answer += Int_theta_recurs(l, lp, m, q, qp, phi, 0, i*pi_v<long double>/n, (i+1)*pi_v<long double>/n);

	return(Answer);
}

Around<complex<long double>> Int_theta_recurs(int l, int lp, int m, long double q, long double qp, long double phi, int level, long double a, long double b)
{
	GKRule QRule;
	QRule.Rule = RuleType::GK16;
	Around<complex<long double>> Answerh(complex<long double>(0,0),complex<long double>(0,0));
	Around<complex<long double>> Answerl(complex<long double>(0,0),complex<long double>(0,0));
	Around<complex<long double>> Answer;
	long double mid = (a+b)/2.l;
	long double dist = b-a;

	Answer = Int_r(l, lp, m, q, qp, phi, mid)*sin(mid);
	Answerh = Answer*QRule.wh(0);
	Answerl = Answer*QRule.wl(0);
	for(int i = 0; i < QRule.length(); i++)
	{
		Answer = Int_r(l, lp, m, q, qp, phi, mid+dist*QRule.Disp(i))*sin(mid+dist*QRule.Disp(i));
		Answerh += Answer*QRule.wh(i+1);
		Answerl += Answer*QRule.wl(i+1);
		Answer = Int_r(l, lp, m, q, qp, phi, mid-dist*QRule.Disp(i))*sin(mid-dist*QRule.Disp(i));
		Answerh += Answer*QRule.wh(i+1);
		Answerl += Answer*QRule.wl(i+1);
	}

	Answer = Around<complex<long double>>(Answerh.Value(),complex<long double>(
		sqrt(pow(Answerh.Value().real()-Answerl.Value().real(),2)+pow(Answerh.Error().real(),2)),
		sqrt(pow(Answerh.Value().imag()-Answerl.Value().imag(),2)+pow(Answerh.Error().imag(),2))));
	if(abs(Answer.RelErr()) > 1e-6 && level < 5)
	{
		return(Int_theta_recurs(l, lp, m, q, qp, phi, level+1, a, mid)+Int_theta_recurs(l, lp, m, q, qp, phi, level+1, mid, b));
	}
	return(Answer);
}

Around<complex<long double>> Int_r(int l, int lp, int m, long double q, long double qp, long double phi, long double theta)
{
	using namespace numbers;
	Around<complex<long double>> Answer(complex<long double>(0,0),complex<long double>(0,0));
	int n = lcm(l, lp);

	for(int i = 0; i < n; i++)
		Answer += Int_r_recurs(l, lp, m, q, qp, phi, theta, 0, i*pi_v<long double>/n, (i+1)*pi_v<long double>/n);

	return(Answer);
}

Around<complex<long double>> Int_r_recurs(int l, int lp, int m, long double q, long double qp, long double phi, long double theta, int level, long double a, long double b)
{
	GKRule QRule;
	QRule.Rule = RuleType::GK16;
	complex<long double> Answerh(0,0);
	complex<long double> Answerl(0,0);
	complex<long double> Answer;
	Around<complex<long double>> Result;
	long double mid = (a+b)/2.l;
	long double dist = b-a;

	Answer = Integrand(l, lp, m, q, qp, phi, theta, mid)*pow(mid,2);
	Answerh = Answer*QRule.wh(0);
	Answerl = Answer*QRule.wl(0);
	for(int i = 0; i < QRule.length(); i++)
	{
		Answer = Integrand(l, lp, m, q, qp, phi, theta, mid+dist*QRule.Disp(i))*pow(mid+dist*QRule.Disp(i),2);
		Answerh += Answer*QRule.wh(i+1);
		Answerl += Answer*QRule.wl(i+1);
		Answer = Integrand(l, lp, m, q, qp, phi, theta, mid-dist*QRule.Disp(i))*pow(mid-dist*QRule.Disp(i),2);
		Answerh += Answer*QRule.wh(i+1);
		Answerl += Answer*QRule.wl(i+1);
	}

	Result = Around<complex<long double>>(Answerh,Answerh-Answerl);
	if(abs(Result.RelErr()) > 1e-6 && level < 5)
	{
		return(Int_r_recurs(l, lp, m, q, qp, phi, theta, level+1, a, mid)+Int_r_recurs(l, lp, m, q, qp, phi, theta, level+1, mid, b));
	}
	return(Result);
}


complex<long double> Integrand(int l, int lp, int m, long double q, long double qp, long double phi, long double theta, long double r)
{
	return(complex<long double>(0,0));
}

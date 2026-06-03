#include<iostream>
#include<iomanip>
#include<complex>
#include<numeric>
#include<cmath>
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

long double f0(long double q, long double qp)
{
	if(q == 0)
		return(1.l/(pow(qp,2)+.0025));
	else if(qp == 0)
		return(1.l/pow(q,2));
	return(log((pow(q+qp,2)+.0025)/(pow(q-qp,2)+.0025))/(4*q*qp));
}

long double f1(long double q, long double qp)
{
	if(q == 0 || qp == 0)
		return(0);
	return((-4*q*qp+(pow(q,2)+pow(qp,2))*log(pow(q+qp,2)/pow(q-qp,2)))/(8*pow(q*qp,2)));
}

long double f2(long double q, long double qp)
{
	if(q == 0 || qp == 0)
		return(0);
	return((-12*q*qp*(pow(q,2)+pow(qp,2))+(3*(pow(q,4)+pow(qp,4))+2*pow(q*qp,2))*log(pow(q+qp,2)/pow(q-qp,2)))/(32*pow(q*qp,3)));
}

long double f3(long double q, long double qp)
{
	if(q == 0 || qp == 0)
		return(0);
	return((-4*q*qp*(15*pow(q,4)+14*pow(q,2)*pow(qp,2)+15*pow(qp,4))+3*(5*pow(q,6)+3*pow(q,4)*pow(qp,2)+3*pow(q,2)*pow(qp,4)+5*pow(qp,6))*log(pow(q+qp,2)/pow(q-qp,2)))/(192*pow(q*qp,4)));
}

long double f4(long double q, long double qp)
{
	if(q == 0 || qp == 0)
		return(0);
	return((-5*q*qp*(pow(q,2)+pow(qp,2))*(21*pow(q,4)-2*pow(q,2)*pow(qp,2)+21*pow(qp,4))+.75l*(35*pow(q,8)+20*pow(q,6)*pow(qp,2)+18*pow(q,4)*pow(qp,4)+20*pow(q,2)*pow(qp,6)+35*pow(qp,8))*log(pow(q+qp,2)/pow(q-qp,2)))/(384*pow(q*qp,5)));
}

long double f5(long double q, long double qp)
{
	if(q == 0 || qp == 0)
		return(0);
	return((-q*qp*(945*pow(q,8)+840*pow(q,6)*pow(qp,2)+814*pow(q,4)*pow(qp,4)+840*pow(q,2)*pow(qp,6)+945*pow(qp,8))+3.75l*(pow(q,2)+pow(qp,2))*(63*pow(q,8)-28*pow(q,6)*pow(qp,2)+58*pow(q,4)*pow(qp,4)-28*pow(q,2)*pow(qp,6)+63*pow(qp,8))*log(pow(q+qp,2)/pow(q-qp,2)))/(3840*pow(q*qp,6)));
}

int main()
{
	Around<long double> real;
	Around<complex<long double>> Ans;
	for(long double q = 0; q < 10; q++)
	{
		for(long double qp = q; qp < 10; qp++)
		{
			Ans = Int_phi(0,0,0,q,qp);
			real = Around<long double>(Ans.Value().real(), Ans.Error().real());
			cout << q << " " << qp << setw(10) << f0(q,qp) << setw(14) << real << setw(14) << (real.Value()/f0(q,qp)-1.l)*100.l << setw(14) << real.RelErr() << endl;
		}
	}
	/*for(int l = 0; l <= 5; l++)
	{
		for(int lp = 0; lp <= 5; lp++)
		{
			for(int m = -min(l,lp); m <= min(l,lp); m++)
			{
				Ans = Int_phi(l, lp, m, 0, 0);
				real = Around<long double>(Ans.Value().real(), Ans.Error().real());
				cout << l << " " << lp << " " << m << " " << flush;
				if(l == lp)
					cout << 1 << setw(14) << real << setw(14) << real.RelErr() << endl;
				else
					cout << 0 << setw(14) << real << setw(14) << real.RelErr() << endl;
			}
		}
	}*/

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
	int n = l+lp+1;

	for(int i = 0; i < n; i++)
		Answer += Int_theta_recurs(l, lp, m, q, qp, phi, 0, acos(2.l*(i+1)/n-1), acos(2.l*i/n-1));

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
	long double dist = (b-a)/2.l;

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
	Answerh = Answerh*dist;
	Answerl = Answerl*dist;

	Answer = Around<complex<long double>>(Answerh.Value(),complex<long double>(
		sqrt(pow(Answerh.Value().real()-Answerl.Value().real(),2)+pow(Answerh.Error().real(),2)),
		sqrt(pow(Answerh.Value().imag()-Answerl.Value().imag(),2)+pow(Answerh.Error().imag(),2))));
	if(abs(Answer.RelErr()) > 1e-5 && level < 5)
	{
		return(Int_theta_recurs(l, lp, m, q, qp, phi, level+1, a, mid)+Int_theta_recurs(l, lp, m, q, qp, phi, level+1, mid, b));
	}

	return(Answer);
}

Around<complex<long double>> Int_r(int l, int lp, int m, long double q, long double qp, long double phi, long double theta)
{
	using namespace numbers;
	long double r0 = pi_v<long double>/fmax(q,qp)+(long double)(max(l,lp));
	long double delta_r = 9.l*pi_v<long double>/(q+qp);
	Around<complex<long double>> Answer;
	Around<complex<long double>> Temp;
	int i = 0;

	if(isinf(delta_r) || isnan(delta_r))
	{
		r0 = 1;
		delta_r = 1;
	}

	Answer = Int_r_recurs(l, lp, m, q, qp, phi, theta, 0, 0, r0);
	do
	{
		Temp = Int_r_recurs(l, lp, m, q, qp, phi, theta, 0, r0, r0+delta_r);
		Answer += Temp;
		i++;
		r0 += delta_r;
	}while(((Temp/Answer).Value().real() > 1e-8 || (Temp/Answer).Value().imag() > 1e-8) || r0-delta_r < 250);

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
	long double dist = (b-a)/2.l;

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
	Answerh = Answerh*dist;
	Answerl = Answerl*dist;

	Result = Around<complex<long double>>(Answerh,Answerh-Answerl);
	if(abs(Result.RelErr()) > 1e-8 && level < 5)
	{
		return(Int_r_recurs(l, lp, m, q, qp, phi, theta, level+1, a, mid)+Int_r_recurs(l, lp, m, q, qp, phi, theta, level+1, mid, b));
	}

	return(Result);
}


complex<long double> Integrand(int l, int lp, int m, long double q, long double qp, long double phi, long double theta, long double r)
{
	return(Y(l, m, theta, phi)*conj(Y(lp, m, theta, phi))*j(l, r*q)*j(lp, r*qp)*exp(-r/20.l)/r);
}

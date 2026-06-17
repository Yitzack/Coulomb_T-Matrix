#include<iostream>
#include<iomanip>
#include<complex>
#include<numeric>
#include<cmath>
#include<chrono>
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

long double f(unsigned int l, unsigned int lp, long double q, long double qp)
{
	using namespace numbers;
	if(q > qp)
	{
		swap(q, qp);
		swap(l, lp);
	}

	if(q == qp)
		return(1./0.);
	else if(l == 0 && lp == 0)
	{
		if(q != 0)
			return(atanh(q/qp)/(q*qp));
		return(pow(qp,-2));
	}
	else if(q == 0)
	{
		if(lp != 0 && l == 0)
			return(sqrt(pi_v<long double>)*tgammal(1+lp/2.)/(pow(qp,2)*tgammal((1+lp)/2.)));
		return(0);
	}
	else
		return(pi_v<long double>*pow(q,l)*tgammal(1.l+l/2.l+lp/2.l)*Hypergeometric2F1((1.l+l-lp)/2.l,1.l+l/2.l+lp/2.l,1.5l+l,pow(q/qp,2))/(2*pow(qp,2+l)*tgammal(.5l-l/2.l+lp/2.l)*tgammal(1.5l+l)));
}

bool This_one;
long double final_r;

int main()
{
	Around<long double> real;
	Around<complex<long double>> Ans;
	long double expected;

	for(int l = 0; l < 10; l++)
	{
		for(int lp = l; lp < 10; lp++)
		{
			for(long double q = 0; q < 10; q++)
			{
				for(long double qp = q+1; qp < 10; qp++)
				{
					final_r = 0;
					Ans = Int_phi(l,lp,0,q,qp);
					real = Around<long double>(Ans.Value().real(), Ans.Error().real());
					expected = f(l,lp,q,qp);
					cout << l << " " << lp << " " << q << " " << qp << setw(12) << expected << setw(14) << real << setw(14) << (real.Value()/expected-1.l)*100.l << setw(14) << real.RelErr() << setw(12) << final_r << setw(12) << flush;
					if(expected == real)
						cout << "right" << endl;
					else
						cout << "wrong" << endl;
				}
			}
		}
	}

	return(0);
}

Around<complex<long double>> Int_phi(int l, int lp, int m, long double q, long double qp)	//Phi integral is analytic.
{
	if((l > 0 && q == 0) || (lp > 0 && qp == 0))
		return(Around<complex<long double>>(complex<long double>(0,0),complex<long double>(0,0)));
	using namespace numbers;
	return(Int_theta(l, lp, m, q, qp, 0.l)*2.l*pi_v<long double>);
}

Around<complex<long double>> Int_theta(int l, int lp, int m, long double q, long double qp, long double phi)
{
	using namespace numbers;
	Around<complex<long double>> Answer(complex<long double>(0,0),complex<long double>(0,0));
	int n = l+lp+1;

	for(int i = 0; i < n; i++)
	{
		if(i == 0)
			This_one = true;
		else
			This_one = false;
		Answer += Int_theta_recurs(l, lp, m, q, qp, phi, 0, acos(2.l*(i+1)/n-1), acos(2.l*i/n-1));
	}

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

	if(This_one && level == 0)
		This_one = true;
	else
		This_one = false;
	Answer = Int_r(l, lp, m, q, qp, phi, mid)*sin(mid);
	This_one = false;
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
	long double delta_r = 3.l*pi_v<long double>/(q+qp);
	static Around<complex<long double>> Answer;
	Around<complex<long double>> Temp;

	if(!This_one)
		return(Answer/(4.l*pi_v<long double>));//*Y(l, m, theta, phi)*conj(Y(lp, m, theta, phi)));

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
		r0 += delta_r;
	}while((((Temp/Answer).Value().real() > 1e-15 || (Temp/Answer).Value().imag() > 1e-15) || r0-delta_r < 500) && !(abs(Temp.Value().real()/jj_on_r(l, lp, q, qp, r0-delta_r, r0)-1.l) < 1e-6 || abs(jj_on_r(l, lp, q, qp, r0)/Answer.Value().real()) < 1e-8));
	Answer += complex<long double>(jj_on_r(l, lp, q, qp, r0),0);
	final_r = r0;

	return(Answer/(4.l*pi_v<long double>));//*Y(l, m, theta, phi)*conj(Y(lp, m, theta, phi)));
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
	return(j(l, r*q)*j(lp, r*qp)/r);
	//return(Y(l, m, theta, phi)*conj(Y(lp, m, theta, phi))*j(l, r*q)*j(lp, r*qp)*exp(-r/20.l)/r);
}

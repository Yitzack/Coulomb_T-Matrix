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
		return(1.l/(pow(q,2)+.0025));
	return(log((pow(q+qp,2)+.0025)/(pow(q-qp,2)+.0025))/(4*q*qp));
}

long double f1(long double q, long double qp)
{
	if(q == 0 || qp == 0)
		return(0);
	return((-4*q*qp+(pow(q,2)+pow(qp,2)+.0025)*log((pow(q+qp,2)+.0025)/(pow(q-qp,2)+.0025)))/(8*pow(q*qp,2)));
}

long double f2(long double q, long double qp)
{
	if(q == 0 || qp == 0)
		return(0);
	long double mu = .05;
	if(q == qp)
		return(-(1248*pow(q,4)+84*pow(q*mu,2)+420*pow(q,2)*(pow(2*q,2)+pow(mu,2))*cos(2*atan((2*q)/mu))+8*pow(q,4)*(-25+12*numbers::egamma_v<long double>+12*log(mu))+12*pow(q*mu,2)*(-25+12*numbers::egamma_v<long double>+12*log(mu))+3*pow(mu,4)*(-25+12*numbers::egamma_v<long double>+12*log(mu))-8*pow(q,4)*(-25+12*numbers::egamma_v<long double>+6*log(pow(mu,2)+pow(2*q,2)))+96*pow(q,3)*(6*mu*atan((2*q)/mu)+q*(-25+12*numbers::egamma_v<long double>+6*log(pow(mu,2)+pow(2*q,2))))+60*pow(q,2)*(pow(2*q,2)+pow(mu,2))*(cos(2*atan((2*q)/mu))*(-25+12*numbers::egamma_v<long double>+6*log(pow(mu,2)+pow(2*q,2)))-12*atan((2*q)/mu)*sin(2*atan((2*q)/mu)))-72*q*pow(pow(mu,2)+pow(2*q,2),1.5l)*sin(3*atan((2*q)/mu))-3*sqrt(1+(pow(2*q,2))/pow(mu,2))*pow(mu,2)*(pow(2*q,2)+pow(mu,2))*(cos(3*atan((2*q)/mu))*(-25+12*numbers::egamma_v<long double>+6*log(pow(mu,2)+pow(2*q,2)))-12*atan((2*q)/mu)*sin(3*atan((2*q)/mu)))-18*q*sqrt(pow(mu,2)+pow(2*q,2))*(pow(2*q,2)+pow(mu,2))*(12*atan((2*q)/mu)*cos(3*atan((2*q)/mu))+(-25+12*numbers::egamma_v<long double>+6*log(pow(mu,2)+pow(2*q,2)))*sin(3*atan((2*q)/mu))))/(192*pow(q,6)));
	return((-312*pow(q,2)*(q-qp)*qp+312*q*(q-qp)*pow(qp,2)-312*pow(q,2)*qp*(q+qp)-312*q*pow(qp,2)*(q+qp)-84*pow(q,2)*(pow(q+qp,2)+pow(mu,2))*cos(2*atan((q+qp)/mu))-252*q*qp*(pow(q+qp,2)+pow(mu,2))*cos(2*atan((q+qp)/mu))-84*pow(qp,2)*(pow(q+qp,2)+pow(mu,2))*cos(2*atan((q+qp)/mu))+84*pow(q,2)*(pow(q-qp,2)+pow(mu,2))*cos(2*atan(abs(q-qp)/mu))-252*q*qp*(pow(q-qp,2)+pow(mu,2))*cos(2*atan(abs(q-qp)/mu))+84*pow(qp,2)*(pow(q-qp,2)+pow(mu,2))*cos(2*atan(abs(q-qp)/mu))+8*pow(q,2)*pow(qp,2)*(-25+12*numbers::egamma_v<long double>+6*log(1+pow(q+qp,2)/pow(mu,2))+12*log(mu))-24*pow(q,2)*qp*(12*mu*atan((q+qp)/mu)+(q+qp)*(-25+12*numbers::egamma_v<long double>+6*log(1+pow(q+qp,2)/pow(mu,2))+12*log(mu)))-24*q*pow(qp,2)*(12*mu*atan((q+qp)/mu)+(q+qp)*(-25+12*numbers::egamma_v<long double>+6*log(1+pow(q+qp,2)/pow(mu,2))+12*log(mu)))-8*pow(q,2)*pow(qp,2)*(-25+12*numbers::egamma_v<long double>+12*log(mu)+6*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2)))-(24*pow(q,2)*(q-qp)*qp*(12*mu*atan(abs(q-qp)/mu)+abs(q-qp)*(-25+12*numbers::egamma_v<long double>+12*log(mu)+6*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2)))))/abs(q-qp)+(24*q*(q-qp)*pow(qp,2)*(12*mu*atan(abs(q-qp)/mu)+abs(q-qp)*(-25+12*numbers::egamma_v<long double>+12*log(mu)+6*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2)))))/abs(q-qp)-12*pow(q,2)*(pow(q+qp,2)+pow(mu,2))*(cos(2*atan((q+qp)/mu))*(-25+12*numbers::egamma_v<long double>+6*log(1+pow(q+qp,2)/pow(mu,2))+12*log(mu))-12*atan((q+qp)/mu)*sin(2*atan((q+qp)/mu)))-36*q*qp*(pow(q+qp,2)+pow(mu,2))*(cos(2*atan((q+qp)/mu))*(-25+12*numbers::egamma_v<long double>+6*log(1+pow(q+qp,2)/pow(mu,2))+12*log(mu))-12*atan((q+qp)/mu)*sin(2*atan((q+qp)/mu)))-12*pow(qp,2)*(pow(q+qp,2)+pow(mu,2))*(cos(2*atan((q+qp)/mu))*(-25+12*numbers::egamma_v<long double>+6*log(1+pow(q+qp,2)/pow(mu,2))+12*log(mu))-12*atan((q+qp)/mu)*sin(2*atan((q+qp)/mu)))+36*q*pow(1+pow(q+qp,2)/pow(mu,2),1.5l)*pow(mu,3)*sin(3*atan((q+qp)/mu))+36*qp*pow(1+pow(q+qp,2)/pow(mu,2),1.5l)*pow(mu,3)*sin(3*atan((q+qp)/mu))+3*sqrt(1+pow(q+qp,2)/pow(mu,2))*pow(mu,2)*(pow(q+qp,2)+pow(mu,2))*(cos(3*atan((q+qp)/mu))*(-25+12*numbers::egamma_v<long double>+6*log(1+pow(q+qp,2)/pow(mu,2))+12*log(mu))-12*atan((q+qp)/mu)*sin(3*atan((q+qp)/mu)))+12*q*sqrt(1+pow(q+qp,2)/pow(mu,2))*mu*(pow(q+qp,2)+pow(mu,2))*(12*atan((q+qp)/mu)*cos(3*atan((q+qp)/mu))+(-25+12*numbers::egamma_v<long double>+6*log(1+pow(q+qp,2)/pow(mu,2))+12*log(mu))*sin(3*atan((q+qp)/mu)))+12*qp*sqrt(1+pow(q+qp,2)/pow(mu,2))*mu*(pow(q+qp,2)+pow(mu,2))*(12*atan((q+qp)/mu)*cos(3*atan((q+qp)/mu))+(-25+12*numbers::egamma_v<long double>+6*log(1+pow(q+qp,2)/pow(mu,2))+12*log(mu))*sin(3*atan((q+qp)/mu)))-3*(q+qp)*sqrt(1+pow(q+qp,2)/pow(mu,2))*mu*(pow(q+qp,2)+pow(mu,2))*(12*atan((q+qp)/mu)*cos(3*atan((q+qp)/mu))+(-25+12*numbers::egamma_v<long double>+6*log(1+pow(q+qp,2)/pow(mu,2))+12*log(mu))*sin(3*atan((q+qp)/mu)))+12*pow(q,2)*(pow(q-qp,2)+pow(mu,2))*(cos(2*atan(abs(q-qp)/mu))*(-25+12*numbers::egamma_v<long double>+12*log(mu)+6*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2)))-12*atan(abs(q-qp)/mu)*sin(2*atan(abs(q-qp)/mu)))-36*q*qp*(pow(q-qp,2)+pow(mu,2))*(cos(2*atan(abs(q-qp)/mu))*(-25+12*numbers::egamma_v<long double>+12*log(mu)+6*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2)))-12*atan(abs(q-qp)/mu)*sin(2*atan(abs(q-qp)/mu)))+12*pow(qp,2)*(pow(q-qp,2)+pow(mu,2))*(cos(2*atan(abs(q-qp)/mu))*(-25+12*numbers::egamma_v<long double>+12*log(mu)+6*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2)))-12*atan(abs(q-qp)/mu)*sin(2*atan(abs(q-qp)/mu)))-(36*q*(q-qp)*pow(pow(q-qp,2)+pow(mu,2),1.5l)*sin(3*atan(abs(q-qp)/mu)))/abs(q-qp)+(36*(q-qp)*qp*pow(pow(q-qp,2)+pow(mu,2),1.5l)*sin(3*atan(abs(q-qp)/mu)))/abs(q-qp)-3*mu*pow(pow(q-qp,2)+pow(mu,2),1.5l)*(cos(3*atan(abs(q-qp)/mu))*(-25+12*numbers::egamma_v<long double>+12*log(mu)+6*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2)))-12*atan(abs(q-qp)/mu)*sin(3*atan(abs(q-qp)/mu)))-(12*q*(q-qp)*pow(pow(q-qp,2)+pow(mu,2),1.5l)*(12*atan(abs(q-qp)/mu)*cos(3*atan(abs(q-qp)/mu))+(-25+12*numbers::egamma_v<long double>+12*log(mu)+6*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2)))*sin(3*atan(abs(q-qp)/mu))))/abs(q-qp)+(12*(q-qp)*qp*pow(pow(q-qp,2)+pow(mu,2),1.5l)*(12*atan(abs(q-qp)/mu)*cos(3*atan(abs(q-qp)/mu))+(-25+12*numbers::egamma_v<long double>+12*log(mu)+6*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2)))*sin(3*atan(abs(q-qp)/mu))))/abs(q-qp)+3*pow(pow(q-qp,2)+pow(mu,2),1.5l)*abs(q-qp)*(12*atan(abs(q-qp)/mu)*cos(3*atan(abs(q-qp)/mu))+(-25+12*numbers::egamma_v<long double>+12*log(mu)+6*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2)))*sin(3*atan(abs(q-qp)/mu))))/(192*pow(q*qp,3)));
}

long double f3(long double q, long double qp)
{
	if(q == 0 || qp == 0)
		return(0);
	long double mu = .05;
	if(q == qp)
		return(-(-28704*pow(q,6)+2736*pow(q,4)*pow(mu,2)+660*pow(q,2)*pow(mu,4)-30096*pow(q,4)*(pow(2*q,2)+pow(mu,2))*cos(2*atan((2*q)/mu))+5940*pow(q,2)*pow(pow(mu,2)+pow(2*q,2),2)*cos(4*atan((2*q)/mu))+48*pow(q,6)*(-49+20*numbers::egamma_v<long double>+20*log(mu))+144*pow(q,4)*pow(mu,2)*(-49+20*numbers::egamma_v<long double>+20*log(mu))+90*pow(q,2)*pow(mu,4)*(-49+20*numbers::egamma_v<long double>+20*log(mu))+15*mu*pow(mu,5)*(-49+20*numbers::egamma_v<long double>+20*log(mu))+48*pow(q,6)*(-49+20*numbers::egamma_v<long double>+10*log(pow(mu,2)+pow(2*q,2)))-576*pow(q,5)*(20*mu*atan((2*q)/mu)+2*q*(-49+20*numbers::egamma_v<long double>+10*log(pow(mu,2)+pow(2*q,2))))-1584*pow(q,4)*(pow(2*q,2)+pow(mu,2))*(cos(2*atan((2*q)/mu))*(-49+20*numbers::egamma_v<long double>+10*log(pow(mu,2)+pow(2*q,2)))-20*atan((2*q)/mu)*sin(2*atan((2*q)/mu)))+20720*pow(q,3)*pow(pow(mu,2)+pow(2*q,2),1.5l)*sin(3*atan((2*q)/mu))+1680*pow(q,3)*sqrt(pow(mu,2)+pow(2*q,2))*(pow(2*q,2)+pow(mu,2))*(20*atan((2*q)/mu)*cos(3*atan((2*q)/mu))+(-49+20*numbers::egamma_v<long double>+10*log(pow(mu,2)+pow(2*q,2)))*sin(3*atan((2*q)/mu)))+810*pow(q,2)*pow(pow(mu,2)+pow(2*q,2),2)*(cos(4*atan((2*q)/mu))*(-49+20*numbers::egamma_v<long double>+10*log(pow(mu,2)+pow(2*q,2)))-20*atan((2*q)/mu)*sin(4*atan((2*q)/mu)))-600*q*pow(pow(mu,2)+pow(2*q,2),2.5l)*sin(5*atan((2*q)/mu))-15*sqrt(1+(pow(2*q,2))/pow(mu,2))*pow(pow(2*q,2)*mu+pow(mu,3),2)*(cos(5*atan((2*q)/mu))*(-49+20*numbers::egamma_v<long double>+10*log(pow(mu,2)+pow(2*q,2)))-20*atan((2*q)/mu)*sin(5*atan((2*q)/mu)))-150*q*sqrt(pow(mu,2)+pow(2*q,2))*pow(pow(mu,2)+pow(2*q,2),2)*(20*atan((2*q)/mu)*cos(5*atan((2*q)/mu))+(-49+20*numbers::egamma_v<long double>+10*log(pow(mu,2)+pow(2*q,2)))*sin(5*atan((2*q)/mu))))/(1920*pow(q,8)));
	return((-(87.l/20.l)*pow(q,3)*(q-qp)*pow(qp,2)-(49*pow(q*qp,3))/20.l+87.l/20.l*pow(q,2)*(q-qp)*pow(qp,3)+87.l/20.l*pow(q,3)*pow(qp,2)*(q+qp)+87.l/20.l*pow(q,2)*pow(qp,3)*(q+qp)+57.l/16.l*pow(q,3)*qp*(pow(q+qp,2)+pow(mu,2))*cos(2*atan((q+qp)/mu))+171.l/20.l*pow(q*qp,2)*(pow(q+qp,2)+pow(mu,2))*cos(2*atan((q+qp)/mu))+57.l/16.l*q*pow(qp,3)*(pow(q+qp,2)+pow(mu,2))*cos(2*atan((q+qp)/mu))-11.l/16.l*pow(q,2)*pow(1+pow(q+qp,2)/pow(mu,2),2)*pow(mu,4)*cos(4*atan((q+qp)/mu))-55.l/32.l*q*qp*pow(1+pow(q+qp,2)/pow(mu,2),2)*pow(mu,4)*cos(4*atan((q+qp)/mu))-11.l/16.l*pow(qp,2)*pow(1+pow(q+qp,2)/pow(mu,2),2)*pow(mu,4)*cos(4*atan((q+qp)/mu))+57.l/16.l*pow(q,3)*qp*(pow(q-qp,2)+pow(mu,2))*cos(2*atan(abs(q-qp)/mu))-171.l/20.l*pow(q*qp,2)*(pow(q-qp,2)+pow(mu,2))*cos(2*atan(abs(q-qp)/mu))+57.l/16.l*q*pow(qp,3)*(pow(q-qp,2)+pow(mu,2))*cos(2*atan(abs(q-qp)/mu))+11.l/16.l*pow(q,2)*pow(pow(q-qp,2)+pow(mu,2),2)*cos(4*atan(abs(q-qp)/mu))-55.l/32.l*q*qp*pow(pow(q-qp,2)+pow(mu,2),2)*cos(4*atan(abs(q-qp)/mu))+11.l/16.l*pow(qp,2)*pow(pow(q-qp,2)+pow(mu,2),2)*cos(4*atan(abs(q-qp)/mu))-1.l/40.l*pow(q*qp,3)*(-49+20*numbers::egamma_v<long double>+20*log(mu)+10*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2)))-(3*pow(q,3)*(q-qp)*pow(qp,2)*(-49*abs(q-qp)+20*numbers::egamma_v<long double>*abs(q-qp)+20*mu*atan(abs(q-qp)/mu)+20*abs(q-qp)*log(mu)+10*abs(q-qp)*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))))/(20*abs(q-qp))+(3*pow(q,2)*(q-qp)*pow(qp,3)*(-49*abs(q-qp)+20*numbers::egamma_v<long double>*abs(q-qp)+20*mu*atan(abs(q-qp)/mu)+20*abs(q-qp)*log(mu)+10*abs(q-qp)*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))))/(20*abs(q-qp))-1.l/40.l*pow(q*qp,3)*(-49+20*numbers::egamma_v<long double>+20*log(mu)+10*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2)))+3.l/20.l*pow(q,3)*pow(qp,2)*(-49*q+20*numbers::egamma_v<long double>*q-49*qp+20*numbers::egamma_v<long double>*qp+20*mu*atan((q+qp)/mu)+20*q*log(mu)+20*qp*log(mu)+10*q*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))+10*qp*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2)))+3.l/20.l*pow(q,2)*pow(qp,3)*(-49*q+20*numbers::egamma_v<long double>*q-49*qp+20*numbers::egamma_v<long double>*qp+20*mu*atan((q+qp)/mu)+20*q*log(mu)+20*qp*log(mu)+10*q*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))+10*qp*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2)))+3.l/16.l*pow(q,3)*qp*(pow(q+qp,2)+pow(mu,2))*(-49*cos(2*atan((q+qp)/mu))+20*numbers::egamma_v<long double>*cos(2*atan((q+qp)/mu))+20*cos(2*atan((q+qp)/mu))*log(mu)+10*cos(2*atan((q+qp)/mu))*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))-20*atan((q+qp)/mu)*sin(2*atan((q+qp)/mu)))+9/20.l*pow(q*qp,2)*(pow(q+qp,2)+pow(mu,2))*(-49*cos(2*atan((q+qp)/mu))+20*numbers::egamma_v<long double>*cos(2*atan((q+qp)/mu))+20*cos(2*atan((q+qp)/mu))*log(mu)+10*cos(2*atan((q+qp)/mu))*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))-20*atan((q+qp)/mu)*sin(2*atan((q+qp)/mu)))+3.l/16.l*q*pow(qp,3)*(pow(q+qp,2)+pow(mu,2))*(-49*cos(2*atan((q+qp)/mu))+20*numbers::egamma_v<long double>*cos(2*atan((q+qp)/mu))+20*cos(2*atan((q+qp)/mu))*log(mu)+10*cos(2*atan((q+qp)/mu))*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))-20*atan((q+qp)/mu)*sin(2*atan((q+qp)/mu)))-37.l/48.l*pow(q,3)*pow(1+pow(q+qp,2)/pow(mu,2),1.5l)*pow(mu,3)*sin(3*atan((q+qp)/mu))-37.l/8.l*pow(q,2)*qp*pow(1+pow(q+qp,2)/pow(mu,2),1.5l)*pow(mu,3)*sin(3*atan((q+qp)/mu))-37.l/8.l*q*pow(qp,2)*pow(1+pow(q+qp,2)/pow(mu,2),1.5l)*pow(mu,3)*sin(3*atan((q+qp)/mu))-37.l/48.l*pow(qp,3)*pow(1+pow(q+qp,2)/pow(mu,2),1.5l)*pow(mu,3)*sin(3*atan((q+qp)/mu))-1.l/16.l*pow(q,3)*sqrt(1+pow(q+qp,2)/pow(mu,2))*mu*(pow(q+qp,2)+pow(mu,2))*(20*atan((q+qp)/mu)*cos(3*atan((q+qp)/mu))-49*sin(3*atan((q+qp)/mu))+20*numbers::egamma_v<long double>*sin(3*atan((q+qp)/mu))+20*log(mu)*sin(3*atan((q+qp)/mu))+10*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))*sin(3*atan((q+qp)/mu)))-3.l/8.l*pow(q,2)*qp*sqrt(1+pow(q+qp,2)/pow(mu,2))*mu*(pow(q+qp,2)+pow(mu,2))*(20*atan((q+qp)/mu)*cos(3*atan((q+qp)/mu))-49*sin(3*atan((q+qp)/mu))+20*numbers::egamma_v<long double>*sin(3*atan((q+qp)/mu))+20*log(mu)*sin(3*atan((q+qp)/mu))+10*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))*sin(3*atan((q+qp)/mu)))-3.l/8.l*q*pow(qp,2)*sqrt(1+pow(q+qp,2)/pow(mu,2))*mu*(pow(q+qp,2)+pow(mu,2))*(20*atan((q+qp)/mu)*cos(3*atan((q+qp)/mu))-49*sin(3*atan((q+qp)/mu))+20*numbers::egamma_v<long double>*sin(3*atan((q+qp)/mu))+20*log(mu)*sin(3*atan((q+qp)/mu))+10*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))*sin(3*atan((q+qp)/mu)))-1.l/16.l*pow(qp,3)*sqrt(1+pow(q+qp,2)/pow(mu,2))*mu*(pow(q+qp,2)+pow(mu,2))*(20*atan((q+qp)/mu)*cos(3*atan((q+qp)/mu))-49*sin(3*atan((q+qp)/mu))+20*numbers::egamma_v<long double>*sin(3*atan((q+qp)/mu))+20*log(mu)*sin(3*atan((q+qp)/mu))+10*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))*sin(3*atan((q+qp)/mu)))-3.l/32.l*pow(q,2)*pow(pow(q+qp,2)+pow(mu,2),2)*(-49*cos(4*atan((q+qp)/mu))+20*numbers::egamma_v<long double>*cos(4*atan((q+qp)/mu))+20*cos(4*atan((q+qp)/mu))*log(mu)+10*cos(4*atan((q+qp)/mu))*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))-20*atan((q+qp)/mu)*sin(4*atan((q+qp)/mu)))-15.l/64.l*q*qp*pow(pow(q+qp,2)+pow(mu,2),2)*(-49*cos(4*atan((q+qp)/mu))+20*numbers::egamma_v<long double>*cos(4*atan((q+qp)/mu))+20*cos(4*atan((q+qp)/mu))*log(mu)+10*cos(4*atan((q+qp)/mu))*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))-20*atan((q+qp)/mu)*sin(4*atan((q+qp)/mu)))-3.l/32.l*pow(qp,2)*pow(pow(q+qp,2)+pow(mu,2),2)*(-49*cos(4*atan((q+qp)/mu))+20*numbers::egamma_v<long double>*cos(4*atan((q+qp)/mu))+20*cos(4*atan((q+qp)/mu))*log(mu)+10*cos(4*atan((q+qp)/mu))*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))-20*atan((q+qp)/mu)*sin(4*atan((q+qp)/mu)))+5.l/32.l*q*pow(1+pow(q+qp,2)/pow(mu,2),2.5l)*pow(mu,5)*sin(5*atan((q+qp)/mu))+5.l/32.l*qp*pow(1+pow(q+qp,2)/pow(mu,2),2.5l)*pow(mu,5)*sin(5*atan((q+qp)/mu))+1.l/128.l*sqrt(1+pow(q+qp,2)/pow(mu,2))*pow(mu,2)*pow(pow(q+qp,2)+pow(mu,2),2)*(-49*cos(5*atan((q+qp)/mu))+20*numbers::egamma_v<long double>*cos(5*atan((q+qp)/mu))+20*cos(5*atan((q+qp)/mu))*log(mu)+10*cos(5*atan((q+qp)/mu))*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))-20*atan((q+qp)/mu)*sin(5*atan((q+qp)/mu)))+3.l/64.l*q*sqrt(1+pow(q+qp,2)/pow(mu,2))*mu*pow(pow(q+qp,2)+pow(mu,2),2)*(20*atan((q+qp)/mu)*cos(5*atan((q+qp)/mu))-49*sin(5*atan((q+qp)/mu))+20*numbers::egamma_v<long double>*sin(5*atan((q+qp)/mu))+20*log(mu)*sin(5*atan((q+qp)/mu))+10*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))*sin(5*atan((q+qp)/mu)))+3.l/64.l*qp*sqrt(1+pow(q+qp,2)/pow(mu,2))*mu*pow(pow(q+qp,2)+pow(mu,2),2)*(20*atan((q+qp)/mu)*cos(5*atan((q+qp)/mu))-49*sin(5*atan((q+qp)/mu))+20*numbers::egamma_v<long double>*sin(5*atan((q+qp)/mu))+20*log(mu)*sin(5*atan((q+qp)/mu))+10*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))*sin(5*atan((q+qp)/mu)))-1.l/128.l*(q+qp)*sqrt(1+pow(q+qp,2)/pow(mu,2))*mu*pow(pow(q+qp,2)+pow(mu,2),2)*(20*atan((q+qp)/mu)*cos(5*atan((q+qp)/mu))-49*sin(5*atan((q+qp)/mu))+20*numbers::egamma_v<long double>*sin(5*atan((q+qp)/mu))+20*log(mu)*sin(5*atan((q+qp)/mu))+10*log((pow(q+qp,2)+pow(mu,2))/pow(mu,2))*sin(5*atan((q+qp)/mu)))+3.l/16.l*pow(q,3)*qp*(pow(q-qp,2)+pow(mu,2))*(-49*cos(2*atan(abs(q-qp)/mu))+20*numbers::egamma_v<long double>*cos(2*atan(abs(q-qp)/mu))+20*cos(2*atan(abs(q-qp)/mu))*log(mu)+10*cos(2*atan(abs(q-qp)/mu))*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))-20*atan(abs(q-qp)/mu)*sin(2*atan(abs(q-qp)/mu)))-9/20.l*pow(q*qp,2)*(pow(q-qp,2)+pow(mu,2))*(-49*cos(2*atan(abs(q-qp)/mu))+20*numbers::egamma_v<long double>*cos(2*atan(abs(q-qp)/mu))+20*cos(2*atan(abs(q-qp)/mu))*log(mu)+10*cos(2*atan(abs(q-qp)/mu))*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))-20*atan(abs(q-qp)/mu)*sin(2*atan(abs(q-qp)/mu)))+3.l/16.l*q*pow(qp,3)*(pow(q-qp,2)+pow(mu,2))*(-49*cos(2*atan(abs(q-qp)/mu))+20*numbers::egamma_v<long double>*cos(2*atan(abs(q-qp)/mu))+20*cos(2*atan(abs(q-qp)/mu))*log(mu)+10*cos(2*atan(abs(q-qp)/mu))*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))-20*atan(abs(q-qp)/mu)*sin(2*atan(abs(q-qp)/mu)))+(37*pow(q,3)*(q-qp)*pow(pow(q-qp,2)+pow(mu,2),1.5l)*sin(3*atan(abs(q-qp)/mu)))/(48*abs(q-qp))-(37*pow(q,2)*(q-qp)*qp*pow(pow(q-qp,2)+pow(mu,2),1.5l)*sin(3*atan(abs(q-qp)/mu)))/(8*abs(q-qp))+(37*q*(q-qp)*pow(qp,2)*pow(pow(q-qp,2)+pow(mu,2),1.5l)*sin(3*atan(abs(q-qp)/mu)))/(8*abs(q-qp))-(37*(q-qp)*pow(qp,3)*pow(pow(q-qp,2)+pow(mu,2),1.5l)*sin(3*atan(abs(q-qp)/mu)))/(48*abs(q-qp))+1/(16*abs(q-qp))*pow(q,3)*(q-qp)*pow(pow(q-qp,2)+pow(mu,2),1.5l)*(20*atan(abs(q-qp)/mu)*cos(3*atan(abs(q-qp)/mu))-49*sin(3*atan(abs(q-qp)/mu))+20*numbers::egamma_v<long double>*sin(3*atan(abs(q-qp)/mu))+20*log(mu)*sin(3*atan(abs(q-qp)/mu))+10*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))*sin(3*atan(abs(q-qp)/mu)))-1/(8*abs(q-qp))*3*pow(q,2)*(q-qp)*qp*pow(pow(q-qp,2)+pow(mu,2),1.5l)*(20*atan(abs(q-qp)/mu)*cos(3*atan(abs(q-qp)/mu))-49*sin(3*atan(abs(q-qp)/mu))+20*numbers::egamma_v<long double>*sin(3*atan(abs(q-qp)/mu))+20*log(mu)*sin(3*atan(abs(q-qp)/mu))+10*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))*sin(3*atan(abs(q-qp)/mu)))+1/(8*abs(q-qp))*3*q*(q-qp)*pow(qp,2)*pow(pow(q-qp,2)+pow(mu,2),1.5l)*(20*atan(abs(q-qp)/mu)*cos(3*atan(abs(q-qp)/mu))-49*sin(3*atan(abs(q-qp)/mu))+20*numbers::egamma_v<long double>*sin(3*atan(abs(q-qp)/mu))+20*log(mu)*sin(3*atan(abs(q-qp)/mu))+10*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))*sin(3*atan(abs(q-qp)/mu)))-1/(16*abs(q-qp))*(q-qp)*pow(qp,3)*pow(pow(q-qp,2)+pow(mu,2),1.5l)*(20*atan(abs(q-qp)/mu)*cos(3*atan(abs(q-qp)/mu))-49*sin(3*atan(abs(q-qp)/mu))+20*numbers::egamma_v<long double>*sin(3*atan(abs(q-qp)/mu))+20*log(mu)*sin(3*atan(abs(q-qp)/mu))+10*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))*sin(3*atan(abs(q-qp)/mu)))+3.l/32.l*pow(q,2)*pow(pow(q-qp,2)+pow(mu,2),2)*(-49*cos(4*atan(abs(q-qp)/mu))+20*numbers::egamma_v<long double>*cos(4*atan(abs(q-qp)/mu))+20*cos(4*atan(abs(q-qp)/mu))*log(mu)+10*cos(4*atan(abs(q-qp)/mu))*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))-20*atan(abs(q-qp)/mu)*sin(4*atan(abs(q-qp)/mu)))-15.l/64.l*q*qp*pow(pow(q-qp,2)+pow(mu,2),2)*(-49*cos(4*atan(abs(q-qp)/mu))+20*numbers::egamma_v<long double>*cos(4*atan(abs(q-qp)/mu))+20*cos(4*atan(abs(q-qp)/mu))*log(mu)+10*cos(4*atan(abs(q-qp)/mu))*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))-20*atan(abs(q-qp)/mu)*sin(4*atan(abs(q-qp)/mu)))+3.l/32.l*pow(qp,2)*pow(pow(q-qp,2)+pow(mu,2),2)*(-49*cos(4*atan(abs(q-qp)/mu))+20*numbers::egamma_v<long double>*cos(4*atan(abs(q-qp)/mu))+20*cos(4*atan(abs(q-qp)/mu))*log(mu)+10*cos(4*atan(abs(q-qp)/mu))*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))-20*atan(abs(q-qp)/mu)*sin(4*atan(abs(q-qp)/mu)))-(5*q*(q-qp)*pow(pow(q-qp,2)+pow(mu,2),2.5l)*sin(5*atan(abs(q-qp)/mu)))/(32*abs(q-qp))+(5*(q-qp)*qp*pow(pow(q-qp,2)+pow(mu,2),2.5l)*sin(5*atan(abs(q-qp)/mu)))/(32*abs(q-qp))-1.l/128.l*mu*pow(pow(q-qp,2)+pow(mu,2),2.5l)*(-49*cos(5*atan(abs(q-qp)/mu))+20*numbers::egamma_v<long double>*cos(5*atan(abs(q-qp)/mu))+20*cos(5*atan(abs(q-qp)/mu))*log(mu)+10*cos(5*atan(abs(q-qp)/mu))*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))-20*atan(abs(q-qp)/mu)*sin(5*atan(abs(q-qp)/mu)))-1/(64*abs(q-qp))*3*q*(q-qp)*pow(pow(q-qp,2)+pow(mu,2),2.5l)*(20*atan(abs(q-qp)/mu)*cos(5*atan(abs(q-qp)/mu))-49*sin(5*atan(abs(q-qp)/mu))+20*numbers::egamma_v<long double>*sin(5*atan(abs(q-qp)/mu))+20*log(mu)*sin(5*atan(abs(q-qp)/mu))+10*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))*sin(5*atan(abs(q-qp)/mu)))+1/(64*abs(q-qp))*3*(q-qp)*qp*pow(pow(q-qp,2)+pow(mu,2),2.5l)*(20*atan(abs(q-qp)/mu)*cos(5*atan(abs(q-qp)/mu))-49*sin(5*atan(abs(q-qp)/mu))+20*numbers::egamma_v<long double>*sin(5*atan(abs(q-qp)/mu))+20*log(mu)*sin(5*atan(abs(q-qp)/mu))+10*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))*sin(5*atan(abs(q-qp)/mu)))+1.l/128.l*pow(pow(q-qp,2)+pow(mu,2),2.5l)*abs(q-qp)*(20*atan(abs(q-qp)/mu)*cos(5*atan(abs(q-qp)/mu))-49*sin(5*atan(abs(q-qp)/mu))+20*numbers::egamma_v<long double>*sin(5*atan(abs(q-qp)/mu))+20*log(mu)*sin(5*atan(abs(q-qp)/mu))+10*log((pow(q-qp,2)+pow(mu,2))/pow(mu,2))*sin(5*atan(abs(q-qp)/mu))))/(pow(q*qp,4)));
}

long double f4(long double q, long double qp)
{
	if(q == 0 || qp == 0)
		return(0);
	return((-5*q*qp*(pow(q,2)+pow(qp,2))*(21*pow(q,4)-2*pow(q,2)*pow(qp,2)+21*pow(qp,4))+.75l*(35*pow(q,8)+20*pow(q,6)*pow(qp,2)+18*pow(q,4)*pow(qp,4)+20*pow(q,2)*pow(qp,6)+35*pow(qp,8))*log((pow(q+qp,2)+.0025)/(pow(q-qp,2)+.0025)))/(384*pow(q*qp,5)));
}

long double f5(long double q, long double qp)
{
	if(q == 0 || qp == 0)
		return(0);
	return((-q*qp*(945*pow(q,8)+840*pow(q,6)*pow(qp,2)+814*pow(q,4)*pow(qp,4)+840*pow(q,2)*pow(qp,6)+945*pow(qp,8))+3.75l*(pow(q,2)+pow(qp,2))*(63*pow(q,8)-28*pow(q,6)*pow(qp,2)+58*pow(q,4)*pow(qp,4)-28*pow(q,2)*pow(qp,6)+63*pow(qp,8))*log((pow(q+qp,2)+.0025)/(pow(q-qp,2)+.0025)))/(3840*pow(q*qp,6)));
}

bool This_one;

int main()
{
	Around<long double> real;
	Around<complex<long double>> Ans;
	for(int l = 0; l < 10; l++)
	{
		for(int lp = l; lp < 10; lp++)
		{
			for(long double q = 0; q < 10; q++)
			{
				for(long double qp = q; qp < 10; qp++)
				{
					Ans = Int_phi(l,lp,0,q,qp);
					real = Around<long double>(Ans.Value().real(), Ans.Error().real());
					if(l == lp && lp == 0)
						cout << l << " " << lp << " " << q << " " << qp << setw(12) << f0(q,qp) << setw(14) << real << setw(14) << (real.Value()/f0(q,qp)-1.l)*100.l << setw(14) << real.RelErr() << endl;
					else if(l == lp && lp == 1)
						cout << l << " " << lp << " " << q << " " << qp << setw(12) << f1(q,qp) << setw(14) << real << setw(14) << (real.Value()/f1(q,qp)-1.l)*100.l << setw(14) << real.RelErr() << endl;
					else if(l == lp && lp == 2)
						cout << l << " " << lp << " " << q << " " << qp << setw(12) << f2(q,qp) << setw(14) << real << setw(14) << (real.Value()/f2(q,qp)-1.l)*100.l << setw(14) << real.RelErr() << endl;
					else if(l == lp && lp == 3)
						cout << l << " " << lp << " " << q << " " << qp << setw(12) << f3(q,qp) << setw(14) << real << setw(14) << (real.Value()/f3(q,qp)-1.l)*100.l << setw(14) << real.RelErr() << endl;
					else if(l == lp && lp >= 4)
						cout << l << " " << lp << " " << q << " " << qp << setw(12) << "not 0" << setw(14) << real << setw(12) << real.RelErr() << endl;
					else if(l != lp)
						cout << l << " " << lp << " " << q << " " << qp << setw(12) << "0" << setw(14) << real << setw(12) << real.RelErr() << endl;
				}
			}
		}
	}

	return(0);
}

Around<complex<long double>> Int_phi(int l, int lp, int m, long double q, long double qp)	//Phi integral is analytic.
{
	if((l > 0 || lp > 0) && (q == 0 || qp == 0))
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
	long double delta_r = 9.l*pi_v<long double>/(q+qp);
	static Around<complex<long double>> Answer;
	Around<complex<long double>> Temp;
	int i = 0;

	if(!This_one)
		return(Answer*Y(l, m, theta, phi)*conj(Y(lp, m, theta, phi)));

	if(isinf(delta_r) || isnan(delta_r))
	{
		r0 = 1;
		delta_r = 1;
	}

	long double Expected;
	Expected = Int_r_recurs(l, lp, m, q, qp, phi, theta, 0, 0, r0).Value().real();
	do
	{
		Temp = Int_r_recurs(l, lp, m, q, qp, phi, theta, 0, r0, r0+delta_r);
		Expected += Temp.Value().real();
		i++;
		r0 += delta_r;
	}while(r0-delta_r < 500);
	Expected += Int_r_recurs(l, lp, m, q, qp, phi, theta, 0, r0, 500).Value().real();

	r0 = pi_v<long double>/fmax(q,qp)+(long double)(max(l,lp));
	if(isinf(r0) || isnan(r0))
		r0 = 1;
	Answer = Int_r_recurs(l, lp, m, q, qp, phi, theta, 0, 0, r0);
	do
	{
		Temp = Int_r_recurs(l, lp, m, q, qp, phi, theta, 0, r0, r0+delta_r);
		Answer += Temp;
		i++;
		r0 += delta_r;
		if(abs(Answer.Value().real()/Expected-1.l) < 1e-8)
		{
			cerr << l << " " << lp << " " << q << " " << qp << " " << Answer << " " << Answer.Value().real()/Expected-1.l << " " << r0-delta_r << " " << i << " " << (r0-delta_r)*(q+qp)/pi_v<long double> << endl;
			break;
		}
	}while(((Temp/Answer).Value().real() > 1e-15 || (Temp/Answer).Value().imag() > 1e-15) || r0-delta_r < 500);

	return(Answer*Y(l, m, theta, phi)*conj(Y(lp, m, theta, phi)));
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
	return(j(l, r*q)*j(lp, r*qp)*exp(-r/20.l)/r);
	//return(Y(l, m, theta, phi)*conj(Y(lp, m, theta, phi))*j(l, r*q)*j(lp, r*qp)*exp(-r/20.l)/r);
}

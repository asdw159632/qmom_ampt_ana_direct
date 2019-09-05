#include <string.h>
#include <complex.h>
//#include <complex>
//#include "pp.h"
#include "pair.h"
//#include "k+k+.h"
double q_xx[3][3][1000];
double qinv_xx[3][3][1000];
double totalq[3][3];

#define erfi

double B(double p,double R,double a)	//Sinyukov
{
	double B=1+2*R/a+3*M_PI*R*R/(4*a*a);
	return B;
}

void readqinv(int j,int par,TH1D* q_inv)
{
	totalq[j][(int)floor(par)]=0;
        for(int i=0;i<q_inv->GetXaxis()->GetNbins();i++){
	        q_xx[j][(int)floor(par)][i]=q_inv->GetBinCenter(i+1);
		qinv_xx[j][(int)floor(par)][i]=q_inv->GetBinContent(i+1);
        }
        totalq[j][(int)floor(par)]=q_inv->GetXaxis()->GetNbins();
}

double getqinv(double x)
{
	int total=totalq[par6][par5];
        for(int i=0;i<total;i++)
        {
                if(i<total-1){
                        if(x<(q_xx[par6][(int)floor(par5)][i]+q_xx[par6][(int)floor(par5)][i+1])/2)
                                return qinv_xx[par6][(int)floor(par5)][i];
                        else
                                continue;
                }
                else
                        return qinv_xx[par6][(int)floor(par5)][i];
        }
        return 0;
}

double erf(double x)
{
#ifdef erfi
	double result;
	double a;
	double b;
	double c;
	x/=sqrt(2);
	a=1;
	b=3;
	c=2/sqrt(M_PI)*x;
	for(int i=1;i<=40;i++)
	{
		c*=(x*x);
		c/=a;
		result=c/b;
		a++;
		b+=2;
	}
	result=result*result+1;
	return result;
#else
	return 1;
#endif
}

Double_t fit1(Double_t *x,Double_t *par)
{
        double xx;
	if(par5==3){
		xx=x[0]/2;
	}
	else
		xx=getqinv(x[0])/2;
        double p1=par[0];
        double p2=par[1]*fm_gev/2;
        double R0=par[2];
        double n=par[3];
        double A=par[4];
	Double_t co;
#ifdef qinv_half
	xx*=2;
	p2*=2
#endif
#ifdef Coulomb
	Double_t R=R0*fm_gev*4/sqrt(M_PI);
	Double_t eta=1/(xx*a);
	Double_t Ac=2*M_PI*eta/(exp(2*M_PI*eta)-1);
	Double_t kl=M_PI/(4*R)*(1+2*R/a);
	Double_t d1=1-kl*kl/(xx*xx)*(1-Ac*B(xx,R,a));
	//Double_t d2=1-kl*kl/(xx*xx)*(1-Ac*B(xx,-R,a));
	double d=(d1);
	if(xx<=kl)
		co=(1-p1+p1*(Ac*B(xx,R,a))*(1+A*exp(-(x[0]*p2*x[0]*p2))*erf(x[0]*p2)))*n;
	else
		co=(1-p1+p1*d*(1+A*exp(-(x[0]*p2*x[0]*p2))*erf(x[0]*p2)))*n;
#else
	co=1+p1*exp(-(x[0]*x[0]*p2*p2))*erf(x[0]*p2);
#endif
	return co;
}

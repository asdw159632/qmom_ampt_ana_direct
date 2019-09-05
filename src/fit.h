#include <string.h>
#include <complex.h>
//#include <complex>
//#include "pp.h"
#include "pair.h"
//#include "k+k+.h"
double _Complex ci=I;
double q_xx[4][1000];
double qinv_xx[4][1000];
double totalq;
#define Sinyukov

#ifdef Sinyukov
double B(double p,double R,double a)	//Sinyukov
{
	double B=1+2*R/a+3*M_PI*R*R/(4*a*a);
	return B;
}
#else

double w(int n)
{
	double r=1;
	for(int i=1;i<=n;i++)r*=i;
	return r;
}

double B(double p,double R,double a)	//Bowler
{
	int nmax=8;
	R/=(4/sqrt(M_PI));
	double result;
	double eta=1/(p*a);
	double_complex An[nmax];
	double_complex Am[nmax];
	An[0]=1;
	Am[0]=1;
	for(int i=1;i<nmax;i++){
		An[i]=(i-1+ci*eta)/(i*i)*An[i-1]*pow(ci,i)*pow(p*R,i);
		Am[i]=(i-1-ci*eta)/(i*i)*Am[i-1]*pow(-ci,i)*pow(p*R,i);
	}
	double_complex B=1;
	for(int i=1;i<nmax;i++){
		for(int j=1;j<nmax;j++){
			B+=An[i]*Am[j]*w(j+i+2)/(2*(j+i+1));}}
	result=real(B);
	return result;
}
#endif

void readqinv(double par5)
{
        const char *osl;        //par5
	totalq=0;
#ifdef byrun
	TFile *result_root=new TFile("../result_root/result-all.root","READ");
#else
	TFile *result_root=new TFile("../result_root/result.root","READ");
#endif
	if(floor(par5)==0)osl="00";
        if(floor(par5)==1)osl="ql";
        if(floor(par5)==2)osl="qs";
        if(floor(par5)==3)osl="qo";
        char *dat_xx=new char[100];
        sprintf(dat_xx,"n_%s_qinv",osl);
	TH1D *q_qinv=(TH1D*)result_root->Get(dat_xx);
        for(int i=0;i<q_qinv->GetXaxis()->GetNbins();i++){
        q_xx[(int)floor(par5)][i]=q_qinv->GetBinCenter(i+1);
	qinv_xx[(int)floor(par5)][i]=q_qinv->GetBinContent(i+1);
        }
        totalq=q_qinv->GetXaxis()->GetNbins();
}

double getqinv(double x)
{
	if(par5!=0)
	{
	        for(int i=0;i<totalq;i++)
	        {
	                if(i<totalq-1){
	                        if(x<(q_xx[(int)floor(par5)][i]+q_xx[(int)floor(par5)][i+1])/2)
        	                        return qinv_xx[(int)floor(par5)][i];
        	                else
        	                        continue;
        	        }
        	        else
        	                return qinv_xx[par5][i];
       		}
		return 0;
	}
	return 0;
}


Double_t fit(Double_t *x,Double_t *par)
{
        double xx;
	if(par[5]==0){
		xx=x[0]/2;
	}
	else{
		xx=getqinv(x[0])/2;
	}
        double p1=par[0];
        double p2=par[1]*fm_gev/2;
        double R0=par[2];
        double n=par[3];
        double A=par[4];
	Double_t co;
#ifdef qinv_half
	p2*=2;
	xx*=2
#endif
#ifdef Coulomb
	Double_t eta=1/(xx*a);
	Double_t R=R0*fm_gev*4/sqrt(M_PI);
	Double_t Ac=2*M_PI*eta/(exp(2*M_PI*eta)-1);
	Double_t kl=M_PI/(4*R)*(1+2*R/a);
	Double_t d1=1-kl*kl/(xx*xx)*(1-Ac*B(xx,R,a));
	//Double_t d2=1-kl*kl/(xx*xx)*(1-Ac*B(xx,-R,a));
	double d=(d1);
	if(xx<=kl)
		co=(1-p1+p1*(Ac*B(xx,R,a))*(1+A*exp(-(x[0]*p2*x[0]*p2))))*n;
	else
		co=(1-p1+p1*d*(1+A*exp(-(x[0]*p2*x[0]*p2))))*n;
#else
	co=1+p1*exp(-(x[0]*x[0]*p2*p2));
#endif
	return co;
}

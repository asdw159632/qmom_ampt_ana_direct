#include <fstream>
//#include <complex>
#include <complex.h>

using namespace std;

//#define Sinyukov
#ifdef Sinyukov
//#define simple
#endif
#define fm_gev 5.0676896
#define nmax 10

double B1(double eta,double epsilon);
double _Complex B2(double eta,double epsilon);
double _Complex gamma(double _Complex c);
double _Complex B3(double eta,double epsilon);
TLorentzVector trace(TLorentzVector q,TLorentzVector K);

double B1(double eta,double epsilon)
{
	double R=epsilon*eta;
	double result=1+2.*R+3.*R*R/(2.);
	return result;
}

double B1(double eta,double epsilon,double nepsilon)
{
	double R=epsilon*eta;
	double Rn=nepsilon*eta;
	double _Complex F=1;
	double _Complex F_s=1;
	double _Complex An=1;
	double _Complex Am=1;
	for(int i=1;i<nmax;i++){
		An*=(i-1+eta*I)/(i*i);
		F+=An*cpow(0-I*R,i);
		Am*=(i-1+eta*I)/(i*i);
		F_s+=Am*cpow(0-I*Rn,i);}
	double result=pow(cabs(F+l*F_s),2);
	return result;
}

double _Complex B4(double eta,double epsilon)
{
//	double R=epsilon*eta;
	double _Complex F=1;
	double _Complex An=1;
	for(int i=1;i<nmax;i++){
		An*=(i-1-eta*I)/(i*i);
		F+=An*cpow(I*epsilon,i);}
	double _Complex result=F;
	return result;
}

double _Complex B2(double eta,double epsilon)
{
	double _Complex a=(0-I)*eta;
	double _Complex b=1.;
	double _Complex c=1.;
	double _Complex chyper=1.;
	double _Complex result=0.;
	if(abs(epsilon)<10.0){
		for (double j=1.;cabs(c)>0.00001;j++){
			c*=I*epsilon*a/(b*j);
			chyper+=c;
			a+=1.;
			b+=1.;
			if (j>=40){
				cerr<<"not coverging!"<<endl;
				break;
			}
		}
		result=chyper;
	}
	else result=B3(eta,epsilon);
//	cout<<result<<" ";
	return result;
}

double _Complex gamma(double _Complex c)
{
	double _Complex cg,cphase;
	int mm,j;
	double x,y,phase,delp,cgmag;
	x=creal(c);
	y=cimag(c);
	#define EULER 0.5772156649015328606
	phase=-EULER*y;
	for(j=1;j<=100000;j++){
		delp=(y/(double)j)-atan(y/(double)j);
		phase=phase+delp;
		if(fabs(delp)<1E-10) goto CGAMMA_ESCAPE;
	}
	printf("oops not accurate enough, increase jmax\n");
	CGAMMA_ESCAPE:
	phase=phase-2.0*M_PI*floor(phase/(2.0*M_PI));
	cphase=cexp(I*phase);
	cgmag=sqrt(M_PI*y/sinh(M_PI*y));
	mm=(int)floor(x+0.5);
	cg=cgmag*cphase;
	if(mm<1){
		for(j=1;j<=-mm+1;j++){
		cg=cg/(1.0+(double)(-j)+I*y);
		}
	}
	if(mm>1) {
		for(j=1;j<=mm-1;j++){
		cg=cg*((double)(j)+I*y);
		}
	}
	return cg;
}


double _Complex B3(double eta,double epsilon)
{
	double _Complex a=(0-I)*eta;
	double _Complex b=1.;
	double _Complex c11=1.;
	double _Complex c21=1.;
	double _Complex c12=1.;
	double _Complex c22=1.;
	double _Complex chype1=1.;
	double _Complex chype2=1.;
	double _Complex bot =1.;
	for (double j=1.;j<=5;j++){
		c11*=(j+a-1.);
		c21*=(j-a);
		c12*=(j-b+a);
		c22*=(j+b-a-1);
		bot*=j;
		chype1+=c11*c12/(bot*cpow(-epsilon*I,j));
		chype2+=c21*c22/(bot*cpow(epsilon*I,j));
	}
	chype1*=cpow(-epsilon*I,0-a)*1/gamma(b-a);
	chype2*=cpow(epsilon*I,a-b)*cexp(epsilon*I)*1/gamma(a);
	double _Complex result=chype1+chype2;
	return result;
}

double G(TVector3 p,TVector3 r,double m1,double m2)
{
	double pp=p.Mag();
	double R=r*p/pp+r.Mag();
	double Rn=-r*p/pp+r.Mag();
	if(pp<2.127e-4)return 0;
//	if(pp==0){R=r.Mag();Rn=r.Mag();}
//	pp=(floor(p.Mag()/dq)+0.5)*dq;
	double m=m1*m2/(m1+m2);
	double a=137.036/m;
	double eta=1/(a*pp);
	double Gi=2*M_PI*eta/(exp(2*M_PI*eta)-1);
	double kl=M_PI/(4*R)*(1+2*R/a);
	//double kln=M_PI/(4*Rn)*(1+2*Rn/a);
	double G;
//	kl=10;
//	if(p*r>30)return 1.;
	if(eta*2*M_PI>15)return 0;
	//if(fabs(p.Mag()*r.Mag())>3.00)return Gi*(1.+l*cos(2.*p*r));
	if(pp<=kl){
#ifdef Sinyukov
#ifdef simple
		G=Gi*B1(eta,R*pp);
#else
//		G=Gi*B1(eta,R*pp,Rn*pp);
		G=Gi*creal((B4(eta,R*pp)*cexp(I*(p*r))+l*B4(eta,Rn*pp)*cexp(0-I*(p*r)))
				*(B4(-eta,-R*pp)*cexp(0-I*(p*r))+l*B4(-eta,-Rn*pp)*cexp(I*(p*r))))*0.5;
//		G=Gi*pow(cabs(B4(eta,R*pp)+l*B4(eta,Rn*pp)),2)*0.5;
#endif
#else
//		G=Gi*pow(cabs(B2(eta,R*pp)+l*B2(eta,Rn*pp)),2)*0.5;
		//G=Gi*creal(cpow(B2(eta,R*pp)*cexp(I*(p*r))
		G=Gi*pow(cabs(B2(eta,R*pp)*cexp(I*(p*r))
				+l*B2(eta,Rn*pp)*cexp(0-I*(p*r))),2)*0.5;
#endif
//		if(G>1.01&&pp>0.06)cout<<"G1 "<<G<<" "<<eta<<" "<<R<<" "<<pp<<endl;
	}
	else{
#ifdef Sinyukov
#ifdef simple
		G=1-kl*kl/(pp*pp)*(1-Gi*B1(eta,R*pp));
#else
//		G=1-kl*kl/(pp*pp)*(1-Gi*B1(eta,R*pp,Rn*pp));
		G=1-kl*kl/(pp*pp)*(1-Gi*creal((B4(eta,R*pp)*cexp(I*(p*r))+l*B4(eta,Rn*pp)*cexp(0-I*(p*r)))
				*(B4(-eta,-R*pp)*cexp(0-I*(p*r))+l*B4(-eta,-Rn*pp)*cexp(I*(p*r))))*0.5);
//		G=1-kl*kl/(pp*pp)*(1-Gi*pow(cabs(B4(eta,R*pp)+l*B4(eta,Rn*pp)),2)*0.5);
#endif
#else
//		G=Gi*pow(cabs(B2(eta,R*pp)+l*B2(eta,Rn*pp)),2)*0.5;
		G=Gi*pow(cabs(B2(eta,R*pp)*cexp(I*(p*r))
				+l*B2(eta,Rn*pp)*cexp(0-I*(p*r))),2)*0.5;
#endif
//		if(G>1.01&&pp>0.06)cout<<"G2 "<<G<<" "<<eta<<" "<<R<<" "<<pp<<endl;
	}
	return G;
}

//psi_ij=G(eta)^(1/2)(1+cos(2*q dot r));
double Psi_ij(TLorentzVector p1,TLorentzVector p2,TLorentzVector r1,TLorentzVector r2)
{
	TLorentzVector q;
	TLorentzVector K;
	K=(p1+p2)*0.5;
	q=(p1-p2)*0.5-trace((p1-p2)*0.5,K);
	q.Boost(-K.Vect()*(1/K.E()));
	TLorentzVector r;
	r=((r1-r2));
	r.Boost(-K.Vect()*(1/K.E()));
	r*=fm_gev;
	//if(fabs(q.P()*r.Rho())>5.0)return 1.;
	if(fabs(r.Vect().Mag())>80)return 1.;
	double m1=p1.Mag();
	double m2=p2.Mag();
#ifdef Coulomb
#ifdef Sinyukov
#ifdef simple
	double Psi=(G(q.Vect(),r.Vect(),m1,m2)+l*G(q.Vect(),-r.Vect(),m1,m2))*
		(1.+l*cos(2.*q.Vect()*r.Vect()))*0.5;
#else
//	double Psi=G(q.Vect(),r.Vect(),m1,m2)*(1.+l*cos(2.*q.Vect()*r.Vect()));
	double Psi=G(q.Vect(),r.Vect(),m1,m2);
#endif
#else
	double Psi=(G(q.Vect(),r.Vect(),m1,m2));
#endif
#else
	double Psi=(1.+l*cos(2.*q.Vect()*r.Vect()));
#endif
//	cout<<" "<<G(q.Vect(),r.Vect(),m1,m2)
//		<<" "<<q.Vect().Mag()<<endl;
//	if(Psi>1.01&&q.Vect().Mag()>0.06)cout<<"Psi "<<Psi<<" "<<q.Vect()*r.Vect()<<" "<<r.Vect().Mag()<<" "<<q.Vect().Mag()<<endl;
	return Psi;
}

TLorentzVector trace(TLorentzVector q,TLorentzVector K)
{
	TLorentzVector trace=((q*K)/K.Mag2())*K;
	return trace;
}

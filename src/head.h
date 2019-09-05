#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <cstdlib>
#include <string.h>
#include <math.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TBenchmark.h>
#include <TInterpreter.h>
#include <stdio.h>
#include <stdlib.h>
#include <TComplex.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH3.h>
#include <TF3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMinuit.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include "TVector3.h"
#include "TVector2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TSystem.h"
#include "hadd.h"

using namespace std;

#define byrun
#define logpt
#define LCMS
int lambda_fixed=0;	// 0 for not fixed; 1 for fixed.
double point[2];
int runnum=0;
int par5=0;

#include "fit.h"
#include "outopt.h"
/*A & B*/
double q[q_step+1];	//q=p1-p2 or q=(p1-p2)/2;
/*C*/
double q_pt[pt_step+1];
/*D*/
double q_phi[phi_step+1];

/*C*/
TGraphErrors *qinv_pt_R=new TGraphErrors(pt_step);
TGraphErrors *side_pt_R=new TGraphErrors(pt_step);
TGraphErrors *beam_pt_R=new TGraphErrors(pt_step);
TGraphErrors *out_pt_R=new TGraphErrors(pt_step);
/*D*/
TGraphErrors *qinv_phi_R=new TGraphErrors(phi_step);
TGraphErrors *side_phi_R=new TGraphErrors(phi_step);
TGraphErrors *beam_phi_R=new TGraphErrors(phi_step);
TGraphErrors *out_phi_R=new TGraphErrors(phi_step);

void Getq(TLorentzVector p1, TLorentzVector p2, TVector3 *q_osl, double *qinv, double *pt,double *phi)
{
	TLorentzVector q;
	TLorentzVector K;
#ifdef qinv_half
	q=(p1-p2)*0.5;
#else
	q=(p1-p2);
#endif
	K=(p1+p2)*0.5;
#ifdef LCMS
	TVector3 Kz;
	Kz.SetXYZ(0,0,K.Pz()/K.E());
	q.Boost(-Kz);
	K.Boost(-Kz);
#endif
	*pt=K.Perp();
	*phi=K.Phi();
	TVector2 q2;
	q2.Set(q.X(),q.Y());
	q2=q2.Rotate(-K.Phi());
	q_osl->SetXYZ(q2.X(),q2.Y(),q.Z());
	*qinv=fabs(q.Mag());
}

void Getr(TLorentzVector p1, TLorentzVector p2, TLorentzVector r1, TLorentzVector r2, double *rs, double *ro, double *rl)
{
	TLorentzVector q;
	TLorentzVector K;
	TLorentzVector r;
	double beta_prep;
	double beta_l;
	double Kt;
#ifdef qinv_half
	q=(p1-p2)*0.5;
#else
	q=(p1-p2);
#endif
	K=(p1+p2)*0.5;
	r=(r1-r2);
#ifdef LCMS
	TVector3 Kz;
	Kz.SetXYZ(0,0,K.Pz()/K.E());
	q.Boost(-Kz);
	K.Boost(-Kz);
	r.Boost(-Kz);
#endif
	Kt=sqrt(pow(K.Px(),2)+pow(K.Py(),2));
	beta_prep=Kt/K.E();
	beta_l=K.Pz()/K.E();

	TVector2 q2;
	TVector2 rr2;
	q2.Set(q.X(),q.Y());
	rr2.Set(r.X(),r.Y());
	q2=q2.Rotate(-K.Phi());
	rr2=rr2.Rotate(-K.Phi());
	*rs=rr2.Y();
	*ro=rr2.X()-beta_prep*r.T();
	*rl=r.Z()-beta_l*r.T();
}

void Getbin(double *q, double *q_pt, double *q_phi)
{
	q[0]=q_min;
	for(int k=1;k<=q_step;k++)
		q[k]=dq*(k)+q_min;
	q_pt[0]=pt_min;
	for(int i=1;i<=pt_step;i++)
#ifdef logpt
		q_pt[i]=floor((pow((double)pt_step,(i-1)/((double)pt_step-1))*dpt+pt_min)*1000)/1000;
#else
		q_pt[i]=dpt*(i)+pt_min;
#endif
	for(int j=0;j<=phi_step;j++)
		q_phi[j]=dphi*j+phi_min;
}

void Fileout(char FileOutput[400][100], double *q_pt, double *q_phi)
{
	int point;
	/*A*/
	sprintf(FileOutput[0],"qinv");
	/*B*/
	sprintf(FileOutput[1],"q_beam");
	sprintf(FileOutput[2],"q_sidewards");
	sprintf(FileOutput[3],"q_outwards");
	/*C*/
	for(int j=1;j<4;j++){
		for(int i=0;i<pt_step;i++){
			point=i+4+j*pt_step-pt_step;
			sprintf(FileOutput[point],"%s_pt%04g"
					,FileOutput[j],floor((q_pt[i]+q_pt[i+1])*500));
		}
	}
	for(int i=0;i<pt_step;i++){
		point=i+4+3*phi_step+3*pt_step;
		sprintf(FileOutput[point],"%s_pt%04g"
				,FileOutput[0],floor((q_pt[i]+q_pt[i+1])*500));
	}
	/*D*/
	for(int j=1;j<4;j++){
		for(int i=0;i<phi_step;i++){
			point=i+4+3*pt_step+j*phi_step-phi_step;
			sprintf(FileOutput[point],"%s_phi%04g"
					,FileOutput[j],floor((q_phi[i]+q_phi[i+1])*500));
		}
	}
	for(int i=0;i<phi_step;i++){
		point=i+4+3*phi_step+4*pt_step;
		sprintf(FileOutput[point],"%s_phi%04g"
				,FileOutput[0],floor((q_phi[i]+q_phi[i+1])*500));
	}
}

void Fileout(char FileOutput[400][100])
{
	int point;
	double q_pt;
	double q_phi;
	/*A*/
	sprintf(FileOutput[0],"qinv");
	/*B*/
	sprintf(FileOutput[1],"q_beam");
	sprintf(FileOutput[2],"q_sidewards");
	sprintf(FileOutput[3],"q_outwards");
	/*C*/
	for(int j=1;j<4;j++){
		for(int i=0;i<pt_step;i++){
			point=i+4+j*pt_step-pt_step;
			q_pt=floor((i+1)*dpt*1000);
			sprintf(FileOutput[point],"%s_pt%04g"
					,FileOutput[j],q_pt);
		}
	}
	for(int i=0;i<pt_step;i++){
		point=i+4+3*pt_step+3*phi_step;
		q_pt=floor((i+1)*dpt*1000);
		sprintf(FileOutput[point],"%s_pt%04g"
				,FileOutput[0],q_pt);
	}
	/*D*/
	for(int j=1;j<4;j++){
		for(int i=0;i<phi_step;i++){
			point=i+4+3*pt_step+j*phi_step-phi_step;
			q_phi=(i+1)*dphi;
			sprintf(FileOutput[point],"%s_phi%04g"
					,FileOutput[j],q_phi);
		}
	}
	for(int i=0;i<phi_step;i++){
		point=i+4+4*pt_step+3*phi_step;
		q_phi=(i+1)*dphi;
		sprintf(FileOutput[point],"%s_phi%04g"
				,FileOutput[0],q_phi);
	}
}

void print(double *q,double *corr,double* err,int *Npair,char *FileOutput)
{
	TCanvas *c1=new TCanvas();
	char dattitle[100];
	char giftitle[100];
	double ex[q_step];
	TGraphErrors *plot=new TGraphErrors(q_step,q,corr,ex,err);
	bool t=1;	//1 for pt;
			//0 for phi;
	t=(strstr(FileOutput,"phi")==NULL);
	MakePlot(plot,t);
//	fit_plot(plot);
	plot->SetTitle(FileOutput);
	plot->Draw(drawopt);
	sprintf(giftitle,"../result_gif/%s.gif",FileOutput);
	c1->Print(giftitle);
	sprintf(dattitle,"../result_dat/%s.dat",FileOutput);
	ofstream result;
	result.open(dattitle);
	for(int i=-1;i<q_step;i++){
		if(i<0){
			result<<line1<<endl;
			result<<line2<<endl;
			continue;
		}
		result<<q[i]<<"	"<<corr[i]<<" "<<err[i]<<" "<<Npair[i]<<endl;
	}
	result.close();
}

void write(TProfile *plot,const char *FileOutput)
{
	plot->SetTitle(FileOutput);
#ifdef byrun
	char rootname[100];
	sprintf(rootname,"../result_root/result_%d.root",runnum);
	TFile result_root(rootname,"UPDATE");
#else
	TFile result_root("../result_root/result.root","UPDATE");
#endif
	plot->Write();
	result_root.Close();
}

void write(TH1D *plot,const char *FileOutput)
{
	plot->SetTitle(FileOutput);
#ifdef byrun
	char rootname[100];
	sprintf(rootname,"../result_root/result_%d.root",runnum);
	TFile result_root(rootname,"UPDATE");
#else
	TFile result_root("../result_root/result.root","UPDATE");
#endif
	plot->Write();
	result_root.Close();
}

void write(TH1D *plot,const char *FileOutput,int Nevent)
{
	TF1 *f1=new TF1("f1","x",dpt_min,dpt_max);
	plot->Divide(f1,2*M_PI*Nevent);
	plot->SetTitle(FileOutput);
#ifdef byrun
	char rootname[100];
	sprintf(rootname,"../result_root/result_%d.root",runnum);
	TFile result_root(rootname,"UPDATE");
#else
	TFile result_root("../result_root/result.root","UPDATE");
#endif
	plot->Write();
	result_root.Close();
}

void write(TH2D *plot,const char *FileOutput)
{
	plot->SetTitle(FileOutput);
#ifdef byrun
	char rootname[100];
	sprintf(rootname,"../result_root/result_%d.root",runnum);
	TFile result_root(rootname,"UPDATE");
#else
	TFile result_root("../result_root/result.root","UPDATE");
#endif
	plot->Write();
	result_root.Close();
}

void T_Rwrite(TH1 *plot,const char *FileOutput)
{
#ifdef byrun
	plot->SetTitle(FileOutput);
	char rootname[100];
	sprintf(rootname,"../result_root/result_%d.root",runnum);
	TFile result_root(rootname,"UPDATE");
#else
	TFile result_root("../result_root/result.root","UPDATE");
#endif
	plot->Write();
	result_root.Close();
}

void T_Rwrite(TH2D *plot,const char *FileOutput)
{
#ifdef byrun
	plot->SetTitle(FileOutput);
	char rootname[100];
	sprintf(rootname,"../result_root/result_%d.root",runnum);
	TFile result_root(rootname,"UPDATE");
#else
	TFile result_root("../result_root/result.root","UPDATE");
#endif
	plot->Write();
	result_root.Close();
}

void T_Rwrite(TProfile *plot,const char *FileOutput)
{
#ifdef byrun
	plot->SetTitle(FileOutput);
	char rootname[100];
	sprintf(rootname,"../result_root/result_%d.root",runnum);
	TFile result_root(rootname,"UPDATE");
#else
	TFile result_root("../result_root/result.root","UPDATE");
#endif
	plot->Write();
	result_root.Close();
}

void print(TProfile *plot,char *FileOutput)
{
	TCanvas *c1=new TCanvas();
	char dattitle[100];
	char giftitle[100];
	TF1 *line=new TF1("Line","1",q_min,q_max);
	MakeLine(line);
	TF1 *fit_plot = new TF1("fit",fit,q_min,q_max,npar);
	MakeFit(fit_plot);
	plot->Fit(fit_plot,fittype);

	TGraphErrors *tet=new TGraphErrors(plot->GetNbinsX());
	for(int i=1;i<=plot->GetNbinsX();i++)
	{
		tet->SetPoint(i,plot->GetBinCenter(i),plot->GetBinContent(i));
		tet->SetPointError(i,0,plot->GetBinError(i));
	}
	tet->RemovePoint(0);
	tet->Fit(fit_plot,fittype);

	MakePlot(plot);
	lambda[par5]=fit_plot->GetParameter(0);
	plot->Draw(drawopt2);
	fit_plot->Draw("same");
	line->Draw("same");
	sprintf(giftitle,"../result_gif/%s.gif",FileOutput);
	c1->Print(giftitle);
	sprintf(dattitle,"../result_dat/%s.dat",FileOutput);
	ofstream result;
	result.open(dattitle);
	for(int i=0;i<q_step;i++){
		if(i<=0){
			result<<line1<<endl;
			result<<"Fit result: lambda="<<fit_plot->GetParameter(0)<<" err_lambda="<<fit_plot->GetParError(0)<<" R="<<fit_plot->GetParameter(1)<<" err_R="<<fit_plot->GetParError(1)<<endl;
			result<<line2<<endl;
			continue;
		}
		result<<plot->GetBinCenter(i)<<"	"<<plot->GetBinContent(i)<<" "<<plot->GetBinError(i)<<" "<<plot->GetBinEntries(i)<<endl;
	}
	result.close();
	ifstream ifexit;
	ifexit.open("../result_dat/fit.dat");
	if(! ifexit){
		ofstream fit_result_first;
		fit_result_first.open("../result_dat/fit.dat");
		fit_result_first<<"type		lambda	+/-	R/fm	+/-"<<endl;
		fit_result_first<<FileOutput<<"	"<<fit_plot->GetParameter(0)<<"	"<<fit_plot->GetParError(0)<<"	"<<fit_plot->GetParameter(1)<<"	"<<fit_plot->GetParError(1)<<endl;
		fit_result_first.close();
	}
	else{
		ofstream fit_result;
		fit_result.open("../result_dat/fit.dat",ios::app);
		fit_result<<FileOutput<<"	"<<fit_plot->GetParameter(0)<<"	"<<fit_plot->GetParError(0)<<"	"<<fit_plot->GetParameter(1)<<"	"<<fit_plot->GetParError(1)<<endl;
		fit_result.close();
	}
	ifexit.close();
	point[0]=fit_plot->GetParameter(1);
	point[1]=fit_plot->GetParError(1);
}

void print(TProfile2D* plot2D,char *FileOutput,int j)
{
	TProfile *plot=plot2D->ProfileX(FileOutput,(int)j+1,(int)j+1);
#ifdef byrun
	TFile *result_root=new TFile("../result_root/result-all.root","UPDATE");
#else
	TFile *result_root=new TFile("../result_root/result.root","UPDATE");
#endif
	plot->SetTitle(FileOutput);
	plot->Write();
	print(plot,FileOutput);
	if(strncmp(FileOutput,"q_sidewards_pt",14)==0){
		side_pt_R->SetPoint(j,(q_pt[j]+q_pt[j+1])/2.,point[0]);
		side_pt_R->SetPointError(j,0,point[1]);
	}
	else if(strncmp(FileOutput,"q_beam_pt",9)==0){
		beam_pt_R->SetPoint(j,(q_pt[j]+q_pt[j+1])/2.,point[0]);
		beam_pt_R->SetPointError(j,0,point[1]);
	}
	else if(strncmp(FileOutput,"q_outwards_pt",13)==0){
		out_pt_R->SetPoint(j,(q_pt[j]+q_pt[j+1])/2.,point[0]);
		out_pt_R->SetPointError(j,0,point[1]);
	}
	else if(strncmp(FileOutput,"q_sidewards_phi",14)==0){
		side_phi_R->SetPoint(j,(q_phi[j]+q_phi[j+1])/2.,point[0]);
		side_phi_R->SetPointError(j,0,point[1]);
	}
	else if(strncmp(FileOutput,"q_beam_phi",9)==0){
		beam_phi_R->SetPoint(j,(q_phi[j]+q_phi[j+1])/2.,point[0]);
		beam_phi_R->SetPointError(j,0,point[1]);
	}
	else if(strncmp(FileOutput,"q_outwards_phi",13)==0){
		out_phi_R->SetPoint(j,(q_phi[j]+q_phi[j+1])/2.,point[0]);
		out_phi_R->SetPointError(j,0,point[1]);
	}
	result_root->Close();
}

void print(TGraphErrors *plot,const char *FileOutput)
{
	TCanvas *c1=new TCanvas();
	char giftitle[100];
	//plot->RemovePoint(0);
	bool t=1;	//1 for pt;
			//0 for phi;
	t=(strstr(FileOutput,"phi")==NULL);
	MakePlot(plot,t);
	plot->SetName(FileOutput);
	plot->SetTitle(FileOutput);
	plot->Draw(drawopt);
	sprintf(giftitle,"../result_gif/%s.gif",FileOutput);
#ifdef byrun
	TFile *result_root=new TFile("../result_root/result-all.root","UPDATE");
#else
	TFile *result_root=new TFile("../result_root/result.root","UPDATE");
#endif
	plot->Write();
	c1->Print(giftitle);
	result_root->Close();
}

void print(TH1D *plot,const char *FileOutput)
{
	TCanvas *c1=new TCanvas();
	char giftitle[100];
	MakePlot(plot);
	plot->SetTitle(FileOutput);
	if(strncmp(FileOutput,"dpt",3)==0)gPad->SetLogy();
	plot->Draw(drawopt2);
	sprintf(giftitle,"../result_gif/%s.gif",FileOutput);
	c1->Print(giftitle);
}
	
void print(TH2D *plot2D,const char*FileOutput)
{
	TGraphErrors *graph1=new TGraphErrors();
	char* graphtitle1=new char[100];
	sprintf(graphtitle1,"%s-Kurtosis",FileOutput);
	TF1 *f1=new TF1("f1","1",q_min,q_max);
	TF1 *f_1=new TF1("f_1","pow(x,1)",q_min,q_max);
	TF1 *f_2=new TF1("f_2","pow(x,2)",q_min,q_max);
	TF1 *f_3=new TF1("f_3","pow(x,3)",q_min,q_max);
	TF1 *f_4=new TF1("f_4","pow(x,4)",q_min,q_max);
	double result_4;
	double result_3;
	double result_2;
	double result_1;
	double result0;
	double result;
	double error;
	double error0;
	double error_1;
	double error_2;
	double error_3;
	double error_4;
	if(strstr(FileOutput,"R")==NULL)
	{
		for(int i=1;i<=floor((plot2D->GetNbinsY()-cut)/pack);i++)
		{
			result=0;
			error=0;
			TH1D *plot=plot2D->ProfileX("plot",(i-1)*pack+1,i*pack);

			TH1D *plot_1=plot2D->ProfileX("plot_1",(i-1)*pack+1,i*pack);
			TH1D *plot_2=plot2D->ProfileX("plot_2",(i-1)*pack+1,i*pack);
			TH1D *plot_3=plot2D->ProfileX("plot_3",(i-1)*pack+1,i*pack);
			TH1D *plot_4=plot2D->ProfileX("plot_4",(i-1)*pack+1,i*pack);

			plot_1->Multiply(f_1,1);
			plot_2->Multiply(f_2,1);
			plot_3->Multiply(f_3,1);
			plot_4->Multiply(f_4,1);
		
			result0=plot->IntegralAndError(1,plot->GetNbinsX(),error0,"width")-f1->Integral(q_min,q_max);
			result_1=(plot_1->IntegralAndError(1,plot_1->GetNbinsX(),error_1,"width")-f_1->Integral(q_min,q_max))/result0;
			error_1=sqrt(pow(error_1/result0,2)+pow(result_1*error0/(result0*result0),2));
			result_2=(plot_2->IntegralAndError(1,plot_2->GetNbinsX(),error_2,"width")-f_2->Integral(q_min,q_max))/result0;
			error_2=sqrt(pow(error_2/result0,2)+pow(result_2*error0/(result0*result0),2));
			result_3=(plot_3->IntegralAndError(1,plot_3->GetNbinsX(),error_3,"width")-f_3->Integral(q_min,q_max))/result0;
			error_3=sqrt(pow(error_3/result0,2)+pow(result_3*error0/(result0*result0),2));
			result_4=(plot_4->IntegralAndError(1,plot_4->GetNbinsX(),error_4,"width")-f_4->Integral(q_min,q_max))/result0;
			error_4=sqrt(pow(error_4/result0,2)+pow(result_4*error0/(result0*result0),2));

			result=(result_4-4*result_3*result_1-3*pow(result_2,2)+12*result_2*pow(result_1,2)-6*pow(result_1,4))/pow(result_2-pow(result_1,2),2)/3;
			error_4=pow(error_4/pow(result_2,2)/3,2);
			error_3=pow(4*error_3*result_1/pow(result_2,2)/3,2);
			error_2=(-6*result_2*error_2+12*error_2*pow(result_1,2))/pow(result_2-pow(result_1,2),2)/3;
			error_2-=result*2*error_2/(result_2-pow(result_1,2));
			error_2*=error_2;
			error_1=(24*result_2*result_1*error_1-24*pow(result_1,3)*error_1)/pow(result_2-pow(result_1,2),2)/3;
			error_1-=result*2*(-2)*result_1*error_1/(result_2-pow(result_1,2));
			error_1*=error_1;
			error=sqrt(error_1+error_2+error_3+error_4);
		
			double last=plot2D->GetYaxis()->GetBinCenter((i-1)*pack+1);
			double next=plot2D->GetYaxis()->GetBinCenter(i*pack);
			graph1->SetPoint(i,(last+next)/2,result);
			//graph->SetPointError(i,0,0);
			graph1->SetPointError(i,0,error);
			delete plot;
			delete plot_1;
			delete plot_2;
			delete plot_3;
			delete plot_4;
		}
		graph1->RemovePoint(0);
		print(graph1,graphtitle1);
 
		TGraphErrors *graph2=new TGraphErrors();
		char* graphtitle2=new char[100];
		sprintf(graphtitle2,"%s-Skewness",FileOutput);
		for(int i=1;i<=floor((plot2D->GetNbinsY()-cut)/pack);i++)
		{
			result=0;
			error=0;
			TH1D *plot=plot2D->ProfileX("plot",(i-1)*pack+1,i*pack);

			TH1D *plot_1=plot2D->ProfileX("plot_1",(i-1)*pack+1,i*pack);
			TH1D *plot_2=plot2D->ProfileX("plot_2",(i-1)*pack+1,i*pack);
			TH1D *plot_3=plot2D->ProfileX("plot_3",(i-1)*pack+1,i*pack);
	
			plot_1->Multiply(f_1,1);
			plot_2->Multiply(f_2,1);
			plot_3->Multiply(f_3,1);
			
			result0=plot->IntegralAndError(1,plot->GetNbinsX(),error0,"width")-f1->Integral(q_min,q_max);
			result_1=(plot_1->IntegralAndError(1,plot_1->GetNbinsX(),error_1,"width")-f_1->Integral(q_min,q_max))/result0;
			error_1=sqrt(pow(error_1/result0,2)+pow(result_1*error0/(result0*result0),2));
			result_2=(plot_2->IntegralAndError(1,plot_2->GetNbinsX(),error_2,"width")-f_2->Integral(q_min,q_max))/result0;
			error_2=sqrt(pow(error_2/result0,2)+pow(result_2*error0/(result0*result0),2));
			result_3=(plot_3->IntegralAndError(1,plot_3->GetNbinsX(),error_3,"width")-f_3->Integral(q_min,q_max))/result0;
			error_3=sqrt(pow(error_3/result0,2)+pow(result_3*error0/(result0*result0),2));
	
			result=(result_3-3*result_2*result_1+2*pow(result_1,3))/pow(result_2-pow(result_1,2),3/2);
			error_3=pow(error_3/pow(result_2,3/2),2);
			error_2=(-3*result_1*error_2)/pow(result_2-pow(result_1,2),3/2);
			error_2-=result*3/2*error_2/(result_2-pow(result_1,2));
			error_2*=error_2;
			error_1=(-3*result_2*error_1+6*pow(result_1,2)*error_1)/pow(result_2-pow(result_1,2),3/2);
			error_1-=result*3/2*(-2)*result_1*error_1/(result_2-pow(result_1,2));
			error_1*=error_1;
			error=sqrt(error_1+error_2+error_3);
		
			double last=plot2D->GetYaxis()->GetBinCenter((i-1)*pack+1);
			double next=plot2D->GetYaxis()->GetBinCenter(i*pack);
			graph2->SetPoint(i,(last+next)/2,result);
			//graph->SetPointError(i,0,0);
			graph2->SetPointError(i,0,error);
			delete plot;
			delete plot_1;
			delete plot_2;
			delete plot_3;
		}
		graph2->RemovePoint(0);
		print(graph2,graphtitle2);

		TGraphErrors *graph3=new TGraphErrors();
		char* graphtitle3=new char[100];
		sprintf(graphtitle3,"%s-R2",FileOutput);
		for(int i=1;i<=floor((plot2D->GetNbinsY()-cut)/pack);i++)
		{
			TH1D *plot=plot2D->ProfileX("plot",(i-1)*pack+1,i*pack);

			TH1D *plot_2=plot2D->ProfileX("plot_2",(i-1)*pack+1,i*pack);
	
			plot_2->Multiply(f_2,1);
			
			result0=plot->IntegralAndError(1,plot->GetNbinsX(),error0,"width")-f1->Integral(q_min,q_max);
			result=(plot_2->IntegralAndError(1,plot_2->GetNbinsX(),error_2,"width")-f_2->Integral(q_min,q_max))/result0;
			result=1/result;
			error=sqrt(pow(error_2/result0,2)+pow(result*error0/(result0*result0),2));
			error=abs(error/(result=result));

			double last=plot2D->GetYaxis()->GetBinCenter((i-1)*pack+1);
			double next=plot2D->GetYaxis()->GetBinCenter(i*pack);
			graph3->SetPoint(i,(last+next)/2,result);
			graph3->SetPointError(i,0,error);
			delete plot;
			delete plot_2;
		}
		graph3->RemovePoint(0);
		print(graph3,graphtitle3);
	}
	else
	{
		TGraphErrors *graph3=new TGraphErrors();
		char* graphtitle3=new char[100];
		sprintf(graphtitle3,"%s",FileOutput);
		for(int i=1;i<=floor((plot2D->GetNbinsY()-cut)/pack);i++)
		{
			TH1D *plot=plot2D->ProfileX("plot",(i-1)*pack+1,i*pack);
			result=pow(plot->GetStdDev(),2);
	
			double last=plot2D->GetYaxis()->GetBinCenter((i-1)*pack+1);
			double next=plot2D->GetYaxis()->GetBinCenter(i*pack);
			graph3->SetPoint(i,(last+next)/2,result);
			graph3->SetPointError(i,0,0);
			delete plot;
		}
		graph3->RemovePoint(0);
		print(graph3,graphtitle3);
	}
}

void T_Rprint(TH1 *plot,const char *FileOutput)
{
	TCanvas *c1=new TCanvas();
	char dattitle[100];
	char giftitle[100];
	MakeT_RPlot(plot);
	plot->SetTitle(FileOutput);
	plot->Draw(drawopt2);
	sprintf(giftitle,"../result_gif/%s.gif",FileOutput);
	c1->Print(giftitle);
	sprintf(dattitle,"../result_dat/%s.dat",FileOutput);
	ofstream result;
	result.open(dattitle);
	for(int i=0;i<plot->GetNbinsX();i++){
		if(i<=0){
			result<<line1<<endl;
			result<<"T_R	N "<<endl;
			continue;
		}
		result<<plot->GetBinCenter(i)<<"	"<<plot->GetBinContent(i)<<endl;
	}
	result.close();
}

void T_Rprint(TProfile *plot,const char *FileOutput)
{
	TCanvas *c1=new TCanvas();
	char dattitle[100];
	char giftitle[100];
	MakeT_RPlot(plot);
	plot->SetTitle(FileOutput);
	plot->Draw(drawopt2);
	sprintf(giftitle,"../result_gif/%s.gif",FileOutput);
	c1->Print(giftitle);
	sprintf(dattitle,"../result_dat/%s.dat",FileOutput);
	ofstream result;
	result.open(dattitle);
	for(int i=0;i<plot->GetNbinsX();i++){
		if(i<=0){
			result<<line1<<endl;
			result<<"T_R	N "<<endl;
			continue;
		}
		result<<plot->GetBinCenter(i)<<"	"<<plot->GetBinContent(i)<<endl;
	}
	result.close();
}

void T_Rprint(TH2D *plot2D,const char *FileOutput)
{
	TH1D *plot=new TH1D();
	for(int i=0;i<q_step;i++){
		char title[100];
		double q=(q_max-q_min)/q_step*(i+0.5);
		sprintf(title,"%s_%04g",FileOutput,q*1000);
		plot=plot2D->ProjectionY(title,i+1,i+1);
		T_Rprint(plot,title);
#ifdef byrun
		TFile *result_root=new TFile("../result_root/result-all.root","UPDATE");
#else
		TFile *result_root=new TFile("../result_root/result.root","UPDATE");
#endif
		plot->SetTitle(title);
		result_root->Close();
	}
}


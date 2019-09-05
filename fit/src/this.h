using namespace std;

#define pt
//#define phi

const char *fittype="ON REMW";
TGraphErrors *q_inv_R[3];
TGraphErrors *q_beam_R[3];
TGraphErrors *q_side_R[3];
TGraphErrors *q_out_R[3];
TMultiGraph *mg_inv=new TMultiGraph();
TMultiGraph *mg_beam=new TMultiGraph();
TMultiGraph *mg_side=new TMultiGraph();
TMultiGraph *mg_out=new TMultiGraph();
const char *nameinv="q_{inv}-R";
const char *namebeam="q_{beam}-R";
const char *nameside="q_{side}-R";
const char *nameout="q_{out}-R";
const char *outputroot="result-fit.root";
int cut=0;
int pack=1;
double lambda[4]={0.6,0.6,0.6,0.6};
double R0=4;
double A=1;
double n=1;

void XY()
{
	mg_inv->SetName(nameinv);
	mg_inv->SetTitle(nameinv);
	/*
#ifdef pt
	mg_inv->GetXaxis()->SetTitle("pt/GeV");
#endif
#ifdef phi
	mg_inv->GetXaxis()->SetTitle("phi/^{#circ}");
#endif
	mg_inv->GetXaxis()->SetRangeUser(xmin,xmax);
	mg_inv->GetYaxis()->SetTitle("R/fm");*/

	mg_beam->SetName(namebeam);
	mg_beam->SetTitle(namebeam);
	/*
#ifdef pt
	mg_beam->GetXaxis()->SetTitle("pt/GeV");
#endif
#ifdef phi
	mg_beam->GetXaxis()->SetTitle("phi/^{#circ}");
#endif
	mg_beam->GetXaxis()->SetRangeUser(xmin,xmax);
	mg_beam->GetYaxis()->SetTitle("R/fm");*/

	mg_side->SetName(nameside);
	mg_side->SetTitle(nameside);
	/*
#ifdef pt
	mg_side->GetXaxis()->SetTitle("pt/GeV");
#endif
#ifdef phi
	mg_side->GetXaxis()->SetTitle("phi/^{#circ}");
#endif
	mg_side->GetXaxis()->SetRangeUser(xmin,xmax);
	mg_side->GetYaxis()->SetTitle("R/fm");*/

	mg_out->SetName(nameout);
	mg_out->SetTitle(nameout);
	/*
#ifdef pt
	mg_out->GetXaxis()->SetTitle("pt/GeV");
#endif
#ifdef phi
	mg_out->GetXaxis()->SetTitle("phi/^{#circ}");
#endif
	mg_out->GetXaxis()->SetRangeUser(xmin,xmax);
	mg_out->GetYaxis()->SetTitle("R/fm");*/
}

void fit(TH1D* plot,TF1 *f)
{
	if(iffix==0){
		f->SetParameter(0,lambda[par5]);
		f->SetParLimits(0,0.1,1.1);}
	if(iffix==1)
		f->FixParameter(0,lambda[par5]);
	f->SetParameter(1,6);
	f->SetParLimits(1,0.1,15);
	f->FixParameter(2,R0);
	f->FixParameter(3,A);
	f->FixParameter(4,n);
}

void fit(TH1D** plot,const char* FileOutput)
{
	TF1 *f[3];
	for(int i=0;i<3;i++)
	{
		double q_min=plot[i]->GetXaxis()->GetBinCenter(1);
		int q_step=plot[i]->GetNbinsX();
		double q_max=plot[i]->GetXaxis()->GetBinCenter(q_step);
		char *fname=new char[100];
		sprintf(fname,"fit_%d",i);
		f[i]=new TF1(fname,fit1,q_min,q_max,5);
		par6=i;
		fit(plot[i],f[i]);
		plot[i]->Fit(f[i],fittype);
		if(iffix==0)
			lambda[par5]=f[i]->GetParameter(0);
		make(plot[i],i);
		make(f[i],i);
		delete fname;
	}
	plot[0]->Draw("E0 X0");
	f[0]->Draw("same");
	plot[1]->Draw("E0 X0 same");
	f[1]->Draw("same");
	plot[2]->Draw("E0 X0 same");
	f[2]->Draw("same");
	output=new TFile(outputroot,"UPDATE");
	legd(plot);
	char *name=new char[100];
	sprintf(name,"%s.gif",FileOutput);
	c1->Print(name);
	c1->SetName(FileOutput);
	c1->Write();
	c1->Clear();
	output->Close();
}

void fit(TH1D** plot,const char* FileOutput,TF1** f)
{
	for(int i=0;i<3;i++){
		double q_min=plot[i]->GetXaxis()->GetBinCenter(1);
		int q_step=plot[i]->GetNbinsX();
		double q_max=plot[i]->GetXaxis()->GetBinCenter(q_step);
		char *fname=new char[100];
		sprintf(fname,"fit_%d",i);
		f[i]=new TF1(fname,fit1,q_min,q_max,5);
		par6=i;
		fit(plot[i],f[i]);
		make(plot[i],i);
		make(f[i],i);
		plot[i]->Fit(fname,fittype);
		if(iffix==0)
			lambda[par5]=f[i]->GetParameter(0);
		delete fname;

		TGraphErrors *tet=new TGraphErrors(plot[i]->GetNbinsX());
		for(int j=1;j<=plot[i]->GetNbinsX();j++)
		{
			tet->SetPoint(j,plot[i]->GetBinCenter(j),plot[i]->GetBinContent(j));
			tet->SetPointError(j,0,plot[i]->GetBinError(j));
		}
		tet->RemovePoint(0);
		tet->Fit(f[i],fittype);
		delete tet;
	}
	plot[0]->Draw("E0 X0");
	f[0]->Draw("same");
	plot[1]->Draw("E0 X0 same");
	f[1]->Draw("same");
	plot[2]->Draw("E0 X0 same");
	f[2]->Draw("same");
	output=new TFile(outputroot,"UPDATE");
	legd(plot);
	c1->SetName(FileOutput);
	c1->Write();
	c1->Clear();
	output->Close();
}

void fit(TH2D** plot2D,TGraphErrors** graph,const char* FileOutput)
{
	TH1D* plot[3];
	TF1* f[3];
	for(int i=1;i<=floor((plot2D[1]->GetNbinsY()-cut)/pack);i++)
	{
		for(int j=0;j<3;j++)
		{
			char *kk=new char[100];
			sprintf(kk,"plot%d",j);
			plot[j]=plot2D[j]->ProfileX(kk,i,i);
			delete kk;
		}
		char* name=new char[100];
#ifdef pt
		sprintf(name,"%s%04g",FileOutput,floor(plot2D[1]->GetYaxis()->GetBinCenter(i)*1000));
#endif
#ifdef phi
		sprintf(name,"%s%04g",FileOutput,floor(plot2D[1]->GetYaxis()->GetBinCenter(i)));
#endif
		fit(plot,name,f);
		delete name;

		for(int j=0;j<3;j++)
		{
			graph[j]->SetPoint(i,plot2D[j]->GetYaxis()->GetBinCenter(i),f[j]->GetParameter(1));
			graph[j]->SetPointError(i,0,f[j]->GetParError(1));
			//delete f[j];
			//delete plot[j];
		}
	}
}

void add()
{
	for(int i=0;i<3;i++)
	{
		make(q_inv_R[i],i);
		make(q_beam_R[i],i);
		make(q_side_R[i],i);
		make(q_out_R[i],i);
		mg_inv->Add(q_inv_R[i]);
		mg_beam->Add(q_beam_R[i]);
		mg_side->Add(q_side_R[i]);
		mg_out->Add(q_out_R[i]);
	}
}

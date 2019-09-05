#include <fstream>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <strstream>
#include <math.h>
#include <string.h>
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
#include <TMath.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH3.h>
#include <TF3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMinuit.h>
#include <TGraph.h>
#include <TGraph2D.h>
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

using namespace std;

#define dpt_dy
#define q_qinv
#define wards_R
#define kurtosis
#define skewness
#define pt
#define phi

	const char *tcn[3];
	int color[3]={2,3,4};
	int style[3]={33,21,28};
	int size[3]={2,1,2};
	TFile *t[3];
	TFile *output;

#ifdef dpt_dy
	TH1D *dy[3];
	TH1D *dpt[3];
#endif

#ifdef q_qinv
	TH1D *ql_qinv[3];
	TH1D *qs_qinv[3];
	TH1D *qo_qinv[3];
#endif

#ifdef wards_R
#ifdef pt
	TGraphErrors *beam_pt_R[3];
	TGraphErrors *side_pt_R[3];
	TGraphErrors *out_pt_R[3];
	TMultiGraph *mg_beam_pt_R=new TMultiGraph("mg_beam_pt_R","mg_beam_pt_R");
	TMultiGraph *mg_side_pt_R=new TMultiGraph("mg_side_pt_R","mg_side_pt_R");
	TMultiGraph *mg_out_pt_R=new TMultiGraph("mg_out_pt_R","mg_out_pt_R");
#endif
#ifdef phi
	TGraphErrors *beam_phi_R[3];
	TGraphErrors *side_phi_R[3];
	TGraphErrors *out_phi_R[3];
	TMultiGraph *mg_beam_phi_R=new TMultiGraph("mg_beam_phi_R","mg_beam_phi_R");
	TMultiGraph *mg_side_phi_R=new TMultiGraph("mg_side_phi_R","mg_side_phi_R");
	TMultiGraph *mg_out_phi_R=new TMultiGraph("mg_out_phi_R","mg_out_phi_R");
#endif
#endif

#ifdef kurtosis
#ifdef pt
	TGraphErrors *beam_pt_k[3];
	TGraphErrors *side_pt_k[3];
	TGraphErrors *out_pt_k[3];
	TMultiGraph *mg_beam_pt_k=new TMultiGraph("mg_beam_pt_k","mg_beam_pt_k");
	TMultiGraph *mg_side_pt_k=new TMultiGraph("mg_side_pt_k","mg_side_pt_k");
	TMultiGraph *mg_out_pt_k=new TMultiGraph("mg_out_pt_k","mg_out_pt_k");
#endif
#ifdef phi
	TGraphErrors *beam_phi_k[3];
	TGraphErrors *side_phi_k[3];
	TGraphErrors *out_phi_k[3];
	TMultiGraph *mg_beam_phi_k=new TMultiGraph("mg_beam_phi_k","mg_beam_phi_k");
	TMultiGraph *mg_side_phi_k=new TMultiGraph("mg_side_phi_k","mg_side_phi_k");
	TMultiGraph *mg_out_phi_k=new TMultiGraph("mg_out_phi_k","mg_out_phi_k");
#endif
#endif

#ifdef skewness
#ifdef pt
	TGraphErrors *beam_pt_s[3];
	TGraphErrors *side_pt_s[3];
	TGraphErrors *out_pt_s[3];
	TMultiGraph *mg_beam_pt_s=new TMultiGraph("mg_beam_pt_s","mg_beam_pt_s");
	TMultiGraph *mg_side_pt_s=new TMultiGraph("mg_side_pt_s","mg_side_pt_s");
	TMultiGraph *mg_out_pt_s=new TMultiGraph("mg_out_pt_s","mg_out_pt_s");
#endif
#ifdef phi
	TGraphErrors *beam_phi_s[3];
	TGraphErrors *side_phi_s[3];
	TGraphErrors *out_phi_s[3];
	TMultiGraph *mg_beam_phi_s=new TMultiGraph("mg_beam_phi_s","mg_beam_phi_s");
	TMultiGraph *mg_side_phi_s=new TMultiGraph("mg_side_phi_s","mg_side_phi_s");
	TMultiGraph *mg_out_phi_s=new TMultiGraph("mg_out_phi_s","mg_out_phi_s");
#endif
#endif

	TCanvas *c1=new TCanvas();
	double legdx1=0.65;
	double legdy1=0.7;
	double legdx2=0.9;
	double legdy2=0.9;

void init()
{
	tcn[0]="../Triangle/result_root/result-all.root";
	tcn[1]="../Chain/result_root/result-all.root";
	tcn[2]="../nofluc/result_root/result-all.root";
}

void XY()
{
#ifdef dpt_dy
	dpt[0]->GetXaxis()->SetTitle("p_{T}(GeV/c)");
	dpt[0]->GetYaxis()->SetTitle("1/#left[(2#pip_{T}N)d^{2}N/dp_{T}dy#right] (GeV/c)^{-2}");
	dpt[0]->GetYaxis()->CenterTitle();
	dpt[0]->SetTitle("");

	dy[0]->GetXaxis()->SetTitle("y");
	dy[0]->GetYaxis()->SetTitle("N");
	dy[0]->GetYaxis()->CenterTitle();
	dy[0]->SetTitle("");
#endif

#ifdef q_qinv
	ql_qinv[0]->GetXaxis()->SetTitle("q_{l}/GeV");
	ql_qinv[0]->GetYaxis()->SetTitle("q_{qinv}/GeV");
	ql_qinv[0]->GetYaxis()->CenterTitle();
	ql_qinv[0]->SetTitle("");

	qs_qinv[0]->GetXaxis()->SetTitle("q_{s}/GeV");
	qs_qinv[0]->GetYaxis()->SetTitle("q_{qinv}/GeV");
	qs_qinv[0]->GetYaxis()->CenterTitle();
	qs_qinv[0]->SetTitle("");

	qo_qinv[0]->GetXaxis()->SetTitle("q_{o}/GeV");
	qo_qinv[0]->GetYaxis()->SetTitle("q_{qinv}/GeV");
	qo_qinv[0]->GetYaxis()->CenterTitle();
	qo_qinv[0]->SetTitle("");
#endif

#ifdef wards_R
#ifdef pt
	//mg_beam_pt_R->GetXaxis()->SetTitle("pt/GeV");
	//mg_beam_pt_R->GetYaxis()->SetTitle("R/fm");
	//mg_beam_pt_R->GetYaxis()->CenterTitle();
	mg_beam_pt_R->SetTitle("beam-pt-R");

	//mg_side_pt_R->GetXaxis()->SetTitle("pt/GeV");
	//mg_side_pt_R->GetYaxis()->SetTitle("R/fm");
	//mg_side_pt_R->GetYaxis()->CenterTitle();
	mg_side_pt_R->SetTitle("sidewards-pt-R");

	//mg_out_pt_R->GetXaxis()->SetTitle("pt/GeV");
	//mg_out_pt_R->GetYaxis()->SetTitle("R/fm");
	//mg_out_pt_R->GetYaxis()->CenterTitle();
	mg_out_pt_R->SetTitle("outwards-pt-R");
#endif
#ifdef phi
	//mg_beam_phi_R->GetXaxis()->SetTitle("phi_{2}/^{#circ}");
	//mg_beam_phi_R->GetYaxis()->SetTitle("R/fm");
	//mg_beam_phi_R->GetYaxis()->CenterTitle();
	mg_beam_phi_R->SetTitle("beam-phi-R");

	//mg_side_phi_R->GetXaxis()->SetTitle("phi_{2}/^{#circ}");
	//mg_side_phi_R->GetYaxis()->SetTitle("R/fm");
	//mg_side_phi_R->GetYaxis()->CenterTitle();
	mg_side_phi_R->SetTitle("sidewards-phi-R");

	//mg_out_phi_R->GetXaxis()->SetTitle("phi_{2}/^{#circ}");
	//mg_out_phi_R->GetYaxis()->SetTitle("R/fm");
	//mg_out_phi_R->GetYaxis()->CenterTitle();
	mg_out_phi_R->SetTitle("outwards-phi-R");
#endif
#endif

#ifdef kurtosis
#ifdef pt
	//mg_beam_pt_k->GetXaxis()->SetTitle("pt/GeV");
	//mg_beam_pt_K->GetYaxis()->SetTitle("#kappa_{beam}(pt)");
	//mg_beam_pt_K->GetYaxis()->CenterTitle();
	mg_beam_pt_k->SetTitle("beam-pt-Kurtosis");

	//mg_side_pt_k->GetXaxis()->SetTitle("pt/GeV");
	//mg_side_pt_k->GetYaxis()->SetTitle("#kappa_{side}(pt)");
	//mg_side_pt_k->GetYaxis()->CenterTitle();
	mg_side_pt_k->SetTitle("sidewards-pt-Kurtosis");

	//mg_out_pt_k->GetXaxis()->SetTitle("pt/GeV");
	//mg_out_pt_k->GetYaxis()->SetTitle("#kappa_{out}(pt)");
	//mg_out_pt_k->GetYaxis()->CenterTitle();
	mg_out_pt_k->SetTitle("outwards-pt-Kurtosis");
#endif
#ifdef phi
	//mg_beam_phi_k->GetXaxis()->SetTitle("phi_{2}/^{#circ}");
	//mg_beam_phi_K->GetYaxis()->SetTitle("#kappa_{beam}(phi_{2})");
	//mg_beam_phi_K->GetYaxis()->CenterTitle();
	mg_beam_phi_k->SetTitle("beam-phi-Kurtosis");

	//mg_side_phi_k->GetXaxis()->SetTitle("phi_{2}/^{#circ}");
	//mg_side_phi_k->GetYaxis()->SetTitle("#kappa_{side}(phi_{2})");
	//mg_side_phi_k->GetYaxis()->CenterTitle();
	mg_side_phi_k->SetTitle("sidewards-phi-Kurtosis");

	//mg_out_phi_k->GetXaxis()->SetTitle("phi_{2}/^{#circ}");
	//mg_out_phi_k->GetYaxis()->SetTitle("#kappa_{out}(phi_{2})");
	//mg_out_phi_k->GetYaxis()->CenterTitle();
	mg_out_phi_k->SetTitle("outwards-phi-Kurtosis");
#endif
#endif

#ifdef skewness
#ifdef pt
	//mg_beam_pt_s->GetXaxis()->SetTitle("pt/GeV");
	//mg_beam_pt_s->GetYaxis()->SetTitle("S_{beam}(pt)");
	//mg_beam_pt_s->GetYaxis()->CenterTitle();
	mg_beam_pt_s->SetTitle("beam-pt-Skewness");

	//mg_side_pt_s->GetXaxis()->SetTitle("pt/GeV");
	//mg_side_pt_s->GetYaxis()->SetTitle("S_{side}(pt)");
	//mg_side_pt_s->GetYaxis()->CenterTitle();
	mg_side_pt_s->SetTitle("sidewards-Skenwness");

	//mg_out_pt_s->GetXaxis()->SetTitle("pt/GeV");
	//mg_out_pt_s->GetYaxis()->SetTitle("S_{out}(pt)");
	//mg_out_pt_s->GetYaxis()->CenterTitle();
	mg_out_pt_s->SetTitle("outwards-pt-Skenwness");
#endif
#ifdef phi
	//mg_beam_phi_s->GetXaxis()->SetTitle("phi_{2}/^{#circ}");
	//mg_beam_phi_s->GetYaxis()->SetTitle("S_{beam}(phi_{2})");
	//mg_beam_phi_s->GetYaxis()->CenterTitle();
	mg_beam_phi_s->SetTitle("beam-phi-Skewness");

	//mg_side_phi_s->GetXaxis()->SetTitle("phi_{2}/^{#circ}");
	//mg_side_phi_s->GetYaxis()->SetTitle("S_{side}(phi_{2})");
	//mg_side_phi_s->GetYaxis()->CenterTitle();
	mg_side_phi_s->SetTitle("sidewards-Skenwness");

	//mg_out_phi_s->GetXaxis()->SetTitle("phi_{2}/^{#circ}");
	//mg_out_phi_s->GetYaxis()->SetTitle("S_{out}(phi_{2})");
	//mg_out_phi_s->GetYaxis()->CenterTitle();
	mg_out_phi_s->SetTitle("outwards-phi-Skenwness");
#endif
#endif
}

void make(TH1D *plot,int i)
{
		plot->SetMarkerSize(size[i]);
		plot->SetMarkerStyle(style[i]);
		plot->SetMarkerColor(color[i]);
}

void make2(TH1D *qinv,int i)
{
		qinv->SetMarkerSize(0.5);
		qinv->SetMarkerStyle(style[i]);
		qinv->SetMarkerColor(color[i]);
		qinv->GetXaxis()->SetNdivisions(505);
		qinv->GetXaxis()->SetLabelSize(0.06);
		qinv->GetXaxis()->SetTitleSize(0.06);
		qinv->GetYaxis()->SetNdivisions(503);
		qinv->GetYaxis()->SetLabelSize(0.06);
		qinv->GetYaxis()->SetTitleSize(0.06);
}

void make(TGraph *plot,int i)
{
		plot->SetMarkerSize(size[i]);
		plot->SetMarkerStyle(style[i]);
		plot->SetMarkerColor(color[i]);
}

void legd(TH1D **plot)
{
        TLegend *legd=new TLegend(legdx1,legdy1,legdx2,legdy2);
        legd->AddEntry(plot[0],"triangle","lep");
        legd->AddEntry(plot[1],"chain","lep");
        legd->AddEntry(plot[2],"nofluc","lep");
        legd->Draw();
}

void legd(TGraphErrors **plot)
{
        TLegend *legd=new TLegend(legdx1,legdy1,legdx2,legdy2);
        legd->AddEntry(plot[0],"triangle","lep");
        legd->AddEntry(plot[1],"chain","lep");
        legd->AddEntry(plot[2],"nofluc","lep");
        legd->Draw();
}

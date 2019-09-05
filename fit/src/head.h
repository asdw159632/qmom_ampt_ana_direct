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

	const char *tcn[3];
	int color[3]={2,3,4};
	int style[3]={33,21,28};
	int size[3]={2,1,2};
	int linestyle[3]={1,10,9};
	TFile *t[3];
	TFile *output;
	TH1D *dy[3];
	TH1D *dpt[3];
	TH1D *ql_qinv[3];
	TH1D *qs_qinv[3];
	TH1D *qo_qinv[3];
	TH1D *q_inv[3];
	TH1D *q_beam[3];
	TH1D *q_side[3];
	TH1D *q_out[3];
	TH2D *q_inv_pt[3];
	TH2D *q_beam_pt[3];
	TH2D *q_side_pt[3];
	TH2D *q_out_pt[3];
	TH2D *q_inv_phi[3];
	TH2D *q_beam_phi[3];
	TH2D *q_side_phi[3];
	TH2D *q_out_phi[3];
	TCanvas *c1=new TCanvas();
	double xmin=0;
	double xmax=0.2;
	double legdx1=0.65;
	double legdy1=0.7;
	double legdx2=0.9;
	double legdy2=0.9;
int par5;
int par6;
int iffix;

void init()
{
	tcn[0]="../Triangle/result_root/result-all.root";
	tcn[1]="../Chain/result_root/result-all.root";
	tcn[2]="../nofluc/result_root/result-all.root";
}

void make(TH1D *plot,int i)
{
		plot->SetMarkerSize(size[i]);
		plot->SetMarkerStyle(style[i]);
		plot->SetMarkerColor(color[i]);
		plot->GetXaxis()->SetTitle("q=(p1-p2)/GeV");
		plot->GetXaxis()->SetNdivisions(505);
		plot->GetXaxis()->SetLabelSize(0.07);
		plot->GetXaxis()->SetTitleSize(0.07);
		plot->GetYaxis()->SetTitle("C(q)");
		plot->GetYaxis()->SetNdivisions(503);
		plot->GetYaxis()->SetLabelSize(0.07);
		plot->GetYaxis()->SetTitleSize(0.07);
		//plot->GetYaxis()->SetRangeUser(0.98,1.31);
}

void make(TGraphErrors *plot,int i)
{
		plot->RemovePoint(0);
		plot->SetMarkerSize(size[i]);
		plot->SetMarkerStyle(style[i]);
		plot->SetMarkerColor(color[i]);
}

void make(TF1 *plot,int i)
{
		plot->SetLineStyle(linestyle[i]);
		plot->SetLineColor(color[i]);
}

void make(TMultiGraph* plot)
{
		/*plot->GetXaxis()->SetNdivisions(505);
		plot->GetXaxis()->SetLabelSize(0.07);
		plot->GetXaxis()->SetTitleSize(0.07);
		plot->GetXaxis()->SetRangeUser(0,1.2);
		plot->GetYaxis()->SetNdivisions(505);
		plot->GetYaxis()->SetLabelSize(0.07);
		plot->GetYaxis()->SetTitleSize(0.07);*/
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
        legd->AddEntry(plot[0],"triangle","p");
        legd->AddEntry(plot[1],"chain","lep");
        legd->AddEntry(plot[2],"nofluc","lep");
        legd->Draw();
}

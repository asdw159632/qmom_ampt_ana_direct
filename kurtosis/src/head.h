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
	TFile *t[3];
	TFile *output;
	TH1D *dy[3];
	TH1D *dpt[3];
	TH1D *ql_qinv[3];
	TH1D *qs_qinv[3];
	TH1D *qo_qinv[3];
	TH2D *q_beam_pt[3];
	TH2D *q_side_pt[3];
	TH2D *q_out_pt[3];
	TH2D *q_beam_phi[3];
	TH2D *q_side_phi[3];
	TH2D *q_out_phi[3];
	TCanvas *c1=new TCanvas();
	double legdx1=0.15;
	double legdy1=0.65;
	double legdx2=0.39;
	double legdy2=0.85;
	double xmin=0;
	double xmax=1.7;
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
}

void make(TGraphErrors *plot,int i)
{
		plot->SetMarkerSize(size[i]);
		plot->SetMarkerStyle(style[i]);
		plot->SetMarkerColor(color[i]);
}

void legd(TGraphErrors **plot)
{
	TLegend *legd=new TLegend(legdx1,legdy1,legdx2,legdy2);
        legd->AddEntry(plot[0],"triangle","lep");
        legd->AddEntry(plot[1],"chain","lep");
        legd->AddEntry(plot[2],"nofluc","lep");
	legd->Draw();
}

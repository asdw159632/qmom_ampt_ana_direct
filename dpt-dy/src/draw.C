#include "head.h"
void draw()
{
	init();
	for(int i=0;i<3;i++)
	{
		t[i]=new TFile(tcn[i]);
		dpt[i]=(TH1D*)t[i]->Get("dpt");
		ql_qinv[i]=(TH1D*)t[i]->Get("n_qs_qinv");
		make(dpt[i],i);
		make2(ql_qinv[i],i);
	}
	ql_qinv[0]->GetXaxis()->SetTitle("q_{s}/GeV");
	ql_qinv[0]->GetYaxis()->SetTitle("q_{qinv}/GeV");
	ql_qinv[0]->GetYaxis()->CenterTitle();
	ql_qinv[0]->SetTitle("");

	gStyle->SetOptTitle(kFALSE);
	gStyle->SetOptStat(0);
	ql_qinv[0]->Draw("E0 X0");
	ql_qinv[1]->Draw("same E0 X0");
	ql_qinv[2]->Draw("same E0 X0");
	auto *legd3=new TLegend(legdx1,legdy1,legdx2,legdy2);
	legd3->AddEntry(dpt[0],"triangle","lep");
	legd3->AddEntry(dpt[1],"chain","lep");
	legd3->AddEntry(dpt[2],"nofluc","lep");
	gPad->SetLogy(0);
	legd3->Draw();
}

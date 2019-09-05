#include "head.h"
#include "fit.h"
#include "this.h"

using namespace std;

int main()
{
	init();
	for(int i=0;i<3;i++)
	{
		t[i]=new TFile(tcn[i],"READ");

		q_inv[i]=(TH1D*)t[i]->Get("q_inv;1");
		q_beam[i]=(TH1D*)t[i]->Get("q_beam;1");
		q_side[i]=(TH1D*)t[i]->Get("q_side;1");
		q_out[i]=(TH1D*)t[i]->Get("q_out;1");

		ql_qinv[i]=(TH1D*)t[i]->Get("n_ql_qinv;1");
		qs_qinv[i]=(TH1D*)t[i]->Get("n_qs_qinv;1");
		qo_qinv[i]=(TH1D*)t[i]->Get("n_qo_qinv;1");

#ifdef pt
		q_inv_pt[i]=(TH2D*)t[i]->Get("q_inv_pt;1");
		q_beam_pt[i]=(TH2D*)t[i]->Get("q_beam_pt;1");
		q_side_pt[i]=(TH2D*)t[i]->Get("q_side_pt;1");
		q_out_pt[i]=(TH2D*)t[i]->Get("q_out_pt;1");

		q_inv_R[i]=new TGraphErrors(floor((q_inv_pt[i]->GetNbinsY()-cut)/pack));
		q_beam_R[i]=new TGraphErrors(floor((q_beam_pt[i]->GetNbinsY()-cut)/pack));
		q_side_R[i]=new TGraphErrors(floor((q_side_pt[i]->GetNbinsY()-cut)/pack));
		q_out_R[i]=new TGraphErrors(floor((q_out_pt[i]->GetNbinsY()-cut)/pack));
#endif
#ifdef phi
		q_inv_phi[i]=(TH2D*)t[i]->Get("q_inv_phi;1");
		q_beam_phi[i]=(TH2D*)t[i]->Get("q_beam_phi;1");
		q_side_phi[i]=(TH2D*)t[i]->Get("q_side_phi;1");
		q_out_phi[i]=(TH2D*)t[i]->Get("q_out_phi;1");

		q_inv_R[i]=new TGraphErrors(floor((q_inv_phi[i]->GetNbinsY()-cut)/pack));
		q_beam_R[i]=new TGraphErrors(floor((q_beam_phi[i]->GetNbinsY()-cut)/pack));
		q_side_R[i]=new TGraphErrors(floor((q_side_phi[i]->GetNbinsY()-cut)/pack));
		q_out_R[i]=new TGraphErrors(floor((q_out_phi[i]->GetNbinsY()-cut)/pack));
#endif
	
		readqinv(i,0,ql_qinv[i]);
		readqinv(i,1,qs_qinv[i]);
		readqinv(i,2,qo_qinv[i]);
	}
	output=new TFile(outputroot,"RECREATE");
	output->Close();
	iffix=0;
	par5=3;
	fit(q_inv,"qinv");
	par5=0;
	fit(q_beam,"q_beam");
	par5=1;
	fit(q_side,"q_sidewards");
	par5=2;
	fit(q_out,"q_outwards");

	iffix=0;
#ifdef pt
	par5=3;
	fit(q_inv_pt,q_inv_R,"q_inv_pt");
	par5=0;
	fit(q_beam_pt,q_beam_R,"q_beam_pt");
	par5=1;
	fit(q_side_pt,q_side_R,"q_side_pt");
	par5=2;
	fit(q_out_pt,q_out_R,"q_out_pt");
	
	add();
	XY();
	
	output=new TFile(outputroot,"UPDATE");
	make(mg_inv);
	mg_inv->Draw("ap");
	legd(q_inv_R);
	c1->Print("qinv_pt_R.gif");
	c1->SetName("qinv_pt_R");
	c1->Write();
	c1->Clear();

	make(mg_beam);
	mg_beam->Draw("ap");
	legd(q_beam_R);
	c1->Print("beam_pt_R.gif");
	c1->SetName("beam_pt_R");
	c1->Write();
	c1->Clear();

	make(mg_side);
	mg_side->Draw("ap");
	legd(q_side_R);
	c1->Print("sidewards_pt_R.gif");
	c1->SetName("side_pt_R");
	c1->Write();
	c1->Clear();

	make(mg_out);
	mg_out->Draw("ap");
	legd(q_out_R);
	c1->Print("outwards_pt_R.gif");
	c1->SetName("out_pt_R");
	c1->Write();
	c1->Clear();
#endif

#ifdef phi
	par5=3;
	fit(q_inv_phi,q_inv_R,"q_inv_phi");
	par5=0;
	fit(q_beam_phi,q_beam_R,"q_beam_phi");
	par5=1;
	fit(q_side_phi,q_side_R,"q_side_phi");
	par5=2;
	fit(q_out_phi,q_out_R,"q_out_phi");
	
	add();
	XY();
	
	output=new TFile(outputroot,"UPDATE");
	make(mg_inv);
	mg_inv->Draw("ap");
	legd(q_inv_R);
	c1->Print("qinv_phi_R.gif");
	c1->SetName("qinv_phi_R");
	c1->Write();
	c1->Clear();

	make(mg_beam);
	mg_beam->Draw("ap");
	legd(q_beam_R);
	c1->Print("beam_phi_R.gif");
	c1->SetName("beam_phi_R");
	c1->Write();
	c1->Clear();

	make(mg_side);
	mg_side->Draw("ap");
	legd(q_side_R);
	c1->Print("sidewards_phi_R.gif");
	c1->SetName("side_phi_R");
	c1->Write();
	c1->Clear();

	make(mg_out);
	mg_out->Draw("ap");
	legd(q_out_R);
	c1->Print("outwards_phi_R.gif");
	c1->SetName("out_phi_R");
	c1->Write();
	c1->Clear();
#endif

	t[0]->Close();
	t[1]->Close();
	t[2]->Close();
	output->Close();
	
	return 0;
}

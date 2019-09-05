#include "head.h"

using namespace std;

int main()
{
	init();
	for(int i=0;i<3;i++)
	{
		t[i]=new TFile(tcn[i],"READ");

#ifdef dpt_dy
		dy[i]=(TH1D*)t[i]->Get("dy");
		dpt[i]=(TH1D*)t[i]->Get("dpt");
		make(dy[i],i);
		make(dpt[i],i);
#endif

#ifdef q_qinv
		ql_qinv[i]=(TH1D*)t[i]->Get("n_ql_qinv");
		qs_qinv[i]=(TH1D*)t[i]->Get("n_qs_qinv");
		qo_qinv[i]=(TH1D*)t[i]->Get("n_qo_qinv");
		make2(ql_qinv[i],i);
		make2(qs_qinv[i],i);
		make2(qo_qinv[i],i);
#endif

#ifdef wards_R
#ifdef pt
		beam_pt_R[i]=(TGraphErrors*)t[i]->Get("beam_pt_R");
		side_pt_R[i]=(TGraphErrors*)t[i]->Get("side_pt_R");
		out_pt_R[i]=(TGraphErrors*)t[i]->Get("out_pt_R");
		make(beam_pt_R[i],i);
		make(side_pt_R[i],i);
		make(out_pt_R[i],i);
		mg_beam_pt_R->Add(beam_pt_R[i]);
		mg_side_pt_R->Add(side_pt_R[i]);
		mg_out_pt_R->Add(out_pt_R[i]);
#endif
#ifdef phi
		beam_phi_R[i]=(TGraphErrors*)t[i]->Get("beam_phi_R");
		side_phi_R[i]=(TGraphErrors*)t[i]->Get("side_phi_R");
		out_phi_R[i]=(TGraphErrors*)t[i]->Get("out_phi_R");
		make(beam_phi_R[i],i);
		make(side_phi_R[i],i);
		make(out_phi_R[i],i);
		mg_beam_phi_R->Add(beam_phi_R[i]);
		mg_side_phi_R->Add(side_phi_R[i]);
		mg_out_phi_R->Add(out_phi_R[i]);
#endif
#endif

#ifdef kurtosis
#ifdef pt
		beam_pt_k[i]=(TGraphErrors*)t[i]->Get("beam-pt-Kurtosis");
		side_pt_k[i]=(TGraphErrors*)t[i]->Get("sidewards-pt-Kurtosis");
		out_pt_k[i]=(TGraphErrors*)t[i]->Get("outwards-pt-Kurtosis");
		make(beam_pt_k[i],i);
		make(side_pt_k[i],i);
		make(out_pt_k[i],i);
		mg_beam_pt_k->Add(beam_pt_k[i]);
		mg_side_pt_k->Add(side_pt_k[i]);
		mg_out_pt_k->Add(out_pt_k[i]);
#endif
#ifdef phi
		beam_phi_k[i]=(TGraphErrors*)t[i]->Get("beam-phi-Kurtosis");
		side_phi_k[i]=(TGraphErrors*)t[i]->Get("sidewards-phi-Kurtosis");
		out_phi_k[i]=(TGraphErrors*)t[i]->Get("outwards-phi-Kurtosis");
		make(beam_phi_k[i],i);
		make(side_phi_k[i],i);
		make(out_phi_k[i],i);
		mg_beam_phi_k->Add(beam_phi_k[i]);
		mg_side_phi_k->Add(side_phi_k[i]);
		mg_out_phi_k->Add(out_phi_k[i]);
#endif
#endif

#ifdef skewness
#ifdef pt
		beam_pt_s[i]=(TGraphErrors*)t[i]->Get("beam-pt-Skewness");
		side_pt_s[i]=(TGraphErrors*)t[i]->Get("sidewards-pt-Skewness");
		out_pt_s[i]=(TGraphErrors*)t[i]->Get("outwards-pt-Skewness");
		make(beam_pt_s[i],i);
		make(side_pt_s[i],i);
		make(out_pt_s[i],i);
		mg_beam_pt_s->Add(beam_pt_s[i]);
		mg_side_pt_s->Add(side_pt_s[i]);
		mg_out_pt_s->Add(out_pt_s[i]);
#endif
#ifdef phi
		beam_phi_s[i]=(TGraphErrors*)t[i]->Get("beam-phi-Skewness");
		side_phi_s[i]=(TGraphErrors*)t[i]->Get("sidewards-phi-Skewness");
		out_phi_s[i]=(TGraphErrors*)t[i]->Get("outwards-phi-Skewness");
		make(beam_phi_s[i],i);
		make(side_phi_s[i],i);
		make(out_phi_s[i],i);
		mg_beam_phi_s->Add(beam_phi_s[i]);
		mg_side_phi_s->Add(side_phi_s[i]);
		mg_out_phi_s->Add(out_phi_s[i]);
#endif
#endif
	}
	XY();
	output=TFile::Open("./result-read.root","RECREATE");

	//gStyle->SetOptTitle(kFALSE);
	gStyle->SetOptStat(0);
#ifdef dpt_dy
	dy[0]->Draw("E0 X0");
	dy[1]->Draw("same E0 X0");
	dy[2]->Draw("same E0 X0");
	legd(dy);
	c1->SetName("dy");
	c1->Write();
	c1->Print("./dy.png");
	c1->Clear();

	dpt[0]->Draw("E0 X0");
	dpt[1]->Draw("same E0 X0");
	dpt[2]->Draw("same E0 X0");
	legd(dpt);
	gPad->SetLogy();
	c1->SetName("dpt");
	c1->Write();
	c1->Print("./dpt.png");
	c1->Clear();
#endif

#ifdef q_qinv
	gPad->SetLogy(0);
	ql_qinv[0]->Draw("PLC PMC");
	ql_qinv[1]->Draw("same PLC PMC");
	ql_qinv[2]->Draw("same PLC PMC");
	legd(ql_qinv);
	c1->Update();
	c1->SetName("ql_qinv");
	c1->Write();
	c1->Print("./ql_qinv.png");
	c1->Clear();

	qs_qinv[0]->Draw("PLC PMC");
	qs_qinv[1]->Draw("same PLC PMC");
	qs_qinv[2]->Draw("same PLC PMC");
	legd(qs_qinv);
	c1->SetName("qs_qinv");
	c1->Write();
	c1->Print("./qs_qinv.png");
	c1->Clear();

	qo_qinv[0]->Draw("PLC PMC");
	qo_qinv[1]->Draw("same PLC PMC");
	qo_qinv[2]->Draw("same PLC PMC");
	legd(qo_qinv);
	c1->SetName("qo_qinv");
	c1->Write();
	c1->Print("./qo_qinv.png");
	c1->Clear();
#endif

#ifdef wards_R
#ifdef pt
	mg_beam_pt_R->Draw("ap");
	legd(beam_pt_R);
	c1->SetName("beam_pt_R");
	c1->Write();
	c1->Print("beam_pt_R.gif");
	c1->Clear();

	mg_side_pt_R->Draw("ap");
	legd(side_pt_R);
	c1->SetName("side_pt_R");
	c1->Write();
	c1->Print("side_pt_R.gif");
	c1->Clear();

	mg_out_pt_R->Draw("ap");
	legd(out_pt_R);
	c1->SetName("out_pt_R");
	c1->Write();
	c1->Print("out_pt_R.gif");
	c1->Clear();
#endif
#ifdef phi
	mg_beam_phi_R->Draw("ap");
	legd(beam_phi_R);
	c1->SetName("beam_phi_R");
	c1->Write();
	c1->Print("beam_phi_R.gif");
	c1->Clear();

	mg_side_phi_R->Draw("ap");
	legd(side_phi_R);
	c1->SetName("side_phi_R");
	c1->Write();
	c1->Print("side_phi_R.gif");
	c1->Clear();

	mg_out_phi_R->Draw("ap");
	legd(out_phi_R);
	c1->SetName("out_phi_R");
	c1->Write();
	c1->Print("out_phi_R.gif");
	c1->Clear();
#endif
#endif

#ifdef kurtosis
#ifdef pt
	mg_beam_pt_k->Draw("ap");
	legd(beam_pt_k);
	c1->SetName("beam-pt-Kurtosis");
	c1->Write();
	c1->Print("beam-pt-Kurtosis.gif");
	c1->Clear();

	mg_side_pt_k->Draw("ap");
	legd(side_pt_k);
	c1->SetName("side-pt-Kurtosis");
	c1->Write();
	c1->Print("side-pt-Kurtosis.gif");
	c1->Clear();

	mg_out_pt_k->Draw("ap");
	legd(out_pt_k);
	c1->SetName("out-pt-Kurtosis");
	c1->Write();
	c1->Print("out-pt-Kurtosis.gif");
	c1->Clear();
#endif
#ifdef phi
	mg_beam_phi_k->Draw("ap");
	legd(beam_phi_k);
	c1->SetName("beam-phi-Kurtosis");
	c1->Write();
	c1->Print("beam-phi-Kurtosis.gif");
	c1->Clear();

	mg_side_phi_k->Draw("ap");
	legd(side_phi_k);
	c1->SetName("side-phi-Kurtosis");
	c1->Write();
	c1->Print("side-phi-Kurtosis.gif");
	c1->Clear();

	mg_out_phi_k->Draw("ap");
	legd(out_phi_k);
	c1->SetName("out-phi-Kurtosis");
	c1->Write();
	c1->Print("out-phi-Kurtosis.gif");
	c1->Clear();
#endif
#endif

#ifdef skewness
#ifdef pt
	mg_beam_pt_s->Draw("ap");
	legd(beam_pt_s);
	c1->SetName("beam-pt-Skewness");
	c1->Write();
	c1->Print("beam-pt-Skewness.gif");
	c1->Clear();

	mg_side_pt_s->Draw("ap");
	legd(side_pt_s);
	c1->SetName("side-pt-Skewness");
	c1->Write();
	c1->Print("side-pt-Skewness.gif");
	c1->Clear();

	mg_out_pt_s->Draw("ap");
	legd(out_pt_s);
	c1->SetName("out-pt-Skewness");
	c1->Write();
	c1->Print("out-pt-Skewness.gif");
	c1->Clear();
#endif
#ifdef phi
	mg_beam_phi_s->Draw("ap");
	legd(beam_phi_s);
	c1->SetName("beam-phi-Skewness");
	c1->Write();
	c1->Print("beam-phi-Skewness.gif");
	c1->Clear();

	mg_side_phi_s->Draw("ap");
	legd(side_phi_s);
	c1->SetName("side-phi-Skewness");
	c1->Write();
	c1->Print("side-phi-Skewness.gif");
	c1->Clear();

	mg_out_phi_s->Draw("ap");
	legd(out_phi_s);
	c1->SetName("out-phi-Skewness");
	c1->Write();
	c1->Print("out-phi-Skewness.gif");
	c1->Clear();
#endif
#endif
	output->Close();
	t[0]->Close();
	t[1]->Close();
	t[2]->Close();
	return 0;
}

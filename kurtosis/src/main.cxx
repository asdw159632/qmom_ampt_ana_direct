#include "head.h"
#include "this.h"

using namespace std;

int main()
{
	init();

	for(int i=0;i<3;i++)
	{
		t[i]=new TFile(tcn[i],"READ");
#ifdef pt
		q_beam_pt[i]=(TH2D*)t[i]->Get("q_beam_pt");
		q_side_pt[i]=(TH2D*)t[i]->Get("q_side_pt");
		q_out_pt[i]=(TH2D*)t[i]->Get("q_out_pt");
		beam[i]=new TGraphErrors(floor((q_beam_pt[i]->GetNbinsY()-cut)/pack));
		side[i]=new TGraphErrors(floor((q_side_pt[i]->GetNbinsY()-cut)/pack));
		out[i]=new TGraphErrors(floor((q_out_pt[i]->GetNbinsY()-cut)/pack));
		print(q_beam_pt[i],beam[i],"beam");
		print(q_side_pt[i],side[i],"side");
		print(q_out_pt[i],out[i],"out");
#endif
#ifdef phi
		q_beam_phi[i]=(TH2D*)t[i]->Get("q_beam_phi");
		q_side_phi[i]=(TH2D*)t[i]->Get("q_side_phi");
		q_out_phi[i]=(TH2D*)t[i]->Get("q_out_phi");
		beam[i]=new TGraphErrors(floor((q_beam_phi[i]->GetNbinsY()-cut)/pack));
		side[i]=new TGraphErrors(floor((q_side_phi[i]->GetNbinsY()-cut)/pack));
		out[i]=new TGraphErrors(floor((q_out_phi[i]->GetNbinsY()-cut)/pack));
		print(q_beam_phi[i],beam[i],"beam");
		print(q_side_phi[i],side[i],"side");
		print(q_out_phi[i],out[i],"out");
#endif
		beam[i]->Draw();
		make(beam[i],i);
		make(side[i],i);
		make(out[i],i);
		mg_l->Add(beam[i]);
		mg_s->Add(side[i]);
		mg_o->Add(out[i]);
	}
	XY();
	output=new TFile("./result-kurtosis.root","RECREATE");
	mg_l->Draw("ap");
	legd(beam);
	c1->SetName("beam");
	c1->Write();
	c1->Print("beam.gif");
	c1->Clear();
	mg_s->Draw("ap");
	legd(side);
	c1->SetName("side");
	c1->Write();
	c1->Print("side.gif");
	c1->Clear();
	mg_o->Draw("ap");
	legd(out);
	c1->SetName("out");
	c1->Write();
	c1->Print("out.gif");
	c1->Clear();
	output->Close();
	t[0]->Close();
	t[1]->Close();
	t[2]->Close();
}

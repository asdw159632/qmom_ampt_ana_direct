using namespace std;

#define pt
#define phi

	TGraphErrors *beam[3];
	TGraphErrors *side[3];
	TGraphErrors *out[3];
	TMultiGraph *mg_l=new TMultiGraph();
	TMultiGraph *mg_s=new TMultiGraph();
	TMultiGraph *mg_o=new TMultiGraph();
	const char *namel="beam";
	const char *names="sidewards";
	const char *nameo="outwards";
	int cut=2;
	int pack=1;

void XY()
{
	mg_l->SetName(namel);
	mg_s->SetName(names);
	mg_o->SetName(nameo);
	mg_l->SetTitle(namel);
	mg_s->SetTitle(names);
	mg_o->SetTitle(nameo);
	/*
#ifdef pt
	mg_l->GetXaxis()->SetTitle("pt/GeV");
	mg_l->GetXaxis()->SetRangeUser(xmin,xmax);
	mg_l->GetYaxis()->SetTitle("S_{beam}(pt)");
	mg_s->GetXaxis()->SetTitle("pt/GeV");
	mg_s->GetXaxis()->SetRangeUser(xmin,xmax);
	mg_s->GetYaxis()->SetTitle("S_{side}(pt)");
	mg_o->GetXaxis()->SetTitle("pt/GeV");
	mg_o->GetXaxis()->SetRangeUser(xmin,xmax);
	mg_o->GetYaxis()->SetTitle("S_{out}(pt)");
#endif
#ifdef phi
	mg_l->GetXaxis()->SetTitle("phi/^{#circ}");
	mg_l->GetXaxis()->SetRangeUser(xmin,xmax);
	mg_l->GetYaxis()->SetTitle("S_{beam}(phi)");
	mg_s->GetXaxis()->SetTitle("phi/^{#circ}");
	mg_s->GetXaxis()->SetRangeUser(xmin,xmax);
	mg_s->GetYaxis()->SetTitle("S_{side}(phi)");
	mg_o->GetXaxis()->SetTitle("phi/^{#circ}");
	mg_o->GetXaxis()->SetRangeUser(xmin,xmax);
	mg_o->GetYaxis()->SetTitle("S_{out}(phi)");
#endif
*/
}

double par(int pow2)
{
	int tt=pow2-1;
	double c2=1;
	double c4=1;
	while(tt>=1)
	{
		c2*=tt/2.;
		tt=tt-2;
	}
	tt=2*pow2-1;
	while(tt>=1)
	{
		c4*=tt/2.;
		tt=tt-2;
	}
	double c=c4/c2/c2;
	return c;
}

void print(TH2D *plot2D,TGraphErrors* graph,const char*FileOutput)
{
	double q_min;
	double q_max;
	int q_step;
	q_min=plot2D->GetXaxis()->GetBinCenter(1);
	q_step=plot2D->GetNbinsX();
	q_max=plot2D->GetXaxis()->GetBinCenter(q_step);
	TF1 *f1=new TF1("f1","1",q_min,q_max);
	TF1 *f_1=new TF1("f_1","pow(x,1)",q_min,q_max);
	TF1 *f_2=new TF1("f_2","pow(x,2)",q_min,q_max);
	TF1 *f_3=new TF1("f_3","pow(x,3)",q_min,q_max);
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
		error_1=sqrt(pow(error_1/result0,2)+pow(result_1*error0/result0,2));
		result_2=(plot_2->IntegralAndError(1,plot_2->GetNbinsX(),error_2,"width")-f_2->Integral(q_min,q_max))/result0;
		error_2=sqrt(pow(error_2/result0,2)+pow(result_2*error0/result0,2));
		result_3=(plot_3->IntegralAndError(1,plot_3->GetNbinsX(),error_3,"width")-f_3->Integral(q_min,q_max))/result0;
		error_3=sqrt(pow(error_3/result0,2)+pow(result_3*error0/result0,2));

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
		graph->SetPoint(i,(last+next)/2,result);
		//graph->SetPointError(i,0,0);
		graph->SetPointError(i,0,error);
		delete plot;
		delete plot_1;
		delete plot_2;
		delete plot_3;
	}
	graph->RemovePoint(0);
	graph->SetName(FileOutput);
	graph->SetTitle(FileOutput);
}

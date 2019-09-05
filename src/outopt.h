#define qinv_half

/*A & B*/
#define q_step 40
#ifdef qinv_half 
double q_max=0.1;
#else
double q_max=0.2;
#endif
double q_min=0;
double qinv=0;
double dq=(q_max-q_min)/q_step;
/*C*/
#define pt_step 5
double pt_max=0.7;
double pt_min=0.2;
double pt=0;
double dpt=(pt_max-pt_min)/pt_step;
/*D*/
#define phi_step 16
double phi_max=M_PI;
double phi_min=-M_PI;
double phi=0;
double dphi=(phi_max-phi_min)/phi_step;

//y
double y_min=-0.5;
double y_max=0.5;

//dpt-dy
int dy_step=120;
double dy_min=-6;
double dy_max=6;
int dpt_step=10;
double dpt_min=0;
double dpt_max=10;

//kurtosis
int cut=0;
int pack=1;

//plot
const char *drawopt="alp";
const char *drawopt2="E0 X0 P0";
//refer to https://root.cern.ch/doc/master/classTAttMarker.html;
double MarkerSize=1;
double MarkerStyle=21;	//22 for fill up triangle;
double MarkerColor=1;	//black;
double LineWidth=1;
double LineStyle=1;	//solid;
double LineColor=1;	//black;

//fit
double lambda[4]={0.6,0.6,0.6,0.6};
double R0=4;
double N=1;
double A=1;
double npar=6;
const char *fittype="ON REMW";
double FitColor=2;
double FitStyle=1;
double FitWidth=1;

//Line
double style=9;
double width=1;
double color=1;

#ifdef qinv_half
const char *line1="!k=(p1-p2)/2";
#else
const char *line1="!k=p1-p2";
#endif
const char *line2="k/GeV/c	C(k)	+/-	PairNum";

void MakePlot(TGraphErrors *plot,bool t)
{
	plot->SetMarkerSize(MarkerSize);
	plot->SetMarkerStyle(MarkerStyle);
	plot->SetMarkerColor(MarkerColor);
	plot->SetLineWidth(LineWidth);
	plot->SetLineStyle(LineStyle);
	plot->SetLineColor(LineColor);
	if(t==1)
		plot->GetXaxis()->SetTitle("pt(GeV/c)");
	else
		plot->GetXaxis()->SetTitle("phi(rad)");
	plot->GetYaxis()->SetTitle("R/fm");
	plot->GetYaxis()->CenterTitle();
}

void MakePlot(TProfile *plot)
{
	plot->SetMarkerSize(MarkerSize);
	plot->SetMarkerStyle(MarkerStyle);
	plot->SetMarkerColor(MarkerColor);
	plot->SetLineWidth(LineWidth);
	plot->SetLineStyle(LineStyle);
	plot->SetLineColor(LineColor);
	plot->GetXaxis()->SetTitle("k(GeV/c)");
	plot->GetYaxis()->SetTitle("Ck");
	plot->GetYaxis()->CenterTitle();
}

void MakePlot(TH1 *plot)
{
	plot->SetMarkerSize(MarkerSize);
	plot->SetMarkerStyle(MarkerStyle);
	plot->SetMarkerColor(MarkerColor);
	plot->SetLineWidth(LineWidth);
	plot->SetLineStyle(LineStyle);
	plot->SetLineColor(LineColor);
}

void MakeFit(TF1 *g)
{
	g->SetParNames("lambda","R","R0","n","A");
	if(lambda_fixed==0){
	        g->SetParameter(0,lambda[par5]);
        	g->SetParLimits(0,0.,1.1);
	}
	else
	        g->FixParameter(0,lambda[par5]);
        g->SetParameter(1,3);
        g->SetParLimits(1,0,20);
        g->FixParameter(2,R0);
        g->FixParameter(3,N);
        g->FixParameter(4,A);
	g->FixParameter(5,par5);
	g->SetLineColor(FitColor);
	g->SetLineWidth(FitWidth);
	g->SetLineStyle(FitStyle);
}

void MakeLine(TF1 *line)
{
	line->SetLineStyle(style);
	line->SetLineWidth(width);
	line->SetLineColor(color);
}

void MakeT_RPlot(TH1 *plot)
{
	plot->SetMarkerSize(MarkerSize);
	plot->SetMarkerStyle(MarkerStyle);
	plot->SetMarkerColor(MarkerColor);
	plot->SetLineWidth(LineWidth);
	plot->SetLineStyle(LineStyle);
	plot->SetLineColor(LineColor);
	plot->GetXaxis()->SetTitle("T_R");
	plot->GetYaxis()->SetTitle("N");
	plot->GetYaxis()->CenterTitle();
}

void MakeT_RPlot(TProfile *plot)
{
	plot->SetMarkerSize(MarkerSize);
	plot->SetMarkerStyle(MarkerStyle);
	plot->SetMarkerColor(MarkerColor);
	plot->SetLineWidth(LineWidth);
	plot->SetLineStyle(LineStyle);
	plot->SetLineColor(LineColor);
	plot->GetXaxis()->SetTitle("q/GeV");
	plot->GetYaxis()->SetTitle("qinv/GeV");
	plot->GetYaxis()->CenterTitle();
}

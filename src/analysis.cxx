#include "head.h"
#include "psi.h"
//#include "psi-pion.h"
//#include "analysis.h"
//#include "particles.h"
#include "AMPT.h"
//#include "./create_list_d.cxx"

using namespace std;

bool condition(Int_t IDAART, Double_t PxAART, Double_t PyAART, Double_t PzAART, Double_t EnergyAART, Double_t X, Double_t Y, Double_t Z, Double_t T,char *pair);
double psi2(AMPT* event);
void participant_plane_correction(Double_t *px, Double_t *py, Double_t *x, Double_t *y, Double_t psi);
int makelist(int argc,char* argv[]);


int main(int argc,char *argv[])
{
  cout<<"//////////////////////////////////////////////"<<endl;
  cout<<endl;
  cout<<"this program start at "<<endl;
  system("date");
  cout<<endl;
  cout<<"//////////////////////////////////////////////"<<endl;

  //TChain *chain = new TChain("AMPT");

  if(argc!=3&&argc!=4) {cout<<"argc error: argc="<<argc<<endl;return 0;}

  const char *FileInput=0;
//  const char *FileOutput=0;
//  char *FilePar=0;

 /* if(argc==1)
  {
    FileInput  = "example.list";
    FileOutput = "example.root";
  }*/

/*  if(argc==3)//执行文件 list文件 任意输出文件名
  {
    FileInput = argv[1];
    FileOutput = argv[2];
  }*/

  if(argc==3)
  { 
    FileInput = "../list/rawdatapath.list";
  }
  Getbin(q,q_pt,q_phi);

  int Noutput=pt_step*4+phi_step*4+4;
  char FileOutput[Noutput][100];
  Fileout(FileOutput,q_pt,q_phi);

#ifdef byrun
for(runnum=1;runnum<2;runnum++){
#endif  
  makelist(runnum,argv);
/*  if(argc==4)
  {
    makelist(3,argv);
    FileInput = "../list/rawdatapath.list";
  }*/
  char FileList[512];
  string temps;
//  char outputlist[10000][512];

  TChain *chain = new TChain("AMPT");
  ifstream inputStream;
  inputStream.open(FileInput);
  if (!(inputStream))
  {
    printf("can not open input-list file\n");
    return 0;
  }
  while(inputStream.good())
  {
    inputStream.getline(FileList,512);

    if(strstr(FileList,".root")==NULL)
    {
      printf("%s is not a root-file adress!!!\n",FileList);
      continue;
    }
    if  ( inputStream.good() )
    {
      TFile *ftmp = new TFile(FileList);
      if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys()))
      {
        printf(" file %s error in opening!!!\n",FileList);
      }
      else
      {
        printf(" read in file %s\n",FileList);
        chain->Add(FileList);
      }
      delete ftmp;
    }
    //if (datnumber>lastdatnumber)lastdatnumber=datnumber;
  }
  //printf("flie number%d\n",fileNumber);
  inputStream.close();

  AMPT* ampt = new AMPT();
  ampt->SetData(chain);

  int mNEvents = (int)chain->GetEntries();
  Double_t impactpar=0;
  Double_t Xavg=0;
  Double_t Yavg=0;
  TVector2 theRho;
  Double_t QPsi_2cos=0;
  Double_t QPsi_2sin=0;
  Double_t Psi_2pp=0;
  Double_t QPsi_3cos=0;
  Double_t QPsi_3sin=0;
  Double_t Psi_3pp=0;
  Int_t n_avg=0;
  Int_t IDAART=0;
  Double_t PxAART=0;
  Double_t PyAART=0;
  Double_t PzAART=0;
  Double_t EnergyAART=0;
  Double_t MassAART=0;
  Double_t XAART=0;
  Double_t YAART=0;
  Double_t ZAART=0;
  Double_t TAART=0;
  TVector2 xy;
  TVector3 mR;
  TVector3 mP;

  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector r1;
  TLorentzVector r2;
  TVector3 q_osl;
  double Psi_12;

  double Rs=0;
  double Ro=0;
  double Rl=0;

  /*A*/
  TProfile *q_inv=new TProfile("q_inv","q_inv",q_step,q);
  /*B*/
  TProfile *q_beam=new TProfile("q_beam","q_beam",q_step,q);
  TProfile *q_side=new TProfile("q_side","q_side",q_step,q);
  TProfile *q_out=new TProfile("q_out","q_out",q_step,q);
  /*C*/
  TProfile2D *q_inv_pt=new TProfile2D("q_inv_pt","q_inv_pt",q_step,q,pt_step,q_pt);
  TProfile2D *q_out_pt=new TProfile2D("q_out_pt","q_out_pt",q_step,q,pt_step,q_pt);
  TProfile2D *q_beam_pt=new TProfile2D("q_beam_pt","q_beam_pt",q_step,q,pt_step,q_pt);
  TProfile2D *q_side_pt=new TProfile2D("q_side_pt","q_side_pt",q_step,q,pt_step,q_pt);

  TH2D *R2_side_pt=new TH2D("R2_side_pt","R2_side_pt",1000,0,10,pt_step,q_pt);
  TH2D *R2_out_pt=new TH2D("R2_out_pt","R2_out_pt",1000,0,10,pt_step,q_pt);
  /*D*/
  TProfile2D *q_inv_phi=new TProfile2D("q_inv_phi","q_inv_phi",q_step,q,phi_step,q_phi);
  TProfile2D *q_out_phi=new TProfile2D("q_out_phi","q_out_phi",q_step,q,phi_step,q_phi);
  TProfile2D *q_beam_phi=new TProfile2D("q_beam_phi","q_beam_phi",q_step,q,phi_step,q_phi);
  TProfile2D *q_side_phi=new TProfile2D("q_side_phi","q_side_phi",q_step,q,phi_step,q_phi);

  TH2D *R2_side_phi=new TH2D("R2_side_phi","R2_side_phi",1000,0,10,phi_step,q_phi);
  TH2D *R2_out_phi=new TH2D("R2_out_phi","R2_out_phi",1000,0,10,phi_step,q_phi);

  /*dpt-dy*/
  TH1D *dpt=new TH1D("dpt","dpt",dpt_step,dpt_min,dpt_max);
  TH1D *dy=new TH1D("dy","dy",dy_step,dy_min,dy_max);

#define n_T_R_plot
#ifdef n_T_R_plot
  /*t-r*/
  TH1D *n_T_R=new TH1D("n_T_R","n_T_R",20,0,100);
  TProfile *n_ql_qinv=new TProfile("n_ql_qinv","n_ql_qinv",q_step*10,q_min,q_max);
  TProfile *n_qs_qinv=new TProfile("n_qs_qinv","n_qs_qinv",q_step*10,q_min,q_max);
  TProfile *n_qo_qinv=new TProfile("n_qo_qinv","n_qo_qinv",q_step*10,q_min,q_max);
#endif

//#define arttime
#ifdef arttime
  Double_t artTime=100.;
#endif
  Double_t temp[10][25000];
  memset(temp,0,sizeof(temp));

  ofstream fileout;
/*  if(mode==1){
  fileout.open(FileOutput);
  fileout<<"line1"<<endl;
  fileout<<"line2"<<endl;
  fileout<<"line3"<<endl;
}*/

  //if(test==1||fileNumber%2==1){

  int TrackCount = 0;
  int EventCount = 0;
//  Double_t psi_2 = 0;
  for(int mNEventCount = 0; mNEventCount < mNEvents; mNEventCount++)
  {

    chain->GetEntry(mNEventCount);
    impactpar = ampt->ImPar();

    //if(impactpar<1)
    if(impactpar<5)
    {
//Boost
      Xavg=0;
      Yavg=0;
      QPsi_2cos=0;
      QPsi_2sin=0;
      Psi_2pp=0;
      QPsi_3cos=0;
      QPsi_3sin=0;
      Psi_3pp=0;
      n_avg=0;
      for(int track=0; track < ampt->Ntrack_HZPC(); track++)
      {
        Xavg += ampt->mXHZPC(track);
        Yavg += ampt->mYHZPC(track);
        n_avg ++;
      }
      Xavg /= n_avg;
      Yavg /= n_avg;
      for(int track=0; track < ampt->Ntrack_HZPC(); track++)
      {
        theRho.Set(ampt->mXHZPC(track)-Xavg,ampt->mYHZPC(track)-Yavg);
        QPsi_2cos += pow(theRho.Mod(),2)*cos(2.*theRho.Phi());
        QPsi_2sin += pow(theRho.Mod(),2)*sin(2.*theRho.Phi());
        QPsi_3cos += pow(theRho.Mod(),2)*cos(3.*theRho.Phi());
        QPsi_3sin += pow(theRho.Mod(),2)*sin(3.*theRho.Phi());
      }
      Psi_2pp=(atan2(QPsi_2sin,QPsi_2cos)+M_PI)/2.;
      Psi_3pp=(atan2(QPsi_3sin,QPsi_3cos)+M_PI)/3.;
      int TrackCount = 0;
//Reading Data
      for(int track=0; track < ampt->Ntrack_AART(); track++)
      {
        IDAART = ampt->mIDAART(track);
        PxAART = ampt->mPxAART(track);
        PyAART = ampt->mPyAART(track);
        PzAART = ampt->mPzAART(track);
        EnergyAART = ampt->mEnergyAART(track);
        MassAART = sqrt(pow(EnergyAART,2)-pow(PxAART,2)-pow(PyAART,2)-pow(PzAART,2));
        if(isnan(MassAART)) MassAART = 0;
        XAART = ampt->mXAART(track);
        YAART = ampt->mYAART(track);
	ZAART = ampt->mZAART(track);
	TAART = ampt->mTimeAART(track);

	xy.Set(XAART-Xavg,YAART-Yavg);
	xy.SetMagPhi(xy.Mod(),xy.Phi()-Psi_2pp);
	XAART=xy.X();
	YAART=xy.Y();
	xy.Set(PxAART,PyAART);
	xy.SetMagPhi(xy.Mod(),xy.Phi()-Psi_2pp);
	PxAART=xy.X();
	PyAART=xy.Y();

	Double_t Pt = sqrt(pow(PxAART,2)+pow(PyAART,2));
	Double_t y = 0.5*log((EnergyAART+PzAART)/(EnergyAART-PzAART));
	/*dpt-dy*/
	if(IDAART==PID1||IDAART==PID2){
	  dy->Fill(y);
	  if(y_min<y && y<y_max)dpt->Fill(Pt);
	}
	else
	  continue;
	if(pt_min<Pt && Pt<pt_max && y_min<y && y<y_max)
	{
//	  participant_plane_correction(&PxAART, &PyAART, &XAART, &YAART, psi_2);
	  mR.SetXYZ(XAART,YAART,ZAART);
	  mP.SetXYZ(PxAART,PyAART,PzAART);
#ifdef arttime
	  mR = mR + (artTime-TAART)/EnergyAART*mP;
	  if(mR.Mag()<0.00001&&mR.Mag()>40)continue;
#endif
//	  if(TAART<20)continue;
	  TrackCount = TrackCount+1;
	  temp[0][TrackCount]=IDAART;
	  temp[1][TrackCount]=PxAART;
	  temp[2][TrackCount]=PyAART;
	  temp[3][TrackCount]=PzAART;
	  temp[4][TrackCount]=EnergyAART;
	  temp[5][TrackCount]=MassAART;
	  temp[6][TrackCount]=mR.X();
	  temp[7][TrackCount]=mR.Y();
	  temp[8][TrackCount]=mR.Z();
#ifdef arttime
	  temp[9][TrackCount]=artTime;
#else
	  temp[9][TrackCount]=TAART;
#endif
	}
	else
	  continue;
      }

      if(TrackCount>0)
      {
        EventCount = EventCount+1;
        for (int itrack=1;itrack<=TrackCount;itrack++){
	  if(temp[0][itrack]!=PID1)continue;
	  p1.SetPxPyPzE(temp[1][itrack],temp[2][itrack],temp[3][itrack],temp[4][itrack]);
	  r1.SetXYZT(temp[6][itrack],temp[7][itrack],temp[8][itrack],temp[9][itrack]);
	  for(int jtrack=itrack+1;jtrack<=TrackCount;jtrack++){
	    if(temp[0][jtrack]!=PID2)continue;
	    p2.SetPxPyPzE(temp[1][jtrack],temp[2][jtrack],temp[3][jtrack],temp[4][jtrack]);
	    r2.SetXYZT(temp[6][jtrack],temp[7][jtrack],temp[8][jtrack],temp[9][jtrack]);
	    Getq(p1,p2,&q_osl,&qinv,&pt,&phi);
	    Getr(p1,p2,r1,r2,&Rs,&Ro,&Rl);
#ifdef test
	    return 0;
#endif
	    Psi_12=Psi_ij(p1,p2,r1,r2);
//	    if(floor(qinv/dq)==0)cout<<Psi_12<<endl;
	  
	    /*A*/
	    q_inv->Fill(qinv,Psi_12);
	    //if(Psi_12>3)cout<<(r1-r2).Vect()*(p1-p2).Vect()*0.5<<" "<<Psi_12<<endl;
	    /*B*/
	    if(fabs(q_osl.Px())<0.03 && fabs(q_osl.Py())<0.03)q_beam->Fill(q_osl.Pz(),Psi_12);
	    if(fabs(q_osl.Px())<0.03 && fabs(q_osl.Pz())<0.03)q_side->Fill(q_osl.Py(),Psi_12);
	    if(fabs(q_osl.Py())<0.03 && fabs(q_osl.Pz())<0.03)q_out->Fill(q_osl.Px(),Psi_12);
	    /*C*/
	    q_inv_pt->Fill(qinv,pt,Psi_12);
	    if(fabs(q_osl.Px())<0.03 && fabs(q_osl.Py())<0.03){
	      q_beam_pt->Fill(q_osl.Pz(),pt,Psi_12);
	      n_ql_qinv->Fill(q_osl.Pz(),qinv);
	    }
	    if(fabs(q_osl.Px())<0.03 && fabs(q_osl.Pz())<0.03){
	      q_side_pt->Fill(q_osl.Py(),pt,Psi_12);
	      R2_side_pt->Fill(Rs,pt);
	      n_qs_qinv->Fill(q_osl.Py(),qinv);
	    }
	    if(fabs(q_osl.Py())<0.03 && fabs(q_osl.Pz())<0.03){
	      q_out_pt->Fill(q_osl.Px(),pt,Psi_12);
	      R2_out_pt->Fill(Ro,pt);
	      n_qo_qinv->Fill(q_osl.Px(),qinv);
	    }
	    /*D*/
	    q_inv_phi->Fill(qinv,phi,Psi_12);
	    if(fabs(q_osl.Px())<0.03 && fabs(q_osl.Py())<0.03){
	      q_beam_phi->Fill(q_osl.Pz(),phi,Psi_12);
	      R2_side_phi->Fill(Rs,phi);
	      //cout<<abs(Rs)<<endl;
	    }
	    if(fabs(q_osl.Px())<0.03 && fabs(q_osl.Pz())<0.03){
	      q_side_phi->Fill(q_osl.Py(),phi,Psi_12);
	      R2_out_phi->Fill(Ro,phi);
	    }
	    if(fabs(q_osl.Py())<0.03 && fabs(q_osl.Pz())<0.03){
	      q_out_phi->Fill(q_osl.Px(),phi,Psi_12);
	    }

#ifdef n_T_R_plot
	    /*t-r*/
	    TLorentzVector qq=(p1-p2);
	    TLorentzVector qlcms=(p1-p2);
	    TVector3 Kz=((p1+p2)*0.5).Vect();
	    double KE=((p1+p2)*0.5).E();
	    qq.Boost(-Kz*(1/KE));
	    Kz.SetX(0);
	    Kz.SetY(0);
	    qlcms.Boost(-Kz*(1/KE));
	    TLorentzVector rr=r1-r2;
	    rr.Boost(-Kz*(1/KE));
	    n_T_R->Fill(rr.Vect().Mag());
#endif
           }
	}
      }
//      if(EventCount==1000)break;
      TrackCount = 0;
      //psi_2 = psi2(ampt);
      //cout<<psi_2<<endl;
    }
    else
      continue;
    TrackCount = 0;
  }
  //cout<<TrackCount<<endl;
  //if(datnumber*2==fileNumber)
  delete ampt;
  delete chain;
  
  cout<<"Reading in "<<EventCount<<" events."<<endl;

  if(EventCount!=0){

    /*A*/
    write(q_inv,FileOutput[0]);

    /*B*/
    write(q_beam,FileOutput[1]);
    write(q_side,FileOutput[2]);
    write(q_out,FileOutput[3]);

    /*C*/
    write(q_inv_pt,"q_inv_pt");
    write(q_beam_pt,"q_beam_pt");
    write(q_side_pt,"q_side_pt");
    write(q_out_pt,"q_out_pt");

    write(R2_side_pt,"R2_side_pt");
    write(R2_out_pt,"R2_out_pt");

    /*D*/
    write(q_inv_phi,"q_inv_phi");
    write(q_beam_phi,"q_beam_phi");
    write(q_side_phi,"q_side_phi");
    write(q_out_phi,"q_out_phi");

    write(R2_side_phi,"R2_side_phi");
    write(R2_out_phi,"R2_out_phi");

    /*dpt-dy*/
    write(dpt,"dpt",EventCount);
    write(dy,"dy");

#ifdef n_T_R_plot
    /*t-r*/
    T_Rwrite(n_T_R,"n_T_R");
    T_Rwrite(n_ql_qinv,"n_ql_qinv");
    T_Rwrite(n_qs_qinv,"n_qs_qinv");
    T_Rwrite(n_qo_qinv,"n_qo_qinv");
#endif
#ifdef byrun
    char rootname[100];
    sprintf(rootname,"../result_root/result_%d.root",runnum);
    if(runnum==1){
      myhadd("../result_root/result-all.root",rootname);
      gSystem->CopyFile("../result_root/result-all.root","../result_root/result-all-old.root");
    }
    else if(runnum>1){
      myhadd("../result_root/result-all.root","../result_root/result-all-old.root",rootname);
      gSystem->CopyFile("../result_root/result-all.root","../result_root/result-all-old.root");
    }
#endif
  }
  /*clear memory*/
  /*A*/
  delete q_inv;
  /*B*/
  delete q_beam;
  delete q_side;
  delete q_out;
  /*C*/
  delete q_inv_pt;
  delete q_out_pt;
  delete q_beam_pt;
  delete q_side_pt;

  delete R2_out_pt;
  delete R2_side_pt;
  /*D*/
  delete q_inv_phi;
  delete q_out_phi;
  delete q_beam_phi;
  delete q_side_phi;

  delete R2_out_phi;
  delete R2_side_phi;

  /*dpt-dy*/
  delete dpt;
  delete dy;

#ifdef n_T_R_plot
  /*t-r*/
  delete n_T_R;
  delete n_ql_qinv;
  delete n_qs_qinv;
  delete n_qo_qinv;
#endif
#ifdef byrun
}

  TFile *all=new TFile("../result_root/result-all.root","READ");
#else
  TFile *all=new TFile("../result_root/result.root","READ");
#endif
  /*A*/
  TProfile *q_inv_r=(TProfile *)all->Get("q_inv");
  /*B*/
  TProfile *q_beam_r=(TProfile*)all->Get("q_beam");
  TProfile *q_side_r=(TProfile*)all->Get("q_side");
  TProfile *q_out_r=(TProfile*)all->Get("q_out");
  /*C*/
  TProfile2D *q_inv_pt_r=(TProfile2D *)all->Get("q_inv_pt");
  TProfile2D *q_beam_pt_r=(TProfile2D*)all->Get("q_beam_pt");
  TProfile2D *q_side_pt_r=(TProfile2D*)all->Get("q_side_pt");
  TProfile2D *q_out_pt_r=(TProfile2D*)all->Get("q_out_pt");
  
  TProfile2D *R2_side_pt_r=(TProfile2D*)all->Get("R2_side_pt");
  TProfile2D *R2_out_pt_r=(TProfile2D*)all->Get("R2_out_pt");
  /*D*/
  TProfile2D *q_inv_phi_r=(TProfile2D *)all->Get("q_inv_phi");
  TProfile2D *q_beam_phi_r=(TProfile2D*)all->Get("q_beam_phi");
  TProfile2D *q_side_phi_r=(TProfile2D*)all->Get("q_side_phi");
  TProfile2D *q_out_phi_r=(TProfile2D*)all->Get("q_out_phi");
  
  TProfile2D *R2_side_phi_r=(TProfile2D*)all->Get("R2_side_phi");
  TProfile2D *R2_out_phi_r=(TProfile2D*)all->Get("R2_out_phi");
  /*dpt-dy*/
  TH1D *dpt_r=(TH1D*)all->Get("dpt");
  TH1D *dy_r=(TH1D*)all->Get("dy");
#ifdef n_T_R_plot
  /*t-r*/
  TH1D *n_T_R_r=(TH1D*)all->Get("n_T_R");
  TProfile *n_ql_qinv_r=(TProfile*)all->Get("n_ql_qinv");
  TProfile *n_qs_qinv_r=(TProfile*)all->Get("n_qs_qinv");
  TProfile *n_qo_qinv_r=(TProfile*)all->Get("n_qo_qinv");
#endif
  
  lambda_fixed=0;

  /*A*/
  par5=0;
  print(q_inv_r,FileOutput[0]);
  for(int i=0;i<pt_step;i++){
      print(q_inv_pt_r,FileOutput[i+3*pt_step+3*phi_step+4],i);
  }
  print(q_inv_pt_r,"qinv-pt");
  for(int i=0;i<phi_step;i++){
      print(q_inv_phi_r,FileOutput[i+3*phi_step+4*pt_step+4],i);
  }
  print(q_inv_phi_r,"qinv-phi");

  /*B*/
  par5=1;
  readqinv(par5);
  print(q_beam_r,FileOutput[1]);
  par5=2;
  readqinv(par5);
  print(q_side_r,FileOutput[2]);
  par5=3;
  readqinv(par5);
  print(q_out_r,FileOutput[3]);

  lambda_fixed=0;

  /*C*/
  for(int i=0;i<pt_step;i++){
      par5=1;
      print(q_beam_pt_r,FileOutput[i+4],i);
      par5=2;
      print(q_side_pt_r,FileOutput[i+pt_step+4],i);
      par5=3;
      print(q_out_pt_r,FileOutput[i+2*pt_step+4],i);
  }
  print(q_beam_pt_r,"beam-pt");
  print(q_out_pt_r,"outwards-pt");
  print(q_side_pt_r,"sidewards-pt");

  print(R2_out_pt_r,"Ro2-pt");
  print(R2_side_pt_r,"Rs2-pt");

  /*D*/
  for(int i=0;i<phi_step;i++){
      par5=1;
      print(q_beam_phi_r,FileOutput[i+3*pt_step+4],i);
      par5=2;
      print(q_side_phi_r,FileOutput[i+3*pt_step+phi_step+4],i);
      par5=3;
      print(q_out_phi_r,FileOutput[i+3*pt_step+2*phi_step+4],i);
  }
  print(q_beam_phi_r,"beam-phi");
  print(q_out_phi_r,"outwards-phi");
  print(q_side_phi_r,"sidewards-phi");

  print(R2_out_phi_r,"Ro2-phi");
  print(R2_side_phi_r,"Rs2-phi");

  /*dpt-dy*/
  print(dpt_r,"dpt_n");
  print(dy_r,"dy_n");

#ifdef n_T_R_plot
  /*t-r*/
  T_Rprint(n_T_R_r,"n_T_R");
  T_Rprint(n_ql_qinv_r,"n_ql_qinv");
  T_Rprint(n_qs_qinv_r,"n_qs_qinv");
  T_Rprint(n_qo_qinv_r,"n_qo_qinv");
#endif
  all->Close();
  
  /*C*/
  print(beam_pt_R,"beam_pt_R");
  print(side_pt_R,"side_pt_R");
  print(out_pt_R,"out_pt_R");
  /*D*/
  print(beam_phi_R,"beam_phi_R");
  print(side_phi_R,"side_phi_R");
  print(out_phi_R,"out_phi_R");

  //delete chain;
  cout<<"//////////////////////////////////////////////"<<endl;
  cout<<endl;
  cout<<"this program finished at "<<endl;
  system("date");
  cout<<endl;
  cout<<"//////////////////////////////////////////////"<<endl;

  return 0;
}

bool condition(Int_t IDAART, Double_t PxAART, Double_t PyAART, Double_t PzAART, Double_t EnergyAART, Double_t X, Double_t Y, Double_t Z, Double_t T,char* pair)
{
//	if(sqrt(X*X+Y*Y+Z*Z)>60) return 0;
	int ID=0;
	Double_t Pt = sqrt(pow(PxAART,2)+pow(PyAART,2));
	Double_t y = 0.5*log((EnergyAART+PzAART)/(EnergyAART-PzAART));
	if(strcmp(pair,"pp")==0)ID=2212;
	if(strcmp(pair,"pi+pi+")==0)ID=211;
	if(strcmp(pair,"k+k+")==0)ID=321;
	//if(IDAART==211 && 0.1<Pt && Pt<3 && -1<y && y<1)
	//if(IDAART==333 && 0.2<Pt && Pt<1.5 && -1<y && y<1)
	//if(IDAART==2212 && 0.2<Pt && Pt<1.5 && -1<y && y<1)
	if(IDAART==ID && pt_min<Pt && Pt<pt_max && y_min<y && y<y_max)
	//if(IDAART==321 && 0.2<Pt && Pt<1.5 && -1<y && y<1)
	//if(IDAART==3122 && 0.2<Pt && Pt<1.5 && -1<y && y<1)
		return 1;
	else 
		return 0;
}

double psi2(AMPT* event)
{
  double PI = TMath::Pi();
  TVector2 participant_nucleon;
  Double_t r2_cos = 0;
  Double_t r2_sin = 0;

  int multi_INUCP = event->Ntrack_INUCP();
  for(int i=0; i<multi_INUCP; i++)
  {
    if(event->mNCOLLP(i)<1) continue;
    participant_nucleon.Set(event->mXINUCP(i),event->mYINUCP(i));
    r2_cos += pow(participant_nucleon.Mod(),2)*cos(2.*participant_nucleon.Phi());
    r2_sin += pow(participant_nucleon.Mod(),2)*sin(2.*participant_nucleon.Phi());
  }
  int multi_INUCT = event->Ntrack_INUCT();
  for(int i=0; i<multi_INUCT; i++)
  {
    if(event->mNCOLLT(i)<1) continue;
    participant_nucleon.Set(event->mXINUCT(i),event->mYINUCT(i));
    r2_cos += pow(participant_nucleon.Mod(),2)*cos(2.*participant_nucleon.Phi());
    r2_sin += pow(participant_nucleon.Mod(),2)*sin(2.*participant_nucleon.Phi());
  }
  Double_t psi = (atan2(r2_sin,r2_cos)+PI)/2.;
  return psi;
}

void participant_plane_correction(Double_t *px, Double_t *py, Double_t *x, Double_t *y, Double_t psi)
{
	Double_t pt = sqrt(pow(*px,2)+pow(*py,2));
	Double_t r = sqrt(pow(*x,2)+pow(*y,2));
	Double_t phi_pt = atan2(*py,*px)-psi;
	Double_t phi = atan2(*y,*x)-psi;
	*px = pt*cos(phi_pt);
	*py = pt*sin(phi_pt);
	*x = r*cos(phi);
	*y = r*sin(phi);
}


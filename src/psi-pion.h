TLorentzVector trace(TLorentzVector q,TLorentzVector K)
{
	TLorentzVector trace=((q*K)/K.Mag2())*K;
	return trace;
}

TComplex ConHyperg(TComplex a, TComplex b, TComplex z)
{
  //here b=1

  TComplex F(1.,0.);//n=0

  TComplex an;

  for(int n=1; n<10; n++)
  {
    an=a;
    for(int in=0; in<n-1; in++)
    {
      an *= (a+in*1.0);
    }
    F += an*TComplex::Power(z,n)/pow(TMath::Factorial(n),2);
  }

  return F;
}

Double_t Psi_ij(TLorentzVector p1Lab,TLorentzVector p2Lab,TLorentzVector r1Lab,TLorentzVector r2Lab)
{
  /*
  double ETotalLab; //total energy of cluster
  TVector3 Pcmtolab; //3d-momenta of cluster
  TVector3 betalabtocm; //velocity of lab relative to center of cluster

  TLorentzVector pcpt1Rcm; //coordinates of particle in center of mass frame
  TLorentzVector pcpt2Rcm;

  double pcpt1Timecm;//time of particle in center of mass frame
  double pcpt2Timecm;
  double maxft; //the last freeze out time

  TVector3 pcpt1R3cm;//3d-coordinates of particle in center of mass frame
  TVector3 pcpt2R3cm;

  TLorentzVector pcpt1Pcm;//momenta of particle in center of mass frame
  TLorentzVector pcpt2Pcm;

  TVector3 pcpt1P3cm;//3d-momenta of particle in center of mass frame
  TVector3 pcpt2P3cm;

  TLorentzVector current_mP; //momenta of cluster

  ETotalLab = p1Lab.Energy() + p2Lab.Energy();
  current_mP = p1Lab+p2Lab;
  Pcmtolab = current_mP.Vect();

  betalabtocm = -1./ETotalLab*Pcmtolab;
  TLorentzRotation LRltoc(betalabtocm);

  pcpt1Rcm = LRltoc*r1Lab;//boost R(4d) from lb to cm frame
  pcpt2Rcm = LRltoc*r2Lab;

  pcpt1Timecm = pcpt1Rcm.T();
  pcpt2Timecm = pcpt2Rcm.T();

  maxft = TMath::Max(pcpt1Timecm,pcpt2Timecm);

  pcpt1R3cm = pcpt1Rcm.Vect();//extract R(3d) from R(4d)
  pcpt2R3cm = pcpt2Rcm.Vect();

  pcpt1Pcm = LRltoc*p1Lab;//boost P(4d) from lb to cm frame
  pcpt2Pcm = LRltoc*p2Lab;

  pcpt1P3cm = pcpt1Pcm.Vect();//extract P(3d) from P(4d)
  pcpt2P3cm = pcpt2Pcm.Vect();

  pcpt1R3cm = pcpt1R3cm + (maxft-pcpt1Timecm)/pcpt1Pcm.Energy()*pcpt1P3cm;//reset R(3) at freeze-out time
  pcpt2R3cm = pcpt2R3cm + (maxft-pcpt2Timecm)/pcpt2Pcm.Energy()*pcpt2P3cm;


  TLorentzVector p1=pcpt1Pcm;
  TLorentzVector p2=pcpt2Pcm;

  TLorentzVector r1(pcpt1R3cm, maxft);
  TLorentzVector r2(pcpt2R3cm, maxft);

  TLorentzVector q=0.5*(p1-p2);
  //TLorentzVector q=(p1-p2);
  TLorentzVector r=(r1-r2);
*/

  TLorentzVector q;
  TLorentzVector K;

  K=(p1Lab+p2Lab)*0.5;
  q=(p1Lab-p2Lab)*0.5-trace((p1Lab-p2Lab)*0.5,K);
  //q=(p1Lab-p2Lab)-trace((p1Lab-p2Lab),K);
  q.Boost(-K.Vect()*(1/K.E()));

  //q = 0.5*q;

  TLorentzVector r;
  //r=((r1Lab-r2Lab)*0.5);
  r=((r1Lab-r2Lab));
  r-=trace(r,K);
  r.Boost(-K.Vect()*(1/K.E()));
  //r*=fm_gev;

 /*
  TLorentzVector q = (p1Lab-p2Lab);
  TLorentzVector r = (r1Lab-r2Lab);
*/

  Double_t hbarc=0.19732696;//GeV fm
  //Double_t aBohr=388;//fm
  Double_t aBohr=137.036/(0.139/2.)*hbarc;//fm, pion


  TComplex CI(0., 1.);
  TComplex C1(1., 0.);

  if(q.P()<2.127e-4) return 0.;//related to 2.*PI*eta>15

  Double_t eta = hbarc/(q.P()*aBohr);
  //Double_t eta = 1./(q.P()*aBohr);
  Double_t ksi_qplus = (q.P()*r.Rho()-q.Vect()*r.Vect())/hbarc;
  Double_t ksi_qminus = (q.P()*r.Rho()+q.Vect()*r.Vect())/hbarc;

  Double_t Ac_eta = 2.*TMath::Pi()*eta/(exp(2.*TMath::Pi()*eta)-1.);

  //if(r.Rho()>aBohr*1e-4) ;
//  if(eta<1e-6) return (1+cos(2.*q.Vect()*r.Vect()));
//  if(eta>1e2)  return 0.;

//  if(r.Rho()>aBohr*5e-2) return 0.;//Ac_eta*(1+cos(2.*q.Vect()*r.Vect()));
  //if(q.P()*r.Rho()/hbarc>5) return 0;

  /*
  if(eta<1e-6) return (1+cos(2.*q.Vect()*r.Vect()));
  if(eta>1e2)  return 0.;
  if(q.P()*r.Rho()/hbarc>5) return 0;
*/
//  if(eta<1e-6) return (1+cos(2.*q.Vect()*r.Vect()));
//  if(eta>1e2) Ac_eta=0;
//  if(q.P()*r.Rho()/hbarc>5) return 0;

  //if(2.*PI*eta<1e-3) return 0.9999*(1+cos(2.*q.Vect()*r.Vect()/hbarc));
  //
  if(2.*TMath::Pi()*eta>15) return 0.;//Ac_eta=0;
  //if(r.Rho()>aBohr*1e-4) return Ac_eta*(1+cos(2.*q.Vect()*r.Vect()/hbarc));

  TComplex F_qplus;
  TComplex F_qminus;
/*
  //if(r.Rho()>aBohr*1e-4)
  if(r.Rho()>aBohr*1e-2)
  //if(r.Rho()>aBohr*5e-2)
  {
    if(r.Rho()>aBohr*4e-2)
    //if(r.Rho()<=aBohr*1e-2)
    {
      F_qplus = 1.;// + (1.-cos(q.Vect()*r.Vect()/q.P()/r.Rho()));//(-CI*eta)*(CI*ksi_qplus/r.Rho()*aBohr);
      F_qminus = 1.;// + (1.+cos(q.Vect()*r.Vect()/q.P()/r.Rho()));//(-CI*eta)*(CI*ksi_qminus/r.Rho()*aBohr);
    }
    else
    {
      F_qplus = 1. + r.Rho()/aBohr*(1.-cos(q.Vect()*r.Vect()/q.P()/r.Rho()));//(-CI*eta)*(CI*ksi_qplus/r.Rho()*aBohr);
      F_qminus = 1. + r.Rho()/aBohr*(1.+cos(q.Vect()*r.Vect()/q.P()/r.Rho()));//(-CI*eta)*(CI*ksi_qminus/r.Rho()*aBohr);
    }

    return 0.5*Ac_eta*(F_qplus*F_qplus+F_qminus*F_qminus+2.*F_qplus*F_qminus*cos(2.*q.Vect()*r.Vect()/hbarc));
  }
  else
  {
    F_qplus = ConHyperg(-CI*eta, C1, ksi_qplus);
    F_qminus = ConHyperg(-CI*eta, C1, ksi_qminus);
  }
*/
  if(r.Rho()*q.P()/hbarc>0.5) return 1.;
  F_qplus = ConHyperg(-CI*eta, C1, ksi_qplus);
  F_qminus = ConHyperg(-CI*eta, C1, ksi_qminus);

  TComplex psi = 0.5*Ac_eta*(F_qplus*TComplex::Conjugate(F_qplus)
         +TComplex::Exp(2.*CI*q.Vect()*r.Vect()/hbarc)*F_qplus*TComplex::Conjugate(F_qminus)
	 +TComplex::Exp(-2.*CI*q.Vect()*r.Vect()/hbarc)*F_qminus*TComplex::Conjugate(F_qplus)
	 +F_qminus*TComplex::Conjugate(F_qminus)
	 );

  //cout<<psi<<endl;
  return psi.Re();
}



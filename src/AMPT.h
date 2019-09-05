#ifndef AMPT_h
#define AMPT_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "TChain.h"
#include "TBranch.h"
#include "TString.h"
#include "TArrayI.h"
#include "TArrayL.h"
#include "TArrayF.h"

const Int_t kMaxtrack = 25000;
typedef struct { 
  Int_t        NEvent;

  Int_t        multiINUCP    ; // Multiplicity
  Int_t        multiINUCT    ; // Multiplicity
  Int_t        multiHZPC    ; // Multiplicity
  Int_t        multiZPC    ; // Multiplicity
  Int_t        multiBART    ; // Multiplicity
  Int_t        multiAART    ; // Multiplicity
  Float_t      impactpar; //impact parameter
  Int_t        NELP     ;//elastic participants projectile
  Int_t        NINP     ;//inelastic participants projectile
  Int_t        NELT     ;//elastic participants target
  Int_t        NINT     ;//inelastic participants target; NINT is a intrinsic fortran function, rename it to NINTHJ
  Int_t        N0       ;//N-N
  Int_t        N01      ;//N-Nwounded
  Int_t        N10      ;//Nwounded-N
  Int_t        N11      ;//Nwounded-Nwounded

  //initial nucleon
  //Projectile
//  Int_t        IDINUC[kMaxtrack];      //KF code of initial nucleon
  Float_t      fXINUCP[kMaxtrack];       //X coordinate of initial nucleon
  Float_t      fYINUCP[kMaxtrack];       //Y coordinate of initial nucleon
  Float_t      fZINUCP[kMaxtrack];       //Z coordinate of initial nucleon
  Int_t        NCOLLP[kMaxtrack];

  //target
  Float_t      fXINUCT[kMaxtrack];       //X coordinate of initial nucleon
  Float_t      fYINUCT[kMaxtrack];       //Y coordinate of initial nucleon
  Float_t      fZINUCT[kMaxtrack];       //Z coordinate of initial nucleon
  Int_t        NCOLLT[kMaxtrack];

  //after hijing before ZPC
  Int_t        IDHZPC[kMaxtrack];      //KF code of particle
  Float_t      fPxHZPC[kMaxtrack];           //X component of the momentum
  Float_t      fPyHZPC[kMaxtrack];           //Y component of the momentum
  Float_t      fPzHZPC[kMaxtrack];           //Z component of the momentum
  Float_t      fEnergyHZPC[kMaxtrack];        //The energy  of this particle
  Float_t      fXHZPC[kMaxtrack];       //X coordinate of the first point
  Float_t      fYHZPC[kMaxtrack];       //Y coordinate of the first point
  Float_t      fZHZPC[kMaxtrack];       //Z coordinate of the first point
  Float_t      fTHZPC[kMaxtrack];       //T coordinate of the first point

  //ZPC
  Int_t        IDZPC[kMaxtrack];      //KF code of particle
  Float_t      fPxZPC[kMaxtrack];           //X component of the momentum
  Float_t      fPyZPC[kMaxtrack];           //Y component of the momentum
  Float_t      fPzZPC[kMaxtrack];           //Z component of the momentum
  Float_t      fEnergyZPC[kMaxtrack];        //The energy  of this particle
  Float_t      fXZPC[kMaxtrack];       //X coordinate of the first point
  Float_t      fYZPC[kMaxtrack];       //Y coordinate of the first point
  Float_t      fZZPC[kMaxtrack];       //Z coordinate of the first point
  Float_t      fTZPC[kMaxtrack];       //T coordinate of the first point

  //before ART
  Int_t        IDBART[kMaxtrack];      //KF code of particle
  Float_t      fPxBART[kMaxtrack];           //X component of the momentum
  Float_t      fPyBART[kMaxtrack];           //Y component of the momentum
  Float_t      fPzBART[kMaxtrack];           //Z component of the momentum
  Float_t      fEnergyBART[kMaxtrack];        //The energy  of this particle
  Float_t      fXBART[kMaxtrack];       //X coordinate of the first point
  Float_t      fYBART[kMaxtrack];       //Y coordinate of the first point
  Float_t      fZBART[kMaxtrack];       //Z coordinate of the first point
  Float_t      fTBART[kMaxtrack];       //T coordinate of the first point

  //after ART
  Int_t        IDAART[kMaxtrack];      //KF code of particle
  Float_t      fPxAART[kMaxtrack];           //X component of the momentum
  Float_t      fPyAART[kMaxtrack];           //Y component of the momentum
  Float_t      fPzAART[kMaxtrack];           //Z component of the momentum
  Float_t      fEnergyAART[kMaxtrack];        //The energy  of this particle
  Float_t      fXAART[kMaxtrack];       //X coordinate of the first point
  Float_t      fYAART[kMaxtrack];       //Y coordinate of the first point
  Float_t      fZAART[kMaxtrack];       //Z coordinate of the first point
  Float_t      fTAART[kMaxtrack];       //T coordinate of the first point
} Cell_t;

class AMPT {
public :
   Cell_t cell;
  
   //event number
   Int_t ImEvent() const{return cell.NEvent;}
   //sys. info.
   Float_t ImPar() const{return cell.impactpar;}
   Int_t Npart() const{return cell.NELP+cell.NINP+cell.NELT+cell.NINT;}
   Int_t Nbin() const{return cell.N0+cell.N01+cell.N10+cell.N11;}

   Int_t Ntrack_INUCP() const{return cell.multiINUCP;}
   Int_t Ntrack_INUCT() const{return cell.multiINUCT;}
   Int_t Ntrack_HZPC() const{return cell.multiHZPC;}
   Int_t Ntrack_ZPC() const{return cell.multiZPC;}
   Int_t Ntrack_BART() const{return cell.multiBART;}
   Int_t Ntrack_AART() const{return cell.multiAART;}

   //track info.
//   Int_t   mIDINUC(Int_t i) const{return cell.IDINUC[i];}
   Float_t mXINUCP(Int_t i) const{return cell.fXINUCP[i];}
   Float_t mYINUCP(Int_t i) const{return cell.fYINUCP[i];}
   Float_t mZINUCP(Int_t i) const{return cell.fZINUCP[i];}
   Int_t   mNCOLLP(Int_t i) const{return cell.NCOLLP[i];}

   Float_t mXINUCT(Int_t i) const{return cell.fXINUCT[i];}
   Float_t mYINUCT(Int_t i) const{return cell.fYINUCT[i];}
   Float_t mZINUCT(Int_t i) const{return cell.fZINUCT[i];}
   Int_t   mNCOLLT(Int_t i) const{return cell.NCOLLT[i];}

   Int_t   mIDHZPC(Int_t i) const{return cell.IDHZPC[i];}
   Float_t mPxHZPC(Int_t i) const{return cell.fPxHZPC[i];}
   Float_t mPyHZPC(Int_t i) const{return cell.fPyHZPC[i];}
   Float_t mPzHZPC(Int_t i) const{return cell.fPzHZPC[i];}
   Float_t mEnergyHZPC(Int_t i) const{return cell.fEnergyHZPC[i];}
   Float_t mXHZPC(Int_t i) const{return cell.fXHZPC[i];}
   Float_t mYHZPC(Int_t i) const{return cell.fYHZPC[i];}
   Float_t mZHZPC(Int_t i) const{return cell.fZHZPC[i];}
   Float_t mTimeHZPC(Int_t i) const{return cell.fTHZPC[i];}

   Int_t   mIDZPC(Int_t i) const{return cell.IDZPC[i];}
   Float_t mPxZPC(Int_t i) const{return cell.fPxZPC[i];}
   Float_t mPyZPC(Int_t i) const{return cell.fPyZPC[i];}
   Float_t mPzZPC(Int_t i) const{return cell.fPzZPC[i];}
   Float_t mEnergyZPC(Int_t i) const{return cell.fEnergyZPC[i];}
   Float_t mXZPC(Int_t i) const{return cell.fXZPC[i];}
   Float_t mYZPC(Int_t i) const{return cell.fYZPC[i];}
   Float_t mZZPC(Int_t i) const{return cell.fZZPC[i];}
   Float_t mTimeZPC(Int_t i) const{return cell.fTZPC[i];}

   Int_t   mIDBART(Int_t i) const{return cell.IDBART[i];}
   Float_t mPxBART(Int_t i) const{return cell.fPxBART[i];}
   Float_t mPyBART(Int_t i) const{return cell.fPyBART[i];}
   Float_t mPzBART(Int_t i) const{return cell.fPzBART[i];}
   Float_t mEnergyBART(Int_t i) const{return cell.fEnergyBART[i];}
   Float_t mXBART(Int_t i) const{return cell.fXBART[i];}
   Float_t mYBART(Int_t i) const{return cell.fYBART[i];}
   Float_t mZBART(Int_t i) const{return cell.fZBART[i];}
   Float_t mTimeBART(Int_t i) const{return cell.fTBART[i];}

   Int_t   mIDAART(Int_t i) const{return cell.IDAART[i];}
   Float_t mPxAART(Int_t i) const{return cell.fPxAART[i];}
   Float_t mPyAART(Int_t i) const{return cell.fPyAART[i];}
   Float_t mPzAART(Int_t i) const{return cell.fPzAART[i];}
   Float_t mEnergyAART(Int_t i) const{return cell.fEnergyAART[i];}
   Float_t mXAART(Int_t i) const{return cell.fXAART[i];}
   Float_t mYAART(Int_t i) const{return cell.fYAART[i];}
   Float_t mZAART(Int_t i) const{return cell.fZAART[i];}
   Float_t mTimeAART(Int_t i) const{return cell.fTAART[i];}

   void SetData(TChain *);

};

#endif


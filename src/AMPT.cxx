#define AMPT_cxx
#include "AMPT.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TChain.h"
#include "TFile.h"
#include <iostream>
using namespace std;

void AMPT::SetData(TChain *chain)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   chain->SetMakeClass(1);

   chain->SetBranchAddress("NEVENT",&cell.NEvent);

   chain->SetBranchAddress("Multi", &cell.multiINUCP);

//   chain->SetBranchAddress("IDINUC", cell.IDINUC);
   chain->SetBranchAddress("XINUCP", cell.fXINUCP);
   chain->SetBranchAddress("YINUCP", cell.fYINUCP);
   chain->SetBranchAddress("ZINUCP", cell.fZINUCP);
   chain->SetBranchAddress("NCOLLP", cell.NCOLLP);

   chain->SetBranchAddress("XINUCT", cell.fXINUCT);
   chain->SetBranchAddress("YINUCT", cell.fYINUCT);
   chain->SetBranchAddress("ZINUCT", cell.fZINUCT);
   chain->SetBranchAddress("NCOLLT", cell.NCOLLT);

   chain->SetBranchAddress("IDHZPC", cell.IDHZPC);
   chain->SetBranchAddress("PxHZPC", cell.fPxHZPC);
   chain->SetBranchAddress("PyHZPC", cell.fPyHZPC);
   chain->SetBranchAddress("PzHZPC", cell.fPzHZPC);
   chain->SetBranchAddress("EnergyHZPC", cell.fEnergyHZPC);
   chain->SetBranchAddress("XHZPC", cell.fXHZPC);
   chain->SetBranchAddress("YHZPC", cell.fYHZPC);
   chain->SetBranchAddress("ZHZPC", cell.fZHZPC);
   chain->SetBranchAddress("TimeHZPC", cell.fTHZPC);

   chain->SetBranchAddress("IDZPC", cell.IDZPC);
   chain->SetBranchAddress("PxZPC", cell.fPxZPC);
   chain->SetBranchAddress("PyZPC", cell.fPyZPC);
   chain->SetBranchAddress("PzZPC", cell.fPzZPC);
   chain->SetBranchAddress("EnergyZPC", cell.fEnergyZPC);
   chain->SetBranchAddress("XZPC", cell.fXZPC);
   chain->SetBranchAddress("YZPC", cell.fYZPC);
   chain->SetBranchAddress("ZZPC", cell.fZZPC);
   chain->SetBranchAddress("TimeZPC", cell.fTZPC);


   chain->SetBranchAddress("IDBART", cell.IDBART);
   chain->SetBranchAddress("PxBART", cell.fPxBART);
   chain->SetBranchAddress("PyBART", cell.fPyBART);
   chain->SetBranchAddress("PzBART", cell.fPzBART);
   chain->SetBranchAddress("EnergyBART", cell.fEnergyBART);
   chain->SetBranchAddress("XBART", cell.fXBART);
   chain->SetBranchAddress("YBART", cell.fYBART);
   chain->SetBranchAddress("ZBART", cell.fZBART);
   chain->SetBranchAddress("TimeBART", cell.fTBART);

   chain->SetBranchAddress("IDAART", cell.IDAART);
   chain->SetBranchAddress("PxAART", cell.fPxAART);
   chain->SetBranchAddress("PyAART", cell.fPyAART);
   chain->SetBranchAddress("PzAART", cell.fPzAART);
   chain->SetBranchAddress("EnergyAART", cell.fEnergyAART);
   chain->SetBranchAddress("XAART", cell.fXAART);
   chain->SetBranchAddress("YAART", cell.fYAART);
   chain->SetBranchAddress("ZAART", cell.fZAART);
   chain->SetBranchAddress("TimeAART", cell.fTAART);
}

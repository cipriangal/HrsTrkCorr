// Compiled code example.  Compiles against libtrkcorr.so

#include "HrsTrkCorr.h"
#include "TRandom.h"
#include <iostream>

using namespace std;

int main(int argc, char **argv) {

   Double_t xsign = -1;  // choice of sign; see below
  
  // instantiate for either Left or Right HRS
   HrsTrkCorr *trkcorr = new HrsTrkCorr(HrsTrkCorr::kLeft);  //or kRight

   trkcorr->Init();   // must do this once in the life of the object

   TRandom *rnd = new TRandom();  

   for (Int_t ievt=0; ievt<500; ievt++) {

     // the world's most trivial monte carlo
     Double_t xbeam = rnd->Uniform(-4,4);        // millimeters
     Double_t tg_ph = rnd->Uniform(-0.02,0.02);  // horizontal tan(angle)
     Double_t tg_th = rnd->Uniform(-0.03,0.03);  // vertical tan(angle)

     trkcorr->Load(xbeam, tg_ph, tg_th);     // must load the class

// which sign ?
// To correct the tracks and bring them closer to sieve holes: xsign = -1
// To move the "ideal" tracks closer to misreconstructed tracks: xsign = +1
     
     tg_ph = tg_ph + xsign * trkcorr->GetDeltaTgPh();
     tg_th = tg_th + xsign * trkcorr->GetDeltaTgTh();
     
     cout << "TRANSPORT angles horizontal tan_phi = "<<tg_ph<<"  vertical tan_theta = "<<tg_th<<endl;
     
     
   }
   return 1;
}




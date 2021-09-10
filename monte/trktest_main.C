// Compiled code example.  Compiles against libtrkcorr.so

#include "../HrsTrkCorr.h"
#include "TRandom.h"
#include <iostream>
#include "TCanvas.h"
#include "TH2F.h"

using namespace std;

int main() {
  //  int main(int argc, char **argv) {

   Double_t xsign = +1;  // choice of sign; see below
  
  // instantiate for either Left or Right HRS
   HrsTrkCorr *trkcorr = new HrsTrkCorr(HrsTrkCorr::kLeft);  //or kRight

   trkcorr->Init();   // must do this once in the life of the object
   trkcorr->SetRtolerance(100);//set how far we can look for holes

   TRandom *rnd = new TRandom();  

//    for (Int_t ievt=0; ievt<500; ievt++) {

//      // the world's most trivial monte carlo
//      Double_t xbeam = rnd->Uniform(-4,4);        // millimeters
//      Double_t tg_ph = rnd->Uniform(-0.02,0.02);  // horizontal tan(angle)
//      Double_t tg_th = rnd->Uniform(-0.03,0.03);  // vertical tan(angle)

//      trkcorr->Load(xbeam, tg_ph, tg_th);     // must load the class

// // which sign ?
// // To correct the tracks and bring them closer to sieve holes: xsign = -1
// // To move the "ideal" tracks closer to misreconstructed tracks: xsign = +1
     
//      cout << "TRANSPORT angles horizontal tan_phi = "<<tg_ph<<"  vertical tan_theta = "<<tg_th<<endl;

//      tg_ph = tg_ph + xsign * trkcorr->GetDeltaTgPh();
//      tg_th = tg_th + xsign * trkcorr->GetDeltaTgTh();
     
//      cout << "TRANSPORT angles horizontal tan_phi = "<<tg_ph<<"  vertical tan_theta = "<<tg_th<<endl;
//      cout << endl;
     
//    }
   const int mxnphi = 100;
   const int mxnth = 100;
   const double phirange = 0.03;
   const double thrange = 0.06;
   Double_t ph_tg;
   Double_t th_tg;
   Double_t phibin = 2*phirange / mxnphi;
   Double_t thbin = 2*thrange / mxnth;   
   TH2F *hdisph[3];
   TH2F *hdisth[3];
   hdisph[0] = new TH2F("hdisph1","Ph shift, X=-3",mxnphi,-1*phirange,phirange,mxnth,-1*thrange,thrange);
   hdisph[1] = new TH2F("hdisph2","Ph shift, X=0",mxnphi,-1*phirange,phirange,mxnth,-1*thrange,thrange);
   hdisph[2] = new TH2F("hdisph3","Ph shift, X=3",mxnphi,-1*phirange,phirange,mxnth,-1*thrange,thrange);
   hdisth[0] = new TH2F("hdisth1","Th shift, x=-3",mxnphi,-1*phirange,phirange,mxnth,-1*thrange,thrange);
   hdisth[1] = new TH2F("hdisth2","Th shift, x=0",mxnphi,-1*phirange,phirange,mxnth,-1*thrange,thrange);
   hdisth[2] = new TH2F("hdisth3","Th shift, x=3",mxnphi,-1*phirange,phirange,mxnth,-1*thrange,thrange);
   
   Double_t disph, disth;
   for (Int_t ix=0; ix<3;ix++) {
     Double_t xb = -3 + 3*ix;
     for (Int_t iph=0; iph<mxnphi;iph++) {
       Double_t ph_tg = -1*phirange + (iph+0.5)*phibin;
       for (Int_t ith=0; ith < mxnth; ith++) {
	 Double_t th_tg = -1*thrange + (ith+0.5)*thbin;
	 
	 trkcorr->Load(xb, ph_tg, th_tg,4,0);     // must load the class
	 // if(abs(ph_tg)<0.01 && abs(th_tg)<0.01)
	 //   trkcorr->Load(xb, ph_tg, th_tg,4,1);     // must load the class
	 // else
	 //   trkcorr->Load(xb, ph_tg, th_tg,4,0);     // must load the class
	 
	 hdisph[ix]->Fill(ph_tg,th_tg,xsign * trkcorr->GetDeltaTgPh());
	 hdisth[ix]->Fill(ph_tg,th_tg,xsign * trkcorr->GetDeltaTgTh());
	 
       }  
     }
   }
   gStyle->SetOptStat(0);
   
   TCanvas *c1 = new TCanvas("c1","c1",1200,800);
   c1->Draw();
   c1->Divide(3,2);
   c1->cd(1); 
   gPad->SetLeftMargin(0.15);
   gPad->SetRightMargin(0.25);
   hdisph[0]->Draw("zcol");
   c1->cd(2);
   gPad->SetLeftMargin(0.15);
   gPad->SetRightMargin(0.25);
   hdisph[1]->Draw("zcol");
   c1->cd(3);
   gPad->SetLeftMargin(0.15);
   gPad->SetRightMargin(0.25);
   hdisph[2]->Draw("zcol");

   // TCanvas *c2 = new TCanvas("c2","c2",1000,400);
   // c2->Draw();
   // c2->Divide(3,1);
   c1->cd(4);
   gPad->SetLeftMargin(0.15);
   gPad->SetRightMargin(0.25);
   hdisth[0]->Draw("zcol");
   c1->cd(5);
   gPad->SetLeftMargin(0.15);
   gPad->SetRightMargin(0.25);
   hdisth[1]->Draw("zcol");
   c1->cd(6);
   gPad->SetLeftMargin(0.15);
   gPad->SetRightMargin(0.25);
   hdisth[2]->Draw("zcol");

   
   return 1;
}




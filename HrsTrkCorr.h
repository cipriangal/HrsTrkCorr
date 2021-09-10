#ifndef __HRSTRKCORR_H
#define __HRSTRKCORR_H

#define MAXCOL  15
#define MAXROW  10
#define MAXXBEAM 3

#include "Rtypes.h"
#include <limits>
#include <stdexcept>

class HrsTrkCorr {

  public:

     enum hrstype_t { kLeft, kRight };   

     Double_t fNotFound, fNotHole, fNotResid;
     Int_t BadIdx;

     HrsTrkCorr( hrstype_t which_hrs );
     ~HrsTrkCorr();

// Must initialize at the start of analysis (once in life of the object).  
    Int_t Init();

    void Print(std::string filename="none");
 
    //set number of sigma for finding neighbour holes
    void SetRtolerance(double val){fRTol=val;}

// Must load TRANSPORT angles each event.  First arg is X_beam (mm).
// 2nd arg tg_ph (horizontal), 3rd arg is tg_th (vertical)

    void Load(Double_t x_beam, Double_t tg_ph, Double_t tg_th, int nNeighbor=1, int debug=0);

// Get sieve hole locations
     Double_t GetHoleTgPh(Int_t ix, Int_t col, Int_t row);
     Double_t GetHoleTgTh(Int_t ix, Int_t col, Int_t row);
  
// X of beam corresponding to residual file
     Double_t GetX(Int_t ix);
  
// Here are corrected TRANSPORT coordinate angles

     Double_t GetCorrTgPh();
     Double_t GetCorrTgTh();

  // The deltas
     Double_t GetDeltaTgPh();
     Double_t GetDeltaTgTh();

  // Residuals at each hole; return of fNotResid means undefined
     Double_t GetResidTgTh(Int_t ixbeam, Int_t col, Int_t row); // theta
     Double_t GetResidTgPh(Int_t ixbeam, Int_t col, Int_t row); // phi
  
  private:
  
     hrstype_t fSpect;
     Double_t fWidTh, fWidPh,fRTol,fXTol;
     Double_t newtgph, newtgth, deltaph, deltath;
     Double_t *tgth_hole, *tgph_hole;
     Int_t *col_has_hole, *hole_found;
     std::string file_holes_m3, file_holes_p3, file_holes_0;
     std::string file_resid_m3, file_resid_p3, file_resid_0;
     FILE *fd;
     Double_t *fXbeam;
     Double_t *thresid, *phresid;
     Int_t CheckIdx(Int_t ix, Int_t col, Int_t row);    
     Int_t did_init;
  
};
#endif//__HRSTRKCORR_H

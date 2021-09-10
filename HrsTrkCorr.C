#include "HrsTrkCorr.h"

#include "TMath.h"
#include <limits>
#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <cstdio>
#include <vector>

//using namespace std;

HrsTrkCorr::HrsTrkCorr(hrstype_t which_spectrometer): did_init(0), fSpect(which_spectrometer) {

  Int_t idx;
  fNotFound = -0.5;  // not found is -999 mrad = -0.999 rad in residual files
  fNotHole = -999;
  fNotResid = -9999;
  BadIdx = -1;
  fWidTh = 0.004;  // assumed "sigma" in tg_th
  fWidPh = 0.002;  // assumed "sigma" in tg_ph
  fRTol = 5;  // this is 5 sigma
  fXTol = 1.0; // X tolerance in mm.
  fXbeam = new Double_t[MAXXBEAM];
  col_has_hole = new Int_t[MAXXBEAM*MAXCOL];
  hole_found = new Int_t [MAXXBEAM*MAXCOL*MAXROW];
  tgth_hole = new Double_t[MAXXBEAM*MAXCOL*MAXROW];
  tgph_hole = new Double_t[MAXXBEAM*MAXCOL*MAXROW];
  thresid = new Double_t[MAXXBEAM*MAXCOL*MAXROW];
  phresid = new Double_t[MAXXBEAM*MAXCOL*MAXROW];

// Initialize the variables.
  for (Int_t ix=0; ix<MAXXBEAM; ix++) fXbeam[ix]=0;
  // don't change these !
  if (MAXXBEAM >= 3) {
    fXbeam[0] =  0;
    fXbeam[1] =  3.8; // mm
    fXbeam[2] = -3.2;
  }

  for (Int_t ix=0; ix<MAXXBEAM; ix++) {
    for (Int_t col=0; col<MAXCOL; col++) {
      idx = ix*MAXCOL + col;
      col_has_hole[idx]=0;
      
      for (Int_t row=0; row<MAXROW; row++) {
        idx = ix*MAXCOL*MAXROW + col*MAXROW + row;
        hole_found[idx]=0;
	tgth_hole[idx]=0;
	tgph_hole[idx]=0;
	thresid[idx]=fNotResid;
	phresid[idx]=fNotResid;
      }
    }
  }

  deltaph=0;
  deltath=0;
  newtgph=0;
  newtgth=0;  

// hole files = where sieve holes are, and
// residuals = measured track residuals
// X_beam is actually horizontal -- so it would be Y in transport.
  if (fSpect == kLeft) {
// X_beam = -3.2 mm
      file_holes_m3 = "./holefiles/holes_LeftHRS_minus3.dat"; // run 2241
// X_beam = +3.8 mm
      file_holes_p3 = "./holefiles/holes_LeftHRS_plus3.dat";  // run 2239
// X_beam = 0
      file_holes_0 = "./holefiles/holes_LeftHRS_0.dat";       // run 2240
      file_resid_m3 = "./holefiles/resid_LeftHRS_minus3.dat"; 
      file_resid_p3 = "./holefiles/resid_LeftHRS_plus3.dat";  
      file_resid_0 = "./holefiles/resid_LeftHRS_0.dat";       
  } else if (fSpect == kRight) {
      file_holes_m3 = "./holefiles/holes_RightHRS_minus3.dat"; // run 21365
      file_holes_p3 = "./holefiles/holes_RightHRS_plus3.dat";  // run 21363
      file_holes_0 = "./holefiles/holes_RightHRS_0.dat";       // run 21364  
      file_resid_m3 = "./holefiles/resid_RightHRS_minus3.dat";
      file_resid_p3 = "./holefiles/resid_RightHRS_plus3.dat";
      file_resid_0 = "./holefiles/resid_RightHRS_0.dat";
  } else {
      std::cout << "HrsTrkCorr::Error: must specify kLeft or kRight in C'tor"<<std::endl;
      std::cout << "try again"<<std::endl;
      exit(0);
  }  
}

HrsTrkCorr::~HrsTrkCorr(){
  delete [] col_has_hole;
  delete [] hole_found;
  delete [] tgth_hole;
  delete [] tgph_hole;
  delete [] thresid;
  delete [] phresid;
}

Int_t HrsTrkCorr::Init() {
// Initialization is here.
// Note, this is separated from construction, per usual C++ practice.
// However, you MUST initialize once in the life of this object.
  FILE *fd;
  Int_t ixbeam,col,row,idx,istrl;
  Float_t ftgth,ftgph;
  Float_t fx0,fx1,fx2,fx3,fx4,fx5;
  std::string filename;
  char strin[80];
  Int_t debug=0;
  if (did_init) return 0;

// read in arrays

  for (ixbeam=0; ixbeam<3; ixbeam++) {

// first read in the sieve holes locations in "angle space"
    switch(ixbeam) {
    case 0:
      filename = file_holes_0;
      break;
    case 1:
      filename = file_holes_m3;
      break;
    case 2:
      filename = file_holes_p3;
      break;
    }
    fd = fopen(filename.c_str(),"r");
    if (fd==NULL) {
      std::cout << "HrsTrkCorr::Error: file "<<filename<<" does not exist"<<std::endl;
      return -1;
    }
    if(debug) std::cout << "Reading file "<<filename<<std::endl;
    while (fgets(strin,80,fd)) {
      if(debug) std::cout << "strin "<<strin<<std::endl;
      istrl = 0;
      if( strstr(strin,"x,y")!=NULL) istrl=1;  // istrl==0 is a line skip
      sscanf(strin,"x,y =  %d %d %f %f ",&row,&col,&ftgth,&ftgph);
      if(debug) std::cout << "row "<<row<<"  col "<<col<<"   tgth "<<ftgth<<"   tgph "<<ftgph<<std::endl;
      if(istrl) {
        if (col>=0 && col<MAXCOL) {
  	   col_has_hole[ixbeam*MAXCOL + col] = 1;
	   if (row>=0 && row<MAXROW) {
	      idx = ixbeam*MAXCOL*MAXROW + col*MAXROW + row;
              hole_found[idx] = 1;
	      tgth_hole[idx] = ftgth;
              tgph_hole[idx] = ftgph;
	   }
	}
      }
    }
    fclose(fd);

// next read in the residuals
    switch(ixbeam) {
    case 0:
      filename = file_resid_0;
      break;
    case 1:
      filename = file_resid_m3;
      break;
    case 2:
      filename = file_resid_p3;
      break;
    }
    fd = fopen(filename.c_str(),"r");
    if (fd==NULL) {
      std::cout << "HrsTrkCorr::Error: file "<<filename<<" does not exist"<<std::endl;
      return -1;
    }
    if(debug) std::cout << "Reading file "<<filename<<std::endl;
    while (fgets(strin,80,fd)) {
      if( strstr(strin,"col=")!=NULL) {
  	  sscanf(strin,"col=%d",&col);
      }
      if( strstr(strin,"tg_th=")!=NULL) {
	sscanf(strin,"tg_th= %f %f %f %f %f %f ",&fx0,&fx1,&fx2,&fx3,&fx4,&fx5);
	// rows 0 to 5
	idx = ixbeam*MAXCOL*MAXROW + col*MAXROW;
        thresid[idx++] = fx0/1000.;  // convert to radians
        thresid[idx++] = fx1/1000.;
        thresid[idx++] = fx2/1000.;
        thresid[idx++] = fx3/1000.;
        thresid[idx++] = fx4/1000.;
        thresid[idx]   = fx5/1000.;
      }
      if( strstr(strin,"tg_ph=")!=NULL) {
	sscanf(strin,"tg_ph= %f %f %f %f %f %f ",&fx0,&fx1,&fx2,&fx3,&fx4,&fx5);
	// rows 0 to 5
	idx = ixbeam*MAXCOL*MAXROW + col*MAXROW;
        phresid[idx++] = fx0/1000.;  // convert to radians
        phresid[idx++] = fx1/1000.;
        phresid[idx++] = fx2/1000.;
        phresid[idx++] = fx3/1000.;
        phresid[idx++] = fx4/1000.;
        phresid[idx]   = fx5/1000.;
      }
      if(debug) std::cout << "strin "<<strin<<std::endl;
    }
    
    fclose(fd);

  }
  
  did_init = 1;

  return did_init;
}

void HrsTrkCorr::Load(Double_t x_beam, Double_t tg_ph, Double_t tg_th, int nNeighbor, int debug) {

// x_beam = horizontal beam position in millimeters.  (this would be called Y in TRANSPORT).
// tg_ph and tg_th are horizontal and vertical TRANSPORT angles (tangents)
// note, the vertical beam position is assumed to be irrelevant.
  
  Double_t thdiff,phdiff;
  Double_t xbdiff0, xbdiff2, xbsum, xfract0, xfract2; 
  Int_t col,row,idx,ix,ix2;
  Double_t rrad,rrmin;
  Int_t locdeb=0;
  if (!did_init) {
    std::cout << "HrsTrkCorr::Load::ERROR:  did not initialize the class"<<std::endl;
    exit(0);
  }

  if(locdeb) std::cout << "------------- enter HrsTrkCorr::Load --------------- "<<std::endl;

  // Assumptions !
  // within +/- fXTol of X=0 use the residuals for x=0
  // Assuming fXBeam[1] = +3.8  and fXBeam[2] = -3.2, so it's not very general.
  // If X_beam within fXTol of fXBeam[1] or X_beam > fXBeam[1]+fXTol
  // use the residuals for fXBeam[1].
  // Similar logic for fXBeam[2].
  // In between, we interpolate in x_beam
  
  ix  = -1;
  ix2 = -1;
  if(x_beam > -1.0*fXTol && x_beam < fXTol) {
    ix=0;
  } else {
    if (fXbeam[1]>fXbeam[2]) {  // fXbeam was defined in constructor
      if (x_beam < fXbeam[2]+fXTol ) ix=2;
      if (x_beam > fXbeam[1]-fXTol ) ix=1;
      if(ix==-1) {  // beam between the positions, will interpolate
        if (x_beam < 0) {
	  ix = 2;    // normally -3.2 mm
	  ix2 = 0;   // normally 0
	} else {
	  ix = 0;    // normally 0
	  ix2 = 1;   // normally +3.8 mm
	}
      }
    } else {
      printf("HrsTrkCorr: fXbeam is screwed up ! \n"); // logically impossible,
              // unless someone plays with the code
      exit(0);
    }    
  }
  if (ix == -1) {
    printf("HrsTrkCorr:: this should be impossible \n");
    exit(0);
  }
  if(locdeb) printf("ix logic %f %d %d \n",x_beam,ix,ix2);

  deltaph=0;
  deltath=0;
  newtgph = tg_ph - deltaph;
  newtgth = tg_th - deltath;

  Int_t colmin=-1;
  Int_t rowmin=-1;
  rrmin = fRTol;
  std::vector<double> dr,dth,dph,vCol,vRow;

// Minimize radius in "angle space" to find hole for X=0
  for (col=0; col<MAXCOL; col++) {
    if (col_has_hole[col]==1) {
        for (row=0; row<MAXROW; row++) {
	   idx =  col*MAXROW + row; // for ixbeam=0
           if ( hole_found[idx]==0 ) continue;
	   thdiff = tg_th - tgth_hole[idx];  
	   phdiff = tg_ph - tgph_hole[idx];
           rrad = TMath::Sqrt( ((thdiff*thdiff)/(fWidTh*fWidTh)) + ((phdiff*phdiff)/(fWidPh*fWidPh)) );
	   
	   vCol.push_back(col);
	   vRow.push_back(row);
	   dr.push_back(rrad);
	   dth.push_back(thdiff);
	   dph.push_back(phdiff);

	   if (rrad < rrmin) {
             rrmin = rrad;
 	     rowmin = row;
             colmin = col;
             if(locdeb) std::cout << " *********************  min "<<rrad<<"   "<<tg_ph<<"   "<<tg_th<<"   "<<tgph_hole[idx]<<"  "<<tgth_hole[idx]<<"  "<<colmin<<"  "<<rowmin<<std::endl;
	   }
	}
      }
  }

  double deltaTh(0),deltaPh(0),sumThWght(0),sumPhWght(0);

// found hole for X=0, now find the delta angles
  if(debug)
    std::cout<<"1by hand "<<colmin<<" "<<rowmin<<" "<<tg_th<<" "<<tg_ph<<std::endl;
  
  if (colmin > -1 && rowmin > -1) {
    if(debug)
      std::cout<<"2by hand "<<colmin<<" "<<rowmin<<std::endl;
    for(int inb=0;inb<nNeighbor;inb++){
      int index = std::distance(dr.begin(), std::min_element(dr.begin(),dr.end()));
      if(debug)
	std::cout<<"index "<<index<<" "<<dr[index]<<std::endl;
      colmin = vCol[index];
      rowmin = vRow[index];
      dr[index]=999;//now it's not the smallest

      switch(ix2) {
	
      case -1:
	
        if ((GetResidTgTh(ix,colmin,rowmin) >fNotFound) &&
	    (GetResidTgPh(ix,colmin,rowmin) >fNotFound)) {  
	  
	  deltaph = GetResidTgPh(ix,colmin,rowmin);
          deltath = GetResidTgTh(ix,colmin,rowmin);
	  
	}
        break;
	
      case 0:
      case 1:
	
	// Use distance to holes to weight the residuals
	
        xbdiff0 = x_beam - fXbeam[ix];
        xbdiff2 = fXbeam[ix2] - x_beam;
        xbsum = xbdiff0 + xbdiff2;
        if (xbsum > 0) {
	  xfract0 = xbdiff0/xbsum;
	} else {
	  xfract0 = 0;
	  printf("HrsTrkCorr: <= 0 sum distance ? %d %d %f %f %f %f %f \n",ix,ix2,fXbeam[ix],fXbeam[ix2],x_beam,xbdiff0,xbdiff2);
	  if(locdeb) sleep(5);
	}
        xfract2 = 1.0 - xfract0;
	
	// only if they are all defined
 	if ((GetResidTgTh(ix,colmin,rowmin) >fNotFound) &&
	    (GetResidTgTh(ix2,colmin,rowmin)>fNotFound) &&
	    (GetResidTgPh(ix,colmin,rowmin) >fNotFound) &&
	    (GetResidTgPh(ix2,colmin,rowmin)>fNotFound)) { 
	  
	  deltaph = xfract0*GetResidTgPh(ix,colmin,rowmin) + xfract2*GetResidTgPh(ix2,colmin,rowmin);
	  deltath = xfract0*GetResidTgTh(ix,colmin,rowmin) + xfract2*GetResidTgTh(ix2,colmin,rowmin);
	}
	
	break;
      }//switch
      if(debug){
	std::cout<<"deltaT "<<deltaTh<<" "<<deltath<<" "<<dth[index]<<std::endl;
	std::cout<<"deltaP "<<deltaPh<<" "<<deltaph<<" "<<dph[index]<<std::endl;
      }
      deltaTh += deltath*dr[index];
      deltaPh += deltaph*dr[index];

      sumThWght += dr[index];
      sumPhWght += dr[index];
      if(debug){
	std::cout<<inb<<" det mins "<<colmin<<" "<<rowmin<<std::endl;
	std::cout<<sumThWght<<" wgt "<<sumPhWght<<std::endl<<std::endl;
      }

    }//neighbor loop
  }//if
  
  if(debug)
    std::cout<<sumThWght<<" Final wgt "<<sumPhWght<<std::endl;
  if(sumPhWght!=0 && sumThWght!=0){
    deltaph = deltaPh/sumPhWght;
    deltath = deltaTh/sumThWght;
  }
  if(debug)
    std::cout<<std::scientific<<deltath<<" deltas "<<deltaph<<std::endl;

  newtgph = tg_ph - deltaph;
  newtgth = tg_th - deltath;

  if(debug){
    std::cout<<newtgth<<" newAngles "<<newtgph<<std::endl<<std::endl<<std::endl;
    // if(sumPhWght!=0 && sumThWght!=0 && abs(deltaph)>0.01)
    //   std::cin.ignore();
  }
  if (locdeb) std::cout << "tg_ph "<<tg_ph<<"  "<<deltaph<<"  "<<newtgph<<"  tg_th "<<tg_th<<"  "<<deltath<<"  "<<newtgth<<std::endl;
  
}


Double_t HrsTrkCorr::GetX(Int_t ix) {
   if (ix < 0 || ix >= MAXXBEAM) return fNotFound;
   return fXbeam[ix];
}

Double_t HrsTrkCorr::GetCorrTgPh() {
  return newtgph;
}

Double_t HrsTrkCorr::GetCorrTgTh() {
  return newtgth;
}

Double_t HrsTrkCorr::GetDeltaTgPh() {
  return deltaph;
} 

Double_t HrsTrkCorr::GetDeltaTgTh() {
  return deltath;
}

Double_t HrsTrkCorr::GetHoleTgTh(Int_t ix, Int_t col, Int_t row) {
  Int_t idx = CheckIdx(ix, col, row);
  if (idx == BadIdx) return fNotHole;
  if (!hole_found[idx]) return fNotHole;
  return tgth_hole[idx];
}

Double_t HrsTrkCorr::GetHoleTgPh(Int_t ix, Int_t col, Int_t row) {
  Int_t idx = CheckIdx(ix, col, row);
  if (idx == BadIdx) return fNotHole;
  if (!hole_found[idx]) return fNotHole;
  return tgph_hole[idx];
}

Double_t HrsTrkCorr::GetResidTgTh(Int_t ix, Int_t col, Int_t row) {
  Int_t idx = CheckIdx(ix, col, row);
  if (idx == BadIdx) return fNotResid;
  if (!hole_found[idx]) return fNotResid;
  return thresid[idx];
}

Double_t HrsTrkCorr::GetResidTgPh(Int_t ix, Int_t col, Int_t row) {
  Int_t idx = CheckIdx(ix, col, row);
  if (idx == BadIdx) return fNotResid;
  if (!hole_found[idx]) return fNotResid;
  return phresid[idx];
}



void HrsTrkCorr::Print(std::string file) {

  char strout[80];
  if (!did_init) {
    std::cout << "HrsTrkCorr::Print::ERROR:  did not initialize the class"<<std::endl;
  }
  Int_t typeprint;  // print to std::cout (if 0), or to a file (if specified)
  if( strstr(file.c_str(),"none")!=NULL) {
    typeprint=0;
  } else {
    typeprint=1;
  }
  fd = fopen(file.c_str(),"w");
  if (fd==NULL) {
      std::cout << "HrsTrkCorr::Print "<<file<<" cannot open"<<std::endl;

  }
  Int_t ixbeam, col, row, idx;
  for (ixbeam=0; ixbeam<MAXXBEAM; ixbeam++) {
    sprintf(strout,"\n\n===== Beam location X[%d] = %f\n",ixbeam,fXbeam[ixbeam]);
    if (typeprint) {
       fputs(strout,fd);
    } else {
      std::cout << strout;
    }
    for (col=0; col<MAXCOL; col++) {
      if(col_has_hole[ixbeam*MAXCOL+col]) {
        sprintf(strout,
	"\nColumn %d  theta_hole  ; theta resid ; phi hole ; phi resid\n",col);
        if (typeprint) {
           fputs(strout,fd);
        } else {
           std::cout << strout;
        }
  	for (row=0; row<MAXROW; row++) {
            idx = ixbeam*MAXCOL*MAXROW + col*MAXROW + row;           
            if(hole_found[idx]) {
	      sprintf(strout,"%6.4f  %8.5f  %6.4f  %8.5f\n",
		    tgth_hole[idx],thresid[idx],tgph_hole[idx],phresid[idx]);
            if (typeprint) {
                fputs(strout,fd);
            } else {
                std::cout << strout;
            }
 
	    }
	}
      }
    }
  }
}

Int_t HrsTrkCorr::CheckIdx(Int_t ix, Int_t col, Int_t row) {
  if ( !did_init ) return BadIdx;
  Int_t idx = ix*MAXCOL*MAXROW + col*MAXROW + row;
  if (idx < 0 || idx >= MAXXBEAM*MAXCOL*MAXROW) return BadIdx;
  return idx;
}
  


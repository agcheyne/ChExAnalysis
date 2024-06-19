#include <TROOT.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TVectorD.h>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

// Function to display progress bar
void displayProgressBar(long long current, long long total) {
    int width = 50; // Width of the progress bar
    float progress = (float)current / total;
    int pos = width * progress;

    std::cout << "[";
    for (int i = 0; i < width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}


//-----------------------------------------------------------------------------------------------------------------------------
// sclose and zclose calculation:
//-----------------------------------------------------------------------------------------------------------------------------

//this section is for the calculation of the closest approach of two tracks
//the two tracks are defined by their coordinates and slopes
//have slides that explain the derivation of the matrix equation (Andrew C)
void calc_sclose_zclose( TVector3 Track1_Coord, TVector3 Track2_Coord, TVector3 Track1_Slope, TVector3 Track2_Slope, double &sclose, double &zclose ){

  double x1 = Track1_Coord.X() - Track1_Coord.Z()*Track1_Slope.X();
  double y1 = Track1_Coord.Y() - Track1_Coord.Z()*Track1_Slope.Y();
  double x2 = Track2_Coord.X() - Track2_Coord.Z()*Track2_Slope.X();
  double y2 = Track2_Coord.Y() - Track2_Coord.Z()*Track2_Slope.Y();

  double xp1 = Track1_Slope.X();
  double yp1 = Track1_Slope.Y();
  double xp2 = Track2_Slope.X();
  double yp2 = Track2_Slope.Y();

  TMatrixD Mclose(2,2);
  TVectorD bclose(2);

  Mclose(0,0) = 1.0 + pow(xp1,2) + pow(yp1,2);
  Mclose(0,1) = -(1.0 + xp1*xp2 + yp1*yp2);
  Mclose(1,0) = Mclose(0,1);
  Mclose(1,1) = 1.0 + pow(xp2,2) + pow(yp2,2);

  bclose(0) = xp1*(x2-x1) + yp1*(y2-y1);
  bclose(1) = xp2*(x1-x2) + yp2*(y1-y2);

  TVectorD zClose = Mclose.Invert() * bclose;

  double z1 = zClose(0);
  double z2 = zClose(1);

  double sClose2 = pow( x1 + xp1*z1 - (x2 + xp2*z2), 2 ) + pow( y1 + yp1*z1 - (y2 + yp2*z2), 2 ) + pow( z1-z2, 2 );

  sclose = sqrt(sClose2);
  zclose = 0.5*(z1 + z2 );
}

//-----------------------------------------------------------------------------------------------------------------------------
// conetest
//-----------------------------------------------------------------------------------------------------------------------------
// this section is for the calculation of the intersection of a track with a box
// the purpose is to check if the track is inside the box
// the track is defined by its coordinates and slopes
// the box is defined by its center, half-lengths, and orientation
// the orientation is defined by the angle theta
// the box is centered at (xcenter, ycenter)
// the box has half-lengths Lx and Ly
// the box is oriented at an angle theta with respect to the x-axis

bool conetest( TVector3 Track1_Coord, TVector3 Track1_Slope, double theta, double zclose, double zback, double Lx=2.0, double Ly=0.6, double xcenter=0.0, double ycenter=0.0 ){
  double xfp, yfp, xpfp, ypfp;

  xfp = Track1_Coord.X() - Track1_Coord.Z() * Track1_Slope.X();
  yfp = Track1_Coord.Y() - Track1_Coord.Z() * Track1_Slope.Y();
  xpfp = Track1_Slope.X();
  ypfp = Track1_Slope.Y();

  double xclose = xfp + xpfp*zclose; // x at the closest approach
  double yclose = yfp + ypfp*zclose; // y at the closest approach

  double xpplus = (xpfp + tan(theta))/(1.-xpfp*tan(theta));
  double xpminus = (xpfp - tan(theta))/(1.+xpfp*tan(theta));
  double ypplus = (ypfp + tan(theta))/(1.-ypfp*tan(theta));
  double ypminus = (ypfp - tan(theta))/(1.+ypfp*tan(theta));

  double xmax = xclose + xpplus * (zback - zclose);
  double xmin = xclose + xpminus * (zback - zclose);
  double ymax = yclose + ypplus * (zback - zclose);
  double ymin = yclose + ypminus * (zback - zclose);

  return ( fabs( xmax - xcenter ) <= Lx/2.0 && fabs( xmin - xcenter ) <= Lx/2.0 && fabs( ymax - ycenter ) <= Ly/2.0 && fabs( ymin - ycenter ) <= Ly/2.0 );
  
}

//-----------------------------------------------------------------------------------------------------------------------------

void GEnRPhcal( Int_t run_no = 9000 ) {

  TChain* C = new TChain("TSkim");
  
//input for the skimmed files
  if( run_no == 9000 ) {
//     C->Add("OUT_DIR/skim_genrp_3000.root");
//     C->Add("OUT_DIR/skim_genrp_4000.root");
    C->Add("OUT_DIR/skim_genrp_5000.root");
    C->Add("OUT_DIR/skim_genrp_6000.root");
    C->Add("OUT_DIR/skim_genrp_7000.root");
    C->Add("OUT_DIR/skim_genrp_8000.root");
    C->Add("OUT_DIR/skim_genrp_9000.root");
  }
  else 
    C->Add(Form("OUT_DIR/skim_genrp_%d.root", run_no));
 //default run number
  TTreeReader Tl(C);

// defining the variables
  TTreeReaderValue<double> scalhel_hel(Tl, "scalhel_hel"); // helicity 

  TTreeReaderValue<double> bb_tr_n(Tl, "bb_tr_n"); // number of tracks

  TTreeReaderValue<double> bb_ps_e(Tl, "bb_ps_e"); // BigBite preshower energy
  TTreeReaderValue<double> bb_sh_e(Tl, "bb_sh_e"); // BigBite shower energy
  TTreeReaderValue<double> sbs_hcal_x(Tl, "sbs_hcal_x"); // SBS HCAL x
  TTreeReaderValue<double> sbs_hcal_y(Tl, "sbs_hcal_y"); // SBS HCAL y
  TTreeReaderValue<double> sbs_hcal_atimeblk(Tl, "sbs_hcal_atimeblk"); // SBS HCAL "ADC time of highest energy block in the largest cluster"
  TTreeReaderValue<double> bb_sh_atimeblk(Tl, "bb_sh_atimeblk"); // BigBite shower "ADC time of highest energy block in the largest cluster"
  TTreeReaderValue<double> e_kine_W2(Tl, "e_kine_W2"); // W2 - squared invariant mass of the electron-nucleon system

  TTreeReaderValue<double> bb_tr_px(Tl, "bb_tr_px"); // BigBite track x-momentum
  TTreeReaderValue<double> bb_tr_py(Tl, "bb_tr_py"); // BigBite track y-momentum
  TTreeReaderValue<double> bb_tr_pz(Tl, "bb_tr_pz"); // BigBite track z-momentum
  TTreeReaderValue<double> bb_tr_vx(Tl, "bb_tr_vx"); // BigBite track x-vertex
  TTreeReaderValue<double> bb_tr_vy(Tl, "bb_tr_vy"); // BigBite track y-vertex
  TTreeReaderValue<double> bb_tr_vz(Tl, "bb_tr_vz"); // BigBite track z-vertex
  TTreeReaderValue<double> bb_tr_p(Tl, "bb_tr_p"); // BigBite track momentum
  TTreeReaderValue<double> bb_tr_x(Tl, "bb_tr_x"); // BigBite track x
  TTreeReaderValue<double> bb_tr_y(Tl, "bb_tr_y"); // BigBite track y
  TTreeReaderValue<double> bb_tr_th(Tl, "bb_tr_th"); // BigBite track theta
  TTreeReaderValue<double> bb_tr_ph(Tl, "bb_tr_ph"); // BigBite track phi

  TTreeReaderValue<double> bb_tr_tg_th(Tl, "bb_tr_tg_th"); //  "Tangent of target theta angle"
  TTreeReaderValue<double> bb_tr_tg_ph(Tl, "bb_tr_tg_ph"); // "Tangent of target phi angle"
  TTreeReaderValue<double> bb_tr_tg_y(Tl, "bb_tr_tg_y"); // "Target y coordinate"


  //-----------------------------------------------------------------------------------------------------------------------------

  const Int_t nphibins = 36;

  TFile *outfile = new TFile(Form("hist/hist_genrp_%i.root",run_no+1),"RECREATE");
  
  Double_t Mp        = 0.93827; // Proton mass (GeV)
  Double_t Eb        = 4.3;  // Beam energy
  Double_t th_bb     = 42.5;  // BigBite angle
  Double_t th_sbs    = 24.7;  // SBS angle
  Double_t pcent     = 2.122;  // Central momentum

// defining histograms
  TH1D* hkin_p        = new TH1D("hkin_p","",100,0.25*pcent,1.25*pcent); //
  TH1D* hkin_th       = new TH1D("hkin_th","",100,-0.3,0.3);
  TH1D* hkin_ph       = new TH1D("hkin_ph","",100,-0.1,0.1);
  TH1D* hkin_yt       = new TH1D("hkin_yt","",100,-0.15,0.15);
  TH1D* hkin_W        = new TH1D("hkin_W","",50,0,2);
  TH2D* hkin2d_thp    = new TH2D("hkin2d_thp","",100, th_bb-6,th_bb+6.,100,0.25*pcent,1.25*pcent);
  TH2D *hbbcal2d_pss  = new TH2D("hbbcal2d_pssh","",100,0.,1.6, 100,0.,1.0); 
// histograms with cuts
  TH1D* hkin_pc       = new TH1D("hkin_pc","",100,0.25*pcent,1.25*pcent);
  TH1D* hkin_thc      = new TH1D("hkin_thc","",100,-0.3,0.3);
  TH1D* hkin_phc      = new TH1D("hkin_phc","",100,-0.1,0.1);
  TH1D* hkin_ytc      = new TH1D("hkin_ytc","",100,-0.15,0.15);
  TH1D* hkin_Wc       = new TH1D("hkin_Wc","",50,0,2);
  TH2D *hbbcal2d_pssc = new TH2D("hbbcal2d_psshc","",100,0.,1.6, 100,0.,1.0); 
  
// histograms for HCAL
  TH1D* hhcal_deltax   = new TH1D("hhcal_deltax","",100,-2.5,2.5);
  TH1D* hhcal_deltay   = new TH1D("hhcal_deltay","",100,-2,2.);
  TH2D* hhcal_deltaxy  = new TH2D("hhcal_deltaxy","",50,-2.,2.,50,-2.5,2.5);
// cuts
  TH1D* hhcal_deltaxc  = new TH1D("hhcal_deltaxc","",100,-2.5,2.5);
  TH1D* hhcal_deltaxcc = new TH1D("hhcal_deltaxcc","",100,-2.5,2.5);
  TH1D* hhcal_deltayc  = new TH1D("hhcal_deltayc","",100,-2.,2.);  
  TH2D* hhcal_deltaxyc = new TH2D("hhcal_deltaxyc","",50,-2.5,2.5,50,-2.5,2.5);

// histograms for polarimeters abd cuts
  TH1D* hpolana_deltaxc   = new TH1D("hpolana_deltaxc","",100,-1.0,1.0);
  TH1D* hpolana_deltaxcc  = new TH1D("hpolana_deltaxcc","",100,-1.0,1.0);
  TH1D* hpolana_deltaxccc = new TH1D("hpolana_deltaxccc","",100,-1.0,1.0);
  TH1D* hpolana_deltayc   = new TH1D("hpolana_deltayc","",100,-1.,1.);
  TH2D* hpolana_deltaxyc  = new TH2D("hpolana_deltaxyc","",50,-1.0,1.0,50,-1.0,1.0);

  TH1D* hpolg_thnp_cx   = new TH1D("hpolg_thnp_cx","",100,0,25);
  TH1D* hpolg_thnp_cxc  = new TH1D("hpolg_thnp_cxc","",100,0,25);
  TH1D* hpolg_thnp_cxcc = new TH1D("hpolg_thnp_cxcc","",100,0,25);
  TH1D* hpolg_phnp_cxp  = new TH1D("hpolg_phnp_cxp","",nphibins,-180,180);
  TH1D* hpolg_phnp_cxm  = new TH1D("hpolg_phnp_cxm","",nphibins,-180,180);
  TH1D* hpolg_phnp_cxpc = new TH1D("hpolg_phnp_cxpc","",nphibins,-180,180);
  TH1D* hpolg_phnp_cxmc = new TH1D("hpolg_phnp_cxmc","",nphibins,-180,180);
  TH1D* hpolg_sclnp     = new TH1D("hpolg_sclnp","",100,0, 0.003);
  TH2D* hpolg_zclnp     = new TH2D("hpolg_zclnp","",50, 4.50, 4.95, 50, 0, 25);

  TH1D* hpolp_thnp_cx   = new TH1D("hpolp_thnp_cx","", 100,0,25);
  TH1D* hpolp_thnp_cxc  = new TH1D("hpolp_thnp_cxc","",100,0,25);
  TH1D* hpolp_phnp_cxp  = new TH1D("hpolp_phnp_cxp","",nphibins,-180,180);
  TH1D* hpolp_phnp_cxm  = new TH1D("hpolp_phnp_cxm","",nphibins,-180,180);
  TH1D* hpolp_sclnp     = new TH1D("hpolp_sclnp","",100,0, 0.003);
  TH2D* hpolp_zclnp     = new TH2D("hpolp_zclnp","",50, 4.50, 4.95,50,0,25);

  //-----------------------------------------------------------------------------------------------------------------------------

// defining the distances
  Double_t hcal_dist    = 9.0;    // distance from the target to the HCAL
  Double_t polana_dist  = 4.724;  // distance from the target to the front of the analyzer
  

  TLorentzVector Tp4(0,0,0,Mp);  // Proton 4-momentum (initial)
  TLorentzVector kp4(0,0,Eb,Eb); // Beam 4-momentum (initial) - only z-component is non-zero
  TLorentzVector Qp4, kpp4, Rp4; // Q2, scattered electron, and recoil proton 4-momenta (undetermined)

  th_bb  *= M_PI/180.; // BigBite angle in radians
  th_sbs *= M_PI/180.; // SBS angle in radians

// defining sbs coordinate system
  TVector3 SBS_zaxis( -sin(th_sbs), 0, cos(th_sbs) ); // SBS z-axis - defined as the nucleon beam direction 
  TVector3 SBS_xaxis(0, -1, 0 );  // SBS x-axis - defined as the negative y-axis
  TVector3 SBS_yaxis = SBS_zaxis.Cross(SBS_xaxis).Unit(); // at right angles to the z-axis and x-axis
  
  TVector3 front_vtx(0.,0.,0);
  TVector3 ana_vtx(0.,0.,0);
  TVector3 rear_vtx(0.,0.,0);
  TVector3 in_track(0.,0.,0);
  TVector3 out_track(0.,0.,0);

  //-----------------------------------------------------------------------------------------------------------------------------

  Long64_t ev = 0;
  //get total events
  Long64_t nEntries = C->GetEntries();

  cout << "Total events: " << nEntries << endl;

  while( Tl.Next() ) {
    //get total events
    ev++;
    //display progress bar
    if (ev % 1000 == 0 || ev == nEntries) {
    displayProgressBar(ev, nEntries);
    }
    
    if( *bb_tr_n <= 0 ) continue;  //if there are no tracks, skip the event
    if( *bb_ps_e < 0.05 ) continue; //if the preshower energy is less than 50 MeV, skip the event
    if( fabs(*sbs_hcal_atimeblk-*bb_sh_atimeblk-42) > 5 ) continue; //if the difference between the SBS HCAL and BigBite shower "ADC time of highest energy block in the largest cluster" is greater than 5, skip the event
    
    Double_t p  = *bb_tr_p; // BigBite track momentum
    Double_t px = *bb_tr_px; // BigBite track x-momentum
    Double_t py = *bb_tr_py; // BigBite track y-momentum
    Double_t pz = *bb_tr_pz; // BigBite track z-momentum

    Double_t vx = *bb_tr_vx; // BigBite track x-vertex
    Double_t vy = *bb_tr_vy; // BigBite track y-vertex
    Double_t vz = *bb_tr_vz;  // BigBite track z-vertex
    TVector3 vertex(vx,vy,vz);  // BigBite track vertex
    
    kpp4.SetPxPyPzE(px,py,pz,p); // Scattered electron 4-momentum set to the BigBite track momentum
    Qp4 = kp4 - kpp4; // Q2 4-momentum = Beam 4-momentum - Scattered electron 4-momentum
    Rp4 = Tp4 + Qp4; // Recoil proton 4-momentum = Proton 4-momentum + Q2 4-momentum

    Rp4.RotateY(th_sbs); // Rotate the recoil proton 4-momentum by the SBS angle in the Hall yz-plane
    
    //angles of the recoil proton constructed from Rp4 vector
    Double_t hcal_th = TMath::ATan(Rp4.Px()/Rp4.Pz()); // Recoil proton polar angle (theta)
    Double_t hcal_ph = TMath::ATan(Rp4.Py()/Rp4.Pz()); // Recoil proton azimuthal angle (phi)
    
    //sbs positions
    Double_t hcal_x = *sbs_hcal_x; // detected x-position of the recoil proton in the SBS HCAL
    Double_t pred_x = -hcal_dist * TMath::Sin( hcal_ph ); // predicted x-position of the recoil proton at the SBS HCAL based on calculated angles

    Double_t hcal_y = *sbs_hcal_y; // similar for y
    Double_t pred_y   = hcal_dist * TMath::Sin(hcal_th);


    
    Double_t delta_x = hcal_x - pred_x; // difference between the detected and predicted x-positions
    Double_t delta_y = hcal_y - pred_y; // similar for y


    //difference is plotted to separated neutrons and protons
    //neutrons are expected to have a smaller x (on 0 ideally)
    //protons are expected to be bent due to the magnetic field
    //the difference in y is expected to be small

    //bigbite kinematics histograms
    hkin_p->Fill(p); // BigBite track momentum
    hkin_th->Fill(*bb_tr_tg_th);    // BigBite track theta
    hkin_ph->Fill(*bb_tr_tg_ph);    // BigBite track phi
    hkin_yt->Fill(*bb_tr_tg_y);     // BigBite track y-coordinate
    hkin_W->Fill(*e_kine_W2);       // W2 - squared invariant mass of the electron-nucleon system

    hbbcal2d_pss->Fill( *bb_sh_e/(*bb_tr_p), *bb_ps_e/(*bb_tr_p) ); //(shower energy over momentum) vs (preshower energy over momentum) - we can veto pions with this
    hhcal_deltax->Fill(delta_x);
    hhcal_deltay->Fill(delta_y);
    hhcal_deltaxy->Fill(delta_y,delta_x);
    
    //some cuts!
    if( *bb_ps_e/(*bb_tr_p) < 0.1 ) continue; //if the preshower energy over momentum is less than 0.1, skip the event (the pion rejection cut)
    if( *e_kine_W2 < 0.7 ) continue;  //if the squared invariant mass of the electron-nucleon system is less than 0.7, skip the event
    if( *e_kine_W2 > 1.4 ) continue;  //if the squared invariant mass of the electron-nucleon system is greater than 1.4, skip the event
    //those last two cuts are to select the elastic peak
    
    //fill the cut histograms
    hkin_pc->Fill(p); 
    hkin_thc->Fill(*bb_tr_tg_th);
    hkin_phc->Fill(*bb_tr_tg_ph);
    hkin_ytc->Fill(*bb_tr_tg_y);
    hkin_Wc->Fill(*e_kine_W2);

    //more cuts
    hbbcal2d_pssc->Fill( *bb_sh_e/(*bb_tr_p), *bb_ps_e/(*bb_tr_p) );
    hhcal_deltaxc->Fill(delta_x);
    hhcal_deltayc->Fill(delta_y);
    hhcal_deltaxyc->Fill(delta_y,delta_x);
    if( fabs(delta_y) < 3*0.23 ) //if the difference in y (hcal_y - pred_y) is less than 3*0.23, fill the delta xcc histogram
      hhcal_deltaxcc->Fill(delta_x); //here, xcc should have only charge-exchange events


//-----------------------------------------------------------------------------------------------------------------------------
//GEMless analysis
//now in the SBS coordinate system

//defining vertex into the SBS coordinate system
TVector3 vertex_SBS( vertex.Dot(SBS_xaxis), vertex.Dot(SBS_yaxis), vertex.Dot(SBS_zaxis) );
front_vtx = vertex_SBS;

Double_t ana_x = 99999;
Double_t ana_y = 99999;
Double_t thsc  = 99999;
Double_t phsc  = 99999;


//working out the predicted position of the recoil proton at the front of the analyzer
Double_t pred_anax = vertex_SBS.X() - (polana_dist + vertex_SBS.Z() ) * TMath::Sin( hcal_ph );
Double_t pred_anay = vertex_SBS.Y() + (polana_dist + vertex_SBS.Z() ) * TMath::Sin( hcal_th );
Double_t pred_anaz = vertex_SBS.Z() + polana_dist;
//?? need to figure out why we add z vertex to the distance

ana_vtx.SetXYZ(pred_anax, pred_anay, pred_anaz ) ;
rear_vtx = ana_vtx; //defining the in and out vertex as the same
     
//intrack
//defining the incoming track as the predicted ana vertex - vertex
in_track.SetXYZ(ana_vtx.X() - vertex.X(), ana_vtx.Y() - vertex.Y(), ana_vtx.Z() - vertex.Z() );
out_track.SetXYZ(hcal_x - ana_vtx.X(), hcal_y - ana_vtx.Y(), hcal_dist - ana_vtx.Z() );
      
in_track  = in_track.Unit();
out_track = out_track.Unit();

TVector3 yaxistemp(0,1,0);
TVector3 xaxistemp = yaxistemp.Cross(in_track).Unit();
yaxistemp = in_track.Cross(xaxistemp).Unit();


thsc = 180./M_PI * acos( out_track.Dot(in_track) ); //calculates the angle between the in and out tracks, in track should be the same as the z the SBS coordinate system
phsc = 180./M_PI * TMath::ATan2( out_track.Dot(yaxistemp), out_track.Dot(xaxistemp) ); //calculates the angle between the x and y components of the out track
    
//-----------------------------------------------------------------------------------------------------------------------------
// np -> pn ChEx in passive (steel) analyser, coordinate system = SBS TRANSPORT
//-----------------------------------------------------------------------------------------------------------------------------
Double_t helicity = *scalhel_hel;
if( run_no >= 4000 ) 
  helicity = -1 * helicity;
    
      double sclose,zclose;
      
      calc_sclose_zclose( front_vtx, rear_vtx, in_track, out_track, sclose, zclose ); //calculates the closest approach of the incoming and outgoing tracks
      //bool conetestnp = conetest( front_vtx, in_track, (thsc/57.3), zclose, polrgem_dist );
      bool conetestnp = true; //?? for now, we don't have the GEMs

      hpolg_sclnp->Fill( sclose );  //fill the sclose histogram
      hpolg_zclnp->Fill( zclose , thsc );  //fill the zclose histogram
      hpolg_thnp_cx->Fill ( thsc ); //fill the theta histogram
      

  //if the theta is greater than 3 degrees, the sclose is less than 0.003, the zclose is less than 0.2, and the phi is less than 180, fill the following histograms
    //  if( thsc >= 3.0 && sclose < 0.003  && fabs(zclose-4.72) < 0.2 && fabs(phsc)<180 ) {  
  	//hpolg_thnp_cxc->Fill ( thsc );

  //defining cuts - removed the sclose cut
        if( thsc >= 6.0  && fabs(zclose-4.72) < 0.2 && fabs(phsc)<180 ) {  
  	hpolg_thnp_cxc->Fill ( thsc );
  	
    //if the helicity is -1, fill the following histograms
    if( helicity == -1 ) { 
  	  hpolg_phnp_cxm->Fill ( phsc ); 
  	  if( conetestnp ) {
  	    hpolg_phnp_cxmc->Fill ( phsc );
  	    hpolg_thnp_cxcc->Fill ( thsc );
  	  }
  	}
    //if the helicity is 1, fill the following histograms
  	else if( helicity == 1 ) {
  	  hpolg_phnp_cxp->Fill ( phsc );
  	  if( conetestnp ) {
  	    hpolg_phnp_cxpc->Fill ( phsc );
  	    hpolg_thnp_cxcc->Fill ( thsc );
  	  }
  	}
      }
      
    }
  
  outfile->Write();
  outfile->Close();

  std::cout << std::endl; // to end the progress bar
}

//-----------------------------------------------------------------------------------------------------------------------------



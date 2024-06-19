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

using namespace std;

//-----------------------------------------------------------------------------------------------------------------------------
// sclose and zclose calculation:
//-----------------------------------------------------------------------------------------------------------------------------

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

bool conetest( TVector3 Track1_Coord, TVector3 Track1_Slope, double theta, double zclose, double zback, double Lx=2.0, double Ly=0.6, double xcenter=0.0, double ycenter=0.0 ){
  double xfp, yfp, xpfp, ypfp;

  xfp = Track1_Coord.X() - Track1_Coord.Z() * Track1_Slope.X();
  yfp = Track1_Coord.Y() - Track1_Coord.Z() * Track1_Slope.Y();
  xpfp = Track1_Slope.X();
  ypfp = Track1_Slope.Y();

  double xclose = xfp + xpfp*zclose;
  double yclose = yfp + ypfp*zclose;

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

void GEnRPAnalysis( Int_t run_no = 9000 ) {

  TChain* C = new TChain("TSkim");

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

  TTreeReader Tl(C);

  TTreeReaderValue<double> scalhel_hel(Tl, "scalhel_hel");

  TTreeReaderValue<double> bb_tr_n(Tl, "bb_tr_n");
  TTreeReaderValue<double> sbs_gemCeF_track_ntrack(Tl, "sbs_gemCeF_track_ntrack");  
  TTreeReaderValue<double> sbs_gemCeR_track_ntrack(Tl, "sbs_gemCeR_track_ntrack");  

  TTreeReaderValue<double> bb_ps_e(Tl, "bb_ps_e");
  TTreeReaderValue<double> bb_sh_e(Tl, "bb_sh_e");
  TTreeReaderValue<double> sbs_hcal_x(Tl, "sbs_hcal_x");
  TTreeReaderValue<double> sbs_hcal_y(Tl, "sbs_hcal_y");
  TTreeReaderValue<double> sbs_hcal_atimeblk(Tl, "sbs_hcal_atimeblk");
  TTreeReaderValue<double> bb_sh_atimeblk(Tl, "bb_sh_atimeblk");

  TTreeReaderValue<double> e_kine_W2(Tl, "e_kine_W2");

  TTreeReaderValue<double> bb_tr_px(Tl, "bb_tr_px");
  TTreeReaderValue<double> bb_tr_py(Tl, "bb_tr_py");
  TTreeReaderValue<double> bb_tr_pz(Tl, "bb_tr_pz");
  TTreeReaderValue<double> bb_tr_vx(Tl, "bb_tr_vx");
  TTreeReaderValue<double> bb_tr_vy(Tl, "bb_tr_vy");
  TTreeReaderValue<double> bb_tr_vz(Tl, "bb_tr_vz");
  TTreeReaderValue<double> bb_tr_p(Tl, "bb_tr_p");
  TTreeReaderValue<double> bb_tr_x(Tl, "bb_tr_x");
  TTreeReaderValue<double> bb_tr_y(Tl, "bb_tr_y");
  TTreeReaderValue<double> bb_tr_th(Tl, "bb_tr_th");
  TTreeReaderValue<double> bb_tr_ph(Tl, "bb_tr_ph");

  TTreeReaderValue<double> bb_tr_tg_th(Tl, "bb_tr_tg_th");
  TTreeReaderValue<double> bb_tr_tg_ph(Tl, "bb_tr_tg_ph");
  TTreeReaderValue<double> bb_tr_tg_y(Tl, "bb_tr_tg_y");

  TTreeReaderValue<double> sbs_gemCeR_track_x(Tl, "sbs_gemCeR_track_x");
  TTreeReaderValue<double> sbs_gemCeR_track_y(Tl, "sbs_gemCeR_track_y");
  TTreeReaderValue<double> sbs_gemCeR_track_xp(Tl, "sbs_gemCeR_track_xp");
  TTreeReaderValue<double> sbs_gemCeR_track_yp(Tl, "sbs_gemCeR_track_yp");
  TTreeReaderValue<double> sbs_gemCeR_track_sclose(Tl, "sbs_gemCeR_track_sclose");
  TTreeReaderValue<double> sbs_gemCeR_track_zclose(Tl, "sbs_gemCeR_track_zclose");
  TTreeReaderValue<double> sbs_gemCeR_track_theta(Tl, "sbs_gemCeR_track_theta");
  TTreeReaderValue<double> sbs_gemCeR_track_phi(Tl, "sbs_gemCeR_track_phi");

  //-----------------------------------------------------------------------------------------------------------------------------

  const Int_t nphibins = 36;

  TFile *outfile = new TFile(Form("hist/hist_genrp_%i.root",run_no),"RECREATE");
  
  Double_t Mp        = 0.93827;
  Double_t Eb        = 4.3;  
  Double_t th_bb     = 42.5; 
  Double_t th_sbs    = 24.7; 
  Double_t pcent     = 2.122;  

  TH1D* hkin_p        = new TH1D("hkin_p","",100,0.25*pcent,1.25*pcent);
  TH1D* hkin_th       = new TH1D("hkin_th","",100,-0.3,0.3);
  TH1D* hkin_ph       = new TH1D("hkin_ph","",100,-0.1,0.1);
  TH1D* hkin_yt       = new TH1D("hkin_yt","",100,-0.15,0.15);
  TH1D* hkin_W        = new TH1D("hkin_W","",50,0,2);
  TH2D* hkin2d_thp    = new TH2D("hkin2d_thp","",100, th_bb-6,th_bb+6.,100,0.25*pcent,1.25*pcent);
  TH2D *hbbcal2d_pss  = new TH2D("hbbcal2d_pssh","",100,0.,1.6, 100,0.,1.0); 

  TH1D* hkin_pc       = new TH1D("hkin_pc","",100,0.25*pcent,1.25*pcent);
  TH1D* hkin_thc      = new TH1D("hkin_thc","",100,-0.3,0.3);
  TH1D* hkin_phc      = new TH1D("hkin_phc","",100,-0.1,0.1);
  TH1D* hkin_ytc      = new TH1D("hkin_ytc","",100,-0.15,0.15);
  TH1D* hkin_Wc       = new TH1D("hkin_Wc","",50,0,2);
  TH2D *hbbcal2d_pssc = new TH2D("hbbcal2d_psshc","",100,0.,1.6, 100,0.,1.0); 
  
  TH1D* hhcal_deltax   = new TH1D("hhcal_deltax","",100,-2.5,2.5);
  TH1D* hhcal_deltay   = new TH1D("hhcal_deltay","",100,-2,2.);
  TH2D* hhcal_deltaxy  = new TH2D("hhcal_deltaxy","",50,-2.,2.,50,-2.5,2.5);

  TH1D* hhcal_deltaxc  = new TH1D("hhcal_deltaxc","",100,-2.5,2.5);
  TH1D* hhcal_deltaxcc = new TH1D("hhcal_deltaxcc","",100,-2.5,2.5);
  TH1D* hhcal_deltayc  = new TH1D("hhcal_deltayc","",100,-2.,2.);  
  TH2D* hhcal_deltaxyc = new TH2D("hhcal_deltaxyc","",50,-2.5,2.5,50,-2.5,2.5);

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

  Double_t hcal_dist    = 9.0;    
  Double_t polana_dist  = 4.724;    
  Double_t polfgem_dist = polana_dist - 0.4;     
  Double_t polrgem_dist = polana_dist + 0.4;     

  TLorentzVector Tp4(0,0,0,Mp); 
  TLorentzVector kp4(0,0,Eb,Eb); 
  TLorentzVector Qp4, kpp4, Rp4; 

  th_bb  *= M_PI/180.;
  th_sbs *= M_PI/180.;

  TVector3 SBS_zaxis( -sin(th_sbs), 0, cos(th_sbs) );
  TVector3 SBS_xaxis(0, -1, 0 );
  TVector3 SBS_yaxis = SBS_zaxis.Cross(SBS_xaxis).Unit();
  
  TVector3 front_vtx(0.,0.,0);
  TVector3 ana_vtx(0.,0.,0);
  TVector3 rear_vtx(0.,0.,0);
  TVector3 in_track(0.,0.,0);
  TVector3 out_track(0.,0.,0);

  //-----------------------------------------------------------------------------------------------------------------------------

  Long64_t ev = 0;

  while( Tl.Next() ) {

    if( *bb_tr_n <= 0 ) continue; 
    if( *bb_ps_e < 0.05 ) continue;
    if( fabs(*sbs_hcal_atimeblk-*bb_sh_atimeblk-42) > 5 ) continue;
    
    Double_t p  = *bb_tr_p;
    Double_t px = *bb_tr_px;
    Double_t py = *bb_tr_py;
    Double_t pz = *bb_tr_pz;

    Double_t vx = *bb_tr_vx;
    Double_t vy = *bb_tr_vy;
    Double_t vz = *bb_tr_vz; 
    TVector3 vertex(vx,vy,vz); 
    
    kpp4.SetPxPyPzE(px,py,pz,p);
    Qp4 = kp4 - kpp4;
    Rp4 = Tp4 + Qp4;

    Rp4.RotateY(th_sbs);
    
    Double_t hcal_th = TMath::ATan(Rp4.Px()/Rp4.Pz());
    Double_t hcal_ph = TMath::ATan(Rp4.Py()/Rp4.Pz());
    
    Double_t hcal_x = *sbs_hcal_x;
    Double_t pred_x = -hcal_dist * TMath::Sin( hcal_ph );
    
    Double_t hcal_y = *sbs_hcal_y;
    Double_t pred_y   = hcal_dist * TMath::Sin(hcal_th);
    
    Double_t delta_x = hcal_x - pred_x;
    Double_t delta_y = hcal_y - pred_y;

    hkin_p->Fill(p);
    hkin_th->Fill(*bb_tr_tg_th);
    hkin_ph->Fill(*bb_tr_tg_ph);
    hkin_yt->Fill(*bb_tr_tg_y);
    hkin_W->Fill(*e_kine_W2);

    hbbcal2d_pss->Fill( *bb_sh_e/(*bb_tr_p), *bb_ps_e/(*bb_tr_p) );
    hhcal_deltax->Fill(delta_x);
    hhcal_deltay->Fill(delta_y);
    hhcal_deltaxy->Fill(delta_y,delta_x);
    
    if( *bb_ps_e/(*bb_tr_p) < 0.1 ) continue;
    if( *e_kine_W2 < 0.7 ) continue; 
    if( *e_kine_W2 > 1.4 ) continue;   
    
    hkin_pc->Fill(p);
    hkin_thc->Fill(*bb_tr_tg_th);
    hkin_phc->Fill(*bb_tr_tg_ph);
    hkin_ytc->Fill(*bb_tr_tg_y);
    hkin_Wc->Fill(*e_kine_W2);

    hbbcal2d_pssc->Fill( *bb_sh_e/(*bb_tr_p), *bb_ps_e/(*bb_tr_p) );
    hhcal_deltaxc->Fill(delta_x);
    hhcal_deltayc->Fill(delta_y);
    hhcal_deltaxyc->Fill(delta_y,delta_x);
    if( fabs(delta_y) < 3*0.23 )
      hhcal_deltaxcc->Fill(delta_x);

    //-----------------------------------------------------------------------------------------------------------------------------
    
    bool SBSGEMF = false;
    bool SBSGEMR = false;
    bool GEMdx   = false;
    bool GEMdxP  = false;
    bool HCALn   = false;
    
    TVector3 vertex_SBS(vertex.Dot(SBS_xaxis), vertex.Dot(SBS_yaxis), 
			vertex.Dot(SBS_zaxis) ); 
    
    Double_t pred_anax = vertex_SBS.X() -( polana_dist + vertex_SBS.Z() ) * TMath::Sin( hcal_ph );
    Double_t pred_anay = vertex_SBS.Y() + (polana_dist + vertex_SBS.Z() ) * TMath::Sin( hcal_th );
    Double_t ana_x = 99999;
    Double_t ana_y = 99999;
    Double_t thsc  = 99999;
    Double_t phsc  = 99999;
    
    if( *sbs_gemCeF_track_ntrack > 0 )
      SBSGEMF = true;
    
    if( *sbs_gemCeR_track_ntrack > 0 )
      SBSGEMR = true;
    
    if ( SBSGEMR ) {

      ana_x = *sbs_gemCeR_track_x + (polana_dist - polfgem_dist)*(*sbs_gemCeR_track_xp); 
      ana_y = *sbs_gemCeR_track_y + (polana_dist - polfgem_dist)*(*sbs_gemCeR_track_yp); 
      
      Double_t delta_anax = ana_x - pred_anax;
      Double_t delta_anay = ana_y - pred_anay;
      
      hpolana_deltaxc->Fill(delta_anax);
      hpolana_deltayc->Fill(delta_anay);
      hpolana_deltaxyc->Fill(delta_anay,delta_anax);

      if( fabs(delta_anay) < 2.0*0.95 ) {
  	hpolana_deltaxcc->Fill(delta_anax);
	if( delta_anax > -0.05 && delta_anax < 0.10 ) {
	  hpolana_deltaxccc->Fill(delta_anax);
	  GEMdx = true;
	}
	if ( delta_anax > -0.48 && delta_anax < -0.06 ) {
	  GEMdxP = true;
	}
      }
    }
    
    if( delta_x > -0.5 )
      HCALn = true;
    
    //-----------------------------------------------------------------------------------------------------------------------------
    // np -> pn ChEx in passive (steel) analyser, coordinate system = SBS TRANSPORT
    //-----------------------------------------------------------------------------------------------------------------------------

    Double_t helicity = *scalhel_hel;
    if( run_no >= 4000 ) 
      helicity = -1 * helicity;
    
    if( SBSGEMR && GEMdx && HCALn && !SBSGEMF ) {
      front_vtx = vertex_SBS;
      ana_vtx.SetXYZ( ana_x, ana_y, polana_dist ) ; 
      rear_vtx.SetXYZ( *sbs_gemCeR_track_x + ((polrgem_dist - polfgem_dist) * *sbs_gemCeR_track_xp), 
  		       *sbs_gemCeR_track_y + ((polrgem_dist - polfgem_dist) * *sbs_gemCeR_track_yp), 
  		       polrgem_dist );
      
      
      in_track.SetXYZ( tan ((ana_vtx.X() - front_vtx.X())/(ana_vtx.Z() - front_vtx.Z())) ,
		       tan ((ana_vtx.Y() - front_vtx.Y())/(ana_vtx.Z() - front_vtx.Z())) ,
		       1. );
      
      out_track.SetXYZ( *sbs_gemCeR_track_xp, 
			*sbs_gemCeR_track_yp, 
			1. );
      
      in_track  = in_track.Unit();
      out_track = out_track.Unit();
      
      TVector3 yaxistemp(0,1,0);
      TVector3 xaxistemp = yaxistemp.Cross(in_track).Unit();
      yaxistemp = in_track.Cross(xaxistemp).Unit();
      
      thsc = 180./M_PI * acos( out_track.Dot(in_track) ); 
      phsc = 180./M_PI * TMath::ATan2( out_track.Dot(yaxistemp), out_track.Dot(xaxistemp) );
      
      double sclose,zclose;
      
      calc_sclose_zclose( front_vtx, rear_vtx, in_track, out_track, sclose, zclose );
      bool conetestnp = conetest( front_vtx, in_track, (thsc/57.3), zclose, polrgem_dist );
      //      conetestnp = true;

      hpolg_sclnp->Fill( sclose ); 
      hpolg_zclnp->Fill( zclose , thsc ); 
      hpolg_thnp_cx->Fill ( thsc );
      
      if( thsc >= 3.0 && sclose < 0.003  && fabs(zclose-4.72) < 0.2 && fabs(phsc)<180 ) { 
  	hpolg_thnp_cxc->Fill ( thsc );
  	if( helicity == -1 ) {
  	  hpolg_phnp_cxm->Fill ( phsc );
  	  if( conetestnp ) {
  	    hpolg_phnp_cxmc->Fill ( phsc );
  	    hpolg_thnp_cxcc->Fill ( thsc );
  	  }
  	}
  	else if( helicity == 1 ) {
  	  hpolg_phnp_cxp->Fill ( phsc );
  	  if( conetestnp ) {
  	    hpolg_phnp_cxpc->Fill ( phsc );
  	    hpolg_thnp_cxcc->Fill ( thsc );
  	  }
  	}
      }
    }

    if( SBSGEMF && SBSGEMR && GEMdxP ) {
      hpolp_sclnp->Fill( *sbs_gemCeR_track_sclose ); 
      hpolp_zclnp->Fill( polana_dist-0.4+*sbs_gemCeR_track_zclose , 57.3*(*sbs_gemCeR_track_theta) ); 
      hpolp_thnp_cx->Fill ( 57.3*(*sbs_gemCeR_track_theta) );
      
      if( *sbs_gemCeR_track_theta >= (3./57.3) && *sbs_gemCeR_track_sclose < 0.003 ) { 
    	hpolp_thnp_cxc->Fill ( 57.3*(*sbs_gemCeR_track_theta) );
    	if( helicity == -1 ) 
    	  hpolp_phnp_cxm->Fill ( 57.3*(*sbs_gemCeR_track_phi) );
    	else if( helicity == 1 ) 
    	  hpolp_phnp_cxp->Fill ( 57.3*(*sbs_gemCeR_track_phi) );
      }
    }
    
    ev++;
  }
  
  outfile->Write();
  outfile->Close();

}

//-----------------------------------------------------------------------------------------------------------------------------

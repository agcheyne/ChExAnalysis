#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>
#include <TVectorD.h>

#include <iostream>
#include <stdio.h>
#include <fstream>

using namespace std;

// options
const Bool_t   PlotKine   = true;
const Bool_t   PlotHCal   = true;
const Bool_t   PlotPolAna = false;
const Bool_t   PlotPolG   = true;
const Bool_t   PlotPolP   = false;

Int_t rebin = 4;

void PlotGEnRPhcal( Int_t run_no = 9001 ) {

  TFile *f = new TFile(Form("hist/hist_genrp_%i.root",run_no),"READ");

  Double_t pcent   = 2.122;  
  Double_t pres    = 0.02;   
  Double_t sh_e    = 0.80;
  Double_t sh_min  = 0.70;
  Double_t sh_max  = 1.10;
  Double_t ps_min  = 0.1;

  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);

  gStyle->SetPadTopMargin(.05);
  gStyle->SetPadLeftMargin(.18);
  gStyle->SetPadRightMargin(.18);
  gStyle->SetPadBottomMargin(.15);

  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleOffset(1.5, "Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetTitleSize(0.055,"X");
  gStyle->SetTitleSize(0.055,"Y");

  gStyle->SetLabelOffset(0.01, "X");
  gStyle->SetLabelOffset(0.01, "Y");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetLabelSize(0.045,"X");
  gStyle->SetLabelSize(0.045,"Y");

  gStyle->SetNdivisions(105,"X");
  gStyle->SetNdivisions(105,"Y");

  gStyle->SetStripDecimals(kFALSE);

  TH1D* hkin_p        = (TH1D*)f->Get("hkin_p");
  TH1D* hkin_th       = (TH1D*)f->Get("hkin_th");
  TH1D* hkin_ph       = (TH1D*)f->Get("hkin_ph");
  TH1D* hkin_yt       = (TH1D*)f->Get("hkin_yt");
  TH1D* hkin_W        = (TH1D*)f->Get("hkin_W");
  
  TH1D* hkin_pc       = (TH1D*)f->Get("hkin_pc");
  TH1D* hkin_thc      = (TH1D*)f->Get("hkin_thc");
  TH1D* hkin_phc      = (TH1D*)f->Get("hkin_phc");
  TH1D* hkin_ytc      = (TH1D*)f->Get("hkin_ytc");
  TH1D* hkin_Wc       = (TH1D*)f->Get("hkin_Wc");
  
  TH2D *hbbcal2d_pss  = (TH2D*)f->Get("hbbcal2d_pssh");
  TH2D *hbbcal2d_pssc  = (TH2D*)f->Get("hbbcal2d_psshc");

  TH1D* hhcal_deltax   = (TH1D*)f->Get("hhcal_deltax");
  TH1D* hhcal_deltaxc  = (TH1D*)f->Get("hhcal_deltaxc");
  TH1D* hhcal_deltaxcc = (TH1D*)f->Get("hhcal_deltaxcc");
  TH1D* hhcal_deltay   = (TH1D*)f->Get("hhcal_deltay");
  TH1D* hhcal_deltayc  = (TH1D*)f->Get("hhcal_deltayc");
  TH2D* hhcal_deltaxy  = (TH2D*)f->Get("hhcal_deltaxy");
  TH2D* hhcal_deltaxyc = (TH2D*)f->Get("hhcal_deltaxyc");
  
  TH1D* hpolana_deltax    = (TH1D*)f->Get("hpolana_deltax");
  TH1D* hpolana_deltay    = (TH1D*)f->Get("hpolana_deltay");
  TH2D* hpolana_deltaxy   = (TH2D*)f->Get("hpolana_deltaxy");
  TH1D* hpolana_deltaxc   = (TH1D*)f->Get("hpolana_deltaxc");
  TH1D* hpolana_deltaxcc  = (TH1D*)f->Get("hpolana_deltaxcc");
  TH1D* hpolana_deltaxccc = (TH1D*)f->Get("hpolana_deltaxccc");
  TH1D* hpolana_deltayc   = (TH1D*)f->Get("hpolana_deltayc");
  TH2D* hpolana_deltaxyc  = (TH2D*)f->Get("hpolana_deltaxyc");

  TH1D* hpolg_thnp_cx   = (TH1D*)f->Get("hpolg_thnp_cx");
  TH1D* hpolg_thnp_cxc  = (TH1D*)f->Get("hpolg_thnp_cxc");
  TH1D* hpolg_thnp_cxcc = (TH1D*)f->Get("hpolg_thnp_cxcc");
  TH1D* hpolg_phnp_cxp  = (TH1D*)f->Get("hpolg_phnp_cxp");
  TH1D* hpolg_phnp_cxm  = (TH1D*)f->Get("hpolg_phnp_cxm");
  TH1D* hpolg_phnp_cxpc = (TH1D*)f->Get("hpolg_phnp_cxpc");
  TH1D* hpolg_phnp_cxmc = (TH1D*)f->Get("hpolg_phnp_cxmc");
  TH1D* hpolg_sclnp     = (TH1D*)f->Get("hpolg_sclnp");
  TH2D* hpolg_zclnp     = (TH2D*)f->Get("hpolg_zclnp");

  TH1D* hpolp_thnp_cx   = (TH1D*)f->Get("hpolp_thnp_cx");
  TH1D* hpolp_thnp_cxc  = (TH1D*)f->Get("hpolp_thnp_cxc");
  TH1D* hpolp_phnp_cxp  = (TH1D*)f->Get("hpolp_phnp_cxp");
  TH1D* hpolp_phnp_cxm  = (TH1D*)f->Get("hpolp_phnp_cxm");
  TH1D* hpolp_sclnp     = (TH1D*)f->Get("hpolp_sclnp");
  TH2D* hpolp_zclnp     = (TH2D*)f->Get("hpolp_zclnp");

  hpolg_phnp_cxp->Rebin(rebin);
  hpolg_phnp_cxm->Rebin(rebin);

  hpolg_phnp_cxpc->Rebin(rebin);
  hpolg_phnp_cxmc->Rebin(rebin);

  hpolp_phnp_cxp->Rebin(rebin);
  hpolp_phnp_cxm->Rebin(rebin);

  Int_t nphibins = hpolg_phnp_cxm->GetXaxis()->GetNbins();
    
  TLatex* tex;
  TLine *line;

  double nnhcal = -14;
  double nphcal = -14;
  double nchex = -14;
  double npgem = -14;
  
  //-----------------------------------------------------------------------------------------------------------------------------
  // Kinematic plots
  //-----------------------------------------------------------------------------------------------------------------------------

  if( PlotKine ) {

    TCanvas* ckine = new TCanvas("ckine","",1200,800);
    ckine->Divide(3,2);

    ckine->cd(1);
    hkin_pc->Draw();
    hkin_pc->GetXaxis()->SetTitle("p_{bbtrack} [GeV/c]");

    ckine->cd(2);
    hkin_thc->Draw();
    hkin_thc->GetXaxis()->SetTitle("#theta_tgt_{bbtrack} [rad]");

    ckine->cd(3);
    hkin_phc->Draw();
    hkin_phc->GetXaxis()->SetTitle("#phi_tgt_{bbtrack} [rad]");

    ckine->cd(4);
    hkin_ytc->Draw();
    hkin_ytc->GetXaxis()->SetTitle("y_tgt_{bbtrack} [m]");

    ckine->cd(5);
    hbbcal2d_pssc->Draw("colz");
    hbbcal2d_pssc->GetXaxis()->SetTitle("E_{SH}/p_{bbtrack}");
    hbbcal2d_pssc->GetYaxis()->SetTitle("E_{PS}/p_{bbtrack}");

    ckine->cd(6);
    hkin_Wc->Draw();
    hkin_Wc->GetXaxis()->SetTitle("W^{2} [GeV^{2}]");

    ckine->Print("kinematics_1.pdf");

  }
  
  //-----------------------------------------------------------------------------------------------------------------------------
  // HCal plots
  //-----------------------------------------------------------------------------------------------------------------------------

  if( PlotHCal ) {

    TCanvas* chcal = new TCanvas("chcal","",1200,800);
    chcal->Divide(2,2);

    chcal->cd(2);
    hhcal_deltaxc->Draw("");
    hhcal_deltaxc->GetXaxis()->SetTitle("HCal (Meas - Predicted) x[m]");

    chcal->cd(1);
    hhcal_deltayc->Draw("");
    hhcal_deltayc->GetXaxis()->SetTitle("HCal (Meas - Predicted) y[m]");
    
    chcal->cd(3);
    hhcal_deltaxyc->Draw("colz");
    hhcal_deltaxyc->GetXaxis()->SetTitle("HCal (Meas - Predicted) y[m]");
    hhcal_deltaxyc->GetYaxis()->SetTitle("HCal (Meas - Predicted) x[m]");

    chcal->cd(4); 
    int pbmin   = 23;
    int pbmax   = 58;
    
    double cmin = -2.0;
    double cmax = 1.0;
    
    hhcal_deltaxcc->GetXaxis()->SetRangeUser(cmin, cmax);
    hhcal_deltaxcc->Draw("");
    hhcal_deltaxcc->GetXaxis()->SetTitle("HCal (meas - pred) x [m]");
    
    float binwidth = hhcal_deltaxcc->GetXaxis()->GetBinWidth(1); 
    
    int pbrange = pbmax - pbmin; 
    float pbins[pbrange], perr[pbrange];
    
    for( int i = pbmin; i < pbmax; i++) { 
      pbins[i - pbmin] = hhcal_deltaxcc->GetBinContent( i ); 
      perr[i - pbmin]  = hhcal_deltaxcc->GetBinError( i ); 
      hhcal_deltaxcc->SetBinContent( i, 0 ); 
      hhcal_deltaxcc->SetBinError( i, 0 ); 
    } 
    
    TF1* bg = new TF1("bg","pol3(0)",cmin,cmax);
    bg->SetLineColor(kBlue);
    hhcal_deltaxcc->Fit("bg","Q","",cmin,cmax);
    
    double par0 = bg->GetParameter(0);
    double par1 = bg->GetParameter(1);
    double par2 = bg->GetParameter(2);
    double par3 = bg->GetParameter(3);
    
    for( int i = pbmin; i < pbmax; i++) { 
      hhcal_deltaxcc->SetBinContent( i, pbins[i - pbmin] ); 
      hhcal_deltaxcc->SetBinError( i, perr[i - pbmin] ); 
    } 

    bg->Draw("same");

    TF1* ppeakbg = new TF1("ppeakbg","pol3(0)+gaus(4)+gaus(7)",cmin,cmax);
    ppeakbg->SetLineColor(kBlack);
    
    ppeakbg->FixParameter(0,par0);
    ppeakbg->FixParameter(1,par1);
    ppeakbg->FixParameter(2,par2);
    ppeakbg->FixParameter(3,par3);
    ppeakbg->SetParameter(4,600);
    ppeakbg->SetParameter(5,-0.9);
    ppeakbg->SetParameter(6,0.15);
    ppeakbg->SetParameter(7,300);
    ppeakbg->SetParameter(8,0);
    ppeakbg->SetParameter(9,0.15);
    
    hhcal_deltaxcc->Fit("ppeakbg","Q","",cmin,cmax);

    float pmean      = ppeakbg->GetParameter(5); 
    float pmean_err  = ppeakbg->GetParError(5); 

    float nmean      = ppeakbg->GetParameter(8); 
    float nmean_err  = ppeakbg->GetParError(8); 
    
    float psigma      = ppeakbg->GetParameter(6); 
    float psigma_err  = ppeakbg->GetParError(6); 
    
    float nsigma      = ppeakbg->GetParameter(9); 
    float nsigma_err  = ppeakbg->GetParError(9); 
    
    float pmin3sig    = pmean - 3*psigma;
    float pplu3sig    = pmean + 3*psigma;

    float nmin3sig    = nmean - 3*nsigma;
    float nplu3sig    = nmean + 3*nsigma;

    float pintegral = (ppeakbg->Integral(pmin3sig,pplu3sig) - bg->Integral(pmin3sig,pplu3sig)) / binwidth;
    float nintegral = (ppeakbg->Integral(nmin3sig,nplu3sig) - bg->Integral(nmin3sig,nplu3sig)) / binwidth;

    TLatex* tex = new TLatex( 0.55, 0.83, Form("#sigma_{n} = %3.2f m",fabs(nsigma)));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.05);
    tex->Draw();

    tex = new TLatex( 0.55, 0.78, Form("#sigma_{p} = %3.2f m",fabs(psigma)));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.05);
    tex->Draw();

    tex = new TLatex( 0.55, 0.73, Form("Nn = %3.0f",fabs(1.0*nintegral)));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.05);
    tex->Draw();

    nphcal = fabs(1.0*pintegral);
    nnhcal = fabs(1.0*nintegral);

    tex = new TLatex( 0.55, 0.68, Form("Np = %3.0f",fabs(1.0*pintegral)));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.05);
    tex->Draw();

    tex = new TLatex( 0.55, 0.63, Form("#Delta x = %3.2f m", fabs(nmean - pmean)));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.05);
    tex->Draw();

    chcal->Print("kinematics5.pdf");
  }

  //-----------------------------------------------------------------------------------------------------------------------------
  // Analyzer plots
  //-----------------------------------------------------------------------------------------------------------------------------

  if( PlotPolAna ) {
      
    TCanvas* cpolana = new TCanvas("cpolana","",1200,800);
    cpolana->Divide(2,2);

    cpolana->cd(2);
    hpolana_deltaxc->Draw("");
    hpolana_deltaxc->GetXaxis()->SetTitle("Analyzer (Meas - Predicted) x[m]");

    cpolana->cd(1);
    hpolana_deltayc->Draw("");
    hpolana_deltayc->GetXaxis()->SetTitle("Analyzer (Meas - Predicted) y[m]");
    
    cpolana->cd(3)->SetLogz(1);
    hpolana_deltaxyc->Draw("colz");
    hpolana_deltaxyc->GetXaxis()->SetTitle("Analyzer (Meas - Predicted) x[m]");
    hpolana_deltaxyc->GetYaxis()->SetTitle("Analyzer (Meas - Predicted) y[m]");

    cpolana->cd(4); 
    int pbmin   = 31;
    int pbmax   = 58;
    
    double cmin = -0.9;
    double cmax = 0.9;
    
    hpolana_deltaxcc->GetXaxis()->SetRangeUser(cmin, cmax);
    hpolana_deltaxcc->Draw("");
    hpolana_deltaxcc->GetXaxis()->SetTitle("Steel Analyzer (meas - pred) x [m]");

    float binwidth = hpolana_deltaxcc->GetXaxis()->GetBinWidth(1); 
    
    int pbrange = pbmax - pbmin; 
    float pbins[pbrange], perr[pbrange];
    
    for( int i = pbmin; i < pbmax; i++) { 
      pbins[i - pbmin] = hpolana_deltaxcc->GetBinContent( i ); 
      perr[i - pbmin]  = hpolana_deltaxcc->GetBinError( i ); 
      hpolana_deltaxcc->SetBinContent( i, 0 ); 
      hpolana_deltaxcc->SetBinError( i, 0 ); 
    } 
    
    TF1* bg = new TF1("bg","pol3(0)",cmin,cmax);
    bg->SetLineColor(kBlue);
    hpolana_deltaxcc->Fit("bg","Q","",cmin,cmax);
    
    double par0 = bg->GetParameter(0);
    double par1 = bg->GetParameter(1);
    double par2 = bg->GetParameter(2);
    double par3 = bg->GetParameter(3);
    
    for( int i = pbmin; i < pbmax; i++) { 
      hpolana_deltaxcc->SetBinContent( i, pbins[i - pbmin] ); 
      hpolana_deltaxcc->SetBinError( i, perr[i - pbmin] ); 
    } 

    bg->Draw("same");

    TF1* ppeakbg = new TF1("ppeakbg","pol3(0)+gaus(4)",cmin,cmax);
    ppeakbg->SetLineColor(kBlack);
    
    ppeakbg->FixParameter(0,par0);
    ppeakbg->FixParameter(1,par1);
    ppeakbg->FixParameter(2,par2);
    ppeakbg->FixParameter(3,par3);
    ppeakbg->SetParameter(4,600);
    ppeakbg->SetParameter(5,-0.28);
    ppeakbg->SetParameter(6,0.06);
    
    hpolana_deltaxcc->Fit("ppeakbg","Q","",cmin,cmax);

    float pmean      = ppeakbg->GetParameter(5); 
    float pmean_err  = ppeakbg->GetParError(5); 

    float psigma      = ppeakbg->GetParameter(6); 
    float psigma_err  = ppeakbg->GetParError(6); 
    
    float pmin3sig    = pmean - 3*psigma;
    float pplu3sig    = pmean + 3*psigma;

    float pintegral = (ppeakbg->Integral(pmin3sig,pplu3sig) - bg->Integral(pmin3sig,pplu3sig)) / binwidth;
    
    hpolana_deltaxccc->SetLineColor(2);
    hpolana_deltaxccc->Draw("same");

    nchex = hpolana_deltaxccc->Integral(0,100);

    tex = new TLatex( 0.45, 0.83, Form("#sigma_{p} = %3.2f m",fabs(psigma)));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.05);
    tex->Draw();

    tex = new TLatex( 0.45, 0.78, Form("#Delta x = %3.2f m", pmean) );
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.05);
    tex->Draw();

    tex = new TLatex( 0.45, 0.73, Form("Np = %3.0f",fabs(1.0*pintegral)));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.05);
    tex->Draw();

    tex = new TLatex( 0.45, 0.68, "NpGEM/NpHCAL"); 
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.05);
    tex->Draw();
	
    tex = new TLatex( 0.45, 0.63, Form("= %3.2f", fabs(1.0*pintegral)/nphcal) );
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.05);
    tex->Draw();

    npgem = fabs(1.0*pintegral);

    cpolana->Print("kinematics9a.pdf");

  }

  //-------------------------------------------------------------------------------


  if( PlotPolG ) {

    TCanvas* cpolnp = new TCanvas("cpolnp","",1200,800);
    cpolnp->Divide(2,2);
    
    cpolnp->cd(1);
    hpolg_sclnp->Draw("");
    hpolp_sclnp->GetXaxis()->SetRangeUser(0, 0.103);
    hpolp_sclnp->GetXaxis()->SetTitle("sclose [m]*");
    hpolp_sclnp->GetYaxis()->SetRangeUser(0, 5000);
   

    tex = new TLatex( 0.35, 0.85, "No GEMs" );
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.075);
    tex->Draw();

    tex = new TLatex( 0.35, 0.75, "" );
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.075);
    tex->Draw();

    tex = new TLatex( 0.35, 0.75, Form("Mean sclose = %3.2f mm", 1000.*hpolg_sclnp->GetMean()) );
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.065);
    tex->Draw();

    cpolnp->cd(2);
 
    hpolg_zclnp->GetXaxis()->SetTitle("zclose test[m]");
    hpolg_zclnp->GetYaxis()->SetTitle("#theta [deg]");
    hpolg_zclnp->GetXaxis()->SetRangeUser(3.72, 5.72 );
    hpolg_zclnp->Draw("colz");

    cpolnp->cd(3);
    hpolg_thnp_cx->Draw("");
    hpolg_thnp_cx->GetXaxis()->SetTitle("#theta [deg]");

    hpolg_thnp_cxcc->Draw("SAME");
    hpolg_thnp_cxcc->SetLineColor(2);

    nchex = hpolg_thnp_cxcc->Integral(0,100);

    tex = new TLatex( 0.45, 0.85, Form("N = %3.0f", nchex) );
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();

    cpolnp->cd(4);

    hpolg_phnp_cxpc->Sumw2();
    hpolg_phnp_cxpc->SetMarkerStyle(21);
    hpolg_phnp_cxpc->SetMarkerColor(6);
    hpolg_phnp_cxpc->SetLineColor(6);
    hpolg_phnp_cxpc->SetMarkerSize(0.5);
    hpolg_phnp_cxpc->Draw("same");
    hpolg_phnp_cxpc->GetXaxis()->SetTitle("#phi [deg]");

    hpolg_phnp_cxmc->Sumw2();
    hpolg_phnp_cxmc->Draw("same");
    hpolg_phnp_cxmc->SetMarkerStyle(21);
    hpolg_phnp_cxmc->SetMarkerColor(1);
    hpolg_phnp_cxmc->SetLineColor(1);
    hpolg_phnp_cxmc->SetMarkerSize(0.5);

    cpolnp->Print("kinematics9b.pdf");  

    //-------------------------------------------------------------------------------

    TCanvas* cpolnp1 = new TCanvas("cpolnp1","",1200,800);

    TH1D* hsumcxc  = new TH1D("hsumcxgc","",nphibins,-180,180);
    TH1D* hdiffcxc = new TH1D("hdiffcxgc","",nphibins,-180,180);
    TH1D* hasymcxc = new TH1D("hasymcxgc","",nphibins,-180,180);
    hsumcxc->Sumw2();
    hdiffcxc->Sumw2();
    hasymcxc->Sumw2();

    double nplusc  = hpolg_phnp_cxpc->Integral(0,100);
    double nminusc = hpolg_phnp_cxmc->Integral(0,100); 
    double ntotc   = nplusc + nminusc;

    hsumcxc->Add( hpolg_phnp_cxpc, hpolg_phnp_cxmc );
    hdiffcxc->Add( hpolg_phnp_cxpc, hpolg_phnp_cxmc, 2*nminusc/ntotc, -2*nplusc/ntotc );
    hasymcxc->Divide( hdiffcxc, hsumcxc, 1, 1);

    hasymcxc->Draw("E1");
    hasymcxc->SetMarkerStyle(21);
    hasymcxc->SetMarkerColor(kBlue);
    hasymcxc->SetLineColor(kBlue);
    hasymcxc->SetMarkerSize(0.5);
    hasymcxc->GetXaxis()->SetTitle("#varphi (degrees)");
    hasymcxc->GetYaxis()->SetTitle("Asymmetry");
    hasymcxc->GetYaxis()->SetRangeUser(-0.05,0.06);

    TF1 *fitcxc = new TF1("fitcxc", "[1]*sin(x/57.29578)-[0]*cos(x/57.29578)",-180,180);
    hasymcxc->Fit(fitcxc,"","",-180,180);

    Double_t Pynpc  = fitcxc->GetParameter(0);
    Double_t ePynpc = fitcxc->GetParError(0);
    Double_t Pxnpc  = fitcxc->GetParameter(1);
    Double_t ePxnpc = fitcxc->GetParError(1);

    tex = new TLatex( 0.4, 0.89, Form("hAyPx^{FP} = %3.3f +/- %3.3f", Pxnpc, ePxnpc));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.05);
    tex->Draw();

    tex = new TLatex( 0.4, 0.82, Form("hAyPy^{FP} = %3.3f +/- %3.3f", Pynpc, ePynpc));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.05);
    tex->Draw();

    tex = new TLatex( 0.4, 0.75, Form("Ay          = %3.3f +/- %3.3f", 
    				      Pxnpc/0.41, ePxnpc/0.41));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.05);
    tex->Draw();

    tex = new TLatex( 0.2, 0.20, Form("NChEx = %3.2f k", ntotc/1000.) );
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.05);
    tex->Draw();

    cpolnp1->Print("kinematics9c.pdf");  

  }

  //-------------------------------------------------------------------------------

  if( PlotPolP ) {

    TCanvas* cpolnpp = new TCanvas("cpolnpp","",1200,800);
    cpolnpp->Divide(2,2);
    
    cpolnpp->cd(1);
    hpolp_sclnp->Draw("");
    hpolp_sclnp->GetXaxis()->SetRangeUser(0, 0.003 );
    hpolp_sclnp->GetXaxis()->SetTitle("sclose [m]");
    hpolp_sclnp->GetYaxis()->SetRangeUser(0, 4550 );
    
    tex = new TLatex( 0.35, 0.85, "Proton Polarimeter" );
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.075);
    tex->Draw();

    cpolnpp->cd(2);
    hpolp_zclnp->Draw("colz");
    hpolp_zclnp->GetXaxis()->SetRangeUser(3.7, 5.7 );
    hpolp_zclnp->GetYaxis()->SetRangeUser(0, 20 );
    hpolp_zclnp->GetXaxis()->SetTitle("zclose [m]");
    hpolp_zclnp->GetYaxis()->SetTitle("#theta [deg]");

    cpolnpp->cd(3);
    hpolp_thnp_cx->Draw("");
    hpolp_thnp_cx->GetXaxis()->SetRangeUser(0, 20. );
    hpolp_thnp_cx->GetXaxis()->SetTitle("#theta [deg]");

    hpolp_thnp_cxc->Draw("SAME");
    hpolp_thnp_cxc->SetLineColor(2);

    double npp1 = hpolp_thnp_cxc->Integral(0,100);

    tex = new TLatex( 0.45, 0.85, Form("N = %3.0f", npp1) );
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();

    cpolnpp->cd(4);

    hpolp_phnp_cxp->Sumw2();
    hpolp_phnp_cxp->SetMarkerStyle(21);
    hpolp_phnp_cxp->SetMarkerColor(6);
    hpolp_phnp_cxp->SetLineColor(6);
    hpolp_phnp_cxp->SetMarkerSize(0.5);
    hpolp_phnp_cxp->Draw("same");
    hpolp_phnp_cxp->GetXaxis()->SetTitle("#phi [deg]");

    hpolp_phnp_cxm->Sumw2();
    hpolp_phnp_cxm->Draw("same");
    hpolp_phnp_cxm->SetMarkerStyle(21);
    hpolp_phnp_cxm->SetMarkerColor(1);
    hpolp_phnp_cxm->SetLineColor(1);
    hpolp_phnp_cxm->SetMarkerSize(0.5);

    cpolnpp->Print("kinematics9d.pdf");  

  //-------------------------------------------------------------------------------

    TCanvas* cpolnp1p = new TCanvas("cpolnp1p","",1200,800);

    TH1D* hsumcx  = new TH1D("hsumcxp","",nphibins,-180,180);
    TH1D* hdiffcx = new TH1D("hdiffcxp","",nphibins,-180,180);
    TH1D* hasymcx = new TH1D("hasymcxp","",nphibins,-180,180);
    hsumcx->Sumw2();
    hdiffcx->Sumw2();
    hasymcx->Sumw2();

    double nplus  = hpolp_phnp_cxp->Integral(0,100);
    double nminus = hpolp_phnp_cxm->Integral(0,100); 
    double ntot   = nplus + nminus;

    hsumcx->Add( hpolp_phnp_cxp, hpolp_phnp_cxm );
    hdiffcx->Add( hpolp_phnp_cxp, hpolp_phnp_cxm, 2*nminus/ntot, -2*nplus/ntot );
    hasymcx->Divide( hdiffcx, hsumcx, 1, 1);

    hasymcx->Draw("E1");
    hasymcx->SetMarkerStyle(21);
    hasymcx->SetMarkerColor(kBlue);
    hasymcx->SetLineColor(kBlue);
    hasymcx->SetMarkerSize(0.5);
    hasymcx->GetXaxis()->SetTitle("#varphi (degrees)");
    hasymcx->GetYaxis()->SetTitle("Proton Asymmetry");
    hasymcx->GetYaxis()->SetRangeUser(-0.05,0.06);

    TF1 *fitcx = new TF1("fitcx", "[1]*sin(x/57.29578)-[0]*cos(x/57.29578)",-180,180);
    hasymcx->Fit(fitcx,"","",-180,180);

    Double_t Pynp  = fitcx->GetParameter(0);
    Double_t ePynp = fitcx->GetParError(0);
    Double_t Pxnp  = fitcx->GetParameter(1);
    Double_t ePxnp = fitcx->GetParError(1);

    tex = new TLatex( 0.4, 0.89, Form("hAyPx^{FP} = %3.3f +/- %3.3f", Pxnp, ePxnp));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.05);
    tex->Draw();

    tex = new TLatex( 0.4, 0.82, Form("hAyPy^{FP} = %3.3f +/- %3.3f", Pynp, ePynp));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.05);
    tex->Draw();

    tex = new TLatex( 0.4, 0.75, Form("Ay          = %3.3f +/- %3.3f", Pxnp/0.41, ePxnp/0.41));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.05);
    tex->Draw();

    cpolnp1p->Print("kinematics9e.pdf");  

  }

  gSystem->Exec(Form("pdfunite  kinem*.pdf pdf/plots_genrp-%d.pdf", run_no));  
  gSystem->Exec("rm  kinema*.pdf");  

  int run = -1;
  float charge = -1.0;
  
  ifstream FileCharge;
  FileCharge.open("log/RunChargeGEnRP.csv");
  string l;
  while( getline(FileCharge, l)) {
    sscanf(l.c_str(), "%d,%f", &run, &charge);
    if( run == run_no )
      break;
  }
  FileCharge.close();

  if( run != -1 && charge != -1 && charge != 0 ) {
    ofstream FileYield;
    FileYield.open("log/RunYieldGEnRP.csv", ios::app);
    FileYield << run_no << "," << nphcal/charge << "," << nnhcal/charge <<  "," << (nphcal/charge)/1.7e6 << "," 
	      << npgem/charge << "," << npgem/nphcal << "," << (npgem/nphcal)*((nphcal/charge)/1.7e6) << "," << nchex/charge << endl;
    FileYield.close();
  }
}

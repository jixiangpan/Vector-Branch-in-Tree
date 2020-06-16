////
////
//// the original file is plot.C from Xin
////
////


//
// test
//
// #include<iostream>
// #include<fstream>
// #include<cmath>
// #include "stdlib.h"
// using namespace std;

// #include<map>
// #include<vector>
// #include<set>

// #include "TH1.h"
// #include "TH2.h"
// #include "TH3.h"
// #include "THStack.h"
// #include "TF1.h"
// #include "TLine.h"
// #include "TMath.h"
// #include "TGraph.h"
// #include "TGraph2D.h"
// #include "TGraphErrors.h"

// #include "TFile.h"
// #include "TTree.h"
// #include "TChain.h"
// #include "TBranch.h"

// #include "TRandom3.h"
// #include "TGaxis.h"
// #include "TStyle.h"

// #include "TCanvas.h"
// #include "TLegend.h"
// #include "TString.h"
// #include "TROOT.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////// Global variables


TGraph *g_dEdx_muon;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

double func_dQdx_from_dEdx_by_ArgoNeut_model(double dEdx, double alpha, double beta)
{
  double result = log(alpha + beta/1.38/0.273*dEdx)/(23.6e-6*beta/1.38/0.273);
  return result;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAIN FUNCTION /////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_mc_obj()
{
  TString roostr = "";
  
  const int electron = 11;
  const int muon     = 13;
  const int proton   = 2212;
  const int pion     = 211;

  int color_muon = kRed;
  int color_proton = kBlack;
  int line_width = 3;
  int line_style = 9;

  const double const_dx = 0.59;

  const double alpha_ArgoNeut = 0.93;
  const double beta_ArgoNeut  = 0.212;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  TFile *file_mc_ave = new TFile("stopping_ave_dQ_dx.root", "read");
  
  TGraph *gh_mc_ave_proton = (TGraph*)file_mc_ave->Get("proton"); gh_mc_ave_proton->SetName("gh_mc_ave_proton");
  TGraph *gh_mc_ave_muon = (TGraph*)file_mc_ave->Get("muon");     gh_mc_ave_muon->SetName("gh_mc_ave_muon");
  //TGraph *gh_mc_ave_pion = (TGraph*)file_mc_ave->Get("pion");
  //TGraph *gh_mc_ave_kaon = (TGraph*)file_mc_ave->Get("kaon");
  //TGraph *gh_mc_ave_electron = (TGraph*)file_mc_ave->Get("electron");
  
  roostr = "canv_gh_basic";
  TCanvas *canv_gh_basic = new TCanvas(roostr, roostr, 900, 650);
  canv_gh_basic->SetLeftMargin(0.15);
  canv_gh_basic->SetRightMargin(0.15);
  canv_gh_basic->SetBottomMargin(0.15);

  roostr = "h2_gh_basic";
  TH2D *h2_gh_basic = new TH2D(roostr, roostr, 100, 0, 100, 200, 0, 200e3);
  h2_gh_basic->Draw();
  h2_gh_basic->SetStats(0);
  h2_gh_basic->SetTitle("");
  h2_gh_basic->GetXaxis()->SetTitleSize(0.05);
  h2_gh_basic->GetXaxis()->SetLabelSize(0.05);
  h2_gh_basic->GetYaxis()->SetTitleSize(0.05);
  h2_gh_basic->GetYaxis()->SetLabelSize(0.05);
  h2_gh_basic->GetXaxis()->SetTitle("Residual range (cm)");
  h2_gh_basic->GetYaxis()->SetTitle("dQ/dx (e/cm)");
  
  //gh_mc_ave_muon->Draw("L same");
  gh_mc_ave_muon->SetLineColor(color_muon);
  gh_mc_ave_muon->SetLineWidth( line_width );  
  //gh_mc_ave_proton->Draw("L same");
  gh_mc_ave_proton->SetLineColor(color_proton);
  gh_mc_ave_proton->SetLineWidth( line_width );

  ////////////
  
  TGraph *g_dEdx_proton = new TGraph();
  ifstream infile4("proton_hybrid_dedx_vs_resrange.txt");
  for (Int_t i=0;i!=648;i++){
    double x(0),y(0);
    infile4 >> x >> y;
    g_dEdx_proton->SetPoint(i,x,y);
  }
  
  //TGraph *g_dEdx_muon = new TGraph();
  g_dEdx_muon = new TGraph();
  ifstream infile5("muon_hybrid_dedx_vs_resrange.txt");
  for (Int_t i=0;i!=732;i++){
    double x(0),y(0);
    infile5 >> x >> y;
    g_dEdx_muon->SetPoint(i,x,y);
  }
  
  TGraph *g_dQdx_model_muon = new TGraph();
  TGraph *g_dQdx_model_proton = new TGraph();
  // default MCC 8 numbers 
  Double_t alpha = 0.93;
  Double_t beta = 0.212;

  int line_g_dEdx_muon = 0;
  double sum_g_dEdx_muon = 0;  
  for (int i=0;i!=g_dEdx_muon->GetN();i++){
    double dEdx = 0;
    double rr = 0;
    double dQdx = 0;
    g_dEdx_muon->GetPoint(i, rr, dEdx);
    dQdx = log(alpha + beta/1.38/0.273*dEdx)/(23.6e-6*beta/1.38/0.273);
    g_dQdx_model_muon->SetPoint(i, rr, dQdx);

    if( rr>=50 && rr<=90 ) {
      line_g_dEdx_muon++;
      sum_g_dEdx_muon += dQdx;
    }
  }
  cout<<" Average dQdx of model in MIP "<<sum_g_dEdx_muon/line_g_dEdx_muon<<endl;
  
  for (int i=0;i!=g_dEdx_proton->GetN();i++){
    double dEdx = 0;
    double rr = 0;
    double dQdx = 0;
    g_dEdx_proton->GetPoint(i, rr, dEdx);
    dQdx = log(alpha + beta/1.38/0.273*dEdx)/(23.6e-6*beta/1.38/0.273);
    g_dQdx_model_proton->SetPoint(i, rr, dQdx);
  }
  
  canv_gh_basic->cd();  
  g_dQdx_model_muon->Draw("L same");
  g_dQdx_model_muon->SetLineColor(color_muon);
  g_dQdx_model_muon->SetLineWidth( line_width );
  g_dQdx_model_muon->SetLineStyle( line_style );  
  g_dQdx_model_proton->Draw("L same");
  g_dQdx_model_proton->SetLineColor(color_proton);
  g_dQdx_model_proton->SetLineWidth( line_width );
  g_dQdx_model_proton->SetLineStyle( line_style );

  ///////////////////////////////////////////////////////////

  TFile *roofile = new TFile("result_mc_truth_100evts.root", "read");
  TTree *tree = (TTree*)roofile->Get("tree_mc_truth");
  
  // Declaration of leaf types
  Int_t           run_no = 0;
  Int_t           subrun_no = 0;
  Int_t           event_no = 0;
  Int_t           muon_no = 0;
  Int_t           electron_no = 0;
  Int_t           proton_no = 0;
  Int_t           pion_no = 0;
  Int_t           p_id = 0;
  vector<double>  *vc_val_dQdx = 0;
  vector<double>  *vc_val_rr = 0;
  Double_t        val_length = 0;
  vector<double>  *vc_val_dQ = 0;
  vector<double>  *vc_val_dE = 0;
  vector<double>  *vc_val_dx = 0;
  vector<double>  *vc_val_x1 = 0;
  vector<double>  *vc_val_y1 = 0;
  vector<double>  *vc_val_z1 = 0;
  vector<double>  *vc_val_x2 = 0;
  vector<double>  *vc_val_y2 = 0;
  vector<double>  *vc_val_z2 = 0;

  // List of branches
  TBranch        *b_run_no;   //!
  TBranch        *b_subrun_no;   //!
  TBranch        *b_event_no;   //!
  TBranch        *b_muon_no;   //!
  TBranch        *b_electron_no;   //!
  TBranch        *b_proton_no;   //!
  TBranch        *b_pion_no;   //!
  TBranch        *b_p_id;   //!
  TBranch        *b_vc_val_dQdx;   //!
  TBranch        *b_vc_val_rr;   //!
  TBranch        *b_val_length;   //!
  TBranch        *b_vc_val_dQ;   //!
  TBranch        *b_vc_val_dE;   //!
  TBranch        *b_vc_val_dx;   //!
  TBranch        *b_vc_val_x1;   //!
  TBranch        *b_vc_val_y1;   //!
  TBranch        *b_vc_val_z1;   //!
  TBranch        *b_vc_val_x2;   //!
  TBranch        *b_vc_val_y2;   //!
  TBranch        *b_vc_val_z2;   //!

  // Set branch addresses and branch pointers 
  tree->SetBranchAddress("run_no", &run_no, &b_run_no);
  tree->SetBranchAddress("subrun_no", &subrun_no, &b_subrun_no);
  tree->SetBranchAddress("event_no", &event_no, &b_event_no);
  tree->SetBranchAddress("muon_no", &muon_no, &b_muon_no);
  tree->SetBranchAddress("electron_no", &electron_no, &b_electron_no);
  tree->SetBranchAddress("proton_no", &proton_no, &b_proton_no);
  tree->SetBranchAddress("pion_no", &pion_no, &b_pion_no);
  tree->SetBranchAddress("p_id", &p_id, &b_p_id);
  tree->SetBranchAddress("vc_val_dQdx", &vc_val_dQdx, &b_vc_val_dQdx);
  tree->SetBranchAddress("vc_val_rr", &vc_val_rr, &b_vc_val_rr);
  tree->SetBranchAddress("val_length", &val_length, &b_val_length);
  tree->SetBranchAddress("vc_val_dQ", &vc_val_dQ, &b_vc_val_dQ);
  tree->SetBranchAddress("vc_val_dE", &vc_val_dE, &b_vc_val_dE);
  tree->SetBranchAddress("vc_val_dx", &vc_val_dx, &b_vc_val_dx);
  tree->SetBranchAddress("vc_val_x1", &vc_val_x1, &b_vc_val_x1);
  tree->SetBranchAddress("vc_val_y1", &vc_val_y1, &b_vc_val_y1);
  tree->SetBranchAddress("vc_val_z1", &vc_val_z1, &b_vc_val_z1);
  tree->SetBranchAddress("vc_val_x2", &vc_val_x2, &b_vc_val_x2);
  tree->SetBranchAddress("vc_val_y2", &vc_val_y2, &b_vc_val_y2);
  tree->SetBranchAddress("vc_val_z2", &vc_val_z2, &b_vc_val_z2);
  
  ///////////////////////////////////////////////////////////

  roostr = "h1_dQdx_model2mc_muon";
  TH1D *h1_dQdx_model2mc_muon = new TH1D(roostr, roostr, 200, 0.9, 1.1);  
  roostr = "h1_dQdx_model2mc_proton";
  TH1D *h1_dQdx_model2mc_proton = new TH1D(roostr, roostr, 200, 0.9, 1.1);

  roostr = "h2_dQdx_rr_2ends_tot";
  TH2D *h2_dQdx_rr_2ends_tot = new TH2D(roostr, roostr, 100,0,100,300,0,300e3);  
  roostr = "h2_dQdx_rr_2ends_muon";
  TH2D *h2_dQdx_rr_2ends_muon = new TH2D(roostr, roostr, 100,0,100,300,0,300e3);  
  roostr = "h2_dQdx_rr_2ends_proton";
  TH2D *h2_dQdx_rr_2ends_proton = new TH2D(roostr, roostr, 100,0,100,300,0,300e3);
  
  roostr = "h1_dx6mm_2ends";
  TH1D *h1_dx6mm_2ends = new TH1D(roostr, roostr, 100, 0, 2);
  roostr = "h1_dQdx_truth_2ends_dx6mm_muon";
  TH1D *h1_dQdx_truth_2ends_dx6mm_muon = new TH1D(roostr, roostr, 300, 0, 300e3);  
  roostr = "h1_dQdx_model_2ends_dx6mm_muon";
  TH1D *h1_dQdx_model_2ends_dx6mm_muon = new TH1D(roostr, roostr, 300, 0, 300e3);
  roostr = "h1_dQdx_truth2model_2ends_dx6mm_muon";
  TH1D *h1_dQdx_truth2model_2ends_dx6mm_muon = new TH1D(roostr, roostr, 200, 0.8, 1.2);
  
  
  roostr = "h2_dQdx_rr_along_tot";
  TH2D *h2_dQdx_rr_along_tot = new TH2D(roostr, roostr, 100,0,100,300,0,300e3);

  roostr = "h1_dx6mm_along";
  TH1D *h1_dx6mm_along = new TH1D(roostr, roostr, 100, 0, 2);
  roostr = "h1_dQdx_truth_along_dx6mm_muon";
  TH1D *h1_dQdx_truth_along_dx6mm_muon = new TH1D(roostr, roostr, 300, 0, 300e3);  
  roostr = "h1_dQdx_model_along_dx6mm_muon";
  TH1D *h1_dQdx_model_along_dx6mm_muon = new TH1D(roostr, roostr, 300, 0, 300e3);
  roostr = "h1_dQdx_truth2model_along_dx6mm_muon";
  TH1D *h1_dQdx_truth2model_along_dx6mm_muon = new TH1D(roostr, roostr, 200, 0.8, 1.2);
  
  roostr = "h2_dQdx_rr_along_muon";
  TH2D *h2_dQdx_rr_along_muon = new TH2D(roostr, roostr, 100,0,100,200,0,200e3);  
  roostr = "h2_dQdx_rr_along_proton";
  TH2D *h2_dQdx_rr_along_proton = new TH2D(roostr, roostr, 100,0,100,200,0,200e3);
  roostr = "h2_dQdx_rr_along_muons_protons";
  TH2D *h2_dQdx_rr_along_muons_protons = new TH2D(roostr, roostr, 100,0,100,200,0,200e3);
  
  
  ///////////////////////////////////////////////////////////

  int entries = tree->GetEntries();

  entries = 2000;
  
  cout<<endl<<" -----> entries "<<entries<<endl<<endl;

  for(int ientry=0; ientry<=entries; ientry++) {
    tree->GetEntry(ientry);

    if(ientry%200==0) cout<<TString::Format(" ---> %6.4f", ientry*1./entries)<<endl;

    // if( 0 ) {
    //   alpha = 0.93;
    //   beta = 0.212;

    //   int size = vc_val_dE->size();
    //   for(int idx=0; idx<size; idx++) {

    // 	double val_dE = vc_val_dE->at(idx);
    // 	double val_dQ = vc_val_dQ->at(idx);
    // 	double val_dx = vc_val_dx->at(idx);
    // 	double val_rr = vc_val_rr->at(idx);

    // 	double mc_val_dEdx = val_dE/val_dx;
    // 	double mc_val_dQdx = val_dQ/val_dx;
	
    // 	double model_val_dQdx = log(alpha + beta/1.38/0.273*mc_val_dEdx)/(23.6e-6*beta/1.38/0.273);

    // 	if( p_id==muon ) {
    // 	  h1_dQdx_model2mc_muon->Fill( mc_val_dQdx/model_val_dQdx );
    // 	}
    // 	if( p_id==proton ) {
    // 	  h1_dQdx_model2mc_proton->Fill( mc_val_dQdx/model_val_dQdx );
    // 	}
	
    //   }// for(int idx=0; idx<size; idx++)      
    // }// if( 1 || (run_no==7017 && subrun_no==1274 && event_no==63725) )

    int size = vc_val_dE->size();
    // cout<<" ---> check "<<size<<endl;
    
    int line_along = size-1;
    double user_dx_along = 0;
    double user_dE_along = 0;
    double user_dQ_along = 0;
    double user_rr_along = 0;

    
    int line_2ends = size-1;
    double user_dx_2ends = 0;
    double user_dE_2ends = 0;
    double user_dQ_2ends = 0;    

    
    double xN_end = vc_val_x2->at(size-1);
    double yN_end = vc_val_y2->at(size-1);
    double zN_end = vc_val_z2->at(size-1);
    
    for( int i=size-1; i>=1; i-- ) {

      double x0 = vc_val_x1->at(line_2ends);
      double y0 = vc_val_y1->at(line_2ends);
      double z0 = vc_val_z1->at(line_2ends);
      
      double x1 = vc_val_x1->at(i);
      double y1 = vc_val_y1->at(i);
      double z1 = vc_val_z1->at(i);

      double x2 = vc_val_x2->at(i);
      double y2 = vc_val_y2->at(i);
      double z2 = vc_val_z2->at(i);            

      double val_dE = vc_val_dE->at(i);
      double val_dQ = vc_val_dQ->at(i);
      
      double dx = sqrt( pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2) );
      user_rr_along += dx;
      
      //////////////////////////////
      
      user_dx_along += dx;
      user_dE_along += val_dE;
      user_dQ_along += val_dQ;
      
      if( user_dx_along > const_dx ) {
	
	double dEdx = user_dE_along/user_dx_along;
	double dQdx = user_dQ_along/user_dx_along;
	
	// double mid_x = (x0 + x2)/2;
	// double mid_y = (y0 + y2)/2;
	// double mid_z = (z0 + z2)/2;
	// double rr = sqrt( pow(mid_x-xN_end,2) + pow(mid_y-yN_end,2) + pow(mid_z-zN_end,2) );

	double rr = user_rr_along;
	int index_mid = (line_along+i)/2;

	for(int j=i; j<index_mid; j++) {
	        
	  double x11 = vc_val_x1->at(j);
	  double y11 = vc_val_y1->at(j);
	  double z11 = vc_val_z1->at(j);

	  double x22 = vc_val_x2->at(j);
	  double y22 = vc_val_y2->at(j);
	  double z22 = vc_val_z2->at(j);            

	  double dxx = sqrt( pow(x11-x22,2) + pow(y11-y22,2) + pow(z11-z22,2) );
	  rr = (rr-dxx);
	}
	
	h2_dQdx_rr_along_tot->Fill(rr, dQdx);

	h1_dx6mm_along->Fill( user_dx_along );
		
	if(p_id==muon) {

	  h1_dQdx_truth_along_dx6mm_muon->Fill( dQdx );
	  
	  double model_dQdx = func_dQdx_from_dEdx_by_ArgoNeut_model( dEdx, alpha_ArgoNeut, beta_ArgoNeut );
	  h1_dQdx_model_along_dx6mm_muon->Fill( model_dQdx );
	  h1_dQdx_truth2model_along_dx6mm_muon->Fill( dQdx/model_dQdx );

	  h2_dQdx_rr_along_muon->Fill( rr, dQdx );
	}

	if(p_id==muon || p_id==proton) {
	  h2_dQdx_rr_along_muons_protons->Fill( rr, dQdx );
	}
	
	
	//cout<<TString::Format(" %4d --> %4d , dx %7.4f", line_along, i, user_dx_along)<<endl;	
	///
	line_along = i-1;	
	user_dx_along = 0;
	user_dE_along = 0;
	user_dQ_along = 0;
      }// if( user_dx_along > const_dx )

      //////////////////////////////
      
      // user_dx_2ends = sqrt( pow(x0-x2,2) + pow(y0-y2,2) + pow(z0-z2,2) );
      // user_dE_2ends += val_dE;
      // user_dQ_2ends += val_dQ;
      
      // if( user_dx_2ends > const_dx ) {	
      // 	double dEdx = user_dE_2ends/user_dx_2ends;
      // 	double dQdx = user_dQ_2ends/user_dx_2ends;

      // 	double mid_x = (x0 + x2)/2;
      // 	double mid_y = (y0 + y2)/2;
      // 	double mid_z = (z0 + z2)/2;

      // 	double rr = sqrt( pow(mid_x-xN_end,2) + pow(mid_y-yN_end,2) + pow(mid_z-zN_end,2) );
	
      // 	h2_dQdx_rr_2ends_tot->Fill(rr, dQdx);

      // 	h1_dx6mm_2ends->Fill( user_dx_2ends );
	
      // 	if(p_id==muon) {
      // 	  h2_dQdx_rr_2ends_muon->Fill(rr, dQdx);
      // 	  h1_dQdx_truth_2ends_dx6mm_muon->Fill( dQdx );
	  
      // 	  double model_dQdx = func_dQdx_from_dEdx_by_ArgoNeut_model( dEdx, alpha_ArgoNeut, beta_ArgoNeut );
      // 	  h1_dQdx_model_2ends_dx6mm_muon->Fill( model_dQdx );
      // 	  h1_dQdx_truth2model_2ends_dx6mm_muon->Fill( dQdx/model_dQdx );
	  
      // 	}
	  
      // 	if(p_id==proton) {
      // 	  h2_dQdx_rr_2ends_proton->Fill(rr, dQdx);
      // 	}
	  
      // 	//cout<<TString::Format(" %4d --> %4d , dx %7.4f", line_2ends, i, user_dx_2ends)<<endl;	
      // 	///
      // 	line_2ends = i-1;	
      // 	user_dx_2ends = 0;
      // 	user_dE_2ends = 0;
      // 	user_dQ_2ends = 0;
      // }// if( user_dx_2ends > const_dx )

      
	
    }
    
    
  }// for(int ientry=0; ientry<=entries; ientry++)


  
  ///////////////////////////////////////////////////////////

  // https://root.cern.ch/root/html/tutorials/tree/hvector.C.html

  ///////////////////////////////////////////////////////////
  /*
  roostr = "canv_h1_dQdx_model2mc_muon";
  TCanvas *canv_h1_dQdx_model2mc_muon = new TCanvas(roostr, roostr, 900, 650);
  canv_h1_dQdx_model2mc_muon->SetLeftMargin(0.15);
  canv_h1_dQdx_model2mc_muon->SetRightMargin(0.15);
  canv_h1_dQdx_model2mc_muon->SetBottomMargin(0.15);

  roostr = "h1_dQdx_model2mc_muon";
  h1_dQdx_model2mc_muon->Draw();
  h1_dQdx_model2mc_muon->SetStats(0);
  h1_dQdx_model2mc_muon->SetTitle("");
  h1_dQdx_model2mc_muon->GetXaxis()->SetTitleSize(0.05);
  h1_dQdx_model2mc_muon->GetXaxis()->SetLabelSize(0.05);
  h1_dQdx_model2mc_muon->GetYaxis()->SetTitleSize(0.05);
  h1_dQdx_model2mc_muon->GetYaxis()->SetLabelSize(0.05);
  h1_dQdx_model2mc_muon->GetXaxis()->SetTitle("dQ/dx MC/model");
  h1_dQdx_model2mc_muon->GetYaxis()->SetTitle("entries");
  h1_dQdx_model2mc_muon->SetLineWidth(2);

  cout<<endl<<" Mean mc/model  "<< h1_dQdx_model2mc_muon->GetMean()<<endl<<endl;
  
  ///////////////////////////////////////////////////////////
 
  roostr = "canv_h1_dQdx_model2mc_proton";
  TCanvas *canv_h1_dQdx_model2mc_proton = new TCanvas(roostr, roostr, 900, 650);
  canv_h1_dQdx_model2mc_proton->SetLeftMargin(0.15);
  canv_h1_dQdx_model2mc_proton->SetRightMargin(0.15);
  canv_h1_dQdx_model2mc_proton->SetBottomMargin(0.15);

  roostr = "h1_dQdx_model2mc_proton";
  h1_dQdx_model2mc_proton->Draw();
  h1_dQdx_model2mc_proton->SetStats(0);
  h1_dQdx_model2mc_proton->SetTitle("");
  h1_dQdx_model2mc_proton->GetXaxis()->SetTitleSize(0.05);
  h1_dQdx_model2mc_proton->GetXaxis()->SetLabelSize(0.05);
  h1_dQdx_model2mc_proton->GetYaxis()->SetTitleSize(0.05);
  h1_dQdx_model2mc_proton->GetYaxis()->SetLabelSize(0.05);
  h1_dQdx_model2mc_proton->GetXaxis()->SetTitle("dQ/dx MC/model");
  h1_dQdx_model2mc_proton->GetYaxis()->SetTitle("entries");
  h1_dQdx_model2mc_proton->SetLineWidth(2);

  cout<<endl<<" Mean mc/model  "<< h1_dQdx_model2mc_proton->GetMean()<<endl<<endl;
  */
  ///////////////////////////////////////////////////////////
 
  roostr = "canv_h2_dQdx_rr_along_tot";
  TCanvas *canv_h2_dQdx_rr_along_tot = new TCanvas(roostr, roostr, 900, 650);
  canv_h2_dQdx_rr_along_tot->SetLeftMargin(0.15);
  canv_h2_dQdx_rr_along_tot->SetRightMargin(0.15);
  canv_h2_dQdx_rr_along_tot->SetBottomMargin(0.15);

  roostr = "h2_dQdx_rr_along_tot";
  h2_dQdx_rr_along_tot->Draw("colz");
  h2_dQdx_rr_along_tot->SetStats(0);
  h2_dQdx_rr_along_tot->SetTitle("");
  h2_dQdx_rr_along_tot->GetXaxis()->SetTitleSize(0.05);
  h2_dQdx_rr_along_tot->GetXaxis()->SetLabelSize(0.05);
  h2_dQdx_rr_along_tot->GetYaxis()->SetTitleSize(0.05);
  h2_dQdx_rr_along_tot->GetYaxis()->SetLabelSize(0.05);
  h2_dQdx_rr_along_tot->GetXaxis()->SetTitle("Residual range (cm)");
  h2_dQdx_rr_along_tot->GetYaxis()->SetTitle("dQ/dx (e/cm)");

  ///////////////////////////////////////////////////////////
 
  // roostr = "canv_h2_dQdx_rr_2ends_tot";
  // TCanvas *canv_h2_dQdx_rr_2ends_tot = new TCanvas(roostr, roostr, 900, 650);
  // canv_h2_dQdx_rr_2ends_tot->SetLeftMargin(0.15);
  // canv_h2_dQdx_rr_2ends_tot->SetRightMargin(0.15);
  // canv_h2_dQdx_rr_2ends_tot->SetBottomMargin(0.15);

  // roostr = "h2_dQdx_rr_2ends_tot";
  // h2_dQdx_rr_2ends_tot->Draw("colz");
  // h2_dQdx_rr_2ends_tot->SetStats(0);
  // h2_dQdx_rr_2ends_tot->SetTitle("");
  // h2_dQdx_rr_2ends_tot->GetXaxis()->SetTitleSize(0.05);
  // h2_dQdx_rr_2ends_tot->GetXaxis()->SetLabelSize(0.05);
  // h2_dQdx_rr_2ends_tot->GetYaxis()->SetTitleSize(0.05);
  // h2_dQdx_rr_2ends_tot->GetYaxis()->SetLabelSize(0.05);
  // h2_dQdx_rr_2ends_tot->GetXaxis()->SetTitle("Residual range (cm)");
  // h2_dQdx_rr_2ends_tot->GetYaxis()->SetTitle("dQ/dx (e/cm)");

  ///////////////////////////////////////////////////////////
 
  // roostr = "canv_h2_dQdx_rr_2ends_muon";
  // TCanvas *canv_h2_dQdx_rr_2ends_muon = new TCanvas(roostr, roostr, 900, 650);
  // canv_h2_dQdx_rr_2ends_muon->SetLeftMargin(0.15);
  // canv_h2_dQdx_rr_2ends_muon->SetRightMargin(0.15);
  // canv_h2_dQdx_rr_2ends_muon->SetBottomMargin(0.15);

  // roostr = "h2_dQdx_rr_2ends_muon";
  // h2_dQdx_rr_2ends_muon->Draw("colz");
  // h2_dQdx_rr_2ends_muon->SetStats(0);
  // h2_dQdx_rr_2ends_muon->SetTitle("");
  // h2_dQdx_rr_2ends_muon->GetXaxis()->SetTitleSize(0.05);
  // h2_dQdx_rr_2ends_muon->GetXaxis()->SetLabelSize(0.05);
  // h2_dQdx_rr_2ends_muon->GetYaxis()->SetTitleSize(0.05);
  // h2_dQdx_rr_2ends_muon->GetYaxis()->SetLabelSize(0.05);
  // h2_dQdx_rr_2ends_muon->GetXaxis()->SetTitle("Residual range (cm)");
  // h2_dQdx_rr_2ends_muon->GetYaxis()->SetTitle("dQ/dx (e/cm)");

  ///////////////////////////////////////////////////////////
 
  // roostr = "canv_h2_dQdx_rr_2ends_proton";
  // TCanvas *canv_h2_dQdx_rr_2ends_proton = new TCanvas(roostr, roostr, 900, 650);
  // canv_h2_dQdx_rr_2ends_proton->SetLeftMargin(0.15);
  // canv_h2_dQdx_rr_2ends_proton->SetRightMargin(0.15);
  // canv_h2_dQdx_rr_2ends_proton->SetBottomMargin(0.15);

  // roostr = "h2_dQdx_rr_2ends_proton";
  // h2_dQdx_rr_2ends_proton->Draw("colz");
  // h2_dQdx_rr_2ends_proton->SetStats(0);
  // h2_dQdx_rr_2ends_proton->SetTitle("");
  // h2_dQdx_rr_2ends_proton->GetXaxis()->SetTitleSize(0.05);
  // h2_dQdx_rr_2ends_proton->GetXaxis()->SetLabelSize(0.05);
  // h2_dQdx_rr_2ends_proton->GetYaxis()->SetTitleSize(0.05);
  // h2_dQdx_rr_2ends_proton->GetYaxis()->SetLabelSize(0.05);
  // h2_dQdx_rr_2ends_proton->GetXaxis()->SetTitle("Residual range (cm)");
  // h2_dQdx_rr_2ends_proton->GetYaxis()->SetTitle("dQ/dx (e/cm)");

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  // roostr = "canv_h1_dx6mm_2ends";
  // TCanvas *canv_h1_dx6mm_2ends = new TCanvas(roostr, roostr, 900, 650);
  // canv_h1_dx6mm_2ends->SetLeftMargin(0.15);
  // canv_h1_dx6mm_2ends->SetRightMargin(0.15);
  // canv_h1_dx6mm_2ends->SetBottomMargin(0.15);

  // roostr = "h1_dx6mm_2ends";
  // h1_dx6mm_2ends->Draw();
  // h1_dx6mm_2ends->SetStats(0);
  // h1_dx6mm_2ends->SetTitle("");
  // h1_dx6mm_2ends->GetXaxis()->SetTitleSize(0.05);
  // h1_dx6mm_2ends->GetXaxis()->SetLabelSize(0.05);
  // h1_dx6mm_2ends->GetYaxis()->SetTitleSize(0.05);
  // h1_dx6mm_2ends->GetYaxis()->SetLabelSize(0.05);
  // h1_dx6mm_2ends->GetXaxis()->SetTitle("dx (cm)");
  // h1_dx6mm_2ends->GetYaxis()->SetTitle("entries");
  // h1_dx6mm_2ends->SetLineWidth(2);

  ///////////////////////////////////////////////////////////
 
  // roostr = "canv_h1_dQdx_truth_2ends_dx6mm_muon";
  // TCanvas *canv_h1_dQdx_truth_2ends_dx6mm_muon = new TCanvas(roostr, roostr, 900, 650);
  // canv_h1_dQdx_truth_2ends_dx6mm_muon->SetLeftMargin(0.15);
  // canv_h1_dQdx_truth_2ends_dx6mm_muon->SetRightMargin(0.15);
  // canv_h1_dQdx_truth_2ends_dx6mm_muon->SetBottomMargin(0.15);

  // roostr = "h1_dQdx_truth_2ends_dx6mm_muon";
  // h1_dQdx_truth_2ends_dx6mm_muon->Draw();
  // h1_dQdx_truth_2ends_dx6mm_muon->SetStats(0);
  // h1_dQdx_truth_2ends_dx6mm_muon->SetTitle("");
  // h1_dQdx_truth_2ends_dx6mm_muon->GetXaxis()->SetTitleSize(0.05);
  // h1_dQdx_truth_2ends_dx6mm_muon->GetXaxis()->SetLabelSize(0.05);
  // h1_dQdx_truth_2ends_dx6mm_muon->GetYaxis()->SetTitleSize(0.05);
  // h1_dQdx_truth_2ends_dx6mm_muon->GetYaxis()->SetLabelSize(0.05);
  // h1_dQdx_truth_2ends_dx6mm_muon->GetXaxis()->SetTitle("dQ/dx (e/cm)");
  // h1_dQdx_truth_2ends_dx6mm_muon->GetYaxis()->SetTitle("entries");
  // h1_dQdx_truth_2ends_dx6mm_muon->SetLineWidth(2);
  // h1_dQdx_truth_2ends_dx6mm_muon->SetLineColor(kRed);
  
  // h1_dQdx_model_2ends_dx6mm_muon->Draw("same");
  // h1_dQdx_model_2ends_dx6mm_muon->SetLineWidth(2);
  // h1_dQdx_model_2ends_dx6mm_muon->SetLineColor(kBlue);
   
  ///////////////////////////////////////////////////////////
 
  // roostr = "canv_h1_dQdx_truth2model_2ends_dx6mm_muon";
  // TCanvas *canv_h1_dQdx_truth2model_2ends_dx6mm_muon = new TCanvas(roostr, roostr, 900, 650);
  // canv_h1_dQdx_truth2model_2ends_dx6mm_muon->SetLeftMargin(0.15);
  // canv_h1_dQdx_truth2model_2ends_dx6mm_muon->SetRightMargin(0.15);
  // canv_h1_dQdx_truth2model_2ends_dx6mm_muon->SetBottomMargin(0.15);

  // roostr = "h1_dQdx_truth2model_2ends_dx6mm_muon";
  // h1_dQdx_truth2model_2ends_dx6mm_muon->Draw();
  // h1_dQdx_truth2model_2ends_dx6mm_muon->SetStats(0);
  // h1_dQdx_truth2model_2ends_dx6mm_muon->SetTitle("");
  // h1_dQdx_truth2model_2ends_dx6mm_muon->GetXaxis()->SetTitleSize(0.05);
  // h1_dQdx_truth2model_2ends_dx6mm_muon->GetXaxis()->SetLabelSize(0.05);
  // h1_dQdx_truth2model_2ends_dx6mm_muon->GetYaxis()->SetTitleSize(0.05);
  // h1_dQdx_truth2model_2ends_dx6mm_muon->GetYaxis()->SetLabelSize(0.05);
  // h1_dQdx_truth2model_2ends_dx6mm_muon->GetXaxis()->SetTitle("MC/model dQdx");
  // h1_dQdx_truth2model_2ends_dx6mm_muon->GetYaxis()->SetTitle("entries");
  // h1_dQdx_truth2model_2ends_dx6mm_muon->SetLineWidth(2);
  // h1_dQdx_truth2model_2ends_dx6mm_muon->SetLineColor(kBlack);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  roostr = "canv_h1_dx6mm_along";
  TCanvas *canv_h1_dx6mm_along = new TCanvas(roostr, roostr, 900, 650);
  canv_h1_dx6mm_along->SetLeftMargin(0.15);
  canv_h1_dx6mm_along->SetRightMargin(0.15);
  canv_h1_dx6mm_along->SetBottomMargin(0.15);

  roostr = "h1_dx6mm_along";
  h1_dx6mm_along->Draw();
  h1_dx6mm_along->SetStats(0);
  h1_dx6mm_along->SetTitle("");
  h1_dx6mm_along->GetXaxis()->SetTitleSize(0.05);
  h1_dx6mm_along->GetXaxis()->SetLabelSize(0.05);
  h1_dx6mm_along->GetYaxis()->SetTitleSize(0.05);
  h1_dx6mm_along->GetYaxis()->SetLabelSize(0.05);
  h1_dx6mm_along->GetXaxis()->SetTitle("dx (cm)");
  h1_dx6mm_along->GetYaxis()->SetTitle("entries");
  h1_dx6mm_along->SetLineWidth(2);

  ///////////////////////////////////////////////////////////
 
  roostr = "canv_h1_dQdx_truth_along_dx6mm_muon";
  TCanvas *canv_h1_dQdx_truth_along_dx6mm_muon = new TCanvas(roostr, roostr, 900, 650);
  canv_h1_dQdx_truth_along_dx6mm_muon->SetLeftMargin(0.15);
  canv_h1_dQdx_truth_along_dx6mm_muon->SetRightMargin(0.15);
  canv_h1_dQdx_truth_along_dx6mm_muon->SetBottomMargin(0.15);

  roostr = "h1_dQdx_truth_along_dx6mm_muon";
  h1_dQdx_truth_along_dx6mm_muon->Draw();
  h1_dQdx_truth_along_dx6mm_muon->SetStats(0);
  h1_dQdx_truth_along_dx6mm_muon->SetTitle("");
  h1_dQdx_truth_along_dx6mm_muon->GetXaxis()->SetTitleSize(0.05);
  h1_dQdx_truth_along_dx6mm_muon->GetXaxis()->SetLabelSize(0.05);
  h1_dQdx_truth_along_dx6mm_muon->GetYaxis()->SetTitleSize(0.05);
  h1_dQdx_truth_along_dx6mm_muon->GetYaxis()->SetLabelSize(0.05);
  h1_dQdx_truth_along_dx6mm_muon->GetXaxis()->SetTitle("dQ/dx (e/cm)");
  h1_dQdx_truth_along_dx6mm_muon->GetYaxis()->SetTitle("entries");
  h1_dQdx_truth_along_dx6mm_muon->SetLineWidth(2);
  h1_dQdx_truth_along_dx6mm_muon->SetLineColor(kRed);
  
  h1_dQdx_model_along_dx6mm_muon->Draw("same");
  h1_dQdx_model_along_dx6mm_muon->SetLineWidth(2);
  h1_dQdx_model_along_dx6mm_muon->SetLineColor(kBlue);
   
  ///////////////////////////////////////////////////////////
 
  roostr = "canv_h1_dQdx_truth2model_along_dx6mm_muon";
  TCanvas *canv_h1_dQdx_truth2model_along_dx6mm_muon = new TCanvas(roostr, roostr, 900, 650);
  canv_h1_dQdx_truth2model_along_dx6mm_muon->SetLeftMargin(0.15);
  canv_h1_dQdx_truth2model_along_dx6mm_muon->SetRightMargin(0.15);
  canv_h1_dQdx_truth2model_along_dx6mm_muon->SetBottomMargin(0.15);

  roostr = "h1_dQdx_truth2model_along_dx6mm_muon";
  h1_dQdx_truth2model_along_dx6mm_muon->Draw();
  h1_dQdx_truth2model_along_dx6mm_muon->SetStats(0);
  h1_dQdx_truth2model_along_dx6mm_muon->SetTitle("");
  h1_dQdx_truth2model_along_dx6mm_muon->GetXaxis()->SetTitleSize(0.05);
  h1_dQdx_truth2model_along_dx6mm_muon->GetXaxis()->SetLabelSize(0.05);
  h1_dQdx_truth2model_along_dx6mm_muon->GetYaxis()->SetTitleSize(0.05);
  h1_dQdx_truth2model_along_dx6mm_muon->GetYaxis()->SetLabelSize(0.05);
  h1_dQdx_truth2model_along_dx6mm_muon->GetXaxis()->SetTitle("MC/model dQdx");
  h1_dQdx_truth2model_along_dx6mm_muon->GetYaxis()->SetTitle("entries");
  h1_dQdx_truth2model_along_dx6mm_muon->SetLineWidth(2);
  h1_dQdx_truth2model_along_dx6mm_muon->SetLineColor(kBlack);

  /////////////////////////////////////////////////////////////////////////////
  
  roostr = "canv_h2_dQdx_rr_along_muon";
  TCanvas *canv_h2_dQdx_rr_along_muon = new TCanvas(roostr, roostr, 900, 650);
  canv_h2_dQdx_rr_along_muon->SetLeftMargin(0.15);
  canv_h2_dQdx_rr_along_muon->SetRightMargin(0.15);
  canv_h2_dQdx_rr_along_muon->SetBottomMargin(0.15);
  roostr = "h2_dQdx_rr_along_muon";
  h2_dQdx_rr_along_muon->Draw("colz");
  h2_dQdx_rr_along_muon->SetStats(0);
  h2_dQdx_rr_along_muon->SetTitle("");
  h2_dQdx_rr_along_muon->GetXaxis()->SetTitleSize(0.05);
  h2_dQdx_rr_along_muon->GetXaxis()->SetLabelSize(0.05);
  h2_dQdx_rr_along_muon->GetYaxis()->SetTitleSize(0.05);
  h2_dQdx_rr_along_muon->GetYaxis()->SetLabelSize(0.05);
  h2_dQdx_rr_along_muon->GetYaxis()->SetTitle("dQ/dx (e/cm)");
  h2_dQdx_rr_along_muon->GetXaxis()->SetTitle("Residual range (cm)");
  h2_dQdx_rr_along_muon->GetYaxis()->SetTitleOffset(1.3);
  h2_dQdx_rr_along_muon->SetLineWidth(2);
  h2_dQdx_rr_along_muon->SetLineColor(kBlack);
  g_dQdx_model_muon->Draw("L same");
   

  roostr = "canv_h2_dQdx_rr_along_muons_protons";
  TCanvas *canv_h2_dQdx_rr_along_muons_protons = new TCanvas(roostr, roostr, 900, 650);
  canv_h2_dQdx_rr_along_muons_protons->SetLeftMargin(0.15);
  canv_h2_dQdx_rr_along_muons_protons->SetRightMargin(0.15);
  canv_h2_dQdx_rr_along_muons_protons->SetBottomMargin(0.15);
  roostr = "h2_dQdx_rr_along_muons_protons";
  h2_dQdx_rr_along_muons_protons->Draw("colz");
  h2_dQdx_rr_along_muons_protons->SetStats(0);
  h2_dQdx_rr_along_muons_protons->SetTitle("");
  h2_dQdx_rr_along_muons_protons->GetXaxis()->SetTitleSize(0.05);
  h2_dQdx_rr_along_muons_protons->GetXaxis()->SetLabelSize(0.05);
  h2_dQdx_rr_along_muons_protons->GetYaxis()->SetTitleSize(0.05);
  h2_dQdx_rr_along_muons_protons->GetYaxis()->SetLabelSize(0.05);
  h2_dQdx_rr_along_muons_protons->GetYaxis()->SetTitle("dQ/dx (e/cm)");
  h2_dQdx_rr_along_muons_protons->GetXaxis()->SetTitle("Residual range (cm)");
  h2_dQdx_rr_along_muons_protons->GetYaxis()->SetTitleOffset(1.3);
  h2_dQdx_rr_along_muons_protons->SetLineWidth(2);
  h2_dQdx_rr_along_muons_protons->SetLineColor(kBlack);
  g_dQdx_model_muon->Draw("L same");
  g_dQdx_model_proton->Draw("L same");
   
  
}

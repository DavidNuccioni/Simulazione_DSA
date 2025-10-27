// VERYTOP
void SetHistoStyle(TH2D* hh);
void SetHistoStyle(TH1D* hh);
void SetGRID(double* xEKN, int NEKN, double EKN1, double EKN2);
void MakeSpectrum(TH1D* hh);

// TOP
void DrawDSASpectrum(){

  bool kSPEC= true;
  
  // ---- ENERGY GRID ----
  const int NEKN= 50; //  N di bins
  double EKNMIN= 8.0; 
  double EKNMAX= 8.e+5;
  double xEKN[NEKN+1]; // punti griglia = n bin +1  
  SetGRID(xEKN, NEKN+1, EKNMIN, EKNMAX); // crea griglia log-uniform

  TH1D* hDSASpectrum_TOT= new TH1D("hDSASpectrum_TOT","DSA Energ Spectrum",NEKN,xEKN);  
  SetHistoStyle(hDSASpectrum_TOT);
  hDSASpectrum_TOT->SetLineColor(kBlack);
  hDSASpectrum_TOT->SetLineWidth(3);

  
  TH1D* hDSASpectrum_US= new TH1D("hDSASpectrum_US","DSA Energy Spectrum UPSTREAM",NEKN,xEKN);
  SetHistoStyle(hDSASpectrum_US);
  hDSASpectrum_US->SetLineColor(kRed+1);


  TH1D* hDSASpectrum_DS= new TH1D("hDSASpectrum_DS","DSA Energ Spectrum DOWNSTREAM",NEKN,xEKN);  
  SetHistoStyle(hDSASpectrum_DS);
  hDSASpectrum_DS->SetLineColor(kAzure+1);
  

  // ---- Pick up data ----
  TFile* hDsaFile= new TFile("../outroot/hDSATimeDep_2025_Epoch_Total.root","READ");
  hDsaFile->cd();


  // ==== Define Tree / Like the one we are reading from file =====
  Double_t pX;  // posizione wrt shock, in unitÃ  di pc
  Double_t pE;  // kinetic energy per nucleon, in GeV/n
  Double_t pD;  // diff coeff of the particle at the end of dsa cm2/s pc2/yr
  Double_t pL;  // diff mean free path pc

  Double_t pT;  // time spent by particle in doing DSA , yr
  Double_t pT1; // start time (SNR age at injection) 
  Double_t pT2; // end time (end of DSA), t.c. pT=pT2-pT1
  Int_t    nT;  // sequential number of epoch [1 - NTIME]
  
  Int_t iSide; // wich side wrt shock: US/DS +1/-1
  Int_t iStat; // status flag: 0,1,2,3

  Int_t  nCross; // number of shock crosses / n of gains
  Int_t  nSteps; // number of steps of that particle / track length in steps 
  Int_t  nTrack; // sequential number of particles IN THAT EPOCHS (we all tracks)
  Int_t  nPart;  // sequential number of particles TOTAL 
    
  
  TTree* PartTree=  (TTree*)hDsaFile->Get("PartTree");

  PartTree->SetBranchAddress("pX",&pX);
  PartTree->SetBranchAddress("pE",&pE);
  PartTree->SetBranchAddress("nT",&nT);
  PartTree->SetBranchAddress("pT",&pT);
  PartTree->SetBranchAddress("pT1",&pT1);
  PartTree->SetBranchAddress("pT2",&pT2);
  PartTree->SetBranchAddress("pD",&pD);
  PartTree->SetBranchAddress("pL",&pL);
  PartTree->SetBranchAddress("iSide",&iSide);
  PartTree->SetBranchAddress("iStat",&iStat);
  PartTree->SetBranchAddress("nCross",&nCross);
  PartTree->SetBranchAddress("nStesp",&nSteps);  
  PartTree->SetBranchAddress("nTrack",&nTrack);
  PartTree->SetBranchAddress("nPart",&nPart);

  

  // ---- scan tree -----
  float pperc=0.;
  float perc;
  Int_t nentries= (Int_t)PartTree->GetEntries();
  

  // ---- event loop ----
  printf("scan on Particle Tree: %d entries\n",nentries);
  for(int pp=0;pp<nentries;pp++){
    perc=pp/(nentries*1.);
    if (perc>=pperc){printf("Processed %5.0f%%\n",pperc*100);pperc+=0.1;}  
    PartTree->GetEntry(pp);

    // -----------
    // NT- volendo possimao introdurre un peso, per il fatto che lo shock quando si ingrandisce
    // spazza via piu particelle, e quindi il termine di sorgente dovrebbe essere proporzionale
    // al raggio dello shock al quadrato Rsh^2. Se l'espansione Ã¨ lineare nel tempo (ejecta dominated phase),
    // allora Rsh proporzionale all'epoca nT. E dunque double Weight = nT*nT;
    // Si veda il paper Kachelriess & Ostapchenko 2011
    // -----------

    // definizione peso
    // double pWeight = nT*nT; // peso, leggi sopra

    
    // fill histo
    // hDSASpectrum_TOT->Fill(pE,Weight);  // con peso
    hDSASpectrum_TOT->Fill(pE); // senza peso

    
    if(iSide== 1)hDSASpectrum_US->Fill(pE); // upstream
    if(iSide==-1)hDSASpectrum_DS->Fill(pE); // downstream
    
  }
  // ---- end of loop scan ----

  // dividi per larghezza del bin, t.c. l'istrogramma dei counts E diventa uno spettro dN/dE 
  // l'asse Y non Ã¨ piÃ¹ "n entries", ma Ã¨ "dN/dE".
  if(kSPEC){
    MakeSpectrum(hDSASpectrum_TOT);
    MakeSpectrum(hDSASpectrum_US);
    MakeSpectrum(hDSASpectrum_DS);
  }

  
  // fit a legge di potenza
  TF1* fPlaw= new TF1("fPlaw","[0]*TMath::Power(x,[1])");
  fPlaw->SetParameter(0,54000);
  fPlaw->SetParameter(1,-2);
  fPlaw->SetLineColor(kGreen+1);

  hDSASpectrum_TOT->Fit(fPlaw,"","",10,10000);
  
  // ---- plot results ----

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  // ---- legend ----
  TLegend* legDsa;
  if(!kSPEC) legDsa = new TLegend(0.45,0.20,0.85,0.38);
  if(kSPEC) legDsa = new TLegend(0.65,0.70,0.90,0.88);
  legDsa->SetBorderSize(0);
  legDsa->SetTextSize(0.040);
  legDsa->SetTextFont(42);
  legDsa->SetFillColor(0);
  legDsa->AddEntry(hDSASpectrum_TOT,"Total","l");
  legDsa->AddEntry(hDSASpectrum_US,"Upstream","l");
  legDsa->AddEntry(hDSASpectrum_DS,"Downstream","l");

  
  TCanvas* ccDSASpectrum= new TCanvas("ccDSASpectrum","DSA Energy Spectrum",1200,800);
  ccDSASpectrum->cd();
  gPad->SetBottomMargin(0.16);
  gPad->SetTopMargin(0.06);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetTicky(1);
  gPad->SetTickx(1);
  
  hDSASpectrum_TOT->Draw("hist");
  legDsa->Draw();
  hDSASpectrum_TOT->Draw("hist same");
  hDSASpectrum_DS->Draw("hist same");
  hDSASpectrum_US->Draw("hist same");


  
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->RedrawAxis();
}

// =================================

// BOTTOM

// create log-uniform kinetic energy grid
 void SetGRID(double* xEKN, int NEKN, double EKN1, double EKN2){ 
   double deltalogz = TMath::Log10(EKN2)-TMath::Log10(EKN1);
   double dx = deltalogz/(NEKN-1);
   double l10 = TMath::Log(10);
   for(int ee=0;ee<NEKN;ee++){
     xEKN[ee] = EKN1*TMath::Exp(l10*ee*dx);
   }
 }

// ========================


// stile per isotgrammi 2D / frames
void SetHistoStyle(TH2D* hh){
  hh->SetTitle(0);
  hh->GetYaxis()->SetNdivisions(504);
  hh->SetLabelFont(42,"X");
  hh->SetLabelFont(42,"Y");
  hh->SetTitleFont(42,"X");
  hh->SetTitleFont(42,"Y");
  hh->GetXaxis()->SetTitle("kinetic energy (GeV)");
  hh->GetYaxis()->SetTitle("n of entries");
  hh->GetXaxis()->SetLabelSize(0.06);
  hh->GetYaxis()->SetLabelSize(0.06);
  hh->GetXaxis()->SetTitleSize(0.06);
  hh->GetYaxis()->SetTitleSize(0.06);
  hh->GetXaxis()->SetTitleOffset(1.375);
  hh->GetYaxis()->SetTitleOffset(1.25);
}


// stile per isotgrammi 1D
void SetHistoStyle(TH1D* hh){
  hh->SetTitle(0);
  hh->GetYaxis()->SetNdivisions(504);
  hh->SetLabelFont(42,"X");
  hh->SetLabelFont(42,"Y");
  hh->SetTitleFont(42,"X");
  hh->SetTitleFont(42,"Y");
  hh->GetXaxis()->SetTitle("kinetic energy (GeV)");
  hh->GetYaxis()->SetTitle("n of entries");
  hh->GetXaxis()->SetLabelSize(0.06);
  hh->GetYaxis()->SetLabelSize(0.06);
  hh->GetXaxis()->SetTitleSize(0.06);
  hh->GetYaxis()->SetTitleSize(0.06);

  hh->GetXaxis()->SetLabelOffset(-0.005);
  hh->GetYaxis()->SetLabelOffset(0.0);

  hh->GetXaxis()->SetTitleOffset(1.2);
  hh->GetYaxis()->SetTitleOffset(1.0);
  hh->SetLineWidth(2);
}


void MakeSpectrum(TH1D* hh){

  int NB= hh->GetXaxis()->GetNbins(); 
  for(int ee=0;ee<NB;ee++){
    double bval = hh->GetBinContent(ee+1);
    double bwid = hh->GetBinWidth(ee+1);
    double nval= bval/bwid; // new value: divided by bin width
    double nerr= sqrt(bval)/bwid; // rescale poissonian error
    hh->SetBinContent(ee+1, nval);    
    hh->SetBinError(ee+1, nerr);    
  }

}

// VERYBOTTOM

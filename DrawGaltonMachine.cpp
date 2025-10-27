void SetStyle(); // plot style, definito alla fine

// MACCHINA DI GALTON, OVVERO RANDOM-WALK 1D

// Descrive la distribuzione di Np palline dopo aver
// attraversato Nl livelli della macchina di Galton.
// Ad ogni livello la pallina urta, puÃ² andare a SX o DX con p=1/2
//
// Il problema Ã¨ equivalente al randomwalk 1D visto, in cui si
// tira una moneta e si va a dx o sx in base a testa o croce


void DrawGaltonMachine(){

  int Np = 10000; // N totale di palline lanciate (o monete)
  int Nl = 25;   // N tot di livelli/urti (o lanci) di ogni pallina

  double X=0;    // Posizione (inizialmente X=0).
  double Sign=0; // segno +/- 1 (direzione, o testa/croce)
  double dX=1;   // Ampiezza passo, unitaria (il passo Ã¨ +/-dX)
  int L = 0;     // N. progressivo del livello, da 1 a Nl
  int P = 0;     // N. progressivo della particella, da 1 a Np

  // definisci stile
  SetStyle();

 
  // tree con tutte le posizioni intermedie/traiettorie
  TTree *TraceTree= new TTree("TraceTree","TraceTreefon_Tree");  
  TraceTree->Branch("X",&X,"X/D");
  TraceTree->Branch("L",&L,"L/I"); // indice livello, tra 1 e Np
  TraceTree->Branch("P",&P,"P/I"); // indice pallina, tra 1 e Np 
  TraceTree->Branch("Nl",&Nl,"Nl/I"); // costante ma meglio averlo nel tree
  TraceTree->Branch("Np",&Np,"Np/I");

  TraceTree->Branch("Sign",&Sign,"Sign/D");
  TraceTree->SetLineColor(kAzure+7);
  TraceTree->SetLineWidth(4);

  // oggetto random Marsenne-Twistor generator
  TRandom3* rand= new TRandom3();
  
  for(int pp=0;pp<Np;pp++){ // loop palline (o monete)

    X= 0;    // posizione iniziale
    L= 0;    // livello iniziale
    P= pp+1; // numero progressivo pallina, da 1 a Np

    TraceTree->Fill(); // fill configurazione iniziale
    
    for(int ll=0;ll<Nl;ll++){ // loop livelli/urti (o lanci)

      // estrazione random del segno +/- (testa/croce)
      double randnum = rand->Uniform(-1.0, 1.0);
      Sign= randnum/TMath::Abs(randnum);

      // incremento del passo
      X += Sign*dX;
      L = ll+1; // numero progressivo livello, da 1 a Nl
 
      TraceTree->Fill();
     
    }
  }

  printf("\n");

  // mostra possibili traiettorie nei vari livelli
  TCanvas* ccXvsL= new TCanvas("ccXvsL","Percorsi possibili",700,700);
  ccXvsL->cd();
  TraceTree->Draw("-L:X","","l");

  // mostra la gaussiana alla fine.
  TCanvas* ccPDF= new TCanvas("ccPDF","distribuzione finale",700,700);
  ccPDF->cd();
  TraceTree->Draw("X","L==Nl","hist");

  /*
  TraceTree->SetLineColor(kRed+1);
  TraceTree->Draw("X","L==10","same");
  */  
}



// ==== SET GRAPHICS STYLE FOR PLOTS ====
void SetStyle(){

  TStyle *myStyle  = new TStyle("MyStyle","My Root Style");   
  int FontType=42;
  
  // General
  myStyle->SetFrameBorderMode(kFALSE);
  
  //Canvas 
  myStyle->SetCanvasBorderMode(kFALSE);
  myStyle->SetCanvasColor(kWhite);
  
  //Pad
  myStyle->SetPadColor(kWhite);
  myStyle->SetPadTickX(kTRUE);
  myStyle->SetPadTickY(kTRUE);

  myStyle->SetPadTopMargin(0.06);
  myStyle->SetPadRightMargin(0.06);
  myStyle->SetPadBottomMargin(0.15);
  myStyle->SetPadLeftMargin(0.15);
  
  //Title
  myStyle->SetOptTitle(kFALSE);
  myStyle->SetTitleColor(kBlack);
  myStyle->SetTitleBorderSize(0);
  myStyle->SetTitleX(0.25);
  myStyle->SetTitleY(0.98);
  myStyle->SetTitleOffset(1.25,"x");
  myStyle->SetTitleOffset(1.25,"y");
  myStyle->SetTitleSize(0.05,"xyz");
  myStyle->SetTitleFont(FontType,"xyz");

  myStyle->SetLabelOffset(0.01,"x");
  myStyle->SetLabelOffset(0.01,"y");

  //Stat
  myStyle->SetOptStat(kFALSE);
  myStyle->SetStatColor(kWhite);
  myStyle->SetStatFont(FontType);
  
  //Legend
  myStyle->SetLegendBorderSize(0);
  
  //Label   
  myStyle->SetLabelFont(FontType,"xyz");
  myStyle->SetLabelSize(0.050,"xyz");
  myStyle->SetLabelColor(kBlack,"xyz");
  
  //Text
  myStyle->SetTextFont(FontType);
  myStyle->SetTextSize(1.2);
  
  //Fit
  myStyle->SetOptFit(kFALSE);
  
  // Histograms & Plots
  myStyle->SetMarkerStyle(20);
  myStyle->SetHistLineWidth(2);
  myStyle->SetErrorX(0.001); // no error bar caps 
  myStyle->SetPalette(1);    
  myStyle->SetPaintTextFormat("5.2f"); // precision in "TEXT" plots
  
  gROOT->SetStyle("MyStyle");
  gROOT->ForceStyle();

  }
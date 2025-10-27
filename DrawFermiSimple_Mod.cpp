// NT2023 very simple simulation of Fermi acceleration
// No random walk. No particle tracing. 
// We just use escape probability P, energy gain Xi, initial energy E0. 

// Every particle starts with energy E0. At each acceleration cycle, the particle has
// probability P of escape the accelerator, and probability 1-P of remaining around and crossing the shock.
// If the particle crosses the shock, its energy increases by a factor Xi. Thus, E1= E0*Xi; E2=E1*Xi, etc.
// If it crosses the shock N time (each with probability 1-P) the final energy is En= E0*Xi^N.


void DrawFermiSimple_Mod(){

  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  // --- initial parameters ----
  double E0 = 1.00;  // initial kinetic energy in GeV
  double Xi = 1.15;  // energy gain Xi=DeltaE/E
  double P  = 0.05;  // escape probability at each cycle

  // Dal meccanismo di Fermi, se P è piccola (P<<1), dovremmo ottenere una distribuzione
  // power-law del tipo dN/dE = E^-gamma, con indice dato da gamma=1+P/Xi (Vedi ad es. tesi Settimo).
  // Lo verifichiamo con la simulazione che segue.

  // --- definizione generatore random ---
  TRandom3* rand= new TRandom3();
  rand->SetSeed(204); // inizializzalizzazione il generatore con numero a caso
  
  int NPART= 50000;    // numero di particelle che vogliamo simulare
  int NCMAX= 1000000;  // numero massimo cicli di accelerazione, per ogni particella;

  
  Double_t E; // energia della particella (inizialmente si ha E=E0, poi verrà inctrementata);
  Int_t    N; // contatore del numero di cicli di accelerazione (cioè di shock crossing)


  // --- definizione istogramma delle energie ---
  // range di energia [E1, E2], ma usiamo i logaritmi. 

  double E1= 1.e+0; // 10 GeV
  double E2= 1.e+9; // 10^9 GeV
  double LogE1 = TMath::Log10( E1 ); // 1 GeV
  double LogE2 = TMath::Log10( E2 ); // 10^9 GeV

  int NofBins = 50; // numero di bin
  TH1F* hLogE = new TH1F("hLogE","Istogramma di Log(E)", NofBins, LogE1, LogE2);
  hLogE->GetXaxis()->SetTitle("log_{10}(E/GeV)");
  hLogE->GetYaxis()->SetTitle("N. entries");
  hLogE->Sumw2(); // calcola e mostra errori poissoniani.
  
  // NT> l'istogramma così definito ha 50 bins ed è equispaziato nella variabile LogE.
  // NT> Per riempire l'istogramma useremo quindi la variabile LogE.

  
  // --- NT2025 nuovo istogramma log-uniformly spaced ----

  // step1: definiamo griglia di energia
  const int NEKN= NofBins; // 50; //  N di bins
  double EKN1= E1;
  double EKN2= E2;
  double xEKN[NEKN+1];     // punti griglia = n bin +1 
  double deltalogz = TMath::Log10(EKN2)-TMath::Log10(EKN1);
  double dx = deltalogz/(NEKN-1);
  double l10 = TMath::Log(10);
  for(int ee=0;ee<NEKN+1;ee++) xEKN[ee] = EKN1*TMath::Exp(l10*ee*dx);
   
  // step2: definiamo l'istogramma di E
  TH1F* hE = new TH1F("hE","Istogramma di E", NEKN, xEKN);
  hE->GetXaxis()->SetTitle("E (GeV)");
  hE->GetYaxis()->SetTitle("N. entries");
  hE->Sumw2(); // calcola e mostra errori poissoniani.
  
  // step2: definiamo l'istogramma di E da cui calcoleremo dN/dE=Spectrum
  TH1F* hSpectrum = new TH1F("hSpectrum","Distribuzione di E, dN/dE", NEKN, xEKN);
  hSpectrum->GetXaxis()->SetTitle("E (GeV)");
  hSpectrum->GetYaxis()->SetTitle("dN/dE");
  hSpectrum->Sumw2(); // calcola e mostra errori poissoniani.

  // -----------------------------------------------------
  
  // --- tree contenente energie finali ----
  TTree* PartTree= new TTree("PartTree","PartTree_Tree"); 
  PartTree->Branch("E",&E,"E/D"); //
  PartTree->Branch("N",&N,"N/I"); // numero cicli
  // NT> Questo tree non lo usiamo veramente. Ma lo teniamo perché può risultare utile.


  
  // --- loop su particelle ----
  for(int pp=0;pp<NPART;pp++){

    
    // definizione condizioni iniziali di ogni particella
    E = E0; // energia della particella (iniziale)
    N = 0;  // numero di cicli di accelerazione subiti (inizialmente 0)
     

    // --- loop sui cicli di accelerazione ---
    for(int kk=0;kk<NCMAX;kk++){
      
      // Ad ogni ciclo di accelerazione kk, la particella ha probabilità P di fuggire dall'acceleratore
      // (e quindi la sua accelerazione termina con l'energia che ha raggiunto)
      // e probabilità 1-P di rimanere nell'accelerazione (e quindi la sua energia continua ad essere incrementanta).
      // Per determinare se la particella fugge o rimane, generiamo un numero random tra 0 e 1 e lo confrontiamo con P:

      
      // esrazione di un numero reale random tra 0 e 1
      double Rand= rand->Uniform(0.0, 1.0);
      
      // Se Rand<P allora c'è fuga, e quindi l'accelerazione termina qui.
      if( Rand < P ) break; // FUGA. Finisce qui.
      
      // Se Rand>P, allora non c'è fuga: c'è attraversamento e il loop continua 

      
      E *= Xi; // incremento energia di un fattore Xi: E_(n) = E_(n-1) * Xi
      N++;     // conteggio cicli di accelerazione
    }


    // registrazione istogramma di LogE
    double LogE= TMath::Log10( E ); 
    hLogE->Fill( LogE );

    hE->Fill( E );
    hSpectrum->Fill( E );
    
    // registrazione entry del tree
    PartTree->Fill();

  }

  // --- draw tree / test ---
  // PartTree->SetLineWidth(2);
  // PartTree->Draw("TMath::Log10(E)","E>0");

  // --- In hSpectrum, dividi per la larghezza del bin. Così da ottenere uno spetro dN/dE  ----
  for(int ee=0;ee<NEKN;ee++){
    double bval = hSpectrum->GetBinContent(ee+1);
    double bwid = hSpectrum->GetBinWidth(ee+1);
    double nval= bval/bwid; // new value: divided by bin width
    double nerr= sqrt(bval)/bwid; // rescale poissonian error
    hSpectrum->SetBinContent(ee+1, nval);    
    hSpectrum->SetBinError(ee+1, nerr);    
  }

  
  // === draw histogram ===

  // definizione stile istogramma
  hLogE->SetLineWidth(1);
  hLogE->SetMarkerStyle(20);
  hLogE->SetMarkerSize(0.8);
  hLogE->SetMarkerColor(kBlue+1);
  hLogE->SetLineColor(kBlue+1);

  hE->SetLineWidth(1);
  hE->SetMarkerStyle(20);
  hE->SetMarkerSize(0.8);
  hE->SetMarkerColor(kRed+1);
  hE->SetLineColor(kRed+1);

  hSpectrum->SetLineWidth(1);
  hSpectrum->SetMarkerStyle(22);
  hSpectrum->SetMarkerSize(0.8);
  hSpectrum->SetMarkerColor(kAzure+1);
  hSpectrum->SetLineColor(kAzure+1);

  
  // definizione di una canvas in cui disegnare l'istogramma dello spettro

  // istogramma di LogE
  TCanvas* ccLogE= new TCanvas("ccLogE","Istrograma di LogE",1000,600);
  ccLogE->cd();
  hLogE->Draw();
  gPad->SetLogy();

  
  // istogramma di E  
  TCanvas* ccE= new TCanvas("ccE","Istogramma di E",1000,600);
  ccE->cd();
  hE->Draw();
  gPad->SetLogx();
  gPad->SetLogy();

  // fit a legge di potenza dell'istogramma di E
  TF1* fPlaw1= new TF1("fPlaw1","[0]*TMath::Power(x,[1])");
  fPlaw1->SetParameter(1,-2); // initialize
  fPlaw1->SetParName(0,"Normalization"); 
  fPlaw1->SetParName(1,"Slope"); 
  fPlaw1->SetLineColor(kRed+1);
  cout<<"Fit dell'istogramma di E:"<<endl;
  hE->Fit(fPlaw1,"","",10,1.e+8);

  // distribuzione dN/dE  
  TCanvas* ccSpectrum= new TCanvas("ccSpectrum","Spettro energie dN/dE",1000,600);
  ccSpectrum->cd();
  hE->Draw();
  gPad->SetLogx();
  gPad->SetLogy();

  // fit a legge di potenza dell'istogramma di E
  TF1* fPlaw2= new TF1("fPlaw2","[0]*TMath::Power(x,[1])");
  fPlaw2->SetParameter(1,-2); // initialize
  fPlaw2->SetParName(0,"Normalization"); 
  fPlaw2->SetParName(1,"Slope"); 
  fPlaw2->SetLineColor(kAzure+1);
  cout<<"Fit dello spettro dN/dE:"<<endl;
  hSpectrum->Fit(fPlaw2,"","",10,1.e+8);


  
}



// void MakeSpectrum(TH1D* hh){
//   int NB= hh->GetXaxis()->GetNbins(); 
//   for(int ee=0;ee<NB;ee++){
//     double bval = hh->GetBinContent(ee+1);
//     double bwid = hh->GetBinWidth(ee+1);
//     double nval= bval/bwid; // new value: divided by bin width
//     double nerr= sqrt(bval)/bwid; // rescale poissonian error
//     hh->SetBinContent(ee+1, nval);    
//     hh->SetBinError(ee+1, nerr);    
//   }
// }


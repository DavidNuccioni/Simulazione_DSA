// NT2023 very simple simulation of Fermi acceleration
// No random walk. No particle tracing. 
// We just use escape probability P, energy gain Xi, initial energy E0. 

// Every particle starts with energy E0. At each acceleration cycle, the particle has
// probability P of escape the accelerator, and probability 1-P of remaining around and crossing the shock.
// If the particle crosses the shock, its energy increases by a factor Xi. Thus, E1= E0*Xi; E2=E1*Xi, etc.
// If it crosses the shock N time (each with probability 1-P) the final energy is En= E0*Xi^N.


void DrawFermiSimple(){

  // --- initial parameters ----
  double E0 = 1.00;  // initial kinetic energy in GeV
  double Xi = 1.15;  // energy gain Xi=DeltaE/E
  double P  = 0.05;  // escape probability at each cycle

  // Dal meccanismo di Fermi, se P Ã¨ piccola (P<<1), dovremmo ottenere una distribuzione
  // power-law del tipo dN/dE = E^-gamma, con indice dato da gamma=1+P/Xi (Vedi ad es. tesi Settimo).
  // Lo verifichiamo con la simulazione che segue.

  // --- definizione generatore random ---
  TRandom3* rand= new TRandom3();
  rand->SetSeed(204); // inizializzalizzazione il generatore con numero a caso
  
  int NPART= 50000;    // numero di particelle che vogliamo simulare
  int NCMAX= 1000000;  // numero massimo cicli di accelerazione, per ogni particella;

  
  Double_t E; // energia della particella (inizialmente si ha E=E0, poi verrÃ  inctrementata);
  Int_t    N; // contatore del numero di cicli di accelerazione (cioÃ¨ di shock crossing)


  // --- definizione istogramma delle energie ---
  // range di energia [E1, E2], ma usiamo i logaritmi. 
  double LogE1 = TMath::Log10( 1.e+0 ); // 1 GeV
  double LogE2 = TMath::Log10( 1.e+9 ); // 10^9 GeV
  int NofBins = 50; // numero di bin
  TH1F* hLogE = new TH1F("hLogE","distribuzione di Log(E)", NofBins, LogE1, LogE2);
  hLogE->GetXaxis()->SetTitle("log_{10}(E/GeV)");
  hLogE->GetYaxis()->SetTitle("dN/dE");
  hLogE->Sumw2(); // calcola e mostra errori poissoniani.
  
  // NT> l'istogramma cosÃ¬ definito ha 50 bins ed Ã¨ equispaziato nella variabile LogE.
  // NT> Per riempire l'istogramma useremo quindi la variabile LogE.

    // --- tree contenente energie finali ----
  TTree* PartTree= new TTree("PartTree","PartTree_Tree"); 
  PartTree->Branch("E",&E,"E/D"); //
  PartTree->Branch("N",&N,"N/I"); // numero cicli
  // NT> Questo tree non lo usiamo veramente. Ma lo teniamo perchÃ© puÃ² risultare utile.


  
  // --- loop su particelle ----
  for(int pp=0;pp<NPART;pp++){

    
    // definizione condizioni iniziali di ogni particella
    E = E0; // energia della particella (iniziale)
    N = 0;  // numero di cicli di accelerazione subiti (inizialmente 0)
     

    // --- loop sui cicli di accelerazione ---
    for(int kk=0;kk<NCMAX;kk++){
      
      // Ad ogni ciclo di accelerazione kk, la particella ha probabilitÃ  P di fuggire dall'acceleratore
      // (e quindi la sua accelerazione termina con l'energia che ha raggiunto)
      // e probabilitÃ  1-P di rimanere nell'accelerazione (e quindi la sua energia continua ad essere incrementanta).
      // Per determinare se la particella fugge o rimane, generiamo un numero random tra 0 e 1 e lo confrontiamo con P:

      
      // esrazione di un numero reale random tra 0 e 1
      double Rand= rand->Uniform(0.0, 1.0);
      
      // Se Rand<P allora c'Ã¨ fuga, e quindi l'accelerazione termina qui.
      if( Rand < P ) break; // FUGA. Finisce qui.
      
      // Se Rand>P, allora non c'Ã¨ fuga: c'Ã¨ attraversamento e il loop continua 

      
      E *= Xi; // incremento energia di un fattore Xi: E_(n) = E_(n-1) * Xi
      N++;     // conteggio cicli di accelerazione
    }


    // registrazione istogramma di LogE
    double LogE= TMath::Log10( E ); 
    hLogE->Fill( LogE );

    // registrazione entry del tree
    PartTree->Fill();

  }


  // --- draw tree / test ---
  // PartTree->SetLineWidth(2);
  // PartTree->Draw("TMath::Log10(E)","E>0");

  
  // === draw histogram ===

  // definizione stile istogramma
  hLogE->SetLineWidth(1);
  hLogE->SetMarkerStyle(20);
  hLogE->SetMarkerSize(0.8);
  hLogE->SetMarkerColor(kBlue+1);
  hLogE->SetLineColor(kBlue+1);

  // definizione di una canvas in cui disegnare l'istogramma dello spettro
  TCanvas* ccLogE= new TCanvas("ccLogE","Spettro energie dN/dE",1000,600);
  ccLogE->cd();
  hLogE->Draw();
  gPad->SetLogy();


  
}

  

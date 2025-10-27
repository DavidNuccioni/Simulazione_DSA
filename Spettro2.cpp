void Spettro2(){

// Variabili della simulazione

int Npart = 50000000;		// Numero di particelle
int NCmax = 1000;		// Numero di attraversamenti massimi
int Ncross;			// Numero di attraversamenti
double E;			// Energia della particella
double E_0 = 1.;		// Energia iniziale della particella
double dE = 1.03;		// Guadagno energetico
double P = 0.026;		// Probabilità di fuga

// Numero random

TRandom3* Nrand= new TRandom3();
Nrand->SetSeed(NCmax);

// Tree dei dati NT- IL TREE NON SERVE IN QUESTA MACRO
TTree* Dati = new TTree("Dati","Dati"); 	
Dati->Branch("E",&E,"E/D"); 
Dati->Branch("Ncross",&Ncross,"Ncross/I");

//Istogramma delle energie (NT- DI LOG-ENERGIA=
double X_min = TMath::Log10( 1.e+0 ); 
double X_max = TMath::Log10( 1.e+9 );  

// NT- VECCHIO ISTOGRAMMA IN LOGE - NON SERVE
 TH1F* dist = new TH1F("dist","", 50, X_min, X_max);
dist->GetXaxis()->SetTitle("log_{10}(E/GeV)");
dist->GetYaxis()->SetTitle("Spettro differenziale dN/dE");
dist->GetYaxis()->SetTitleOffset(0.85);
dist->GetXaxis()->SetTitleOffset(1.00);
dist->GetXaxis()->SetTitleSize(0.04);
dist->GetYaxis()->SetTitleSize(0.04);
dist->Sumw2();
dist->SetLineWidth(1);
dist->SetMarkerStyle(20);
dist->SetMarkerSize(0.8);
dist->SetMarkerColor(kBlack);
dist->SetLineColor(kBlack);


// NT- NUOVO ISTOGRAMMA DELLO SPETTRO DI ENERGIA  


  // step1: definiamo griglia di energia log-uniforme
  const int NEKN= 50; //  N di bins
  double EKN1= 1.e+0;
  double EKN2= 1.e+8;
  double xEKN[NEKN+1];  // griglia di energia. N. punti griglia = n bins +1 

  double deltalogz = TMath::Log10(EKN2)-TMath::Log10(EKN1);
  double dx = deltalogz/(NEKN-1);
  double l10 = TMath::Log(10);
  for(int ee=0;ee<NEKN+1;ee++) xEKN[ee] = EKN1*TMath::Exp(l10*ee*dx);

  // step2: definiamo l'istogramma di E da cui calcoleremo dN/dE=Spectrum
  TH1F* hSpectrum = new TH1F("hSpectrum","Distribuzione di E, dN/dE", NEKN, xEKN);
  hSpectrum->GetXaxis()->SetTitle("E (GeV)");
  hSpectrum->GetYaxis()->SetTitle("dN/dE");
  hSpectrum->Sumw2(); // calcola e mostra errori poissoniani.


  // NT- la divisione per la larghezza del bin va fatta DOPO.  

  // loop eventi
for(int i=0; i<Npart; i++){		// Sceglie una particella e ne inizializza le variabili

	E = E_0;
	Ncross = 0;
	
	for(int j=0; NCmax; j++){	// Inizia il meccanismo di accelerazione 
	
		double Prand = Nrand->Uniform(0.0, 1.0);
		if(Prand<P){break;} 
		else{ E = E*dE; Ncross++;}
		}
		
	double LogE= TMath::Log10( E ); // Il processo è concluso e si salvano i dati
    	dist->Fill( LogE);

	hSpectrum->Fill( E ); // NT- riempiamo energie
	
	Dati->Fill();
    	}

// Scalatura grafico SI VA MESSA QUI, DOPO IL LOOP. MA USA IL GRAFICO IN CUI L'ASSE X è L'ENERGIA


int NB= hSpectrum->GetXaxis()->GetNbins(); 
for(int ee=0;ee<NB;ee++){
  double bval = hSpectrum->GetBinContent(ee+1);
  double bwid = hSpectrum->GetBinWidth(ee+1);
  double nval= bval/bwid; // new value: divided by bin width
  double nerr= sqrt(bval)/bwid; // rescale poissonian error
  hSpectrum->SetBinContent(ee+1, nval);    
  hSpectrum->SetBinError(ee+1, nerr);    
  }

// Creazione del grafico

// NT- stile plain
 gROOT->SetStyle("Plain");
 gStyle->SetOptStat(0);

 
TCanvas* Spettro= new TCanvas("Spettro","Spettro energie dN/dE",1000,600);
Spettro->cd();
//dist->Draw();
 hSpectrum->Draw();
 gPad->SetLogy();
 gPad->SetLogx(); // NT- ANCHE LOGX ORA

 Dati->SetLineColor(kBlack);// NT- perché cambi il colore nel tree Dati? non lo usi
gPad->SetRightMargin(0.05);
gPad->SetTopMargin(0.05);
gPad->SetLeftMargin(0.12);
gPad->SetBottomMargin(0.12);



  // NT- fit a legge di potenza dell'istogramma di E
  TF1* fPlaw1= new TF1("fPlaw1","[0]*TMath::Power(x,[1])");
  fPlaw1->SetParameter(1,-2); // initialize
  fPlaw1->SetParName(0,"Normalization"); 
  fPlaw1->SetParName(1,"Slope"); 
  fPlaw1->SetLineColor(kRed+1);
  cout<<"Fit dell'istogramma di E:"<<endl;
  hSpectrum->Fit(fPlaw1,"","",10,1.e+8);

 
return;
}

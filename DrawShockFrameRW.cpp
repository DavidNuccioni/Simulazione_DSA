void DrawShockFrameRW(){

  // RWSIM.
  // Random walk simulation 1D in the Shock reference frame (SF= shock frame). 


  // ------- constants -----------------------

  const double pc2cm    = 3.08568025e+18; // parsec to cm
  const double year2sec = 31556926.0;     // year to second
  const double c_lig    = 29979245800.0;  // speed of light cm/s

  const double _t_min = 1.*year2sec;     // minimum time considered: 1 year
  const double _t_max = 10000.*year2sec; // maximum time allowed: per una SNR, ordine di 10-15k years

  const double _rcomp    = 4.0;      // compression ratio u_us/u_ds [NB 3.5->2.2]
  const double _n_gas    = 2.0;      // cas density upstream (non ci serve) cm-3
  const double _u_shock  = 0.5e+8;   // speed of the shock cm/s
  const double _b_fie    = 1.;       // magnetic field in microgauss (muG)  

  // NB the b_field is used to calculate the diffusion coefficient D
  // We will assume D proportional to R/B (rigidity/B_field), ie to the gyroradius.
  
  // fluid speed in both sides of the shock: upstream and downstream
  double _u_us= _u_shock;   
  double _u_ds= _u_shock/_rcomp;

  
  // ---------------------------------------


  int Z = 1;             // atomic  number (charge) 
  int A = 1;             // mass number (n of nucleons)
  double Ekn0  = 1.0;   // initial kinetic energy, GeV/n
  double Time0 =_t_min;  // initial time, sec

  int NTRACKS = 1;         // n particles to track
  int NSTEPMAX= 4000000;   // max number of steps

  // ------- particle track kinematic --------------
  // queste variabili verranno aggiornate ad ogni step

  double  X = 0.;     // posizione 1D
  double dX = 0.;     // incremento
  double Ekn = Ekn0;  // energia cinetica
  double Time= Time0; // tempo

  int IsDS = 0;     // IsDownStream? =1 se la particella Ã¨ downstream, 0 se Ã¨ upstream
  int UDS  = 0;     // In che lato Ã¨ rispetto allo shock [-1: downstream. +1: upstream] (insomma, come sopra)
  int IsCross = 0;  // =1 se la particella ha attraversato lo shock.
  int NCross = 0;   // Numero di attraversamenti fatti fino ad ora

  double XMin=0; // valore minimo di X (quindi sarÃ  negativo, US) registrato alla fine della simulazione
  double XMax=0; // valore massimo di X (quindi sarÃ  positivo, DS) registrato alla fine della simulazione
      
  int StatusFlag = 0; // mi dice perchÃ© l'accelerazione Ã¨ terminata.
  int NSTEPS= 0;      // numero di steps 
  int NTRACK= 0;      // numero di tracce (o di particelle)

  TTree *TraceTree= new TTree("TraceTree","TraceTreefon_Tree");  
  TraceTree->Branch("X",&X,"X/D");
  TraceTree->Branch("dX",&dX,"dX/D");
  TraceTree->Branch("Ekn",&Ekn,"Ekn/D");
  TraceTree->Branch("Time",&Time,"Time/D");
  TraceTree->Branch("IsCross",&IsCross,"IsCross/I");
  TraceTree->Branch("NCross",&NCross,"NCross/I");
  TraceTree->Branch("UDS",&UDS,"UDS/I");

  TraceTree->Branch("NSTEPS",&NSTEPS,"NSTEPS/I");
  TraceTree->Branch("NTRACK",&NTRACK,"NTRACK/I");
  TraceTree->Branch("StatusFlag",&StatusFlag,"StatusFlag/I");

  // -------------------------------------------

  // oggetto random  
  TRandom3* rand= new TRandom3();

  // loop su NTRACK, ma abbiamo solo una traccia mi sa
  for(int nn=0; nn<NTRACKS;nn++){

    NTRACK++;

    rand->SetSeed(nn+1);

    StatusFlag = 0;
    NSTEPS= 0;

  // init track
  X = 0.;
  dX = 0.;
  Ekn = Ekn0;
  Time= Time0;
  IsDS = 0;
  IsCross = 0;
  NCross = 0;

  
  // loop sui vari step della particella
  for(int ii=0;ii<NSTEPMAX;ii++){

    NSTEPS++;

    // check where particle is: downstream or upstream?
    if(X>0){IsDS = 0;}else{IsDS = 1;}

    double Etot = Ekn*A; // energia cinetica totale
    double Ptot = Etot;  // impulso
    double Rig  = Ptot/Z; // rigiditÃ 

    // calcolo del coefficiente di diffusione D~ Rig/BField
    double D0= (3.29772e+22)*Rig/_b_fie; // bohm diff [cm2/s]
    double L0= 3.*D0/c_lig; // diffusion mean free path (o diffusion lenth)

    // NB per la diffusion length abbiamo usato approssimazione ultra-relativistica
    
    // double L0old= (330.e+10)*Rig/_b_fie; // ok if B is [muG]
    // double Test= L0 - (330.e+10)*Rig/_b_fie; // ok if B is [muG]
    // cout<< "L0new: "<<L0 <<"   L0old: "<< L0old << "   Diff: "<<L0 - L0old << endl;

    L0 /= TMath::Sqrt(3.); // correggo L0 perchÃ© siamo in 1D! (la formula di prima Ã¨ per il caso 3D)


    // timestep - ur approx
    double dT = L0/c_lig; // tempo impiegato per muoversi di uno step pari a L0
    Time += dT; // incremento il tempo.

    
    // ==== simulazione diffusione : random step ====

    // Ho due opzioni (UserÃ² la prima delle due):
    // Prima opzione: genero un numero +1/-1 che Ã¨ come il lancio di una moneta, com probabilitÃ  50%/50%
    // Seconda opzione: genero una variabile continua gaussiana avente media 0 e larghezza 1.

    // prima opzione:
    double randnum = rand->Uniform(-1.0, 1.0); // generazione numero random.
    double sign= randnum/TMath::Abs(randnum); // estrazione del segno (+/-) [segno di un numero = numero/|numero|]

    // seconda opzione (non usata)
    // double sign = rand->Gaus(0.0, 1.0);     // SMEARED SIGN (numero random gaussiano)

    // definizione dello step casuale dXdiff
    double dXdiff = sign*L0;


    // === simulazione di convection/advection === 

    double Vconv; // same sign us and ds! ? 
    if(IsDS==1){ Vconv = -_u_ds; } 
    if(IsDS==0){ Vconv = -_u_us; }
    double dXconv = Vconv*dT; // aggiungiamo uno step nella direzione del fluido.

    // === calcolo dello step totale = diffusione + convezione ====
    dX = dXdiff + dXconv;
    X = X + dX; // aggiornamento posizione

    
    // === aggiornamento valori minimi e massimi di X nel corso della simulazione ===
    if(X<XMin) XMin=X; // aggiornamento XMin
    if(X>XMax) XMax=X; // aggiornamento XMax
    

    // cout<< "[STEP "<< ii << "]  [TRACK " << nn << "] "<<"  dXdiff: "<<dXdiff<< "   dXconv: "<< dXconv <<"  X : "<< X << endl;

    // === check crossing ===
    // dobbiamo capire se c'Ã¨ stato un attraversamento dello shock. 
    // cioÃ¨, ora che abbiamo la nuova posizione Xnew= Xold+dX, ci chiediamo se Xnew si trova nella stessa parte dello shock di Xold.
    
    int oldIsDS = IsDS;   // old IsDS
    if(X>0){IsDS = 0;}else{IsDS = 1;} // new IsDS
    if(oldIsDS != IsDS){IsCross = 1;}else{IsCross=0;} // changed?
    if(IsDS==1){UDS=-1;}else{UDS=+1;} // UDS non aggiunge nulla rispetto a IsDS, ma torna utile per i plots.

    // se c'Ã¨ stato attraversamento, allora c'Ã¨ energy gain
    if(IsCross == 1){
      NCross++;
      Ekn += (_u_shock/c_lig)*Ekn;  // Correct in the ShockFrame?
    }

    // === fine simulazione di questa traccia ===

    // la simulazione termina quando StatusFlag viene settata diversa da 0.
    // CioÃ¨ quando Ã¨ passato un tempo molto grande (Time> _t_max) [criterio fisicamente sensato naturale]. 
    // Oppure quando il numero di step ha raggiunto il valore massimo (ii >= NSTEPMAX) [meno sensato ma utile per la CPU]
    
    if(Time > 250.*year2sec){ StatusFlag=5.;} //  (ragginto  END OF EJ PHASE

    if(Time > _t_max){ StatusFlag = 1; } // end ST stage
    if(ii>=NSTEPMAX-1){ StatusFlag = 2;} // max steps reached 

    // fill tree
    TraceTree->Fill();

    if(StatusFlag!=0){ break; }

  }
  // fine del loop dei vari steps


  
  printf("\n");

    cout<< "ENERGY:  "<< Ekn << endl;
    cout<< "CROSSES: "<< NCross << endl;
    cout<< "T/TMAX:  "<< Time/_t_max << endl;
    cout<< "NSTEPS:  "<< ((double)NSTEPS/(double)NSTEPMAX)*100. <<" perc"<<endl;
    cout<< "Is DS?: "<< IsDS << endl;
    cout<< "STATUS:  "<< StatusFlag << endl;  
    cout<< endl;
}
  // fine del loop delle varie tracce/particelle (solo UNA).


  
  // ------- DRAW / FACCIAMO DEI PLOTS ------------
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TCanvas* ccShockFrameTrajectory= new TCanvas("ccShockFrameTrajectory","ccShockFrameTrajectory",500,650);
  ccShockFrameTrajectory->Divide(1,3,0.01,0.01);

  // range temporale per i grafici
  double TMIN = _t_min/year2sec; // estremo inferiore
  double TMAX = Time/year2sec; // come estremo TMAX prendiamo "Time", il tempo effettivo della simulazione. 
  //TMAX= _t_max/100./year2sec; //TMAX Ã¨ spesso inferiore di _t_max, cioÃ¨ la simulazione puÃ² finire se NSteps diventa grande...

  // === frame del grafico X vs time ===
  double XMIN = 1.3*XMin/pc2cm; // estremi per il plot
  double XMAX = 1.3*XMax/pc2cm;
  TH2F* hXFrame= new TH2F("hXFrame","hXFrame",200,TMIN,TMAX, 200,XMIN, XMAX);
  hXFrame->SetTitle(0);
  hXFrame->GetXaxis()->SetTitle("Time (yrs)");
  hXFrame->GetYaxis()->SetNdivisions(506);
  hXFrame->GetYaxis()->SetTickLength(0.02);
  hXFrame->SetLabelFont(52,"X");
  hXFrame->SetLabelFont(52,"Y");
  hXFrame->SetTitleFont(52,"X");
  hXFrame->SetTitleFont(52,"Y");
  hXFrame->GetYaxis()->SetTitle("X (pc)");
  hXFrame->GetYaxis()->SetTitleOffset(0.75);
  hXFrame->GetXaxis()->SetTitleOffset(1.00);
  hXFrame->GetYaxis()->SetLabelOffset(0.005);
  hXFrame->GetXaxis()->SetLabelOffset(+0.005);
  hXFrame->GetXaxis()->SetLabelSize(0.06);
  hXFrame->GetYaxis()->SetLabelSize(0.06);
  hXFrame->GetXaxis()->SetTitleSize(0.06);
  hXFrame->GetYaxis()->SetTitleSize(0.06);

  // frame Ekn (l'asse y Ã¨ 'energia cinetica e lo aggiustiamo in base a Ekn0 e Ekn raggiunta)
  double EMIN = Ekn0/2.;
  double EMAX = Ekn*2.;
  TH2F* hEknFrame= new TH2F("hEknFrame","hEknFrame",200,TMIN,TMAX, 200, EMIN, EMAX);
  hEknFrame->SetTitle(0);
  hEknFrame->GetXaxis()->SetTitle("Time (yrs)");
  hEknFrame->GetYaxis()->SetNdivisions(506);
  hEknFrame->GetYaxis()->SetTickLength(0.02);
  hEknFrame->SetLabelFont(52,"X");
  hEknFrame->SetLabelFont(52,"Y");
  hEknFrame->SetTitleFont(52,"X");
  hEknFrame->SetTitleFont(52,"Y");
  hEknFrame->GetYaxis()->SetTitle("E (GeV/n)");
  hEknFrame->GetYaxis()->SetTitleOffset(0.75);
  hEknFrame->GetXaxis()->SetTitleOffset(1.00);
  hEknFrame->GetYaxis()->SetLabelOffset(0.005);
  hEknFrame->GetXaxis()->SetLabelOffset(+0.005);
  hEknFrame->GetXaxis()->SetLabelSize(0.06);
  hEknFrame->GetYaxis()->SetLabelSize(0.06);
  hEknFrame->GetXaxis()->SetTitleSize(0.06);
  hEknFrame->GetYaxis()->SetTitleSize(0.06);

  // frame Cross
  double NCMIN = 1.; 
  double NCMAX = NCross*1.5; // max number of NCross per il grafico
  TH2F* hCrossFrame= new TH2F("hCrossFrame","hCrossFrame",200,TMIN, TMAX, 200, NCMIN, NCMAX);
  hCrossFrame->SetTitle(0);
  hCrossFrame->GetXaxis()->SetTitle("Time (yrs)");
  hCrossFrame->GetYaxis()->SetNdivisions(506);
  hCrossFrame->GetYaxis()->SetTickLength(0.02);
  hCrossFrame->SetLabelFont(52,"X");
  hCrossFrame->SetLabelFont(52,"Y");
  hCrossFrame->SetTitleFont(52,"X");
  hCrossFrame->SetTitleFont(52,"Y");
  hCrossFrame->GetYaxis()->SetTitle("N Crosses");
  hCrossFrame->GetYaxis()->SetTitleOffset(0.75);
  hCrossFrame->GetXaxis()->SetTitleOffset(1.00);
  hCrossFrame->GetYaxis()->SetLabelOffset(0.005);
  hCrossFrame->GetXaxis()->SetLabelOffset(+0.005);
  hCrossFrame->GetXaxis()->SetLabelSize(0.06);
  hCrossFrame->GetYaxis()->SetLabelSize(0.06);
  hCrossFrame->GetXaxis()->SetTitleSize(0.06);
  hCrossFrame->GetYaxis()->SetTitleSize(0.06);


  ccShockFrameTrajectory->cd(1);
  TraceTree->SetLineColor(kBlack);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  hXFrame->Draw();
  gPad->SetLogx();

  // NT- la posizione X Ã¨ in cm. Divido per 3.08568025e+18 per esprimerlo in parsec.
  // Il tempo Time Ã¨ in secondi. Divido per 31556926 per esprimerlo in anni. 
  TraceTree->Draw("X/(3.08568025e+18):Time/31556926.","NTRACK==1 && (abs(X)>=0.0)","l same");

  TLine* line= new TLine(TMIN,0.0,TMAX,0.0); // linea blu orizzontale che indica lo shock a X=0.
  line->SetLineColor(kBlue+1);
  line->SetLineStyle(9);
  line->Draw();

  TraceTree->SetLineColor(kYellow);
  TraceTree->Draw("UDS:Time/31556926.","NTRACK==1","l same"); // disegna linea gialla in corrispondenza dell'attraversamento
  // infatti UDS passa dal valore -1 al valore +1 quando la particella passa da DS a US
  // gPad->SetGridy();
  TraceTree->SetLineColor(kBlack);
  TraceTree->Draw("X/(3.08568025e+18):Time/31556926.","NTRACK==1 && (abs(X)>=0.0)","l same");



  ccShockFrameTrajectory->cd(2);
  TraceTree->SetLineColor(kBlack);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);

  hEknFrame->Draw();
  gPad->SetLogx();
  gPad->SetLogy();

  TraceTree->Draw("Ekn:Time/31556926.","NTRACK==1","l same");
  TraceTree->SetLineColor(kYellow);
  TraceTree->Draw("1000*UDS:Time/31556926.","NTRACK==1","l same");
  //gPad->SetGridy();
  TraceTree->SetLineColor(kBlack);
  TraceTree->Draw("Ekn:Time/31556926.","NTRACK==1","l same"); 
  TraceTree->Draw("1.01*Ekn:Time/31556926.","NTRACK==1","l same");


  ccShockFrameTrajectory->cd(3);
  TraceTree->SetLineColor(kBlack);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);

  hCrossFrame->Draw();
  gPad->SetLogx();
  gPad->SetLogy();

  TraceTree->Draw("NCross:Time/31556926.","NTRACK==1","l same");
  TraceTree->Draw("1.01*NCross:Time/31556926.","NTRACK==1","l same");
  //gPad->SetGridy();
}


  // we obtain that <X> prop to sqrt(N)
  // <X> = K*sqrt(N) with K=sqrt(2/pi)!!
  

  // TraceTree->Draw("X:N");
  // cout<<"NSteps: " <<NSTEPMAX<< endl;
  // cout<<"dX:     " <<dX<< endl;
  // cout<<"dX*N:   "<< dX*NSTEPMAX<< endl;
  // cout<<"sqrt:   "<< sqrt(dX*NSTEPMAX)<< endl;
  // cout<<"--->X:  "<< X << endl;
    
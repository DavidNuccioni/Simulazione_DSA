#include "TROOT.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
using namespace std;

//void DrawShockFrameRW(int NParticles, int NEpochs, int FirstEpoch=0);
void DrawShockFrameRW(int NParticles, int NEpochs);
void DrawShockFrameRW(int NPARTICLES, int NEpochs, int WhichEpoch);


int main(int argc, char** argv){

  // --- make a loop over all epochs ----
  int NParticles= 1000;  
  int NEpochs   = 150; 
  // DrawShockFrameRW(NParticles, NEpochs);

  // ---- run over one selected epoch ----
  int ThisEpoch = 0; //  from 1 to NEpoch. Put zero to loop over all epochs
  if( argc==2 ) ThisEpoch= (int)atoi(argv[1]);
  DrawShockFrameRW(NParticles, NEpochs, ThisEpoch);

  return 0;
}




void DrawShockFrameRW(int NParticles, int NEpochs){
  DrawShockFrameRW(NParticles, NEpochs, 0);
}
 
void DrawShockFrameRW(int NPARTICLES, int NEpochs, int WhichEpoch){
  
  
  // Random walk simulation 1D in the Shock reference frame
  // Vedere ad es. Kachelriess et al. 2011 ApJ "ANTIMATTER PRODUCTION IN SUPERNOVA REMNANTS"
  // https://iopscience.iop.org/article/10.1088/0004-637X/733/2/119/pdf

  // ogni particella termina l'accelerazione con una flag iStat
  // iStat= 1 AGE - Reached ST stage
  // iStat= 2 ESCAPE - Particle escape 
  // iStat= 3 CPU - NSTEPMAX reached 
  
  int NSTEPMAX = 1e+8; // 700 000 000;  // Use N>500K for E0=15. Use 700K for E0=10. Use N>1M for E0=1.
  bool kSMEAR  =  false;      // GAUSSIAN SMEARING IN STEP
  double Mp    = 0.93827229;  // nucleon mass

  
  // ------- constants -----------------------
  const double pc2cm = 3.08568025e+18;
  const double year2sec= 31556926.;
  const double c_lig= 29979245800.;   // cm/s
  const double _t_max = 15000.*year2sec; // SNR/ST age [Truelove & McKee 1994]

  const double _rcomp    = 4; //3.5;    // u_us/u_ds [NB 3.5->2.2]
  const double _u_shock  = 1.e+8; //0.5e+8;   // cm/2
  const double _b_fie    = 1.; // [muG]  
  double _u_us= _u_shock;  // nello shock frame: uus=Ush 
  double _u_ds= _u_shock/_rcomp; // uds= Ush/r
  
  // ---------------------------------------

  int Z = 1;
  int A = 1; 
  double Ekn0=10.; // GeV/n INJECTION ENERGY
  double Time0=0; //_t_max - _t_max/10.; // FIRST EPOCH

  // ------- particle kinematic --------------

  // working quantities: cm and s units
  double  X = 0.;
  double dX = 0.;
  double Ekn = Ekn0;

  double Tdsa= 0;      // DSA duration [sec]
  double Time1=Time0; // start time (SNR age at injection) [sec]
  double Time2=Time0; // end time (SNR age at escape) [sec]



  int IsDS = 0;
  // int UDS  = 0;
  int IsCross = 0;
  int StatusFlag  = 0; // particle status flag

  ///  int NofCross = 0; // sequential cross number USIAMO nCross DOPO
  ///  int NofSteps= 0; // sequential step number USIAMO nSteps DOPO

  double D0= 0; // diffusion coeff cm2/s
  double L0= 0; // diffusion mean free path cm

  
  // ==== NT2025 ====  Escape e criteri per fermare la simulazione
  // 1 Il tempo raggiungee TMAX= T_ST, cioÃ¨ l'etÃ  max in cui la SNR accelera.
  // 2 Escape: la particella Ã¨ cosÃ¬ lontana che nel tempo rimasto (T-TMAX) non tornerÃ  allo shock.
  // Siccome L ~ sqrt(D*T), il criterio Ã¨ che siamo a |X|>>sqrt[D*(T-TMAX)], meglio se DS!!
  // 3 Si raggiunge il n di steps massimi della simulazione
   

  // ------ time grid -----
  const int NTIMES = NEpochs; // n of epochs
  double _T1 = 0.0;
  double _T2  = _t_max; // max time
  double xTime[NTIMES];
  for(int tt=0;tt<NTIMES;tt++){
    xTime[tt] = _T1 + tt*(_T2-_T1)/(NTIMES);
  }

  

  // ------ CR Particle File and Tree ------------------
  TFile* PartFile= new TFile(Form("hDSATimeDep_2025_Epoch_%d.root",WhichEpoch), "RECREATE" );

  // ==== variabili per tree NT2025 =====
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
    
  
  TTree* PartTree= new TTree("PartTree","Part_Tree");
  PartTree->Branch("pX",&pX,"pX/D");
  PartTree->Branch("pE",&pE,"pE/D");
  PartTree->Branch("nT",&nT,"nT/I");
  PartTree->Branch("pT",&pT,"pT/D");
  PartTree->Branch("pT1",&pT1,"pT1/D");
  PartTree->Branch("pT2",&pT2,"pT2/D");
  PartTree->Branch("pD",&pD,"pD/D");
  PartTree->Branch("pL",&pL,"pL/D");

  PartTree->Branch("iSide",&iSide,"iSide/I");
  PartTree->Branch("iStat",&iStat,"iStat/I");

  PartTree->Branch("nCross",&nCross,"nCross/I");
  PartTree->Branch("nStesp",&nSteps,"nSteps/I");
  
  PartTree->Branch("nTrack",&nTrack,"nTrack/I");
  PartTree->Branch("nPart",&nPart,"nPart/I");

  // -------------------------------------------

  TRandom3* rand= new TRandom3();
  rand->SetSeed(1); 
  // NB se abbiamo N tracks e M time0, usare sia N che M per seed

  nT   = 0;
  nPart= 0; // sequential n of particle / regardless the epoch 


  if(WhichEpoch>0){ nPart = NPARTICLES*(WhichEpoch-1); }

  // --- time loop / over difeent epochs of SNR life ----
  for(int tt=0;tt<NTIMES;tt++){

    // --- se l'epoca Ã¨ specificata, si runnna solo quella ---
    if( (WhichEpoch > 0) && (WhichEpoch != tt+1) ) continue;

    Time1 = xTime[tt]; // start time
    Time2 = Time1;       // to be incremented by Tdsa later
    
    nTrack=0; // sequential n of track = n of particle in that epoch

    nT= tt+1; // conta epoche
    
    // --- particle / track  loop ---    
    for(int pp=0; pp<NPARTICLES;pp++){
      // rand->SetSeed(pp+1);
      
      nTrack++; // cioÃ¨ nPart=pp / PER ORA
      nPart++;  // part incrementa a tutte le epoche e non si annulla
      /// if(WhichEpoch>1){ nPart = pp+1 + NPARTICLES*WhichEpoch; } // NT: NO, NON QUI!
      
      iStat = 0;
      nSteps= 0;

      // --- init track ---
      X    = 0.;
      dX   = 0.;
      Ekn  = Ekn0; // ekn of injection
      IsDS = 0;
      IsCross = 0;
      
      nCross  = 0;
      StatusFlag= 0;

      Tdsa = 0; // time spent in dsa / acc time of the track
      
      
      D0 = 0;
      L0 = 0;
       
    // --- steps inside the track loop ----
    for(int ii=0;ii<NSTEPMAX;ii++){

      nSteps++; 

      // determine side 
      if(X<0){IsDS = 1;} else{IsDS = 0;}
      if(X<0){iSide = -1;} else{iSide = 1;}
      
      // calc rigidity HE. use UR formula Etot ~ Ptot ~ Ekn*A ~ Rig*Z, for Ekn>100
      double Rig= (A/Z)*Ekn;

      // calc rigidity LE. use standard relativistic formula for Ekn<100
      if(Ekn<100) Rig= (A/Z)*sqrt( Ekn*Ekn + 2.*Ekn*Mp ); // proton mass
      
      // calc diffusion coefficient cm2/s
      D0= 3.29772e+22*Rig/_b_fie; // bohm diff coeff D0 [cm2/s] 

      // diffusion length cm
      L0= 3.*D0/c_lig; // mean free path
      // if(Ekn<100) L0 /=  (Z/A)*(Z/A)*Rig/sqrt( (Z/A)*(Z/A)*(Rig*Rig) + Mp*Mp ); // use speed beta*c_lig 

      //L0 /= TMath::Sqrt(3.); // ??

      
      // *** NB from above seems that L0 is propto E instead going as sqrt(E)   ***
      // *** but the time-step is also dT=L0/c. So, we recover L0 ~ sqrt(D0*dT) ***      
      // double L0 = TMath::Sqrt(D0*dT);      
      // double L0old= (330.e+10)*Rig/_b_fie; // ok if B is [muG]
      // double Test= L0 - (330.e+10)*Rig/_b_fie; // ok if B is [muG]
      // cout<< "L0new: "<<L0 <<"   L0old: "<< L0old << "   Diff: "<<L0 - L0old << endl;
      
 
      // === timestep - ur approx ===
      double dT = L0/c_lig;  
      Tdsa += dT;

    
      // === random diffusion step ===
      double randnum = rand->Uniform(-1.0, 1.0);
      double sign= randnum/TMath::Abs(randnum); // SIMPLE RANDOM SIGN
      if(kSMEAR) sign *= rand->Gaus(0.0, 1.0); // GAUSS SMEARED SIGN
      double dXdiff = sign*L0;

    
      // === convection step ====
      double Vconv; // same sign us and ds! ? [sÃ¬, nel shock frame sono entrambe verso sinistra]
      if(IsDS==1){ Vconv = -_u_ds; } 
      if(IsDS==0){ Vconv = -_u_us; }
      double dXconv = Vconv*dT;
      
    
      // === increment position ===
      dX = dXdiff + dXconv;
      X = X + dX;
      

      // === check if crossing occurs / there is a corss if the side has changed ===
      int oldIsDS = IsDS;   // old IsDS
      if(X>0){IsDS = 0;}else{IsDS = 1;} // new IsDS
      if(oldIsDS != IsDS){IsCross = 1;}else{IsCross=0;} // CROSS CONDITION
      //// if(IsDS==1){UDS=1;}else{UDS=-1;} // NT2011 DS/US=+/- !!! NT2025 non la usiamo
      
      if(IsCross == 1){
	nCross++;
	Ekn += (_u_shock/c_lig)*Ekn;  // CORRECT IN SF? NT2025 rivedere se l'incremento Ã¨ corretto [pare di sÃ¬]
      }


       
    // === end of track simulation: three conditions ====

      // remaining time
      Time2 = Time1 + Tdsa;
      double TimeLeft = _t_max - Time2; // dove Time2=Time1+Tdsa =SNR age. CioÃ¨ TimeLeft = _t_max -Time1-Tdsa.
      
      // criterio 1: TIME. DSA ends because the SNR reaches of ST age.
      if(TimeLeft <= 0){ StatusFlag = 1; } 
      // if(Tdsa > _t_max){ StatusFlag = 1; } 

      // criterio 2: ESCAPE. DSA ends bc the particle went too far from the shock, it will never reach it again.
      //if(IsDS &&Time<_t_max &&fabs(X)>2.*sqrt(fabs(_t_max-Time)*D0)){ StatusFlag = 2;} 
      if(IsDS && TimeLeft>0 && fabs(X)> 3.*sqrt(TimeLeft*D0) ){ StatusFlag = 2; } // ESCAPE

      // criterio 3: CPU. Maximum number of simulation steps reached. 
      if(ii>=NSTEPMAX-1){ StatusFlag = 3; } // Cpu. Max number of simulation steps reached 


      // --- end of dsa: record some quantity and break ---
      if(StatusFlag!=0){
	pD = D0*(year2sec)/(pc2cm*pc2cm); // cm2/s -> pc2/yr 
	pL = L0/pc2cm; // pc
	break;
      }

    
    }    
    // ----- end of step loop -----
    
    // ===============================

  
    // ========= Printout on screen ===========
    
    printf("\n");  
    cout<< "PART N.:  "<< nPart << endl;
    cout<< "TRACK N.:  "<< nTrack << endl;
    cout<< "Epoch N." << nT <<endl;
    cout<< "ENERGY:  "<< Ekn << endl;
    //cout<< "CROSSES: "<< nCross << endl;
    cout<<"Epoch T:" <<pT1;
    cout<< "T1/TMAX:  "<< 100.*Time1/_t_max <<" perc" << endl;
    cout<< "T2/TMAX:  "<< 100.*Time2/_t_max <<" perc" << endl;
    cout<< "nSteps:  "<< ((double)nSteps/(double)NSTEPMAX)*100. <<" perc" << endl;
    cout<< "Is DS?: "<< IsDS << endl;
    cout<< "STATUS:  "<< StatusFlag << endl;  
    cout<< endl;


    
    // variabili per tree NT2025
    pX  = X/pc2cm;  // posizione in parsec
    pE  = Ekn;  // kinetic energy per nucleon GeV/n
    pT  = Tdsa/year2sec;   // time spent in accelerator years
    pT1 = Time1/year2sec;
    pT2 = Time2/year2sec;

    //nCross = nCross; // number of shock crosses / n of gains
    iStat = StatusFlag; // status flag: 0,1,2,3
    
  // giÃ  registrati
  // iSide = -UDS; // side wrt shock at the end. US/DS +1/-1
  // nStep; // number of steps for that particle
  // nPart; // sequential number of particles
    
    
    PartTree->Fill();
    
    }
    // end of track/particle loop


  }
  // end of time/epoch loop
  
  PartFile->Write();
  PartFile->Close();
  

}



  // ------- DRAW ------------
  /*
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TCanvas* ccShockFrameTrajectory= new TCanvas("ccShockFrameTrajectory","ccShockFrameTrajectory",500,650);
  ccShockFrameTrajectory->Divide(1,3,0.01,0.01);

  double TMIN = _t_min/year2sec/1.3;
  double TMAX = _t_max/3./year2sec;

  // frame X 
  TH2F* hXFrame= new TH2F("hXFrame","hXFrame",200,TMIN,TMAX, 200,-0.025,0.055);
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

  // frame Ekn
  TH2F* hEknFrame= new TH2F("hEknFrame","hEknFrame",200,TMIN,TMAX, 200,6.0,1000.0);
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
  TH2F* hCrossFrame= new TH2F("hCrossFrame","hCrossFrame",200,TMIN,TMAX, 200,50.0,3000.0);
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

  TraceTree->Draw("X/(3.08568025e+18):Time/31556926.","nTrack==1 && (abs(X)>=0.0)","l same");

  TLine* line= new TLine(TMIN,0.0,TMAX,0.0);
  // TLine* line= new TLine(7.4,0.0,9.6,0.0);
  line->SetLineColor(kBlue+1);
  line->SetLineStyle(9);
  line->Draw();

  TraceTree->SetLineColor(kYellow);
  TraceTree->Draw("UDS:Time/31556926.","nTrack==1","l same");
  // gPad->SetGridy();
  TraceTree->SetLineColor(kBlack);
  TraceTree->Draw("X/(3.08568025e+18):Time/31556926.","nTrack==1 && (abs(X)>=0.0)","l same");



  ccShockFrameTrajectory->cd(2);
  TraceTree->SetLineColor(kBlack);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);

  hEknFrame->Draw();
  gPad->SetLogx();
  gPad->SetLogy();

  TraceTree->Draw("Ekn:Time/31556926.","nTrack==1","l same");
  TraceTree->SetLineColor(kYellow);
  TraceTree->Draw("1000*UDS:Time/31556926.","nTrack==1","l same");
  //gPad->SetGridy();
  TraceTree->SetLineColor(kBlack);
  TraceTree->Draw("Ekn:Time/31556926.","nTrack==1","l same"); 
  TraceTree->Draw("1.01*Ekn:Time/31556926.","nTrack==1","l same");


  ccShockFrameTrajectory->cd(3);
  TraceTree->SetLineColor(kBlack);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);

  hCrossFrame->Draw();
  gPad->SetLogx();
  gPad->SetLogy();

  TraceTree->Draw("nCross:Time/31556926.","nTrack==1","l same");
  TraceTree->Draw("1.01*nCross:Time/31556926.","nTrack==1","l same");
  //gPad->SetGridy();
  */

  // we obtain that <X> prop to sqrt(N)
  // <X> = K*sqrt(N) with K=sqrt(2/pi)!!
  
    

/*
  
// conversion rig->ekn
double GetEKNvsRIG(double Rig, double Z, double A){
  double Mp= 0.93827229;
  double ekn= sqrt( (Z/A)*(Z/A)*(Rig*Rig) + Mp*Mp ) - Mp;
  //printf("---> EKN: %3.3f \n",ekn);
  return ekn;
}

// conversion ekn->rig
double GetRIGvsEKN(double Ekn, double Z, double A){
  double Mp= 0.93827229;
  double rig= (A/Z)*sqrt( Ekn*Ekn + 2.*Ekn*Mp );
  //printf("---> RIG: %4.4f \n",rig);
  return rig;

}

// jacobian jrig = dekn/drig vs rig
double GetJvsRIG(double Rig, double Z, double A){ // bug fixed
  double Mp= 0.93827229;
  double jrig = (Z/A)*(Z/A)*Rig/sqrt( (Z/A)*(Z/A)*(Rig*Rig) + Mp*Mp );
  //printf("---> JRIG: %4.4f \n",jrig);
  return  jrig;
}

// jacobian jekn = drig/dekn vs ekn
double GetJvsEKN(double Ekn, double Z, double A){ // sostituito 2.->1. NT: PERCHE'??
  double Mp= 0.93827229;
  double jekn = (A/Z)*(1.*Mp + 1.*Ekn)/sqrt(Ekn*Ekn + 2.*Mp*Ekn);
  //printf("---> JEKN: %4.4f \n",jekn);
  return jekn;
}
*/
// main06.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It studies event properties of LEP1 events.

#include "Pythia8/Pythia.h"
#include "TH1D.h"// T for root
#include "TProfile.h"
#include "TLorentzVector.h" //four vectors
#include "TFile.h"
#include "TROOT.h"

using namespace Pythia8; // don't need to write pythia 8 for pythia packages

bool shouldkeep(Pythia & pythia,Particle & par); // def a bool fn shouldkeep

int main(int argc, char** argv) {
  Pythia pythia; // pythia class object
  const double pi = 3.14159265358979323;

  int nloop = 10000;
  TString base="ILLZZqq";
  if (argc >= 2 ) nloop=atoi(argv[1]); // seen these in c, CML arguments
  if (argc >= 3 ) base=argv[2];
  bool doW = false;
  if (base == "w" or base == "W	") doW = true;
  std::cout<<"main06ext starting for "<<nloop<<" events, with doW="<<doW<<std::endl;

  // Allow no substructure in e+- beams: normal for corrected LEP data.
  pythia.readString("PDF:lepton = off");// parton density fn off
  // Process selection.
  if (doW) {
    pythia.readString("WeakSingleBoson:ffbar2W = on");

    // Switch off all W decays and then switch back on those to quarks.
    pythia.readString("24:onMode = off");
    pythia.readString("24:onIfAny = 3 4");

    // Switch of decays for pi0 and J/psi
    pythia.readString("111:mayDecay = off");
    pythia.readString("443:mayDecay = off");

    // e-nu initialization at W mass.
    pythia.readString("Beams:idA =  12");
    pythia.readString("Beams:idB = -11");
    double mW = pythia.particleData.m0(24);
    pythia.settings.parm("Beams:eCM", mW);
  } else { //Z
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");// Weak bosons -ff ->Z
    // Switch off all Z0 decays and then switch back on those to quarks.
    pythia.readString("23:onMode = off");
    pythia.readString("23:onIfAny = 1 2 3 4 5"); //codes related to quarks pdg codes(23=z)

    // Switch of decays for pi0 and J/psi
    //pythia.readString("111:mayDecay = off"); // turns them off
    //pythia.readString("443:mayDecay = off");

    // LEP1 initialization at Z0 mass.
    pythia.readString("Beams:idA =  11"); //electron and positron
    pythia.readString("Beams:idB = -11");
    double mZ = pythia.particleData.m0(23);//get mass from pythia
    pythia.settings.parm("Beams:eCM", mZ); //energy of Z boson exact collison we want
  }

  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0"); // turns pseudo random off which will repeat themselves

  pythia.init();

  // Histograms
  TH1D * h_npart=new TH1D("npart","Pythia number of final-state particles",98,1.5,99.5);
  TH1D * h_nCharge = new TH1D("nCharge","Pythia number of final-state, charged particles", 30,1.,61.);
  TH1D * h_sumE = new TH1D("sumE ","Pythia sum E "+base+" ; Energy, GeV",120,0.,120.);
  TH1D * h_thrust = new TH1D("thrust","Pythia thrust"+base+"; 1-T",100,0.,0.5);
  //defined 2 new Histograms maybe they need to be 2d??


  TH1D * h_nDurham = new TH1D("Durham jet multiplicity","Durham jet multiplicity", 40, -0.5, 39.5);
  TH1D * h_eDifDurham = new TH1D("Durham e_{i} - e_{i+1}","Durham e_{i} - e_{i+1}", 100, -5.,45.);
  // create a new histogram for the masses on this line
  TH1D * h_mass = new TH1D("mass wide range","mass wide range", 95, 0.,95.); // 100 bins 5-45? // should be mz
  //TH1D * h_massmz = new TH1D("mass ","mass", 250, 89,92.); // 100 bins 5-45? // should be mz
  TH1D *hist_mass[3]; // Array for each of the different mass options
	TH1D *hist_fmass;
hist_fmass = new TH1D("mass peak","mass", 10, 85,95.);
//x axis number of kaons and y is the mean mass 0-100
  TProfile * h_gamfrac = new TProfile("gamfrac","gamfrac",98,1.5,99.5,0.,1.);
  TProfile * h_jpsifrac = new TProfile("jpsifrac","jpsifrac",98,1.5,99.5,0.,1.);
  TProfile * h_kfrac = new TProfile("Kfrac","Kfrac",98,1.5,99.5,0.,1.);
  TProfile * h_pi0frac = new TProfile("pi0frac","pi0frac",98,1.5,99.5,0.,1.);
  TProfile * h_pipmfrac = new TProfile("pipmfrac","pipmfrac",98,1.5,99.5,0.,1.);
  TProfile * h_protonfrac = new TProfile("protonfrac","protonfrac",98,1.5,99.5,0.,1.);
  TProfile * h_neutronfrac = new TProfile("neutronfrac","neutronfrac",98,1.5,99.5,0.,1.);
  TProfile * h_lepfrac = new TProfile("lepfrac","lepfrac",98,1.5,99.5,0.,1.);

  Thrust thr;//event quantity - TYPE OF JETS
  //array with each jet finder
  //ClusterJet durham2("Durham",2, 2, false, false); // Another type of jet
  //ClusterJet durham1("Durham",2, 1, false, false);
  //ClusterJet durham0("Durham",2, 0, false, false);
  //ClusterJet finder[3] = {ClusterJet("Durham",2,2),ClusterJet("Durham",2,2)}]
  //craetiing each class of jet
ClusterJet* finder[3] ;
for (int j = 0; j < 3; j++) {
  TString name = "mass";
  name += j;
  finder[j] = new ClusterJet("Durham",2,j);
  hist_mass[j] = new TH1D(name,name, 45, 50.,95.);
}
  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nloop; ++iEvent) {
    if (!pythia.next()) continue;

    int nparticles = 0;
    int ncharged = 0;
    double sumE = 0;

    int ngam = 0;
    int njpsi = 0;
    int nK = 0;
    int npi0 = 0;
    int npipm = 0;
    int nproton = 0;
    int nneutron = 0;
    int nlep = 0;
    int nneutrino = 0;
    int nother = 0;

    std::vector<TLorentzVector> partons;
    std::vector<TLorentzVector> particles;

    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle & par = pythia.event[i];
      //      std::cout<<"Particle = "<<par.id()<<"; Status = "<<par.status()<<"; isFinal = "<<par.isFinal()<<";pT = "<<par.pT()<<"; m = "<<par.m()<<std::endl;ize and masses custom- maximise number of 4 jet events

      if (par.status() == -23) {
	// Initial partons.
	//	std::cout<<"Status is -23; particle is: "<<par.id()<<std::endl;

	TLorentzVector thisParton(par.px(),par.py(),par.pz(),par.e());
	partons.push_back(thisParton);
      }

      bool keep = shouldkeep(pythia,par);
      if (keep) {
	sumE += par.e();
	nparticles++;
	switch (abs(par.id())) {
	case 11:
	case 13:
	case 15:
	  nlep++; break;
	case 12:
	case 14:
	case 16:
	  nneutrino++; break;
	case 22:
	  ngam++; break;
	case 443:
	  njpsi++; break;
	case 111:
	  npi0++; break;
	case 211:
	  npipm++; break;
	case 130:
	case 310:
	case 311:
	case 321:
	case 313:
	case 323:
	  nK++; break;
	case 2112:
	  nneutron++; break;
	case 2212:
	  nproton++; break;
	default:
	  nother++;
	}
	TLorentzVector thisParticle(par.px(),par.py(),par.pz(),par.e());
	particles.push_back(thisParticle);
      }
      if (par.isFinal() && par.isCharged()) ++ncharged;
    }

    if (partons.size() == 2) { // Pythia was asked to generate W/Z->qq. did it?
	h_npart->Fill( nparticles );
	h_nCharge->Fill( ncharged);
	h_sumE->Fill(sumE);

	h_pi0frac->Fill(nparticles,double(npi0)/double(nparticles));
	h_pipmfrac->Fill(nparticles,double(npipm)/double(nparticles));
	h_kfrac->Fill(nparticles,double(nK)/double(nparticles));
	h_jpsifrac->Fill(nparticles,double(njpsi)/double(nparticles));
	h_gamfrac->Fill(nparticles,double(ngam)/double(nparticles));
	h_protonfrac->Fill(nparticles,double(nproton)/double(nparticles));
	h_neutronfrac->Fill(nparticles,double(nneutron)/double(nparticles));
	h_lepfrac->Fill(nparticles,double(nlep+nneutrino)/double(nparticles));

	// Find and histogram thrust (an  property descrinig how collimated this event was)
	if (thr.analyze( pythia.event )) {
	  h_thrust->Fill( 1 - thr.thrust() );
	}

//size and masses custom- maximise number of 4 jet events

  if (finder[2]->analyze( pythia.event, 0.01, 0.)) { //0.01 width 2nd one is mass check manual
	finder[0]->analyze( pythia.event, 0.01, 0.);
	finder[1]->analyze( pythia.event, 0.01, 0.);
    if (iEvent < 3) finder[2]->list();
    h_nDurham->Fill( finder[2]->size() );
    //for (int j = 0; j < finder[2]->size() - 1; ++j)
      //h_eDifDurham-> Fill( finder[]->p(j).e() - finder[2]->p(j+1).e() ); //four vector-check manual
      // first check finder[]->size == 2
      if (finder[2]->size() == 2 && finder[1]->size() == 2) {
    // add both 4 vectors
    for (int j = 0; j < 3; j++) {
      Vec4 total_4vector = finder[j]->p(0) + finder[j]->p(1);
  // Get the mass from the total 4_vector
  double mass_jets = total_4vector.mCalc();
  // Fill it into a histogram
  hist_mass[j]-> Fill(mass_jets);

	if (j==2)
{
hist_fmass-> Fill(mass_jets);
}

    }





      }
    }




  }
// End of event loop. Statistics. Output histograms.
}

pythia.stat();
 // End of event loop. Print statistics. Output the histograms to a root file

	// Add date and time on a new line
  TString opname="pythia_"+base+".root";
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(opname);
  if (!f) {
    f = new TFile(opname,"RECREATE");
  }
  hist_mass[0]->Write();
  hist_mass[1]->Write();
  hist_mass[2]->Write();
  hist_fmass->Write();
  h_npart->Write();
  h_nCharge->Write();
  h_sumE->Write();
  h_thrust->Write();


  h_gamfrac->Write();
  h_jpsifrac->Write();
  h_kfrac->Write();
  h_pi0frac->Write();
  h_pipmfrac->Write();
  h_protonfrac->Write();
  h_neutronfrac->Write();
  h_lepfrac->Write();
  h_nDurham->Write();
  h_eDifDurham->Write();

  h_mass->Write();



  f->Write();
  f->Close();

  return 0;
}

bool shouldkeep(Pythia & pythia,Particle & par) {

  bool keep = false;

  if (par.isFinal()) {
    keep = true;
  }
  return keep;
}

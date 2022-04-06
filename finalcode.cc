
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
  TString base="Zqq";
  if (argc >= 2 ) nloop=atoi(argv[1]); // seen these in c, CML arguments
  if (argc >= 3 ) base=argv[2];
  bool doW = false;
  if (base == "w" or base == "W	") doW = true;
  std::cout<<"main06ext starting for "<<nloop<<" events, with doW="<<doW<<std::endl;

  // Allow no substructure in e+- beams: normal for corrected LEP data.
  pythia.readString("PDF:lepton = off");// parton density fn off
  // Process selection.


  //z decay
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");// Weak bosons -ff ->Z
    // Switch off all Z0 decays and then switch back on those to quarks.
    pythia.readString("23:onMode = off");

    pythia.readString("23:onIfAny = 4"); //codes related to quarks pdg codes(23=z)

    // Switch of decays for pi0 and J/psi
    //pythia.readString("310:mayDecay = off"); // less kaons on
    //pythia.readString("443:mayDecay = off");
    pythia.readString("111:mayDecay = off");

    // LEP1 initialization at Z0 mass.
    pythia.readString("Beams:idA =  11"); //electron and positron
    pythia.readString("Beams:idB = -11");
    double mZ = pythia.particleData.m0(23);//get mass from pythia
    pythia.settings.parm("Beams:eCM", mZ); //energy of Z boson exact collison we want


  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0"); // turns pseudo random off which will repeat themselves

  pythia.init();

  // Histograms



  TH1D * h_nDurham = new TH1D("Durham jet multiplicity","Durham jet multiplicity",  95, 0.,95.);

  TH1D *hist_mass[3];


  TProfile * protonnum  = new TProfile("proton effect","Number of protons and effect on avg mass",8,-0.5,7.5,82,96 );
  TProfile * pionum  = new TProfile("charged pion effect","Number of charged pions and effect on avg mass",36,0. ,36.,75,100 );
  TProfile * kaon0  = new TProfile("kaon effect","Number of neutral kaon and effect on avg mass",8,-0.5,7.5,82,96 );
  TProfile * neutronnum  = new TProfile("neutron effect","Number of neutron and effect on avg mass",7,0.,7.,82,96 );
  double bins[10] = {0,0.5,1,1.5,2,4,7,10,20,40};
  TProfile * kaonpm  = new TProfile("charged kaon effect","Number of charged kaons and effect on avg mass",11,-0.5,10.5,82,96 );
  TProfile * kaonpmMom = new TProfile("Charged Kaon momentum","Momentum of single charged kaon events and avg mass",9,bins,82,96);


  TH1D *h_neutrinonum = new TH1D("Number of neutrinos ","Number of neutrinos ",20,0,20);
  TH1D *h_elecnum = new TH1D("Number of Leptons ","Number of Leptons ",20,0,20);
  TH1D *protons= new TH1D("Number of protons ","Number of protons ",7,-.5,6.5);
  TH1D *neutrons= new TH1D("Number of neutrons ","Number of neutrons ",7,-.5,6.5);
  TH1D *kaons0= new TH1D("Number of neutral kaons ","Number of neutral kaons;x;y",16,-.5,15.5);
  TH1D *kaonspm= new TH1D("Number of charged kaons ","Number of charged kaons ",16,-.5,15.5);
  TH1D *pionspm= new TH1D("Number of charged pions ","Number of charged pions ",51,-.5,50.5);
  TH1D *pions0= new TH1D("Number of neutral pions ","Number of neutral pions ",51,-.5,50.5);
  TH1D *kaonspmmmntm= new TH1D("Momentums of single kaons ","momentum of single kaons ",40,.0,40.0);
  TH1D *totalmom= new TH1D("Momentums of all kaons ","momentum of all kaons ",40,.0,40.0);






  Thrust thr;//event quantity - TYPE OF JETS

ClusterJet* finder[3] ;
for (int j = 0; j < 3; j++) {
  TString name = "mass";
  name += j;
  finder[j] = new ClusterJet("Durham",2,j);
  hist_mass[j] = new TH1D(name,name, 95, 0.,95.);
  hist_mass[j]->GetXaxis()->SetTitle("Mass(GeV/c^{2})");
  hist_mass[j]->GetYaxis()->SetTitle("Number of events");
  hist_mass[j]->GetXaxis()->SetLabelSize(.05);
  hist_mass[j]->GetXaxis()->SetTitleSize(.045);
  hist_mass[j]->GetYaxis()->SetTitleOffset(1.05);
  hist_mass[j]->GetYaxis()->SetTitleSize(.05);
  hist_mass[j]->GetYaxis()->SetLabelSize(.045);

}

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nloop; ++iEvent) {
    if (!pythia.next()) continue;

    int nparticles = 0;
    int ncharged = 0;
    double sumE = 0;
    double kaonpmMomentum=0;
    int ngam = 0;
    int njpsi = 0;
    int nK0 = 0;
    int nKpm = 0;
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
    //std::cout<<"Particle = "<<par.id()<<"; Status = "<<par.status()<<"; isFinal = "<<par.isFinal()<<";pT = "<<par.pT()<<"; m = "<<par.m()<<std::endl;
//size and masses custom- maximise number of 4 jet events

      if (par.status() == -23) {
	// Initial partons.
	//	std::cout<<"Status is -23; particle is: "<<par.id()<<std::endl;

	TLorentzVector thisParton(par.px(),par.py(),par.pz(),par.e());
	partons.push_back(thisParton);
      }
// pdg codes look up
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
	case 130: //ko
	case 310://ko
	//case 311: //K0
	//case 313: //K*
	  nK0++; break;
  case 321:// kpm
  case 323: //k*+
  kaonpmMomentum = par.pAbs();
  totalmom ->Fill(kaonpmMomentum);
  nKpm++; break;
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


    if (nneutrino > 0 ) continue;



    protons ->Fill(nproton);
    pionspm ->Fill(npipm);
    pions0 ->Fill(npi0);
    neutrons ->Fill(nneutron);
    h_neutrinonum->Fill(nneutrino);
    h_elecnum->Fill(nlep);

    kaons0 ->Fill(nK0);
    kaonspm ->Fill(nKpm);


  if (finder[2]->analyze( pythia.event, 0.01, 0.)) //0.01 width 2nd one is mass check manual
  {
	   finder[0]->analyze( pythia.event, 0.01, 0.);
	   finder[1]->analyze( pythia.event, 0.01, 0.);
    if (iEvent < 3) finder[2]->list();
      h_nDurham->Fill( finder[2]->size() );

      if (finder[2]->size() == 2 && finder[1]->size() == 2) {
    // add both 4 vectors
    for (int j = 0; j < 3; j++)
    {
      Vec4 total_4vector = finder[j]->p(0) + finder[j]->p(1);
      // Get the mass from the total 4_vector
      double mass_jets = total_4vector.mCalc();
      // Fill it into a histogram
      hist_mass[j]-> Fill(mass_jets);
      if (j==1)
      {

        protonnum ->Fill(nproton,mass_jets);

        neutronnum ->Fill(nneutron,mass_jets);
        pionum ->Fill(npipm,mass_jets);
        kaon0 ->Fill(nK0,mass_jets);
        kaonpm ->Fill(nKpm,mass_jets);
        if (nKpm == 1)
        {
        kaonpmMom->Fill(kaonpmMomentum,mass_jets);//profile plot
        kaonspmmmntm ->Fill(kaonpmMomentum);//hist
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


  protonnum ->Write();
  kaon0 -> Write();
  kaonpm ->Write();
  neutronnum -> Write();
  pionum ->Write();
  kaonpmMom->Write();
  h_neutrinonum ->Write();
  h_elecnum ->Write();
  protons ->Write();
  kaons0 ->Write();
  kaonspm ->Write();
  pionspm ->Write();
  pions0 ->Write();
  neutrons ->Write();
  kaonspmmmntm ->Write();
  totalmom ->Write();


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

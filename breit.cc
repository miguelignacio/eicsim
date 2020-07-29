// main92.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// Modified by Rene Brun and Axel Naumann to put the Pythia::event
// into a TTree.


// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for saving Pythia events as trees in a file.
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"

#include "Pythia8Plugins/FastJet3.h"
#include "fastjetbase.h"
using namespace Pythia8;
#include <cmath>
#include <vector>
#include <math.h>

#include <TROOT.h>
#include "TVector2.h"
#include "TVector3.h"
#include <TLorentzVector.h>

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

#include "fastjet/contrib/CentauroPlugin.hh"
//#include "CentauroPlugin.hh"
//#include "fastjet/EECambridgePlugin.hh"

//#include <iostream>
//#include <fstream>
//using namespace std;

fastjet::PseudoJet operator-(const fastjet::PseudoJet &p)
{
  return fastjet::PseudoJet(-p.px(), -p.py(), -p.pz(), p.E());
}

//-------------------------------------------------------------------
// Lorentz transformations (Use For Breit Frame)
//-------------------------------------------------------------------
// boost the vector p with the boost given by (bx,by,bz)
double Vec4Dot(Vec4 a , Vec4 b){                                                                                                                                                                                     double s = a[0] * b[0]  - ( a[1] * b[1] + a[2] * b[2] +  a[3] * b[3] );                                                                                                                                            return s;
}

fastjet::PseudoJet boost(const fastjet::PseudoJet p,
			 const double &bx, const double &by, const double &bz){
  double b2 = bx*bx + by*by + bz*bz;
  //assert(b2 < 1.0);
  if(b2 >= 1.0)
    {
      b2 = 0.999;
    }

  double gamma = 1.0/sqrt(1.0 - b2);
  double bp = bx*p.px() + by*p.py() + bz*p.pz();
  double gamma2 = (b2 > 0.0 ? (gamma - 1.0)/b2 : 0.0);

  return fastjet::PseudoJet(p.px() + gamma2*bp*bx + gamma*bx*p.E(),
			    p.py() + gamma2*bp*by + gamma*by*p.E(),
			    p.pz() + gamma2*bp*bz + gamma*bz*p.E(),
			    gamma*(p.E() + bp));
}

// boost the vector p with the boost given by b (i.e. (bx/bE,by/bE,bz/bE))
fastjet::PseudoJet boost(const fastjet::PseudoJet p, const fastjet::PseudoJet b){
  return boost(p, b.px()/b.E(), b.py()/b.E(), b.pz()/b.E());
}

// rotation around the x axis
fastjet::PseudoJet rotateX(const fastjet::PseudoJet p, const double& psi){
  double sp = sin(psi);
  double cp = cos(psi);
  return fastjet::PseudoJet(p.px(),
			    cp*p.py()-sp*p.pz(),
			    sp*p.py()+cp*p.pz(),
			    p.E());
}

// rotation around the y axis
fastjet::PseudoJet rotateY(const fastjet::PseudoJet p, const double& psi){
  double sp = sin(psi);
  double cp = cos(psi);
  return fastjet::PseudoJet(cp*p.px()+sp*p.pz(),
			    p.py(),
			    cp*p.pz()-sp*p.px(),
			    p.E());
}

// rotation around the z axis
fastjet::PseudoJet rotateZ(const fastjet::PseudoJet p, const double& psi){
  double sp = sin(psi);
  double cp = cos(psi);
  return fastjet::PseudoJet(cp*p.px()-sp*p.py(),
			    sp*p.px()+cp*p.py(),
			    p.pz(),
			    p.E());
}

int main() {

  int nEvent    = 5e3;
  int    power   = -1;     // -1 = anti-kT; 0 = C/A; 1 = kT.
  double R       =  1.0;  // Jet size.
  double pTMin   = 0.0;
  double etaMax  = 4.0;    // Pseudorapidity range of detector.

  double eProton   = 100.0;
  double eElectron = 10;
  double Q2min     = 100.0;

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;
  
  pythia.readString("Random:setSeed=on");
  pythia.readString("Random:seed=0");

  pythia.readString("Beams:idB=11");


  pythia.readString("Beams:idA=2212");

  pythia.settings.parm("Beams:eB", eElectron);
  pythia.settings.parm("Beams:eA", eProton);

  pythia.readString("Beams:frameType=2");
  pythia.readString("Init:showChangedSettings=on");
  pythia.readString("Main:timesAllowErrors=10000");

  //neutral current
   pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  //charged current
  // pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");

  pythia.settings.parm("PhaseSpace:Q2Min",Q2min);
  pythia.readString("SpaceShower:pTmaxMatch=2");
  pythia.readString("SpaceShower:dipoleRecoil=on");
  pythia.init();

  
  fastjet::contrib::CentauroPlugin * centauro_plugin = new fastjet::contrib::CentauroPlugin(1.0);
  fastjet::JetDefinition jetDefCentauro(centauro_plugin);

  //jetDefCentauro.set_recombination_scheme(fastjet::E_scheme);
  jetDefCentauro.set_recombination_scheme(fastjet::WTA_modp_scheme);
  
  std::vector <fastjet::PseudoJet> fjInputs;
  std::vector <fastjet::PseudoJet> fjInputs_breit;

  // Set up the ROOT TFile and TTree.
  TFile *file = TFile::Open("breit.root","recreate");
  TTree *T = new TTree("T","jet Tree");


  //Tree variables:
  UInt_t ntrials, evid, quark_id;
  Float_t xsec, x, y, Q2, W2;
  Float_t quark_e, quark_eta, quark_pt, quark_p, quark_phi, quark_theta;


  UInt_t nconstituents;
  UInt_t n_charged;

  float jet_pt;
  float jet_qt;
  float jet_qt_up;
  float jet_qt_down;
  
  float jet_phi;
  float jet_rap;
  float jet_eta;
  float jet_theta;
  float jet_p;
  float jet_e;
  float jet_z;
  float jet_z_up;
  float jet_z_down;        

  float jet_lab_eta;
  float jet_lab_e;
  float jet_lab_pt;

  std::vector<float> h_z;
  std::vector<float> h_j;
  std::vector<float> h_pid;
  std::vector<float> h_eta;
  std::vector<float> h_rap;
  std::vector<float> h_pt;
  std::vector<double> h_charge;
  
  T->Branch("ntrials", &ntrials, "ntrials/I");
  T->Branch("evid", &evid, "evid/I");
  T->Branch("xsec", &xsec, "xsec/F");
  T->Branch("x", &x, "x/F");
  T->Branch("y", &y, "y/F");
  T->Branch("Q2", &Q2, "Q2/F");
  T->Branch("W2", &W2, "W2/F");

  T->Branch("quark_id", &quark_id, "quark_id/I");
  T->Branch("quark_e", &quark_e, "quark_e/F");
  T->Branch("quark_eta", &quark_eta, "quark_eta/F");
  T->Branch("quark_pt", &quark_pt,"quark_pt/F");
  T->Branch("quark_p",&quark_p,"quark_p/F");
  T->Branch("quark_phi", &quark_phi,"quark_phi/F");
  T->Branch("quark_theta", &quark_theta,"quark_theta/F");

  T->Branch("n_total",&nconstituents, "nconstituents/I");
  T->Branch("n_charged", &n_charged, "n_charged/I");
  T->Branch("jet_pt", &jet_pt, "jet_pt/F");
  T->Branch("jet_qt", &jet_qt, "jet_qt/F");
  T->Branch("jet_qt_down", &jet_qt_down, "jet_qt_down/F");
  T->Branch("jet_qt_up", &jet_qt_up, "jet_qt_up/F");
  
  T->Branch("jet_phi", &jet_phi, "jet_phi/F");
  T->Branch("jet_rap",&jet_rap, "jet_rap/F");
  T->Branch("jet_eta", &jet_eta, "jet_eta/F");
  T->Branch("jet_theta", &jet_theta, "jet_theta/F");
  T->Branch("jet_p", &jet_p, "jet_p/F");
  T->Branch("jet_e", &jet_e, "jet_e/F");
  T->Branch("jet_z", &jet_z, "jet_z/F");
  T->Branch("jet_z_up", &jet_z_up, "jet_z_up/F");
  T->Branch("jet_z_down", &jet_z_down, "jet_z_down/F");         
  
  T->Branch("jet_lab_pt", &jet_lab_pt, "jet_lab_pt/F");
  T->Branch("jet_lab_eta", &jet_lab_eta, "jet_lab_eta/F");
  T->Branch("jet_lab_e" ,  &jet_lab_e, "jet_lab_e/F");
  //hadron variables
  T->Branch("z", &h_z);
  T->Branch("jt", &h_j);
  T->Branch("pid", &h_pid);
  T->Branch("eta", &h_eta);
  T->Branch("rap", &h_rap);
  T->Branch("pt", &h_pt);
  T->Branch("charge",&h_charge);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    fjInputs.clear();
    fjInputs_breit.clear();
    //empty vectors
    quark_id = 0;
    quark_e = 0;
    quark_pt = 0;
    quark_eta = 0;
    quark_phi = 0;
    quark_p = 0;
    quark_theta = 0;


    //general event info
    evid = iEvent;
    xsec = pythia.info.sigmaGen();
    ntrials = pythia.info.nTried();

    // four-momenta of proton, electron, virtual photon/Z^0/W^+-.

    // get electron index:
    int e_index =6 ;
    float maxpt = 0;
    for (int i = 0; i < event.size(); i++)
      {
	if (event[i].isFinal() && event[i].id()==11)
	  {
	    double pt = TMath::Sqrt(event[i].px()*event[i].px() + event[i].py()*event[i].py()) ;
	    if(pt>maxpt){
	        e_index = i;
		maxpt = pt;
	    }
        }
      }
    //std::cout << " e index " << e_index << std::endl;
    //std::cout <<  event[e_index].px() << " " << event[e_index].py() << " " << event[e_index].pz() << std::endl;
    //std::cout <<  event[6].px() << " " << event[6].py() << " " << event[6].pz() << std::endl;
    //e_index = 6;
    
    Vec4 pProton = event[1].p();
    Vec4 peIn = event[4].p();
    Vec4 peOut = event[e_index].p();
    Vec4 pPhoton = peIn - peOut;
    Vec4 ptot_had = pPhoton + pProton;

    // Q2, W2, Bjorken x, y, nu.
    Q2 = -pPhoton.m2Calc();
    W2 = (pProton + pPhoton).m2Calc();
    x = Q2 / (2. * pProton * pPhoton);
    y = (pProton * pPhoton) / (pProton * peIn);

    //define boost and rotations to get to the Breit-Frame
    fastjet::PseudoJet proton(pProton.px(), pProton.py(), pProton.pz(), pProton.e());
    fastjet::PseudoJet gamma( pPhoton.px(), pPhoton.py(), pPhoton.pz(), pPhoton.e());
    fastjet::PseudoJet boost_vector = -(gamma + 2.0*x*proton);
    
    fastjet::PseudoJet temp = boost(proton,boost_vector);

    Double_t phi_p = temp.phi();
    Double_t theta_p = TMath::ATan2(temp.perp(),temp.pz());

    fastjet::PseudoJet boosted_gamma = boost(gamma, boost_vector);
    boosted_gamma = rotateZ(boosted_gamma, -phi_p);
    boosted_gamma = rotateY(boosted_gamma, -theta_p);

    fastjet::PseudoJet boosted_proton = boost(proton,boost_vector);
    boosted_proton = rotateZ(boosted_proton,-phi_p);
    boosted_proton = rotateY(boosted_proton, -theta_p);
    //checkout gamma after boosting and rotations:
    // std::cout << " boosted_gamma" << boosted_gamma.px() << " " << boosted_gamma.py() << " " << boosted_gamma.pz() << " " << boosted_gamma.e() << std::endl;
    //std::cout << " boosted_proton " << boosted_proton.px() << " " << boosted_proton.py() << " " << boosted_proton.pz() << " " << boosted_proton.e() << std::endl;  
    //std::cout << "Q2 " << Q2  << " Q " << TMath::Sqrt(Q2) << std::endl;

    //fastjet::PseudoJet boostingback = rotateY(boosted_proton,+theta_p);
    //boostingback  = rotateZ(boostingback,+phi_p);
    //boostingback = boost(boostingback, -boost_vector);
    //std::cout << " boostedback " << boostingback.px() << " " << boostingback.py() << " " << boostingback.pz() << " " << boostingback.e() << std::endl;
    //std::cout << " proton " << proton.px() << " " << proton.py() << " " << proton.pz() << " " << proton.e() << std::endl;    
    //Q is the pz of the photon in the breit frame: check.
    //std::cout <<  " ////////" << std::endl;
    // get struck quark index:
    int q;
    for (int i = 0; i < event.size(); i++)
      {
	if (event[i].status() == -23 && event[i].id() != 11)
	  {
	    q = i;
	    break;
	  }
      }

    fastjet::PseudoJet quark(event[q].px(), event[q].py(), event[q].pz(),event[q].e());
    quark = boost(quark, boost_vector);
    quark = rotateZ(quark, -phi_p);
    quark = rotateY(quark, -theta_p);

    fastjet::PseudoJet initial_quark(event[3].px(), event[3].py(), event[3].pz(),event[3].e());
    initial_quark = boost(initial_quark, boost_vector);
    initial_quark = rotateZ(initial_quark, -phi_p);
    initial_quark = rotateY(initial_quark, -theta_p);

    quark_id = event[q].id();
    quark_pt = quark.pt();
    quark_eta = quark.eta();
    quark_phi = quark.phi_std();
    quark_p = sqrt(quark.modp2());
    quark_theta = acos( event[q].pz() /sqrt(quark.modp2()));

    //std::cout << " initial_quark_pt " << initial_quark.pt() << std::endl;
    //std::cout << initial_quark.px() << " " << initial_quark.py() << " " << initial_quark.pz() << " " << initial_quark.e() << std::endl;

    //std::cout << " struck quark_pt " << quark.pt() << std::endl;
    //std::cout << quark.px() << " " << quark.py() << " " << quark.pz() << " " << quark.e() << std::endl;

    Vec4   pTemp;
    double mTemp;
    int nAnalyze = 0;
    //loop over particles in the event and store them as input for FastJet
    //std::cout << " partiles in event " << event.size() << std::endl;
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {
        if(i ==e_index) continue; //remove the highest pT electron
	// Require visible/charged particles inside detector.
	if ( !event[i].isVisible()  ) continue;
	
	// Create a PseudoJet from the complete Pythia particle.
	fastjet::PseudoJet particle_lab(event[i].px(), event[i].py(), event[i].pz(),event[i].e());
	if (particle_lab.pt()<0.100) continue;
        if(abs(particle_lab.eta())>4.0)continue;
	// Define Breit Fram Kinematics
	//std::cout << event[i].id() << std::endl;
	//	std::cout << event[i].px() << "    " << event[i].py() << "    " << event[i].pz() << "    " << event[i].e() << std::endl;	
	
	fastjet::PseudoJet boosted_particle = boost(particle_lab,boost_vector);
	boosted_particle = rotateZ(boosted_particle, -phi_p);
	boosted_particle = rotateY(boosted_particle, -theta_p);

	//std::cout << boosted_particle.px() <<"    " << boosted_particle.py() << "    " << boosted_particle.pz() <<"     " << boosted_particle.e() << std::endl;
	particle_lab.set_user_info(new MyUserInfo(event[i].id(),i,event[i].charge()));
	boosted_particle.set_user_info(new MyUserInfo(event[i].id(),i,event[i].charge()));

	
	fjInputs.push_back(particle_lab);
	fjInputs_breit.push_back(boosted_particle);

	++nAnalyze;
	
	  
      } //end loop over particles
    //std::cout << "###########################################" << std::endl;

    // Run Fastjet algorithm and sort jets in pT order.
    vector <fastjet::PseudoJet> inclusiveJets, sortedJets;

    //do Breit-frame rotationally invariant
    //fastjet::ClusterSequence clustSeq(fjInputs_breit, jetDefBF);
    //inclusiveJets = clustSeq.inclusive_jets();
    //sortedJets    = sorted_by_E(inclusiveJets);

    //do Breit-frame Centauro
    fastjet::ClusterSequence clustSeq(fjInputs_breit, jetDefCentauro);
    inclusiveJets = clustSeq.inclusive_jets();
    sortedJets = sorted_by_E(inclusiveJets);
    
    // do lab frame:
    //fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
    //inclusiveJets  = clustSeq.inclusive_jets();
    //sortedJets     = sorted_by_pt(inclusiveJets);
    
    
    //loop over jets
    //std::cout << " //////////////////////////// " << std::endl;
    for (unsigned ijet= 0; ijet < sortedJets.size();ijet++) {
      fastjet::PseudoJet jet = sortedJets[ijet];

      //boost back to lab frame:
      fastjet::PseudoJet jet_lab = rotateY(jet,+theta_p);
      jet_lab  = rotateZ(jet_lab,+phi_p);
      jet_lab = boost(jet_lab, -boost_vector);
      jet_lab_e = jet_lab.e();
      jet_lab_eta = jet_lab.eta();
      jet_lab_pt = jet_lab.perp();
      vector<fastjet::PseudoJet> constituents = jet.constituents();
      double ncharged_constituents = fastjet::SelectorIsCharged().count(constituents);
      
      nconstituents = constituents.size();
      n_charged = ncharged_constituents;
      //std::cout << jet.perp() << std::endl;
      double jetpt = jet.perp();
      double jetrap = jet.rap();
      double jetphi = jet.phi_std();
      double theta = jet.pz()/sqrt(jet.modp2()); 
      jet_theta=acos(theta);
      
      jet_pt=jet.perp();
      jet_phi= jet.phi_std();
      
      double pxj, pyj, pzj, p_jet, e_jet;

      pxj = jet.px();
      pyj = jet.py();
      pzj = jet.pz();
      p_jet = sqrt(jet.modp2());
      e_jet  = jet.e();

      Vec4 n   (0, 0, +1, +1); //was -
      Vec4 nbar(0, 0, -1, +1); //was +  
      Vec4 pJ ( jet.px(), jet.py() ,jet.pz() , jet.e() ) ;
      Vec4 p (boosted_proton.px(), boosted_proton.py(), boosted_proton.pz(), boosted_proton.e());
      Vec4 q (boosted_gamma.px(), boosted_gamma.py(), boosted_gamma.pz(), boosted_gamma.e());
      double y_jet = - log( Vec4Dot(nbar, pJ) / Vec4Dot(n, pJ) ) / 2. ;
      double eJ =  jet.e() ;
      double z = Vec4Dot(n, pJ) / sqrt(Q2);
      double z_jet = Vec4Dot(p, pJ) / Vec4Dot(p, q);
      double z_e =  2 * eJ / sqrt(Q2);
      double qT = jet.pt() / z ;                       
      
          
      jet_qt=qT;

      jet_p=p_jet;
      jet_e=eJ;
      jet_z=z;
      
      
      jet_rap=y_jet;
      jet_eta= jet.eta();

            
      //Plot what would happen if you move up the JES scale or down.
      double JES = 1.02;
      fastjet::PseudoJet jet_up = jet_lab*JES;
      jet_up = boost(jet_up, boost_vector);
      jet_up = rotateZ(jet_up, -phi_p);
      jet_up = rotateY(jet_up, -theta_p);
      Vec4 pJ_up ( jet_up.px(), jet_up.py() ,jet_up.pz() , jet_up.e() );
      double z_up = Vec4Dot(n, pJ_up) / sqrt(Q2);
      double qT_up = jet_up.pt() / z_up ;
      
      fastjet::PseudoJet jet_down = jet_lab*(1.0/JES);
      jet_down = boost(jet_down, boost_vector);
      jet_down = rotateZ(jet_down, -phi_p);
      jet_down = rotateY(jet_down, -theta_p);
      Vec4 pJ_down ( jet_down.px(), jet_down.py() ,jet_down.pz() , jet_down.e() );
      double z_down = Vec4Dot(n, pJ_down) / sqrt(Q2);
      double qT_down = jet_down.pt() / z_down ;
      
      jet_z_up=z_up;
      jet_z_down= z_down;
      jet_qt_up = qT_up;
      jet_qt_down = qT_down;
      // std::cout << " zup " << z_up << " z " << z << " zdo " << z_down << std::endl;
      //std::cout << " qTup " << qT_up << " qt " << qT << " qTdo " << qT_down << std::endl;




      //loop over constituents
      for (unsigned n = 0; n < constituents.size(); n++){
	fastjet::PseudoJet _p = constituents[n]; //.user_info<FJUtils::PythiaUserInfo>().getParticle();
        //cout << "constituent size" <<constituents.size() << endl;
	//if(_p.user_info<MyUserInfo>().charge()!=0){

	  double pxh, pyh, pzh, cross;
	  pxh = _p.px();
	  pyh = _p.py();
	  pzh = _p.pz();
	  	  
	  cross = sqrt( pow((pyj*pzh-pyh*pzj),2.0) + pow((pxj*pzh-pzj*pxh),2.0) + pow((pxj*pyh-pyj*pxh),2.0) );
	  h_z.push_back((pxj*pxh + pyj*pyh + pzj*pzh) / (p_jet*p_jet));
	  h_j.push_back(cross / p_jet);
	  h_pid.push_back(_p.user_info<MyUserInfo>().pdg_id());
	  h_charge.push_back(_p.user_info<MyUserInfo>().charge());
	  h_rap.push_back(_p.rap());
	  h_eta.push_back(_p.eta());
	  h_pt.push_back(_p.pt());
	  
	  //}
      }//end loop over constituents

      T->Fill();
      //empty vectors
      jet_pt =0;
      jet_qt =0;
      jet_qt_up =0;
      jet_qt_down=0;
      jet_phi =0;
      jet_rap =0;
      jet_eta =0;
      jet_theta =0;
      jet_p =0;
      jet_e=0;
      jet_z=0;
      jet_z_up=0;
      jet_z_down=0;
      //empty vectors
      h_z.clear();
      h_j.clear();
      h_pid.clear();
      h_eta.clear();
      h_rap.clear();
      h_pt.clear();
      h_charge.clear(); 
      jet_lab_pt = 0;
      jet_lab_eta =0;
      jet_lab_e =0;
      nconstituents=0;
      n_charged=0;


      
    }//end loop over jets


    // T->Fill(); //fill ttree
  }  // End of event loop.

  pythia.stat();
  // Statistics. Histograms.
  //  Write tree.
  T->Print();
  T->Write();
  delete file;

  // Done.
  return 0;
}

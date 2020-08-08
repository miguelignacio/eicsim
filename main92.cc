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

  
  int nEvent    = 1e6;
  int    power   = -1;     // -1 = anti-kT; 0 = C/A; 1 = kT.
  double R       =  1.0;  // Jet size.
  double pTMin   = 5.0;
  double etaMax  = 4.0;    // Pseudorapidity range of detector.

  double eProton   = 100.0;
  double eElectron = 5;
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
  //pythia.readString("SpaceShower:pTmaxMatch=2");
  pythia.readString("PDF:lepton=off");
  pythia.readString("TimeShower:QEDshowerByL=off");
  pythia.readString("SpaceShower:dipoleRecoil=on");
  pythia.init();

  fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme ;   //WTA_modp_sche
  fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, R,power);
  fastjet::JetDefinition jetDefBF( fastjet::ee_genkt_algorithm, R , power, recomb_scheme );
  std::vector <fastjet::PseudoJet> fjInputs;
  std::vector <fastjet::PseudoJet> fjInputs_breit;

  // Set up the ROOT TFile and TTree.
  TFile *file = TFile::Open("pytree.root","recreate");
  TTree *T = new TTree("T","jet Tree");


  //Tree variables:
  UInt_t ntrials, evid, quark_id;
  Float_t xsec, x, y, Q2, W2;
  Float_t e_pt, e_phi, e_rap, e_eta, e_theta, e_p;
  Float_t quark_e, quark_eta, quark_pt, quark_p, quark_phi, quark_theta;


  std::vector<UInt_t> nconstituents;
  std::vector<UInt_t> n_charged;

  std::vector<float> jet_pt;
  std::vector<float> jet_qt;
  std::vector<float> jet_qt_breit;
  std::vector<float> jet_phi;
  std::vector<float> jet_rap;
  std::vector<float> jet_eta;
  std::vector<float> jet_theta;
  std::vector<float> jet_p;
  std::vector<float> jet_e;
  std::vector<float> jet_z;
  std::vector<float> jet_z_jet;
  std::vector<float> jet_z_e;
  std::vector<float> jet_charge;
  
  std::vector<float> dphi_e_jet;
  std::vector<float> dR_q_jet;
  std::vector<float> dphi_q_jet;
  std::vector<float> deta_q_jet;

  std::vector<float> h_z;
  std::vector<float> h_j;
  std::vector<float> h_pid;
  std::vector<float> h_eta;
  std::vector<float> h_rap;
  std::vector<float> h_pt;
  std::vector<float> h_p;
  std::vector<float> h_theta;

  std::vector<float> in_jet_pt;
  std::vector<float> in_jet_eta;
  std::vector<float> in_jet_dphi;
  std::vector<float> in_jet_qt;
  
  std::vector<double> h_charge;

  T->Branch("ntrials", &ntrials, "ntrials/I");
  T->Branch("evid", &evid, "evid/I");
  T->Branch("xsec", &xsec, "xsec/F");
  T->Branch("x", &x, "x/F");
  T->Branch("y", &y, "y/F");
  T->Branch("Q2", &Q2, "Q2/F");
  T->Branch("W2", &W2, "W2/F");

  //jet variables
  T->Branch("e_pt", &e_pt, "e_pt/F");
  T->Branch("e_phi", &e_phi, "e_phi/F");
  T->Branch("e_rap",&e_rap, "e_rap/F");
  T->Branch("e_eta", &e_eta, "e_eta/F");
  T->Branch("e_p", &e_p, "e_p/F");
  T->Branch("e_theta", &e_theta, "e_theta/F");

  T->Branch("quark_id", &quark_id, "quark_id/I");
  T->Branch("quark_e", &quark_e, "quark_e/F");
  T->Branch("quark_eta", &quark_eta, "quark_eta/F");
  T->Branch("quark_pt", &quark_pt,"quark_pt/F");
  T->Branch("quark_p",&quark_p,"quark_p/F");
  T->Branch("quark_phi", &quark_phi,"quark_phi/F");
  T->Branch("quark_theta", &quark_theta,"quark_theta/F");

  T->Branch("n_total",&nconstituents);
  T->Branch("n_charged", &n_charged);
  T->Branch("jet_pt", &jet_pt);
  T->Branch("jet_qt", &jet_qt);
  T->Branch("jet_qt_breit", &jet_qt_breit);

  T->Branch("jet_phi", &jet_phi);
  T->Branch("jet_rap",&jet_rap);
  T->Branch("jet_eta", &jet_eta);
  T->Branch("jet_p", &jet_p);
  T->Branch("jet_e", &jet_e);
  T->Branch("jet_z", &jet_z);
  T->Branch("jet_z_jet", &jet_z_jet);
  T->Branch("jet_z_e", &jet_z_e);
  T->Branch("jet_charge", &jet_charge);
  
  T->Branch("jet_theta", &jet_theta);
  T->Branch("dphi_e_jet", &dphi_e_jet);
  T->Branch("dR_q_jet",&dR_q_jet);
  T->Branch("dphi_q_jet",&dphi_q_jet);
  T->Branch("deta_q_jet",&deta_q_jet);

  //hadron variables
  T->Branch("z", &h_z);
  T->Branch("jt", &h_j);
  T->Branch("pid", &h_pid);
  T->Branch("eta", &h_eta);
  T->Branch("rap", &h_rap);
  T->Branch("pt", &h_pt);
  T->Branch("p", &h_p);
  T->Branch("theta", &h_theta);

  T->Branch("charge",&h_charge);
  T->Branch("in_jet_pt", &in_jet_pt);
  T->Branch("in_jet_eta", &in_jet_eta);
  T->Branch("in_jet_dphi",&in_jet_dphi);
  T->Branch("in_jet_qt", &in_jet_qt);
  
  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    fjInputs.clear();
    fjInputs_breit.clear();
    //empty vectors
    h_z.clear();
    h_j.clear();
    h_pid.clear();
    h_eta.clear();
    h_rap.clear();
    h_pt.clear();
    h_p.clear();
    h_theta.clear();

    h_charge.clear();
    in_jet_pt.clear();
    in_jet_eta.clear();
    in_jet_dphi.clear();
    in_jet_qt.clear();
      
    jet_pt.clear();
    jet_qt.clear();
    jet_qt_breit.clear();


    jet_phi.clear();
    jet_rap.clear();
    jet_eta.clear();
    jet_theta.clear();
    jet_p.clear();
    jet_e.clear();
    jet_z.clear();
    jet_z_e.clear();
    jet_z_jet.clear();
    jet_charge.clear();
    
    dphi_e_jet.clear();
    dR_q_jet.clear();
    dphi_q_jet.clear();
    deta_q_jet.clear();

    e_pt = 0;
    e_phi = 0;
    e_rap = 0;
    e_eta = 0;
    e_theta = 0;
    e_p = 0;

    quark_id = 0;
    quark_e = 0;
    quark_pt = 0;
    quark_eta = 0;
    quark_phi = 0;
    quark_p = 0;
    quark_theta = 0;

    nconstituents.clear();
    n_charged.clear();

    //general event info
    evid = iEvent;
    xsec = pythia.info.sigmaGen();
    ntrials = pythia.info.nTried();

    // four-momenta of proton, electron, virtual photon/Z^0/W^+-.
    Vec4 pProton = event[1].p();
    Vec4 peIn = event[4].p();
    Vec4 peOut = event[6].p();
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
    //std::cout << " boosted_gamma" << boosted_gamma.px() << " " << boosted_gamma.py() << " " << boosted_gamma.pz() << " " << boosted_gamma.e() << std::endl;
    //std::cout << " boosted_proton " << boosted_proton.px() << " " << boosted_proton.py() << " " << boosted_proton.pz() << " " << boosted_proton.e() << std::endl;  
     //std::cout << "Q2 " << Q2  << " Q " << TMath::Sqrt(Q2) << std::endl;
    //Q is the pz of the photon in the breit frame: check.

    // get struck quark index
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
    //initial_quark = boost(initial_quark, boost_vector);
    //initial_quark = rotateZ(initial_quark, -phi_p);
    //initial_quark = rotateY(initial_quark, -theta_p);

    quark_id = event[q].id();
    quark_pt = quark.pt();
    quark_eta = quark.eta();
    quark_phi = quark.phi_std();
    quark_p = sqrt(quark.modp2());
    quark_theta = acos( event[q].pz() /sqrt(quark.modp2()));

    // std::cout << " initial_quark_pt " << initial_quark.pt() << std::endl;
    //std::cout << initial_quark.px() << " " << initial_quark.py() << " " << initial_quark.pz() << " " << initial_quark.e() << std::endl;

    //std::cout << " struck quark_pt " << quark.pt() << std::endl;
    //std::cout << quark.px() << " " << quark.py() << " " << quark.pz() << " " << quark.e() << std::endl;

    Vec4   pTemp;
    double mTemp;
    int nAnalyze = 0;
    //loop over particles in the event and store them as input for FastJet
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {

	// Require visible/charged particles inside detector.
	if ( !event[i].isVisible()  ) continue;
	if ( event[i].mother1()==6 ) continue; //remove scattered lepton
	if(event[i].id()==11) continue;
	//if (abs(event[i].eta()) > etaMax) continue;

	// Create a PseudoJet from the complete Pythia particle.
	fastjet::PseudoJet particle_lab(event[i].px(), event[i].py(), event[i].pz(),event[i].e());
	if (particle_lab.pt()<0.100) continue;
        if(abs(particle_lab.eta())>4.0)continue;
	// Define Breit Fram Kinematics

	
	fastjet::PseudoJet boosted_particle = boost(particle_lab,boost_vector);
	boosted_particle = rotateZ(boosted_particle, -phi_p);
	boosted_particle = rotateY(boosted_particle, -theta_p);

	particle_lab.set_user_info(new MyUserInfo(event[i].id(),i,event[i].charge()));
	boosted_particle.set_user_info(new MyUserInfo(event[i].id(),i,event[i].charge()));

	
	fjInputs.push_back(particle_lab);
	fjInputs_breit.push_back(boosted_particle);

	++nAnalyze;
      } //end loop over particles


    fastjet::PseudoJet electron(event[6].px(), event[6].py(), event[6].pz(),event[6].e());
    e_pt = electron.pt();
    e_phi = electron.phi_std();
    e_rap= electron.rap();
    e_eta=electron.eta();
    e_theta=  acos( event[6].pz() /sqrt(electron.modp2())) ;
    e_p = sqrt(electron.modp2());


    // Run Fastjet algorithm and sort jets in pT order.
    vector <fastjet::PseudoJet> inclusiveJets, sortedJets;

    //do Breit-frame
    //fastjet::ClusterSequence clustSeq(fjInputs_breit, jetDefBF);
    // inclusiveJets = clustSeq.inclusive_jets();
    // sortedJets    = sorted_by_E(inclusiveJets);
    // do lab frame:
    fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
    inclusiveJets  = clustSeq.inclusive_jets();
    sortedJets     = sorted_by_pt(inclusiveJets);
    
    
    //loop over jets
    //std::cout << " //////////////////////////// " << std::endl;
    for (unsigned ijet= 0; ijet < sortedJets.size();ijet++) {
      fastjet::PseudoJet jet = sortedJets[ijet];

      // std::cout << " ept" << e_pt << " jet " << jet.perp() << std::endl;
      //std::cout << "jet eta " << jet.rap() << " jet e " << jet.e() << std::endl;
      if(jet.perp()<5.0) continue;
      vector<fastjet::PseudoJet> constituents = jet.constituents();
      fastjet::PseudoJet hardest = fastjet::SelectorNHardest(1)(constituents)[0];
      vector<fastjet::PseudoJet> neutral_hadrons  =
	( fastjet::SelectorIsHadron() && fastjet::SelectorIsNeutral())(constituents);
      double neutral_hadrons_pt = join(neutral_hadrons).perp();
      int ncharged_constituents = fastjet::SelectorIsCharged().count(constituents);
      double nphotons_constituents = fastjet::SelectorId(22).count(constituents);
      // if(constituents.size()<2 and fastjet::SelectorId(11).count(constituents)==1) continue; //remove jets with 1 particle = electron


      nconstituents.push_back(constituents.size());
      n_charged.push_back(ncharged_constituents);
      double jetpt = jet.perp();
      double jetrap = jet.rap();
      double jetphi = jet.phi_std();

      jet_pt.push_back( jet.perp());
      jet_phi.push_back( jet.phi_std());
      //jet_rap.push_back( jet.rap());
      //jet_eta.push_back( jet.eta());


      dphi_e_jet.push_back(TMath::Abs(TVector2::Phi_mpi_pi(jet.phi_std()-e_phi)));


      Float_t deta = jet.eta()-quark_eta;
      Float_t dphi = TMath::Abs(TVector2::Phi_mpi_pi(jet.phi_std()-quark_phi));
      Float_t dR = TMath::Sqrt(dphi*dphi + deta*deta);

      dR_q_jet.push_back(dR);
      dphi_q_jet.push_back( TMath::Abs(TVector2::Phi_mpi_pi(jet.phi_std()-quark_phi)));
      deta_q_jet.push_back( jet.eta()-quark_eta);

      double pxj, pyj, pzj, p_jet, e_jet;

      pxj = jet.px();
      pyj = jet.py();
      pzj = jet.pz();
      p_jet = sqrt(jet.modp2());
      e_jet  = jet.e();

      // std::cout << "jet " << pxj << " " << pyj << " " << pzj << " " << e_jet << std::endl;
      //double z_jet = 
      //double z_breit = z_e;//2*e_jet/sqrt(Q2);
      //std::cout << " z breit " << z_breit << std::endl;

      //TVector2 pj_breit(pxj, pyj);
      //TVector2 quark_breit(initial_quark.px(), initial_quark.py());
      //TVector2 q_breit = pj_breit/z_breit + quark_breit;
      //jet_qt_breit.push_back(q_breit.Mod());
      // std::cout << " qt breit " << q_breit.Mod() << std::endl;

      TVector2 jet2(jet.px(),jet.py());
      TVector2 electron2(electron.px(), electron.py());
      TVector2 qT = jet2 + electron2;

      jet_qt_breit.push_back(qT.Mod());
      jet_qt.push_back(qT.Mod());

      double theta = jet.pz()/p_jet;
      jet_p.push_back(0);
      jet_e.push_back(0);
      jet_z.push_back(0);
      jet_z_e.push_back(0);
      jet_z_jet.push_back(0);
      
      jet_rap.push_back(jet.rap());
      jet_eta.push_back(jet.eta());
      jet_theta.push_back(acos(theta));


     
      double Q = 0;

      //loop over constituents
      if(constituents.size()<2 and fastjet::SelectorId(11).count(constituents)==1) continue;
      for (unsigned n = 0; n < constituents.size(); n++){
	fastjet::PseudoJet _p = constituents[n]; //.user_info<FJUtils::PythiaUserInfo>().getParticle();
	if(_p.user_info<MyUserInfo>().charge()!=0){

	  double pxh, pyh, pzh, cross;
	  pxh = _p.px();
	  pyh = _p.py();
	  pzh = _p.pz();

	  //int q = _p.user_info<MyUserInfo>().charge();
	  double zh = _p.pt()/jet.perp();
	  Q += pow(zh,0.3)*_p.user_info<MyUserInfo>().charge(); 
	  cross = sqrt( pow((pyj*pzh-pyh*pzj),2.0) + pow((pxj*pzh-pzj*pxh),2.0) + pow((pxj*pyh-pyj*pxh),2.0) );
	  h_z.push_back((pxj*pxh + pyj*pyh + pzj*pzh) / (p_jet*p_jet));
	  h_j.push_back(cross / p_jet);
	  h_pid.push_back(_p.user_info<MyUserInfo>().pdg_id());
	  h_charge.push_back(_p.user_info<MyUserInfo>().charge());
	  h_rap.push_back(_p.rap());
	  h_eta.push_back(_p.eta());
	  h_pt.push_back(_p.pt());
	  h_p.push_back( sqrt(_p.modp2() ));
	  h_theta.push_back( acos( _p.pz() /sqrt(_p.modp2()))) ;

	  in_jet_pt.push_back(jet.perp());
	  in_jet_eta.push_back(jet.eta());
	  in_jet_dphi.push_back(TMath::Abs(TVector2::Phi_mpi_pi(jet.phi_std()-e_phi)));
	  in_jet_qt.push_back(qT.Mod());
	}
      }//end loop over constituents

      jet_charge.push_back(Q);

    }//end loop over jets


    T->Fill(); //fill ttree
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

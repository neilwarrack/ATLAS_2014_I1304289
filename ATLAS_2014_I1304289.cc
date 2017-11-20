// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
// #include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Tools/JetUtils.hh" 
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/WFinder.hh"

namespace Rivet {


  /// Parton level differential top quark pair production cross sections in pp collisions at 
  /// centre-of-mass energy of 7 TeV, measured in the lepton + jets channel (decays via Taus are included).
  class ATLAS_2014_I1304289 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2014_I1304289);

    
    /// @name Analysis methods (?!)
    //@{

    /// Book histograms and initialise projections before the run
    void init() {


      // declare partonic tops with leptonic and hadronic decay mode, includes by default decays via taus. 
      declare(PartonicTops(PartonicTops::E_MU),     "LeptonicPartonTops");
      declare(PartonicTops(PartonicTops::HADRONIC), "HadronicPartonTops");

      // declare jets (using the 'anti-k_t' algorithm, radius parameter R=0.4)
      const FinalState fs;
      FastJets fj04 (fs, FastJets::ANTIKT, 0.4); 
      declare(fj04, "AntiKt04");
      //      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "AntiKt04");
      //      Cut cuts = Cuts::pT > 35*GeV;
      // WFinder wboson(fs, cuts);
      
      // Identify muons
      IdentifiedFinalState muonfs ( Cuts::abseta < 2.5 && Cuts::pT > 25*GeV ) ;
      muonfs.acceptId(PID::MUON); // what does this line do?
      declare(muonfs, "Muon");
      
      
      // exclude electrons in gap between barrel and endcap calorimeters
      IdentifiedFinalState electronfs((Cuts::abseta > 1.52 && Cuts::abseta < 2.47 && Cuts::pT > 15*GeV)
				      || (Cuts::abseta < 1.37 && Cuts::pT > 15*GeV));
      electronfs.acceptId(PID::ELECTRON);
      declare(electronfs, "Electron");
      
      // declare mising energy projection
      addProjection(MissingMomentum(fs), "MissingMomenta");



      // Book histos
      _hSL_hadronicTopPt   = bookHisto1D(1,1,1);
      _hSL_ttbarMass       = bookHisto1D(2,1,1);
      _hSL_topPtTtbarSys   = bookHisto1D(3,1,1);
      _hSL_topAbsYTtbarSys = bookHisto1D(4,1,1);
    }


    void analyze(const Event& event) {

      //      const WFinder& wboson = apply<WFinder>(event, "WBoson");
      //if (wboson.empty()) vetoEvent ;
     
      // Find tops
      const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();
      const Particles hadronicpartontops = apply<ParticleFinder>(event, "HadronicPartonTops").particlesByPt();

      //find muons & electrons
      const Particles muons     = apply<IdentifiedFinalState>(event, "Muon"    ).particlesByPt();
      const Particles electrons = apply<IdentifiedFinalState>(event, "Electron").particlesByPt();      

      // Veto event if it does not "pass either a single-electron or single-muon trigger"
      /*
      const bool isSemiLeptonic = (leptonicpartontops.size() == 1 && hadronicpartontops.size() == 1 );
      if ( !isSemiLeptonic ) vetoEvent;
      */
      const bool hasSingleLepton = ((muons.size() + electrons.size()) == 1);
      if ( !hasSingleLepton ) vetoEvent;

      // Use only jets with p_T>25GeV and abs(eta) < 2.5
      Jets jets = apply<FastJets>(event, "AntiKt04").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);

      // Missing energy cut
      const MissingMomentum& misMom = applyProjection<MissingMomentum>(event, "MissingMomenta");
      double Pmiss = misMom.missingMomentum().pT(); 
      if (Pmiss<=30*GeV) vetoEvent;
      //      double missingp_TPhi = misMom.azimuthalAngle();
      FourMomentum mis4mom;
      mis4mom = misMom.visibleMomentum();
      double missingp_TPhi = mis4mom.azimuthalAngle();     

      double wmt, leptonPhi, leptonPT;

      if (muons.size() == 1){
	leptonPhi = muons[0].phi();
	leptonPT = muons[0].pT();
      } else { // if (electrons.size() == 1) 
	leptonPhi = electrons[0].phi();
	leptonPT = electrons[0].pT();
      }
      wmt = sqrt( 2*leptonPT * Pmiss * (1-cos(leptonPhi - missingp_TPhi)));      
      const bool wmassAboveThreashold = ( wmt > 35*GeV );
      if (wmassAboveThreashold) vetoEvent;

      // Sort electrons and jets
      Particles isolated_electrons;
      Particles isolated_muons;
      Jets isolated_jets;
      
      // Discard jets within cone of R=.2 of an electron with p_T>15GeV
      for (const Jet j : jets) {
	for  (const Particle& electron : electrons) {
	  if (deltaR(electron, j) >= 0.2) continue ;
	}
	isolated_jets.push_back(j);
      }
      
      // Discard electrons within cone of R=.4 of an isolated jet
      for  (const Particle& electron : electrons) {
	for (const Jet isol_j : isolated_jets) {
	  if (deltaR(electron, isol_j) >= 0.4) continue;
	}
	isolated_electrons.push_back(electron);     
      }
      
      // Keep only events with >4 isolated jets (of which one must be b-tagged). 
      if (isolated_jets.size() < 4) vetoEvent;
      if (!any(jets, hasBTag())) vetoEvent;

      // remove muons within R=.4 of isolated jets
      for  (const Particle& muon : muons) {
	for (const Jet isol_j : isolated_jets) {
	  if (deltaR(muon, isol_j) >= 0.4) continue;
	}
	isolated_muons.push_back(muon);     
      }
      
      // Parton level at full phase space
      // Fill top quarks defined in the parton level, full phase space
      const FourMomentum tLP4 = leptonicpartontops[0];
      const FourMomentum tHP4 = hadronicpartontops[0];
      const FourMomentum ttbarP4 = tLP4 + tHP4;
      
      const double weight = event.weight();
      _hSL_hadronicTopPt->fill(tHP4.pT(), weight);
      _hSL_ttbarMass->fill(ttbarP4.mass(), weight);
      _hSL_topPtTtbarSys->fill(ttbarP4.pT(), weight);
      _hSL_topAbsYTtbarSys->fill(ttbarP4.absrap(), weight);
    }
    
    
    /// Normalise histograms
    void finalize() {
      scale({_hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys, _hSL_topAbsYTtbarSys}, crossSection()/picobarn/sumOfWeights());
    }
    
    /// @name Histograms
    Histo1DPtr _hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys, _hSL_topAbsYTtbarSys;
    
  };
  
  
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1304289);
  
  
}

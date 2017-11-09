// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
// #include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Tools/JetUtils.hh" // is this correct ??
#include "Rivet/Projections/MissingMomentum.hh"
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

      // ----- from ALEPH_2016_I1492968 -------
      // declare mising energy projection
      addProjection(MissingMomentum(fs), "MissingMomenta");



      // Book histos
      _hSL_hadronicTopPt   = bookHisto1D("d01-x01-y01");
      _hSL_ttbarMass       = bookHisto1D("d02-x01-y01");
      _hSL_topPtTtbarSys   = bookHisto1D("d03-x01-y01");
      _hSL_topAbsYTtbarSys = bookHisto1D("d05-x01-y01");
    }


    void analyze(const Event& event) {

      // Find tops
      const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();
      const Particles hadronicpartontops = apply<ParticleFinder>(event, "HadronicPartonTops").particlesByPt();

      // Veto all non-semileptonic events     
      const bool isSemiLeptonic = (leptonicpartontops.size() == 1 && hadronicpartontops.size() == 1 );
      if ( !isSemiLeptonic ) vetoEvent;
      

      // Keep jets with p_T>25GeV and abs(eta) < 2.5
      // Jets jets = apply<FastJets>(event, "AntiKt04").jetsByPt(Cuts::pT > 25*GeV, Cuts::abseta < 2.5); // ERROR: doesn't expect 2 arguments!
      Jets jets = apply<FastJets>(event, "AntiKt04").jetsByPt(Cuts::pT > 25*GeV);

      // Keep only events with >4 jets (of which one must be b-tagged). 
      if (jets.size() < 4) vetoEvent;
      if (!any(jets, hasBTag())) vetoEvent;
     

      // Missing energy cut
      // ------- (FROM ALEPH_2016_I1492968.cc) -----------------
      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MissingMomenta");
      double Pmiss = met.missingMomentum().pT();   // ALEPH analysis uses p(), not pT()
      if (Pmiss>30*GeV) vetoEvent;
      // -------------------------------------------------------


    

 
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
      normalize({_hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys, _hSL_topAbsYTtbarSys});
    }
    
    /// @name Histograms
    Histo1DPtr _hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys, _hSL_topAbsYTtbarSys;

};
  
  
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1304289);
  
  
}

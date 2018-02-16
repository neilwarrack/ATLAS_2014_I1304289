// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PartonicTops.hh"



namespace Rivet {



  /// Parton level differential top quark pair production cross sections in pp collisions at 
  /// centre-of-mass energy of 7 TeV, measured in the lepton + jets channel (decays via Taus are included).
  class ATLAS_2014_I1304289 : public Analysis {
  public:


    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2014_I1304289);

    
    void init() {

      // declare tops with leptonic and hadronic decay mode, includes by default decays via taus. 
      declare(PartonicTops(PartonicTops::E_MU_TAU), "LeptonicPartonTops");
      declare(PartonicTops(PartonicTops::HADRONIC), "HadronicPartonTops");

      // Book histos
      _hSL_hadronicTopPt   = bookHisto1D(1,1,1);
      _hSL_ttbarMass       = bookHisto1D(2,1,1);
      _hSL_topPtTtbarSys   = bookHisto1D(3,1,1);
      _hSL_topAbsYTtbarSys = bookHisto1D(4,1,1);
    }


    void analyze(const Event& event) {

      // Find tops & veto if not semileptonic.
      const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();
      const Particles hadronicpartontops = apply<ParticleFinder>(event, "HadronicPartonTops").particlesByPt();
      const bool isSemileptonic  = (leptonicpartontops.size() == 1 && hadronicpartontops.size() == 1 );


      
      if ( leptonicpartontops.size() + hadronicpartontops.size() == 1 ) {   cout<<"1 "; }
      if ( leptonicpartontops.size() + hadronicpartontops.size() == 0 ) {   cout<<"2 "; }
      if ( leptonicpartontops.size() + hadronicpartontops.size() >= 3 ) {   cout<<"3 "; }
      //      if ( !isSemileptonic && !isDileptonic && !isFullyHadronic ) {   cout<<"4 "; }
      if (leptonicpartontops.size() == 1 ) {cout << "5 "; }
     
      if ( !isSemileptonic) { cout << "8" << endl; vetoEvent ;}     



      
      // Fill top quarks defined in the parton level, full phase space
      if ( isSemileptonic ){
      const FourMomentum t1P4 = leptonicpartontops[0] ;
      const FourMomentum t2P4 = hadronicpartontops[0] ;
      const FourMomentum ttbarP4 = t1P4 + t2P4 ;
      const double weight = event.weight() ;

      cout << "9" ; // SUCCESS!!

      //_hSL_hadronicTopPt->fill(t1P4.pT(), weight) ;
      _hSL_hadronicTopPt->fill(t2P4.pT(), weight) ;
      _hSL_ttbarMass->fill(ttbarP4.mass(), weight) ;
      _hSL_topPtTtbarSys->fill(ttbarP4.pT(), weight) ;
      _hSL_topAbsYTtbarSys->fill(ttbarP4.absrap(), weight) ;      
      }
 cout << endl;
    
    }

    
    /// Normalise histograms
    void finalize() {

      //const float BR = 0.438 ; // branching ratio of ttbar -> l+jets channel
      //double scale_factor = crossSection()*picobarn/sumOfWeights() ;
      //double scale_factor01 = crossSection()*picobarn/sumOfWeights()/1000 ;

      //scale({_hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys}, scale_factor) ;
      //scale({_hSL_topAbsYTtbarSys}, scale_factor01) ;

      // Unit-normalized
      normalize({_hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys}, 1000); //< Factor from x/y unit mismatch in HepData
      normalize(_hSL_topAbsYTtbarSys);
      // xsec-normalized
      //scale( {_hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys, _hSL_topAbsYTtbarSys}, 
      //     crossSection()/picobarn/sumOfWeights()/0.438);

      //MSG_INFO(crossSection()) ;    
    }
    
    /// @name Histograms
    Histo1DPtr _hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys, _hSL_topAbsYTtbarSys ;
    
  } ;
  
  
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1304289) ;
  
  
}

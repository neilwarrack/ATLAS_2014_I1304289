#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PartonicTops.hh"



namespace Rivet {


  /// Parton-level top-pair normalized differential distributions at 7 TeV. 
  class ATLAS_2014_I1304289 : public Analysis {
  public:


    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2014_I1304289) ;

    
    void init() {

      // declare tops with leptonic and hadronic decay mode, includes by default decays via taus. 
      declare(PartonicTops(PartonicTops::E_MU), "LeptonicPartonTops") ;
      declare(PartonicTops(PartonicTops::HADRONIC), "HadronicPartonTops") ;

      // Book histos
      _hSL_hadronicTopPt   = bookHisto1D(1,1,1) ;
      _hSL_ttbarMass       = bookHisto1D(2,1,1) ;
      _hSL_topPtTtbarSys   = bookHisto1D(3,1,1) ;
      _hSL_topAbsYTtbarSys = bookHisto1D(4,1,1) ;
    }


    void analyze(const Event& event) {

      // Find tops & veto if not semileptonic.
      const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt() ;
      const Particles hadronicpartontops = apply<ParticleFinder>(event, "HadronicPartonTops").particlesByPt() ;
      const bool semileptonic  = (leptonicpartontops.size() == 1 && hadronicpartontops.size() == 1 ) ;
      if ( !semileptonic) vetoEvent ;     
      
      // Fill top quarks defined in the parton level, full phase space
      const FourMomentum t1P4 = leptonicpartontops[0] ;
      const FourMomentum t2P4 = hadronicpartontops[0] ;
      const FourMomentum ttbarP4 = t1P4 + t2P4 ;
      const double weight = event.weight() ;


      _hSL_hadronicTopPt->fill(t2P4.pT(), weight) ;
      _hSL_ttbarMass->fill(ttbarP4.mass(), weight) ;
      _hSL_topPtTtbarSys->fill(ttbarP4.pT(), weight) ;
      _hSL_topAbsYTtbarSys->fill(ttbarP4.absrap(), weight) ;      
    }

    
    /// Normalise histograms
    void finalize() {


      // Unit-normalized
      normalize({_hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys}, 1000) ; //< Factor from x/y unit mismatch in HepData
      normalize(_hSL_topAbsYTtbarSys) ;

    }

    
    /// @name Histograms
    Histo1DPtr _hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys, _hSL_topAbsYTtbarSys ;
    
  } ;
  
  
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1304289) ;
  
  
}

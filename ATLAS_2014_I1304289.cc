// -*- C++ -*-
#include "Rivet/Analysis.hh"
//#include "Rivet/Projections/FinalState.hh"
// #include "Rivet/Projections/VetoedFinalState.hh"
//#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
//#include "Rivet/Tools/JetUtils.hh" 
//#include "Rivet/Projections/MissingMomentum.hh"
//#include "Rivet/Projections/IdentifiedFinalState.hh"
//#include "Rivet/Projections/WFinder.hh"



namespace Rivet {



  /// Parton level differential top quark pair production cross sections in pp collisions at 
  /// centre-of-mass energy of 7 TeV, measured in the lepton + jets channel (decays via Taus are included).
  class ATLAS_2014_I1304289 : public Analysis {
  public:



    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2014_I1304289);

    
    /// Book histograms and initialise projections before the run
    void init() {



      // declare tops with leptonic and hadronic decay mode, includes by default decays via taus. 
      declare(PartonicTops(PartonicTops::E_MU_TAU), "LeptonicPartonTops");
      declare(PartonicTops(PartonicTops::HADRONIC), "HadronicPartonTops");

      
      // declare jets (using the 'anti-k_t' algorithm, radius parameter R=0.4)
      const FinalState fs;
      FastJets fj04 (fs, FastJets::ANTIKT, 0.4); 
      declare(fj04, "AntiKt04");
      //      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "AntiKt04");
      //      Cut cuts = Cuts::pT > 35*GeV;
      // WFinder wboson(fs, cuts);
      
      
      // Identify muons
      FinalState muonfs ( Cuts::abspid == PID::MUON && Cuts::abseta < 2.5 ) ;
      declare(muonfs, "Muon");
      
      
      // exclude electrons in gap between barrel and endcap calorimeters
      FinalState electronfs (Cuts::abspid == PID::ELECTRON && (Cuts::abseta < 1.37 || (Cuts::abseta > 1.52 && Cuts::abseta < 2.47)));
      declare(electronfs, "Electron");
      
      
      // declare mising energy projection
      addProjection(MissingMomentum(fs), "MissingMomenta") ;

    
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
      const bool isDileptonic    = (leptonicpartontops.size() == 2 && hadronicpartontops.size() == 0 );
      const bool isFullyHadronic = (leptonicpartontops.size() == 0 && hadronicpartontops.size() == 2 );

      //if ( !isSemileptonic && !isDileptonic && !isFullyHadronic ) {   cout<<"1"<<endl; vetoEvent ; }
      //if ( !isSemileptonic && !isDileptonic) {   cout<<"1"<<endl; vetoEvent ; }
     
      
      // Fill top quarks defined in the parton level, full phase space
      if ( isSemileptonic ){
      const FourMomentum t1P4 = leptonicpartontops[0] ;
      const FourMomentum t2P4 = hadronicpartontops[0] ;
      const FourMomentum ttbarP4 = t1P4 + t2P4 ;
      const double weight = event.weight() ;

      //_hSL_hadronicTopPt->fill(t1P4.pT(), weight) ;
      _hSL_hadronicTopPt->fill(t2P4.pT(), weight) ;
      _hSL_ttbarMass->fill(ttbarP4.mass(), weight) ;
      _hSL_topPtTtbarSys->fill(ttbarP4.pT(), weight) ;
      _hSL_topAbsYTtbarSys->fill(ttbarP4.absrap(), weight) ;      
      cout <<"9"<<endl;
      }
     
      /*
      if ( isFullyHadronic){
      const FourMomentum t1P4 = hadronicpartontops[0] ;
      const FourMomentum t2P4 = hadronicpartontops[1] ;
      const FourMomentum ttbarP4 = t1P4 + t2P4 ;
      const double weight = event.weight() ;

      _hSL_hadronicTopPt->fill(t1P4.pT(), weight) ;
      //_hSL_hadronicTopPt->fill(t2P4.pT(), weight) ;
      _hSL_ttbarMass->fill(ttbarP4.mass(), weight) ;
      _hSL_topPtTtbarSys->fill(ttbarP4.pT(), weight) ;
      _hSL_topAbsYTtbarSys->fill(ttbarP4.absrap(), weight) ;      
      cout<<"9"<<endl;
      }
           
      if ( isDileptonic ){
      const FourMomentum t1P4 = leptonicpartontops[0] ;
      const FourMomentum t2P4 = leptonicpartontops[1] ;
      const FourMomentum ttbarP4 = t1P4 + t2P4 ;
      const double weight = event.weight() ;

      _hSL_hadronicTopPt->fill(t1P4.pT(), weight) ; // the leading leptonic top
      //_hSL_hadronicTopPt->fill(t2P4.pT(), weight) ;
      _hSL_ttbarMass->fill(ttbarP4.mass(), weight) ;
      _hSL_topPtTtbarSys->fill(ttbarP4.pT(), weight) ;
      _hSL_topAbsYTtbarSys->fill(ttbarP4.absrap(), weight) ;      
      cout <<"9"<<endl;
      }
      */
      /*
      if ( !isSemileptonic && !isDileptonic && !isFullyHadronic && (leptonicpartontops.size() == 1)){
      const FourMomentum t1P4 = leptonicpartontops[0] ;
      const double weight = event.weight() ;
      _hSL_hadronicTopPt->fill(t1P4.pT(), weight) ;
      cout <<"5"<<endl;
      }

      if ( !isSemileptonic && !isDileptonic && !isFullyHadronic && (hadronicpartontops.size() == 1)){
      const FourMomentum t1P4 = hadronicpartontops[0] ;
      const double weight = event.weight() ;
      _hSL_hadronicTopPt->fill(t1P4.pT(), weight) ;
      cout <<"6"<<endl;
      }
      
      */
      
      if ( leptonicpartontops.size() + hadronicpartontops.size() == 1 ) {   cout<<"1"; }
      if ( leptonicpartontops.size() + hadronicpartontops.size() == 0 ) {   cout<<"2"; }
      if ( leptonicpartontops.size() + hadronicpartontops.size() >= 3 ) {   cout<<"3"; }
      if ( !isSemileptonic && !isDileptonic && !isFullyHadronic ) {   cout<<"4"; }
      if (leptonicpartontops.size() == 1 ) {cout << "5"; }
      cout << endl;
      
      //      MSG_INFO(100) ;
    
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

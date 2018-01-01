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

    
    /// Book histograms and initialise projections before the run
    void init() {



      // declare tops with leptonic and hadronic decay mode, includes by default decays via taus. 
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
      //    cout << "test" << endl;
      // Find tops & veto if not semileptonic.
      const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();
      const Particles hadronicpartontops = apply<ParticleFinder>(event, "HadronicPartonTops").particlesByPt();
      const bool isSemiLeptonic = (leptonicpartontops.size() == 1 && hadronicpartontops.size() == 1 );
      //      if ( !isSemiLeptonic ) { MSG_INFO(0.1) ; vetoEvent ; }
      if ( !isSemiLeptonic ) {   cout<<"1"<<endl; vetoEvent ; }
     

      // Find leptons & jets; build vectors.
      const Particles muons     = apply<IdentifiedFinalState>(event, "Muon"    ).particlesByPt();
      const Particles electrons = apply<IdentifiedFinalState>(event, "Electron").particlesByPt();      
      const Jets jets = apply<FastJets>(event, "AntiKt04").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);



      // Isolate electrons and jets.
      // Discard jets within cone of R=0.2 of electron
      Jets isolated_jets;
      bool tooClose = false ;
      for (const Jet j : jets) {
	for (const Particle& electron : electrons) {
	  if (deltaR(electron, j) <= 0.2) tooClose = true ; 
	}
	if ( !tooClose ) isolated_jets.push_back(j) ;
	tooClose = false ; // reset flag
      }

      /// Discard electrons within cone of R=.4 of an isolated jet
      Particles isolated_electrons;
      for (const Particle& electron : electrons) {
	for (const Jet ije : isolated_jets) {
	  if (deltaR(ije, electron) <= 0.4) tooClose = true ;
	}
	if ( !tooClose ) isolated_electrons.push_back(electron);     
	tooClose = false ; // reset flag     
      }

      /// remove muons within R=.4 of isolated jets
      Particles isolated_muons;
      for (const Particle& muon : muons) {
	for (const Jet ijm : isolated_jets) {
	  if (deltaR(ijm, muon) <= 0.4) tooClose = true ;
	}
	if ( !tooClose ) isolated_muons.push_back(muon);     
	tooClose = false ; // reset flat    
      }


      // Require event to contain exactly one isolated lepton which fired the trigger
      const bool hasSingleIsolatedLepton = ((isolated_muons.size() + isolated_electrons.size()) == 1);
      //      if ( !hasSingleIsolatedLepton ) { MSG_INFO(0.5) ; vetoEvent ; }
      if ( !hasSingleIsolatedLepton ) { cout<<"2"<<endl ; vetoEvent ; }



      // Missing energy cut.
      const MissingMomentum& misMom = applyProjection<MissingMomentum>(event, "MissingMomenta");
      const double Pmiss = misMom.missingMomentum().pT();
      //      if (Pmiss<=30*GeV){ MSG_INFO(0.2) ; vetoEvent ; }
      if (Pmiss<=30*GeV){  cout<<"3"<<endl ; vetoEvent ; }


      /// Keep only events with >4 isolated jets (of which one must be b-tagged). 
      //if ( isolated_jets.size() < 4 ) { MSG_INFO(0.3) ; vetoEvent ; }
      if ( isolated_jets.size() < 4 ) { cout<<"4"<<endl;         vetoEvent ; }
      //if ( !any( jets, hasBTag() ) )  { MSG_INFO(0.4) ; vetoEvent ; }
      if ( !any( jets, hasBTag() ) )  { cout << "5"<<endl ; vetoEvent ; }

      
      //// store variables to compute W boson transverse mass
      double leptonPhi = 0.0 ;
      double leptonPT  = 0.0 ;
      if (isolated_muons.size() == 1){
	//// Muon trigger threshold
	if (isolated_muons[0].pT() < 18*GeV) {
	  //	  MSG_INFO(0.61) ; vetoEvent ;
	  cout<<"6"<<endl ; vetoEvent ;
	} else {
	  leptonPhi = isolated_muons[0].phi() ;
	  leptonPT  = isolated_muons[0].pT() ;
	}
      } else { 
	//// Electron trigger threshold
	if (isolated_electrons[0].pT() < 22*GeV) {
	  //	  MSG_INFO(0.62); vetoEvent ;
	  cout<<"7"<<endl; vetoEvent ;	
	} else {
	  leptonPhi = isolated_electrons[0].phi() ;
	  leptonPT  = isolated_electrons[0].pT() ;
	}
      }
      
      
      // Transverse W boson mass cut
      FourMomentum mis4mom ;
      mis4mom = misMom.visibleMomentum() ;
      double missingp_TPhi = mis4mom.phi() ;
      double wmt = sqrt( 2*leptonPT * Pmiss * (1-cos(leptonPhi - missingp_TPhi))) ;      
      const bool wmassAboveThreashold = ( wmt > 35*GeV ) ;
      // if (!wmassAboveThreashold) { MSG_INFO(0.7) ; vetoEvent ; }
      if (!wmassAboveThreashold) { cout<<"8"<<endl ; vetoEvent ; }
    

      // Fill top quarks defined in the parton level, full phase space

      const FourMomentum tLP4 = leptonicpartontops[0] ;
      const FourMomentum tHP4 = hadronicpartontops[0] ;
      const FourMomentum ttbarP4 = tLP4 + tHP4 ;
      const double weight = event.weight() ;

     
      _hSL_hadronicTopPt->fill(tHP4.pT(), weight) ;
      _hSL_ttbarMass->fill(ttbarP4.mass(), weight) ;
      _hSL_topPtTtbarSys->fill(ttbarP4.pT(), weight) ;
      _hSL_topAbsYTtbarSys->fill(ttbarP4.absrap(), weight) ;

      //      MSG_INFO(100) ;
      cout<<"9"<<endl ;
    }

    
    /// Normalise histograms
    void finalize() {

      const float BR = 0.438 ; // branching ratio of ttbar -> l+jets channel
      double scale_factorTeV = BR*crossSection()*picobarn*TeV/sumOfWeights() ;
      double scale_factorGeV = BR*crossSection()*picobarn/sumOfWeights() ;

      scale({_hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys}, scale_factorTeV) ;
      scale({_hSL_topAbsYTtbarSys}, scale_factorGeV) ;
   
      MSG_INFO(crossSection()) ;    
    }
    
    /// @name Histograms
    Histo1DPtr _hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys, _hSL_topAbsYTtbarSys ;
    
  } ;
  
  
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1304289) ;
  
  
}

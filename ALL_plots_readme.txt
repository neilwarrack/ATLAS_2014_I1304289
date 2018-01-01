# *************
# ***plots12***
# *************
$ run-pythia -c Top:all=on -e 7000 -n 200000 -c PhaseSpace:bias2Selection=on -c PhaseSpace:bias2SelectionPow=2 -L P8Hook_TTbarSemiLep.so -u TTbarSemiLep -o fifo.hepmc

(output cross-section: 1.676136e+02 pb)

$ rivet -a ATLAS_2014_I1304289 -x 177.3 fifo.hepmc


      const float BR = 0.438 ; // branching ratio of ttbar -> l+jets channel
      double scale_factorTeV = BR*crossSection()*picobarn*TeV/sumOfWeights() ;
      double scale_factorGeV = BR*crossSection()*picobarn/sumOfWeights() ;

      scale({_hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys}, scale_factorTeV) ;
      scale({_hSL_topAbsYTtbarSys}, scale_factorGeV) ;





# *************
# ***plots13***
# *************
run-pythia -c Top:all=on -e 7000 -n 200000 -c PhaseSpace:bias2Selection=on -c PhaseSpace:bias2SelectionPow=2 -L P8Hook_TTbarSemiLep.so -u TTbarSemiLep -o fifo.hepmc &

rivet -a ATLAS_2014_I1304289 fifo.hepmc


      const float BR = 0.438 ; // branching ratio of ttbar -> l+jets channel
      double scale_factorTeV = BR*crossSection()*picobarn*TeV/sumOfWeights() ;
      double scale_factorGeV = BR*crossSection()*picobarn/sumOfWeights() ;

      scale({_hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys}, scale_factorTeV) ;
      scale({_hSL_topAbsYTtbarSys}, scale_factorGeV) ;

# *************
# ***plots14***
# *************
run-pythia -c Top:all=on -e 7000 -n 5000 -L P8Hook_TTbarSemiLep.so -u TTbarSemiLep -o fifo.hepmc &

rivet -a ATLAS_2014_I1304289

 const float BR = 0.438 ; // branching ratio of ttbar -> l+jets channel
      double scale_factorTeV = BR*crossSection()*picobarn*TeV/sumOfWeights() ;
      double scale_factorGeV = BR*crossSection()*picobarn/sumOfWeights() ;

      scale({_hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys}, scale_factorTeV) ;
      scale({_hSL_topAbsYTtbarSys}, scale_factorGeV) ;

# *************
# ***plots15***
# *************
run-pythia -c Top:all=on -e 7000 -n 200000 -L P8Hook_TTbarSemiLep.so -u TTbarSemiLep -o fifo.hepmc &
rivet -a ATLAS_2014_I1304289 -x 177.3 fifo.hepmc > ana.log


      const float BR = 0.438 ; // branching ratio of ttbar -> l+jets channel
      double scale_factorTeV = BR*crossSection()*picobarn*TeV/sumOfWeights() ;
      double scale_factorGeV = BR*crossSection()*picobarn/sumOfWeights() ;

      scale({_hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys}, scale_factorTeV) ;
      scale({_hSL_topAbsYTtbarSys}, scale_factorGeV) ;

# *************
# ***plots16***
# *************
$ run-pythia -c Top:all=on -e 7000 -n 200000 -c PhaseSpace:bias2Selection=on -c PhaseSpace:bias2SelectionPow=2 -o fifo.hepmc


Cross Section (nb) = 1.157e-06
Event 200000 (1:13:21 elapsed)
Finished event loop at 2017-12-31 11:35:55
Cross-section = 1.671228e+02 pb


$ rivet -a ATLAS_2014_I1304289 -x 177.3 fifo.hepmc


      const float BR = 0.438 ; // branching ratio of ttbar -> l+jets channel
      double scale_factorTeV = BR*crossSection()*picobarn*TeV/sumOfWeights() ;
      double scale_factorGeV = BR*crossSection()*picobarn/sumOfWeights() ;

      scale({_hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys}, scale_factorTeV) ;
      scale({_hSL_topAbsYTtbarSys}, scale_factorGeV) ;

# *************
# ***plots17***
# *************
$ run-pythia -c Top:all=on -e 7000 -n 200000 -c PhaseSpace:bias2Selection=on -c PhaseSpace:bias2SelectionPow=2 -o fifo.hepmc


$ rivet -a ATLAS_2014_I1304289 fifo.hepmc > veto.log



Cross Section (nb) = 1.157e-06
Event 200000 (1:13:52 elapsed)
Finished event loop at 2017-12-31 15:37:34
Cross-section = 1.671228e+02 pb


      const float BR = 0.438 ; // branching ratio of ttbar -> l+jets channel
      double scale_factorTeV = BR*crossSection()*picobarn*TeV/sumOfWeights() ;
      double scale_factorGeV = BR*crossSection()*picobarn/sumOfWeights() ;

      scale({_hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys}, scale_factorTeV) ;
      scale({_hSL_topAbsYTtbarSys}, scale_factorGeV) ;

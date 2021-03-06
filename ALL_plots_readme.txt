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
# *****************
# ****plots 18*****
# *****************

Cross Section (nb) = 5.897e-08
Event 200000 (1:13:04 elapsed)
Finished event loop at 2018-01-02 13:10:40
Cross-section = 1.671899e+02 pb
Rivet run completed at 2018-01-02 13:10:40, time elapsed = 1:13:04
Histograms written to /home/ppe/n/nwarrack/public_ppe/Rivet/RivetAnalyses/MyRivetAnalyses/ATLAS_2014_I1304289/Rivet.yoda
[1]+  Done                    run-pythia -c Top:all=on -e 7000 -n 200000 -c PhaseSpace:bias2Selection=on -c PhaseSpace:bias2SelectionPow=3 -o fifo.hepmc

rivet -a ATLAS_2014_I1304289 -x 177.3 fifo.hepmc > veto.log

 const float BR = 0.438 ; // branching ratio of ttbar -> l+jets channel
      double scale_factorTeV = BR*crossSection()*picobarn*TeV/sumOfWeights() ;
      double scale_factorGeV = BR*crossSection()*picobarn/sumOfWeights() ;

      scale({_hSL_hadronicTopPt, _hSL_ttbarMass, _hSL_topPtTtbarSys}, scale_factorTeV) ;
      scale({_hSL_topAbsYTtbarSys}, scale_factorGeV) ;




199995 lines read from file
veto code 1: 140240
veto code 2: 45467
veto code 3: 2034
veto code 4: 4382
veto code 5: 46
veto code 6: 0
veto code 7: 1139
veto code 8: 587
passed: 6089

# *****************
# ****plots19******
# *****************

Cross Section (nb) = 5.857e-08
Event 200000 (1:13:48 elapsed)
Finished event loop at 2018-01-02 16:07:02
Cross-section = 1.673150e+02 pb
Rivet run completed at 2018-01-02 16:07:02, time elapsed = 1:13:48
Histograms written to /home/ppe/n/nwarrack/public_ppe/Rivet/RivetAnalyses/MyRivetAnalyses/ATLAS_2014_I1304289/Rivet.yoda
[1]+  Done                    run-pythia -c Top:all=on -e 7000 -n 200000 -c PhaseSpace:bias2Selection=on -c PhaseSpace:bias2SelectionPow=3 -L P8Hook_TTbarSemiLep.so -u TTbarSemiLep -o fifo.hepmc

rivet -a ATLAS_2014_I1304289 -x 177.3 fifo.hepmc > veto.log


199995 lines read from file
veto code 1: 62754
veto code 2: 104464
veto code 3: 4797
veto code 4: 10002
veto code 5: 123
veto code 6: 0
veto code 7: 2736
veto code 8: 1407
passed: 13701
# ***************
# ****plots20****
# ***************

Cross Section (nb) = 1.528e-06
Event 20000 (0:08:06 elapsed)
Finished event loop at 2018-01-08 17:06:16
Cross-section = 2.322857e+02 pb
Rivet run completed at 2018-01-08 17:06:16, time elapsed = 0:08:06
Histograms written to /home/ppe/n/nwarrack/public_ppe/Rivet/RivetAnalyses/MyRivetAnalyses/ATLAS_2014_I1304289/Rivet.yoda
[1]+  Done                    run-pythia -n 20000 -c Beams:eCM=7000 -c Top:all=on -c Beams:frameType=1 -c PhaseSpace:bias2Selection=on -o fifo.hepmc -c PhaseSpace:bias2SelectionPow=2 -L P8Hook_TTbarSemiLep.so -u TTbarSemiLep


rivet -a ATLAS_2014_I1304289 fifo.hepmc > ana.log


veto code 1: 6667
veto code 2: 10377
veto code 3: 473
veto code 4: 880
veto code 5: 14
veto code 6: 0
veto code 7: 258
veto code 8: 134
passed: 1181
# *************
# ***plots21***
# *************

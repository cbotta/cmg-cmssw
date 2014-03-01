import operator 
import itertools
import copy
from math import *

from ROOT import TLorentzVector, TVectorD
from ROOT import TriggerBitChecker

from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from CMGTools.RootTools.fwlite.Event import Event
from CMGTools.RootTools.statistics.Counter import Counter, Counters
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.RootTools.physicsobjects.Lepton import Lepton
from CMGTools.RootTools.physicsobjects.Photon import Photon
from CMGTools.RootTools.physicsobjects.Electron import Electron
from CMGTools.RootTools.physicsobjects.Muon import Muon
from CMGTools.RootTools.physicsobjects.Jet import Jet

from CMGTools.RootTools.utils.DeltaR import * 
from CMGTools.TTHAnalysis.leptonMVA import LeptonMVA

from CMGTools.RootTools.analyzers.TreeAnalyzerNumpy import TreeAnalyzerNumpy
from CMGTools.TTHAnalysis.analyzers.ntuple import *
import os

def var( tree, varName, type=float ):
    tree.var(varName, type)

def fill( tree, varName, value ):
    tree.fill( varName, value )

class ttHLepFRAnalyzer( TreeAnalyzerNumpy ):
    def __init__(self, cfg_ana, cfg_comp, looperName ):
        super(ttHLepFRAnalyzer,self).__init__(cfg_ana,cfg_comp,looperName) 
        self.leptonMVA = LeptonMVA("%s/src/CMGTools/TTHAnalysis/data/leptonMVA/%%s_BDTG.weights.xml" % os.environ['CMSSW_BASE'], self.cfg_comp.isMC)
        self.triggerCheckers = []
        import ROOT
        for T in self.cfg_ana.triggers:
            trigVec = ROOT.vector(ROOT.string)()
            trigVec.push_back("HLT_%s_v*" % T)
            self.triggerCheckers.append( (T.replace("_eta2p1",""), TriggerBitChecker(trigVec)) )

    def declareHandles(self):
        super(ttHLepFRAnalyzer, self).declareHandles()
        #self.handles['met'] = AutoHandle( 'cmgPFMET', 'std::vector<cmg::BaseMET>' )
        self.handles['met'] = AutoHandle( 'cmgPFMETRaw', 'std::vector<cmg::BaseMET>' )
        self.handles['nopumet'] = AutoHandle( 'nopuMet', 'std::vector<reco::PFMET>' )
        self.handles['TriggerResults'] = AutoHandle( ('TriggerResults','','HLT'), 'edm::TriggerResults' )

    def declareVariables(self):
        tr = self.tree

        isMC = self.cfg_comp.isMC 

        tr = self.tree
        var( tr, 'run', int)
        var( tr, 'lumi', int)
        var( tr, 'evt', int)
        var( tr, 'nVert', int)
        var( tr, 'nLepGood', int)
        var( tr, 'nLepFRLoose', int)
       
        var( tr, 'nJet', int)
        var( tr, 'nJetFwd', int)
        var( tr, 'nBJetLoose', int)
        var( tr, 'nBJetMedium', int)
        var( tr, 'nBJetTight', int)

        for T in self.cfg_ana.triggers:
            var( tr, "Trig_Event_%s" % T.replace("_eta2p1",""), int )

        bookLepton(tr,"Probe", isMC)
        bookJet(tr,"Jet", isMC)
        var(tr, 'dphi')
        var(tr, 'dphi_tp')
        var(tr, 'dr_tp')
        var(tr, 'mtw_probe')
        var(tr, 'mtw2_probe')
        var(tr, 'met')
        var(tr, 'metPhi')

    def beginLoop(self):
        super(ttHLepFRAnalyzer,self).beginLoop()
        self.counters.addCounter('pairs')
        count = self.counters.counter('pairs')
        count.register('all events')
        count.register('one lepton')

    def process(self, iEvent, event):
        self.readCollections( iEvent )
        event.met = self.handles['met'].product()[0]
        event.cleanFRNoIdJetsAll.sort(key = lambda l : (l.pt()* l.rawFactor()), reverse = True)

        self.tree.reset() 
        tr = self.tree
        fill( tr, 'run', event.run) 
        fill( tr, 'lumi',event.lumi)
        fill( tr, 'evt', event.eventId)    
        fill( tr, 'nVert', len(event.goodVertices) )
        fill(tr, 'nLepGood', len(event.selectedLeptons))
        fill(tr, 'nLepFRLoose', len(event.looseFRLeptons))

        fill( tr, 'nJet',        len(event.cleanFRNoIdJetsAll) )
        fill(tr, 'met', event.met.pt())

        fill(tr, 'metPhi', event.met.phi())

        self.counters.counter('pairs').inc('all events')

        triggerResults = self.handles['TriggerResults'].product() 
        for T,C in self.triggerCheckers:
            fill( tr, "Trig_Event_%s" % T, C.check(iEvent.object(), triggerResults))

        def mtw(x1,x2):
            return sqrt(2*x1.pt()*x2.pt()*(1-cos(x1.phi()-x2.phi())))
        def mtw2(lep,met):
            pmu = ROOT.TLorentzVector()
            pmet = ROOT.TLorentzVector()
            pmet.SetPtEtaPhiM(lep.pt(), 0., lep.phi(), lep.mass());
            pmu.SetPtEtaPhiM(met.pt(), 0., met.phi(), 0.);           
            return (pmet+pmu).M()

        for lep in event.looseFRLeptons:
             lep.mvaValue = -99.0
             #self.leptonMVA.addMVA(lep)
             if abs(lep.pdgId()) == 11: continue
             fillLepton(tr, "Probe", lep)
             fill(tr, 'mtw_probe',  mtw(lep, event.met))
             fill(tr, 'mtw2_probe', mtw2(lep, event.met))
             for jet in event.cleanFRNoIdJetsAll:
                 fillJet(tr,"Jet",jet)
                 dphi = deltaPhi(jet.phi(),lep.phi())
                 fill(tr, 'dphi', dphi)
                 break
             tr.tree.Fill()
             break
         
        return True

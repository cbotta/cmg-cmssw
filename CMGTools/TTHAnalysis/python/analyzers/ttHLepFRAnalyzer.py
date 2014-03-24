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
from CMGTools.RootTools.physicsobjects.JetReCalibrator import JetReCalibrator, Type1METCorrection
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
        if self.cfg_comp.isMC:
            self.type1MET = Type1METCorrection("START53_V20","AK5PF",  False)
        else:
            self.type1MET = Type1METCorrection("GR_P_V42_AN4","AK5PF", True) 

    def declareHandles(self):
        super(ttHLepFRAnalyzer, self).declareHandles()
        self.handles['met'] = AutoHandle( 'cmgPFMET', 'std::vector<cmg::BaseMET>' )
        self.handles['metRaw'] = AutoHandle( 'cmgPFMETRaw', 'std::vector<cmg::BaseMET>' )
        self.handles['muons'] = AutoHandle('cmgMuonSel',"std::vector<cmg::Muon>")            
        self.handles['nopumet'] = AutoHandle( 'nopuMet', 'std::vector<reco::PFMET>' )
        self.handles['TriggerResults'] = AutoHandle( ('TriggerResults','','HLT'), 'edm::TriggerResults' )
        self.handles['rho'] = AutoHandle( ('kt6PFJets','rho',''), 'double' )

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
        var(tr, 'met')
        var(tr, 'metPhi')
        var(tr, 'mtw_probe_raw')
        var(tr, 'met_raw')
        var(tr, 'metPhi_raw')
        var(tr, 'mtw_probe_old')
        var(tr, 'met_old')
        var(tr, 'metPhi_old')
        var(tr, 'mtw_probe_t1')
        var(tr, 'met_t1')
        var(tr, 'metPhi_t1')

    def beginLoop(self):
        super(ttHLepFRAnalyzer,self).beginLoop()
        self.counters.addCounter('pairs')
        count = self.counters.counter('pairs')
        count.register('all events')
        count.register('one lepton')

    def process(self, iEvent, event):
        self.readCollections( iEvent )
        event.met = self.handles['met'].product()[0]
        event.metRaw = self.handles['metRaw'].product()[0]
        event.metOld = event.met.__class__(event.met)
        event.metT1  = event.met.__class__(event.metRaw)
        if hasattr(event, 'deltaMetFromJEC'):
            import ROOT
            px,py = event.met.px()+event.deltaMetFromJEC[0], event.met.py()+event.deltaMetFromJEC[1]
            px0, py0 = event.metRaw.px(), event.metRaw.py()
            dpx0, dpy0 = self.type1MET.getMETCorrection( event.allJetsUsedForMET, float(self.handles['rho'].product()[0]),  self.handles['muons'].product())
            px0 += dpx0; py0 += dpy0
            #print "run ",event.run," lumi ", event.lumi," event ", event.eventId, ": MET correction: correct %+7.3f %+7.3f    by hand %+7.3f %+7.3f    diff %+7.3f %+7.3f (phi %+5.3f) " % (event.met.px()-event.metRaw.px(), event.met.py()-event.metRaw.py(), px0-event.metRaw.px(), py0-event.metRaw.py(), event.met.px()-px0,event.met.py()-py0, atan2(event.met.py()-py0, event.met.px()-px0))
            #print "run ",event.run," lumi ", event.lumi," event ", event.eventId, ": old value: ",event.met.pt()," new value: ",hypot(px,py),"  raw met: ",event.metRaw.pt()," type 1 by hand: ",hypot(px0,py0)
            event.met.setP4(  ROOT.reco.Particle.LorentzVector(px,py, 0, hypot(px,py)))
            event.metT1.setP4(ROOT.reco.Particle.LorentzVector(px0,py0, 0, hypot(px0,py0)))
            

        event.cleanFRNoIdJetsAll.sort(key = lambda l : (l.pt()), reverse = True)

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
        fill(tr, 'met_raw', event.metRaw.pt())
        fill(tr, 'metPhi_raw', event.metRaw.phi())
        fill(tr, 'met_t1', event.metT1.pt())
        fill(tr, 'metPhi_t1', event.metT1.phi())
        fill(tr, 'met_old', event.metOld.pt())
        fill(tr, 'metPhi_old', event.metOld.phi())


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
             fill(tr, 'mtw_probe_raw',  mtw(lep, event.metRaw))
             fill(tr, 'mtw_probe_old',  mtw(lep, event.metOld))
             fill(tr, 'mtw_probe_t1',  mtw(lep, event.metT1))
             #fill(tr, 'mtw2_probe', mtw2(lep, event.met))
             for jet in event.cleanFRNoIdJetsAll:
                 fillJet(tr,"Jet",jet)
                 dphi = deltaPhi(jet.phi(),lep.phi())
                 fill(tr, 'dphi', dphi)
                 break
             tr.tree.Fill()
             break
         
        return True

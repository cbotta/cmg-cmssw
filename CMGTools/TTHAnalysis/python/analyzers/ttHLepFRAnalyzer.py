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
        self.tree.reset() 
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

        var( tr, 'tagType', int )

        bookLepton(tr,"TagLepton", isMC)
        bookJet(tr,"Jet", isMC)
        bookJet(tr,"TagJet", isMC)
        bookJet(tr,"SecondJet", isMC)
        var( tr, 'hasSecondB', int)
        #bookParticle(tr,"TagPhoton", isMC)

        for I in (5,8,12,17,24,40):
            var(tr,"Trig_Probe_Mu%d"%I, int)
            var(tr,"Trig_Tag_Mu%d"%I, int)
        var(tr,"Trig_Probe_RelIso1p0Mu5", int)
        var(tr,"Trig_Tag_RelIso1p0Mu5",   int)
        var(tr,"Trig_Probe_1ElT", int)
        var(tr,"Trig_Probe_1ElL", int)
        var(tr,"Trig_Pair_2Mu", int)
        var(tr,"Trig_Pair_MuEG", int)

        var(tr, 'dphi')
        var(tr, 'dphi_tp')
        var(tr, 'dr_tp')
        var(tr, 'mll')
        var(tr, 'mtw_probe')
        var(tr, 'mtw_tag')
        var(tr, 'mtwNoPU_probe')
        var(tr, 'mtwNoPU_tag')
        var(tr, 'met')
        var(tr, 'metPhi')
        var(tr, 'metNoPU')
        var(tr, 'metPhiNoPU')
        var(tr, 'minMllAFAS')
        var(tr, 'minMllAFOS')
        var(tr, 'minMllSFOS')
        var(tr, 'bestMllSFOS')

    def beginLoop(self):
        super(ttHLepFRAnalyzer,self).beginLoop()
        self.counters.addCounter('pairs')
        count = self.counters.counter('pairs')
        count.register('all events')
        count.register('one lepton')
        count.register('loose jet')
        count.register('high pt jet')
        count.register('bjet tag')
        count.register('mu tag')

    def process(self, iEvent, event):
        self.readCollections( iEvent )
        event.met = self.handles['met'].product()[0]
        event.metNoPU = self.handles['nopumet'].product()[0]
        event.cleanFRNoIdJetsAll.sort(key = lambda l : (l.pt()* l.rawFactor()), reverse = True)

        tr = self.tree
        fill( tr, 'run', event.run) 
        fill( tr, 'lumi',event.lumi)
        fill( tr, 'evt', event.eventId)    
        fill( tr, 'nVert', len(event.goodVertices) )
        fill(tr, 'nLepGood', len(event.selectedLeptons))
        fill(tr, 'nLepFRLoose', len(event.looseFRLeptons))

        fill( tr, 'nJet',        len(event.cleanFRNoIdJetsAll) )
        fill( tr, 'nJetFwd',     len(event.cleanFRNoIdJetsFwd) )
        fill( tr, 'nBJetLoose',  sum([j.btag('combinedSecondaryVertexBJetTags') > 0.274 for j in event.cleanFRNoIdJetsAll]) )
        fill( tr, 'nBJetMedium', sum([j.btag('combinedSecondaryVertexBJetTags') > 0.679 for j in event.cleanFRNoIdJetsAll]) )
        fill( tr, 'nBJetTight',  sum([j.btag('combinedSecondaryVertexBJetTags') > 0.898 for j in event.cleanFRNoIdJetsAll]) )

        fill(tr, 'met', event.met.pt())
        fill(tr, 'metPhi', event.met.phi())
        fill(tr, 'metNoPU', event.metNoPU.pt())
        fill(tr, 'metPhiNoPU', event.metNoPU.phi())

        self.counters.counter('pairs').inc('all events')

        # di-lepton veto: simple, for now
        if len(event.selectedLeptons) > 9999.0: ############################## originally was 1 
            return False
        self.counters.counter('pairs').inc('one lepton')

        triggerResults = self.handles['TriggerResults'].product() 
        for T,C in self.triggerCheckers:
            fill( tr, "Trig_Event_%s" % T, C.check(iEvent.object(), triggerResults))

        def mtw(x1,x2):
            return sqrt(2*x1.pt()*x2.pt()*(1-cos(x1.phi()-x2.phi())))

        def passTrigg(lep,triggers):
            for t in triggers:
                if lep.sourcePtr().triggerObjectMatchesByPath('HLT_%s_v*' % t, 1, 0).size() > 0: return True
            return False
        def passFilter(lep,filters):
            for t in filters:
                if lep.sourcePtr().triggerObjectMatchesByFilter(t).size() > 0: return True
            return False

        def passSingleMu(lep):
            return passTrigg(lep, "Mu12 Mu24_eta2p1 Mu30_eta2p1 Mu40 Mu40_eta2p1 Mu17 Mu8".split())

        def passSingleIsoMu(lep):
            return passTrigg(lep, "IsoMu24 IsoMu24_eta2p1".split())

        def passSingleElL(lep):
            return passTrigg(lep, "Ele17_CaloIdL_CaloIsoVL Ele8_CaloIdL_CaloIsoVL".split())

        def passSingleElT(lep):
            return passTrigg(lep, "Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL".split())
        
        def passDoubleMu(lep1,lep2):
            mu17_1 = passFilter(lep1, 'hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17 hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17'.split())
            mu17_2 = passFilter(lep2, 'hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17 hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17'.split())
            mu8_1 = passFilter(lep1, 'hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8 hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8'.split())
            mu8_2 = passFilter(lep2, 'hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8 hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8'.split())
            mu8tk_1 = passFilter(lep1, ['hltDiMuonGlbFiltered17TrkFiltered8'])
            mu8tk_2 = passFilter(lep2, ['hltDiMuonGlbFiltered17TrkFiltered8'])
            if lep1.sourcePtr().triggerObjectMatchesByPath('HLT_Mu17_Mu8_v*',1,0).size() and (mu17_1+mu17_2 >= 1) and (mu8_1+mu8_2 >= 2): return True
            if lep2.sourcePtr().triggerObjectMatchesByPath('HLT_Mu17_Mu8_v*',1,0).size() and (mu17_1+mu17_2 >= 1) and (mu8_1+mu8_2 >= 2): return True
            if lep1.sourcePtr().triggerObjectMatchesByPath('HLT_Mu17_TkMu8_v*',1,0).size() and ((mu17_1 and mu8tk_2) or (mu17_2 and mu8tk_1)): return True
            if lep2.sourcePtr().triggerObjectMatchesByPath('HLT_Mu17_TkMu8_v*',1,0).size() and ((mu17_1 and mu8tk_2) or (mu17_2 and mu8tk_1)): return True
            return False

        def passMuEG(mu,el):
            if (mu.sourcePtr().triggerObjectMatchesByFilter("hltL1Mu12EG7L3MuFiltered17").size() > 0 and
                el.sourcePtr().triggerObjectMatchesByFilter("hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter").size()):
                    return True
            if (mu.sourcePtr().triggerObjectMatchesByFilter("hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8").size() > 0 and
                el.sourcePtr().triggerObjectMatchesByFilter("hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter").size()):
                    return True
            return False
        
        
################################# originally this line was not here
        for lep in event.looseFRLeptons:
             lep.mvaValue = -99.0
             #self.leptonMVA.addMVA(lep)
             if abs(lep.pdgId()) == 11: continue
             fillLepton(tr, "Probe", lep)
             fill(tr, 'mtw_probe', mtw(lep, event.met))
             for jet in event.cleanFRNoIdJetsAll:
                 fillJet(tr,"Jet",jet)
                 dphi = deltaPhi(jet.phi(),lep.phi())
                 fill(tr, 'dphi', dphi)
                 break
             tr.tree.Fill()
             break
         
         

################################# originally this was not commented out
#
#
#
#         # now fill probes
#         for lep in event.selectedLeptons: 
#             fillLepton(tr, "Probe", lep)
#             fill(tr, 'mtw_probe', mtw(lep, event.met))
#             fill(tr, 'mtwNoPU_probe', mtw(lep, event.metNoPU))
#             fill(tr, 'Trig_Probe_Mu5',  passTrigg(lep,['Mu5']))
#             fill(tr, 'Trig_Probe_Mu8',  passTrigg(lep,['Mu8']))
#             fill(tr, 'Trig_Probe_Mu12', passTrigg(lep,['Mu12']))
#             fill(tr, 'Trig_Probe_Mu17', passTrigg(lep,['Mu17']))
#             fill(tr, 'Trig_Probe_Mu24', passTrigg(lep,['Mu24_eta2p1']))
#             fill(tr, 'Trig_Probe_Mu40', passTrigg(lep,['Mu40_eta2p1']))
#             fill(tr, 'Trig_Probe_RelIso1p0Mu5',  passTrigg(lep,['RelIso1p0Mu5']))
#             fill(tr, 'Trig_Probe_1ElL', passSingleElL(lep))
#             fill(tr, 'Trig_Probe_1ElT', passSingleElT(lep))
#             fill(tr, 'minMllAFAS', min([999.]+[(lep.p4()+l2.p4()).M() for l2 in event.looseLeptons ]))
#             fill(tr, 'minMllAFOS', min([999.]+[(lep.p4()+l2.p4()).M() for l2 in event.looseLeptons if l2.charge() == -lep.charge()]))
#             fill(tr, 'minMllSFOS', min([999.]+[(lep.p4()+l2.p4()).M() for l2 in event.looseLeptons if l2.pdgId()  == -lep.pdgId() ]))
#             bestMll = -1
#             for l2 in event.looseLeptons:
#                 if l2.pdgId()  == -lep.pdgId():
#                    mll = (lep.p4()+l2.p4()).M()
#                    if (bestMll == -1 or abs(mll-91.2) < abs(bestMll-91.2)): 
#                         bestMll = mll
#             fill(tr, 'bestMllSFOS', bestMll)
#             # first search for a jet tag
#             fill(tr,'tagType',1) # jet
#             for jet in event.cleanJets:
#                 dr   = deltaR(jet.eta(),jet.phi(),lep.eta(),lep.phi())
#                 dphi = deltaPhi(jet.phi(),lep.phi())
#                 if (dr < 0.8): continue
#                 self.counters.counter('pairs').inc('loose jet')
#                 if jet.pt() < self.cfg_ana.tagJetSel['minPt']: continue
#                 self.counters.counter('pairs').inc('high pt jet')
#                 if jet.btag('combinedSecondaryVertexBJetTags') < self.cfg_ana.tagJetSel['minBTag']: continue
#                 self.counters.counter('pairs').inc('bjet tag')
#                 fillJet(tr, "TagJet", jet)
#                 fill(tr, 'dr_tp',   dr)
#                 fill(tr, 'dphi_tp', dphi)
#                 extraJet = None
#                 for j2 in event.cleanJets:
#                     if deltaR(j2.eta(),j2.phi(),jet.eta(),jet.phi()) > 0.7:
#                         if extraJet == None or extraJet.btag('combinedSecondaryVertexBJetTags') < j2.btag('combinedSecondaryVertexBJetTags'):
#                             extraJet = j2
#                 if extraJet != None:
#                     fillJet(tr, "SecondJet", extraJet) 
#                 fill( tr, 'hasSecondB', extraJet != None)
#                 tr.tree.Fill()
#                 break
#             # then search for a lepton tag
#             fill(tr,'tagType',13) # muon
#             for tag in event.looseLeptons:
#                 tag.mvaValue = -99.0
#                 if abs(tag.pdgId()) != 13: continue
#                 if tag.pt() < 10: continue
#                 if tag.sourcePtr().userFloat("isPFMuon") < 0.5: continue
#                 if tag.sip3D() < self.cfg_ana.tagLeptonSel['minSip3D']: continue
#                 if tag.relIso() < self.cfg_ana.tagLeptonSel['minRelIso']: continue
#                 dr   = deltaR(tag.eta(),tag.phi(),lep.eta(),lep.phi())
#                 dphi = deltaPhi(tag.phi(),lep.phi())
#                 if (dr < 0.8): continue
#                 fillLepton(tr, "TagLepton", tag)
#                 fill(tr, 'mtw_tag', mtw(tag, event.met))
#                 fill(tr, 'mtwNoPU_tag', mtw(tag, event.metNoPU))
#                 fill(tr, 'dr_tp',   dr)
#                 fill(tr, 'dphi_tp', dphi)
#                 fill(tr, 'mll', (tag.p4()+lep.p4()).M())
#                 fill(tr, 'Trig_Tag_RelIso1p0Mu5',  passTrigg(tag,['RelIso1p0Mu5']))
#                 fill(tr, 'Trig_Tag_Mu5',  passTrigg(tag,['Mu5']))
#                 fill(tr, 'Trig_Tag_Mu8',  passTrigg(tag,['Mu8']))
#                 fill(tr, 'Trig_Tag_Mu12', passTrigg(tag,['Mu12']))
#                 fill(tr, 'Trig_Tag_Mu17', passTrigg(tag,['Mu17']))
#                 fill(tr, 'Trig_Tag_Mu24', passTrigg(tag,['Mu24_eta2p1']))
#                 fill(tr, 'Trig_Tag_Mu40', passTrigg(tag,['Mu40_eta2p1']))
#                 fill(tr, 'Trig_Pair_2Mu',  abs(lep.pdgId()) == 13 and passDoubleMu(tag,lep))
#                 fill(tr, 'Trig_Pair_MuEG', abs(lep.pdgId()) == 11 and passMuEG(tag,lep))
#                 self.counters.counter('pairs').inc('mu tag')
#                 extraJet = None
#                 for j2 in event.cleanJets:
#                     if deltaR(j2.eta(),j2.phi(),tag.eta(),tag.phi()) > 0.7:
#                         if extraJet == None or extraJet.btag('combinedSecondaryVertexBJetTags') < j2.btag('combinedSecondaryVertexBJetTags'):
#                             extraJet = j2
#                 if extraJet != None:
#                     fillJet(tr, "SecondJet", extraJet) 
#                 fill( tr, 'hasSecondB', extraJet != None)
#                 tr.tree.Fill()
#                 break
 
        return True

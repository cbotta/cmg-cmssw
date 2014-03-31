from math import *
from os.path import basename

import sys
sys.argv.append('-b-')
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv.remove('-b-')

from CMGTools.TTHAnalysis.plotter.mcPlots import *

class FakeRateSimple:
    def __init__(self,plotFileName, denDir, numDir, options):
        self._plotFileName = plotFileName
        self._denDir = denDir
        self._numDir = numDir
        self._plots = PlotFile(plotFileName,options)
        self._numFile = ROOT.TFile.Open(self._numDir+"/"+basename(self._plotFileName.replace(".txt",".root")))
        self._denFile = ROOT.TFile.Open(self._denDir+"/"+basename(self._plotFileName.replace(".txt",".root")))
        self._options = options
    def makePlotsBySource(self,mca):
        for p in self._plots.plots():
            asig = mca.listSignals()[0]
            abkg = mca.listBackgrounds()[0]
            data = [self._numFile.Get(p.name + "_data"), self._denFile.Get(p.name + "_data")]
            if not data[0]: continue
            sig  = [self._numFile.Get(p.name + "_"+asig).Clone("snum"), self._denFile.Get(p.name + "_"+asig).Clone("sden")]
            bkg  = [self._numFile.Get(p.name + "_"+abkg).Clone("bnum"), self._denFile.Get(p.name + "_"+abkg).Clone("bden")]
            for i in 0,1:
                sig[i].Reset(); bkg[i].Reset();
            for proc in mca.listSignals():
                if self._numFile.Get(p.name + "_" + proc):
                    sig[0].Add(self._numFile.Get(p.name + "_" + proc))
                    sig[1].Add(self._denFile.Get(p.name + "_" + proc))
            for proc in mca.listBackgrounds():
                if self._numFile.Get(p.name + "_" + proc):
                    bkg[0].Add(self._numFile.Get(p.name + "_" + proc))
                    bkg[1].Add(self._denFile.Get(p.name + "_" + proc))
            mc = [ sig[0].Clone("mcnum"), sig[1].Clone("mcden") ]
            for i in 0,1 : mc[i].Add(bkg[i])
            if "TH1" in data[0].ClassName():
                color = { 'data':1, 'qcd':ROOT.kOrange+10, 'ewk':ROOT.kCyan+2, 'mc':ROOT.kOrange+4 }
                for l,h in ('data',data),('qcd',sig),('ewk',bkg),('mc',mc):
                    text = "    Fake rate vs %s for %s\n" % (p.name,l)
                    text += "%3s   %8s  %8s    %6s +/- %6s\n" % ("bin", "xmin ", "xmax ", "value", "error ")
                    text += "%3s   %8s  %8s    %6s-----%6s\n" % ("---", "------", "------", "------", "------")
                    for b in xrange(1,h[0].GetNbinsX()+1):
                        n,d = h[0].GetBinContent(b), h[1].GetBinContent(b)
                        f = n/float(d) if d > 0 else 0; 
                        if l == "data":
                            df = sqrt(f*(1-f)/d) if d > 0 else 0
                        else:
                            # get average weight of events (at numerator)
                            wavg = (h[0].GetBinError(b)**2) /h[0].GetBinError(b) if h[0].GetBinError(b) else 1
                            df = sqrt(f*(1-f)/(d/wavg)) if wavg > 0 and (d/wavg) > 0  and f > 0 and f <  1 else 0
                        text += "%3d   % 8.3f  % 8.3f    %.4f +/- %.4f\n" % (b, h[0].GetXaxis().GetBinLowEdge(b),h[0].GetXaxis().GetBinUpEdge(b), f,df)
                        h[0].SetBinContent(b, f) 
                        h[0].SetBinError(b, df)
                    c1 = ROOT.TCanvas("FR_"+p.name+"_"+l, p.name, 600, 400)
                    h[0].GetYaxis().SetTitle("Fake rate");
                    h[0].SetLineColor(color[l])
                    h[0].SetMarkerColor(color[l])
                    if l != "ewk": h[0].GetYaxis().SetRangeUser(0.0,0.4 if self._options.maxRatioRange[1] > 1 else self._options.maxRatioRange[1]);
                    else:          h[0].GetYaxis().SetRangeUser(0.8,1.0);
                    h[0].SetLineWidth(2)
                    h[0].Draw("E1")
                    ROOT.gStyle.SetErrorX(0.5)
                    #doTinyCmsPrelim(hasExpo = False, textSize = 0.035)
                    for ext in self._options.printPlots.split(","):
                        if ext == "txt": 
                            dump = open("%s/FR_%s_%s.%s" % (self._denDir, p.name, l, ext), "w")
                            dump.write(text)
                        else:
                            c1.Print("%s/FR_%s_%s.%s" % (self._denDir, p.name, l, ext))
                c1 = ROOT.TCanvas("FR_"+p.name+"_stack", p.name, 600, 400)
                sig[0].Draw("E1")
                mc[0].Draw("E1 SAME")
                bkg[0].Draw("E1 SAME")
                data[0].Draw("E1 SAME")
                for ext in self._options.printPlots.split(","):
                    if ext == "txt": continue
                    c1.Print("%s/FR_%s_%s.%s" % (self._denDir, p.name, "stack", ext))
            elif "TH2" in data[0].ClassName():
                for l,h in ('data',data),('qcd',sig),('ewk',bkg),('mc',mc):
                    text = "    Fake rate vs %s for %s; xvar = %s, yvar = %s\n" % (p.name,l, p.getOption('XTitle','x'), p.getOption('YTitle','x'))
                    text += "%3s %3s   %8s  %8s   %8s  %8s    %6s +/- %6s\n" % ("bx", "by", "xmin ", "xmax ", "ymin ", "ymax ", "value", "error ")
                    text += "%3s %3s   %8s  %8s   %8s  %8s    %6s-----%6s\n" % ("---","---",  "------", "------", "------", "------", "------", "------")
                    for bx in xrange(1,h[0].GetNbinsX()+1):
                      for by in xrange(1,h[0].GetNbinsX()+1):
                        n,d = h[0].GetBinContent(bx,by), h[1].GetBinContent(bx,by)
                        f = n/float(d) if d > 0 else 0; 
                        if l == "data":
                            df = sqrt(f*(1-f)/d) if d > 0 else 0
                        else:
                            # get average weight of events (at numerator)
                            wavg = (h[0].GetBinError(bx,by)**2) /h[0].GetBinError(bx,by) if h[0].GetBinError(bx,by) else 1
                            df = sqrt(f*(1-f)/(d/wavg)) if wavg > 0 and (d/wavg) > 0 and f > 0 and f < 1 else 0
                        text += "%3d %3d   % 8.3f  % 8.3f   % 8.3f  % 8.3f    %.4f +/- %.4f\n" % (bx,by, h[0].GetXaxis().GetBinLowEdge(bx),h[0].GetXaxis().GetBinUpEdge(bx), h[0].GetYaxis().GetBinLowEdge(by),h[0].GetYaxis().GetBinUpEdge(by), f,df)
                        h[0].SetBinContent(bx,by, f) 
                        h[0].SetBinError(bx,by, df)
                    c1 = ROOT.TCanvas("FR_"+p.name+"_"+l, p.name, 900, 800)
                    c1.SetRightMargin(0.20)
                    ROOT.gStyle.SetErrorX(0.5)
                    ROOT.gStyle.SetPaintTextFormat(".3f")
                    ROOT.gStyle.SetTextFont(62)
                    h[0].GetZaxis().SetTitle("Fake rate");
                    if l != "ewk": h[0].GetZaxis().SetRangeUser(0.0,0.4);
                    else:          h[0].GetZaxis().SetRangeUser(0.8,1.0);
                    h[0].Draw("COLZ TEXT90E")
                    h[0].SetMarkerSize(1.5)
                    #doTinyCmsPrelim(hasExpo = False, textSize = 0.035)
                    for ext in self._options.printPlots.split(","):
                        if ext == "txt": 
                            dump = open("%s/FR_%s_%s.%s" % (self._denDir, p.name, l, ext), "w")
                            dump.write(text)
                        else:
                            c1.Print("%s/FR_%s_%s.%s" % (self._denDir, p.name, l, ext))
            else:
                raise RuntimeError, "No idea how to handle a " + data[0].ClassName()

class FakeRateMET1Bin:
    def __init__(self,plotFileName, denDir, numDir, options):
        self._plotFileName = plotFileName
        self._denDir = denDir
        self._numDir = numDir
        self._plots = PlotFile(plotFileName,options)
        self._numFile = ROOT.TFile.Open(self._numDir+"/"+basename(self._plotFileName.replace(".txt",".root")))
        self._denFile = ROOT.TFile.Open(self._denDir+"/"+basename(self._plotFileName.replace(".txt",".root")))
        self._options = options
    def integral(self,h,xmin,xmax):
        n, n2 = 0, 0
        for b in xrange(1,h.GetNbinsX()+1):
            if (h.GetXaxis().GetBinCenter(b) > xmin and h.GetXaxis().GetBinCenter(b) < xmax):
                n  += h.GetBinContent(b)
                n2 += h.GetBinError(b)**2
        return [n, sqrt(n2)]
    def frFromRange(self,h,xmin,xmax):
        n = self.integral(h[0],xmin,xmax)
        d = self.integral(h[1],xmin,xmax)
        wavg = n[1]**2/n[0]
        f = n[0]/float(d[0]) if d[0] > 0 else 0
        df = sqrt(f*(1-f)/(d[0]/wavg)) if wavg > 0 and (d[0]/wavg) > 0  and f > 0 and f <  1 else 0
        return (f,df) 
    def rslp(self,hdata,hewk,low,high):
        data_s = self.integral(hdata[1],low[0],low[1])
        data_l = self.integral(hdata[1],high[0],high[1])
        ewk_s = self.integral(hewk[1],low[0],low[1])
        ewk_l = self.integral(hewk[1],high[0],high[1])
        return (data_l[0]/data_s[0]) / (ewk_l[0]/ewk_s[0])
    def makePlotsBySource(self,mca,pname="met"):
        for p in self._plots.plots():
            if p.name != pname: continue
            asig = mca.listSignals()[0]
            abkg = mca.listBackgrounds()[0]
            data = [self._numFile.Get(p.name + "_data"), self._denFile.Get(p.name + "_data")]
            if not data[0]: continue
            sig  = [self._numFile.Get(p.name + "_"+asig).Clone("snum"), self._denFile.Get(p.name + "_"+asig).Clone("sden")]
            bkg  = [self._numFile.Get(p.name + "_"+abkg).Clone("bnum"), self._denFile.Get(p.name + "_"+abkg).Clone("bden")]
            for i in 0,1:
                sig[i].Reset(); bkg[i].Reset();
            for proc in mca.listSignals():
                if self._numFile.Get(p.name + "_" + proc):
                    sig[0].Add(self._numFile.Get(p.name + "_" + proc))
                    sig[1].Add(self._denFile.Get(p.name + "_" + proc))
            for proc in mca.listBackgrounds():
                if self._numFile.Get(p.name + "_" + proc):
                    bkg[0].Add(self._numFile.Get(p.name + "_" + proc))
                    bkg[1].Add(self._denFile.Get(p.name + "_" + proc))
            mc = [ sig[0].Clone("mcnum"), sig[1].Clone("mcden") ]
            for i in 0,1 : mc[i].Add(bkg[i])
            met_s = (  0., 20.); met_l = ( 45., 80. )
            f_s   = self.frFromRange(data, met_s[0], met_s[1])
            f_l   = self.frFromRange(data, met_l[0], met_l[1])
            r_slp = self.rslp(data,bkg,met_s,met_l)
            f_qcd = (f_s[0] - f_l[0]*r_slp)/(1-r_slp)
            if self._options.globalRebin:
                for i in 0,1:
                    for s in data,sig,bkg,mc:
                        s[i].Rebin(self._options.globalRebin)
            if "TH1" in data[0].ClassName():
                color = { 'data':1, 'qcd':ROOT.kOrange+10, 'ewk':ROOT.kCyan+2, 'mc':ROOT.kOrange+4 }
                for l,h in ('data',data),('qcd',sig),('ewk',bkg),('mc',mc):
                    text = "    Fake rate vs %s for %s\n" % (p.name,l)
                    text += "%3s   %8s  %8s    %6s +/- %6s\n" % ("bin", "xmin ", "xmax ", "value", "error ")
                    text += "%3s   %8s  %8s    %6s-----%6s\n" % ("---", "------", "------", "------", "------")
                    for b in xrange(1,h[0].GetNbinsX()+1):
                        n,d = h[0].GetBinContent(b), h[1].GetBinContent(b)
                        f = n/float(d) if d > 0 else 0; 
                        if l == "data":
                            df = sqrt(f*(1-f)/d) if d > 0 else 0
                        else:
                            # get average weight of events (at numerator)
                            wavg = (h[0].GetBinError(b)**2) /h[0].GetBinContent(b) if h[0].GetBinError(b) else 1
                            df = sqrt(f*(1-f)/(d/wavg)) if wavg > 0 and (d/wavg) > 0  and f > 0 and f <  1 else 0
                        text += "%3d   % 8.3f  % 8.3f    %.4f +/- %.4f\n" % (b, h[0].GetXaxis().GetBinLowEdge(b),h[0].GetXaxis().GetBinUpEdge(b), f,df)
                        h[0].SetBinContent(b, f) 
                        h[0].SetBinError(b, df)
                    c1 = ROOT.TCanvas("FR_"+p.name+"_"+l, p.name, 600, 400)
                    h[0].GetYaxis().SetTitle("Fake rate");
                    h[0].SetLineColor(color[l])
                    h[0].SetMarkerColor(color[l])
                    if l != "ewk": h[0].GetYaxis().SetRangeUser(0.0,0.4 if self._options.maxRatioRange[1] > 1 else self._options.maxRatioRange[1]);
                    else:          h[0].GetYaxis().SetRangeUser(0.8,1.0);
                    h[0].SetLineWidth(2)
                    h[0].Draw("E1")
                    ROOT.gStyle.SetErrorX(0.5)
                    #doTinyCmsPrelim(hasExpo = False, textSize = 0.035)
                    for ext in self._options.printPlots.split(","):
                        if ext == "txt": 
                            dump = open("%s/FR_%s_%s.%s" % (self._denDir, p.name, l, ext), "w")
                            dump.write(text)
                        else:
                            c1.Print("%s/FR_%s_%s.%s" % (self._denDir, p.name, l, ext))
                c1 = ROOT.TCanvas("FR_"+p.name+"_stack", p.name, 600, 400)
                sig[0].Draw("E1")
                mc[0].Draw("E1 SAME")
                bkg[0].Draw("E1 SAME")
                data[0].Draw("E1 SAME")
                lsub = ROOT.TLine(sig[0].GetXaxis().GetXmin(), f_qcd, sig[0].GetXaxis().GetXmax(), f_qcd);
                lsub.SetLineWidth(3) 
                lsub.SetLineColor(ROOT.kGreen+2) 
                lsub.Draw("SAME")
                for ext in self._options.printPlots.split(","):
                    if ext == "txt": continue
                    c1.Print("%s/FR_%s_%s.%s" % (self._denDir, p.name, "stack", ext))
            elif "TH2" in data[0].ClassName():
                for l,h in ('data',data),('qcd',sig),('ewk',bkg),('mc',mc):
                    text = "    Fake rate vs %s for %s; xvar = %s, yvar = %s\n" % (p.name,l, p.getOption('XTitle','x'), p.getOption('YTitle','x'))
                    text += "%3s %3s   %8s  %8s   %8s  %8s    %6s +/- %6s\n" % ("bx", "by", "xmin ", "xmax ", "ymin ", "ymax ", "value", "error ")
                    text += "%3s %3s   %8s  %8s   %8s  %8s    %6s-----%6s\n" % ("---","---",  "------", "------", "------", "------", "------", "------")
                    for bx in xrange(1,h[0].GetNbinsX()+1):
                      for by in xrange(1,h[0].GetNbinsX()+1):
                        n,d = h[0].GetBinContent(bx,by), h[1].GetBinContent(bx,by)
                        f = n/float(d) if d > 0 else 0; 
                        if l == "data":
                            df = sqrt(f*(1-f)/d) if d > 0 else 0
                        else:
                            # get average weight of events (at numerator)
                            wavg = (h[0].GetBinError(bx,by)**2) /h[0].GetBinError(bx,by) if h[0].GetBinError(bx,by) else 1
                            df = sqrt(f*(1-f)/(d/wavg)) if wavg > 0 and (d/wavg) > 0 and f > 0 and f < 1 else 0
                        text += "%3d %3d   % 8.3f  % 8.3f   % 8.3f  % 8.3f    %.4f +/- %.4f\n" % (bx,by, h[0].GetXaxis().GetBinLowEdge(bx),h[0].GetXaxis().GetBinUpEdge(bx), h[0].GetYaxis().GetBinLowEdge(by),h[0].GetYaxis().GetBinUpEdge(by), f,df)
                        h[0].SetBinContent(bx,by, f) 
                        h[0].SetBinError(bx,by, df)
                    c1 = ROOT.TCanvas("FR_"+p.name+"_"+l, p.name, 900, 800)
                    c1.SetRightMargin(0.20)
                    ROOT.gStyle.SetErrorX(0.5)
                    ROOT.gStyle.SetPaintTextFormat(".3f")
                    ROOT.gStyle.SetTextFont(62)
                    h[0].GetZaxis().SetTitle("Fake rate");
                    if l != "ewk": h[0].GetZaxis().SetRangeUser(0.0,0.4);
                    else:          h[0].GetZaxis().SetRangeUser(0.8,1.0);
                    h[0].Draw("COLZ TEXT90E")
                    h[0].SetMarkerSize(1.5)
                    #doTinyCmsPrelim(hasExpo = False, textSize = 0.035)
                    for ext in self._options.printPlots.split(","):
                        if ext == "txt": 
                            dump = open("%s/FR_%s_%s.%s" % (self._denDir, p.name, l, ext), "w")
                            dump.write(text)
                        else:
                            c1.Print("%s/FR_%s_%s.%s" % (self._denDir, p.name, l, ext))
            else:
                raise RuntimeError, "No idea how to handle a " + data[0].ClassName()


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] mc.txt cuts.txt plots.txt")
    addPlotMakerOptions(parser)
    (options, args) = parser.parse_args()
    ROOT.gROOT.ProcessLine(".x tdrstyle.cc")
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    mca  = MCAnalysis(args[0],options)
    if len(args) == 5 and args[4] == "MET1Bin":
        FR = FakeRateMET1Bin(args[1], args[2], args[3], options) 
        FR.makePlotsBySource(mca)
    else:
        FR = FakeRateSimple(args[1], args[2], args[3], options) 
        FR.makePlotsBySource(mca)


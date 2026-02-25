import ROOT
import math

print("Loading in plotHelper.py")

variables = {"obj": {}, "evt": {}}
variables["obj"]["pt"] =  {"nbins": 30, "xmin": 0,     "xmax": 200,    "accessor": ".Perp()",  "label": "p_{T} [GeV]"}
variables["obj"]["eta"] = {"nbins": 20, "xmin": -3,    "xmax": 3,       "accessor": ".Eta()",   "label": "#eta"}
variables["obj"]["phi"] = {"nbins": 20, "xmin": -3.5,  "xmax": 3.5,     "accessor": ".Phi()",   "label": "#phi"}
variables["obj"]["E"] =   {"nbins": 30, "xmin": 0,     "xmax": 200,    "accessor": ".E()",      "label": "Energy [GeV]"}
variables["evt"]["n"] =   {"nbins": 20, "xmin": 0,     "xmax": 20,      "accessor": "",          "label": "number per event"}

colors =         [ROOT.TColor.GetColor("#FFA900"),   # Sunny Orange
                ROOT.TColor.GetColor("#3629AC"),    # Dark Blue
                ROOT.TColor.GetColor("#2FC494"),    # Seafoam
                ROOT.TColor.GetColor("#F65866"),    # Pink
                ROOT.TColor.GetColor("#0E81C4"),    # Light Blue
                ]

def getTLV(obj):
    obj_p = obj.getMomentum()
    obj_e = obj.getEnergy()
    obj_tlv = ROOT.TLorentzVector()
    obj_tlv.SetPxPyPzE(obj_p[0], obj_p[1], obj_p[2], obj_e)
    return obj_tlv

def getPt(trk, b_field = 5):
    return 3e-4*abs(b_field/trk.getOmega())

def getP(trk, b_field = 5):
    return getPt(trk)*math.sqrt(1+trk.getTanLambda()**2)

def getTrackTLV(trk, m = .000511, b_field = 5):
    pt = getPt(trk, b_field)
    p = getP(trk, b_field)
    px = pt*math.cos(trk.getPhi())
    py = pt*math.sin(trk.getPhi())
    pz = pt*trk.getTanLambda()
    E = math.sqrt(p**2 + m**2) 
    trk_tlv = ROOT.TLorentzVector()
    trk_tlv.SetPxPyPzE(px, py, pz, E)
    return trk_tlv

def getClusterTLV(cluster):
    pos = cluster.getPosition()
    E   = cluster.getEnergy()
    r = math.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2)
    if r == 0: return ROOT.TLorentzVector()
    px = E * pos[0] / r
    py = E * pos[1] / r
    pz = E * pos[2] / r
    tlv = ROOT.TLorentzVector()
    tlv.SetPxPyPzE(px, py, pz, E)
    return tlv

def fillObjHists(hists, objtype, obj):
    for var in variables["obj"]:
        hists[objtype+"_"+var].Fill(eval("obj"+variables["obj"][var]["accessor"]))

def getMaximum(hists, isLog=False):
    max_factor = 1.1
    if isLog: max_factor = 10
    maximum = max(hists[i].GetMaximum() for i in range(len(hists)))
    return maximum*max_factor

def colorHists(hists):
    i=0
    for h in hists:
        h.SetLineColor(colors[i])
        h.SetMarkerColor(colors[i])
        h.SetMarkerStyle(20) # Added to make error bars look like points
        i += 1
        i %= len(colors)
    return

def plotHistograms(h_map, save_name, xlabel="", ylabel="", interactive=False, logy=False, atltext="", with_bib=False):
    can = ROOT.TCanvas("can", "can", 800, 600)
    can.SetLeftMargin(0.12)
    can.SetBottomMargin(0.12)
    can.SetTopMargin(0.12)
    can.SetRightMargin(0.05)
    
    h_keys = list(h_map.keys())
    h_values = list(h_map.values())

    if len(h_keys)<1:
        print("No histograms found in h_map. Drawing blank canvas.")
        return

    maxy = 1.5*getMaximum(h_values)
    miny = 0
    if logy:
        maxy*=1e4
        miny = 1e-1
        can.SetLogy(1)

    colorHists(h_values)
    h_map[h_keys[0]].SetTitle("")
    h_map[h_keys[0]].GetXaxis().SetTitle(xlabel)
    h_map[h_keys[0]].GetYaxis().SetTitle(ylabel)
    h_map[h_keys[0]].GetXaxis().SetTitleSize(0.045)
    h_map[h_keys[0]].GetYaxis().SetTitleSize(0.045)
    h_map[h_keys[0]].GetXaxis().SetLabelSize(0.04)
    h_map[h_keys[0]].GetYaxis().SetLabelSize(0.04)
    h_values[0].SetMinimum(miny)
    h_values[0].SetMaximum(maxy)
    
    # Draw with "E" for error bars and "hist" for the line
    h_map[h_keys[0]].Draw("E hist")
    
    for k in h_keys[1:]:
        h_map[k].Draw("E hist same")

    # Legend implementation
    leg = ROOT.TLegend(0.65, 0.65, 0.93, 0.85)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    for k in h_keys:
        leg.AddEntry(h_map[k], k, "lp") # lp means Line and Point (for error bars)
    leg.Draw()

    if atltext != "":
        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextSize(0.04)
        if isinstance(atltext, list):
            text.SetTextFont(62)
            text.DrawLatex(0.15, 0.88, atltext[0])
            text.SetTextFont(42)
            text.DrawLatex(0.15, 0.83, atltext[1])
            if len(atltext) > 2:
                text.DrawLatex(0.15, 0.78, atltext[2])
            if len(atltext) > 3:
                text.DrawLatex(0.67, 0.88, atltext[-1])
        else:
            text.DrawLatex(0.15, 0.85, atltext)

    if interactive: 
        input("Press Enter to continue...")

    ROOT.gStyle.SetOptStat(0)
    can.SaveAs(save_name)
    can.Close()
    return

def plotEfficiencies(eff_map, save_name, xlabel="", ylabel="", xrange="", with_bib=False):
    can = ROOT.TCanvas("can", "can", 800, 600)
    can.SetLeftMargin(0.12)
    can.SetBottomMargin(0.12)
    can.SetTopMargin(0.08)
    can.SetRightMargin(0.05)
    
    if len(eff_map) < 1:
        return

    for i, k in enumerate(eff_map):
        eff_map[k].SetMarkerColor(colors[i % len(colors)])
        eff_map[k].SetLineColor(colors[i % len(colors)])
        eff_map[k].SetMarkerStyle(20 + i)
        eff_map[k].SetMarkerSize(1.2)
        eff_map[k].SetLineWidth(2)

    for i, k in enumerate(eff_map):
        if i == 0:
            eff_map[k].Draw("ape")
            eff_map[k].SetTitle("")
            eff_map[k].GetXaxis().SetTitle(xlabel)
            eff_map[k].GetYaxis().SetTitle(ylabel)
            eff_map[k].GetXaxis().SetTitleSize(0.045)
            eff_map[k].GetYaxis().SetTitleSize(0.045)
            eff_map[k].GetXaxis().SetLabelSize(0.04)
            eff_map[k].GetYaxis().SetLabelSize(0.04)
            eff_map[k].SetMinimum(0)
            eff_map[k].SetMaximum(1.4)
            eff_map[k].GetXaxis().SetTickLength(0)
            eff_map[k].GetYaxis().SetTickLength(0)
        else:
            eff_map[k].Draw("pe same")
            
    line = ROOT.TLine()
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kGray+1)
    line.SetLineWidth(1)
    xmin = eff_map[list(eff_map.keys())[0]].GetXaxis().GetXmin()
    xmax = eff_map[list(eff_map.keys())[0]].GetXaxis().GetXmax()
    line.DrawLine(xmin, 1.0, xmax, 1.0)
    
    leg = ROOT.TLegend(0.65, 0.72, 0.93, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    for k in eff_map:
        clean_label = k.replace("pT slice ", "")
        leg.AddEntry(eff_map[k], clean_label, "lp")
    leg.Draw()

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextAlign(11)
    text.SetTextSize(0.04)
    text.SetTextFont(62)
    concept_text = ROOT.TLatex()
    concept_text.SetNDC()
    concept_text.SetTextAlign(31)
    concept_text.SetTextSize(0.035)
    concept_text.SetTextFont(42)
    concept_text.DrawLatex(0.94, 0.92, "MAIA Detector Concept")
    text.DrawLatex(0.15, 0.92, "Muon Collider")
    bib_text = "Simulation, with BIB" if with_bib else "Simulation, no BIB"
    text.SetTextFont(42)
    text.DrawLatex(0.15, 0.88, bib_text)
    text.DrawLatex(0.15, 0.84, "|#eta| < 2.4")
    ROOT.gStyle.SetOptStat(0)
    can.SaveAs(save_name)
    can.Close()
    return

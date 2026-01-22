import ROOT
import math

print("Loading in plotHelper.py")

variables = {"obj": {}, "evt": {}}
variables["obj"]["pt"] =  {"nbins": 30, "xmin": 0,     "xmax": 3000,    "accessor": ".Perp()",  "label": "p_{T} [GeV]"}
variables["obj"]["eta"] = {"nbins": 20, "xmin": -3,    "xmax": 3,       "accessor": ".Eta()",   "label": "#eta"}
variables["obj"]["phi"] = {"nbins": 20, "xmin": -3.5,  "xmax": 3.5,     "accessor": ".Phi()",   "label": "#phi"}
variables["obj"]["E"] =   {"nbins": 30, "xmin": 0,     "xmax": 5000,    "accessor": ".E()",     "label": "Energy [GeV]"}
variables["evt"]["n"] =   {"nbins": 20, "xmin": 0,     "xmax": 20,      "accessor": "",         "label": "number per event"}

colors =        [ROOT.TColor.GetColor("#FFA900"),   # Sunny Orange
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

# Get pT from track object
# Taken from here
# https://bib-pubdb1.desy.de/record/81214/files/LC-DET-2006-004%5B1%5D.pdf
def getPt(trk, b_field = 5):
    return 3e-4*abs(b_field/trk.getOmega())

def getP(trk, b_field = 5):
    return getPt(trk)*math.sqrt(1+trk.getTanLambda()**2)

def getTrackTLV(trk, m = .106, b_field = 5):
    pt = getPt(trk, b_field)
    p = getP(trk, b_field)
    px = pt*math.cos(trk.getPhi())
    py = pt*math.sin(trk.getPhi())
    pz = pt*trk.getTanLambda()
    E = math.sqrt(p**2 + m**2) # Assume this is a muon
    trk_tlv = ROOT.TLorentzVector()
    trk_tlv.SetPxPyPzE(px, py, pz, E)
    return trk_tlv

# for clusters
def getClusterTLV(cluster):
    p = cluster.getPosition()
    E = cluster.getEnergy()
    cluster_tlv = ROOT.TLorentzVector()
    cluster_tlv.SetPxPyPzE(p[0], p[1], p[2], E)
    return cluster_tlv

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
        i += 1
        i %= len(colors)
    return

def plotHistograms(h_map, save_name, xlabel="", ylabel="", interactive=False, logy=False, atltext="", with_bib=False):
    can = ROOT.TCanvas("can", "can", 800, 600)
    #can.SetGrid(1, 1)  # Add grid
    can.SetLeftMargin(0.12)
    can.SetBottomMargin(0.12)
    can.SetTopMargin(0.12)  # Increase top margin for header text
    can.SetRightMargin(0.05)
    
    h_keys = list(h_map.keys())
    h_values = list(h_map.values())

    if len(h_keys)<1:
        print("No histograms found in h_map. Drawing blank canvas.")
        return

    # Get maxes/mins of all hists
    maxy = 1.5*getMaximum(h_values)
    miny = 0
    if logy:
        maxy*=1e4
        miny = 1e-1
        can.SetLogy(1)

    # Draw histograms
    colorHists(h_values)
    h_map[h_keys[0]].SetTitle("")  # Remove title
    h_map[h_keys[0]].GetXaxis().SetTitle(xlabel)
    h_map[h_keys[0]].GetYaxis().SetTitle(ylabel)
    h_map[h_keys[0]].GetXaxis().SetTitleSize(0.045)
    h_map[h_keys[0]].GetYaxis().SetTitleSize(0.045)
    h_map[h_keys[0]].GetXaxis().SetLabelSize(0.04)
    h_map[h_keys[0]].GetYaxis().SetLabelSize(0.04)
    h_values[0].SetMinimum(miny)
    h_values[0].SetMaximum(maxy)
    h_map[h_keys[0]].Draw("hist")
    
    for k in h_keys[1:]:
        h_map[k].Draw("hist same")

    # Create better legend
    #leg = ROOT.TLegend(0.65, 0.15, 0.93, 0.45)
    #leg.SetBorderSize(1)
    #leg.SetFillColor(0)
    #leg.SetFillStyle(1001)
    #leg.SetTextSize(0.035)
    #leg.SetHeader("pT slices", "C")
    
    #for k in h_keys:
        # Clean up legend labels
        #clean_label = k.replace("electronGun_pT_", "").replace("_", "-") + " GeV"
        #leg.AddEntry(h_map[k], clean_label, "l")
    #leg.Draw()
    if atltext != "":
        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextSize(0.04)
        if isinstance(atltext, list):
            # First line (Muon Collider) - top left, bold
            text.SetTextFont(62)  # Bold font
            text.DrawLatex(0.15, 0.88, atltext[0])  # "Muon Collider"
            text.SetTextFont(42)  # Regular font
            text.DrawLatex(0.15, 0.83, atltext[1])  # "Simulation, no BIB" - moved down
        
        # Third line (|eta| < 2.4) - under Simulation text
            if len(atltext) > 2:
                text.DrawLatex(0.15, 0.78, atltext[2])  # "|eta| < 2.4" - moved to left side
        
        # Last item (MAIA Detector Concept) - top right
            if len(atltext) > 3:
                text.DrawLatex(0.67, 0.88, atltext[-1])  # Keep this at the top right
        else:
            text.DrawLatex(0.15, 0.85, atltext)
    if interactive: 
        input("Press Enter to continue...")  # Updated for Python 3

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

    # Use your existing color scheme
    for i, k in enumerate(eff_map):
        eff_map[k].SetMarkerColor(colors[i % len(colors)])
        eff_map[k].SetLineColor(colors[i % len(colors)])
        eff_map[k].SetMarkerStyle(20 + i)  # Different markers
        eff_map[k].SetMarkerSize(1.2)
        eff_map[k].SetLineWidth(2)

    for i, k in enumerate(eff_map):
        if i == 0:
            eff_map[k].Draw("ape")
            eff_map[k].SetTitle("")  # Remove title for cleaner look
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
    # Create legend with header like in your desired style
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
    text.SetTextAlign(11)  # Left aligned
    text.SetTextSize(0.04)
    text.SetTextFont(62)  # Bold font like in image
    concept_text = ROOT.TLatex()
    concept_text.SetNDC()
    concept_text.SetTextAlign(31)  # Right aligned
    concept_text.SetTextSize(0.035)
    concept_text.SetTextFont(42)
    concept_text.DrawLatex(0.94, 0.92, "MAIA Detector Concept")
    text.DrawLatex(0.15, 0.92, "Muon Collider")
    bib_text = "Simulation, with BIB" if with_bib else "Simulation, no BIB"
    text.SetTextFont(42)  # Regular font
    text.DrawLatex(0.15, 0.88, bib_text)
    
    
    # Add eta cut information (hardcoded since cut at |eta| < 2)
    text.DrawLatex(0.15, 0.84, "|#eta| < 2.4")
    ROOT.gStyle.SetOptStat(0)  # Remove stats box
    can.SaveAs(save_name)
    can.Close()
    return

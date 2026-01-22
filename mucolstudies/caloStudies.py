import pyLCIO
import glob
import ctypes
import math
import ROOT

exec(open("./plotHelper.py").read())

# ############## SETUP #############################

ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1
obj_type = "el" 
magnetic_field = 5.00
#calibration_mip = 0.0001575
#calibration_mip_to_reco = 0.00641222630095
#sampling_scaling = calibration_mip_to_reco/calibration_mip
calibration_mip = 1
calibration_mip_to_reco = 1
sampling_scaling = 1
append = "fragment_study" # Updated append name for both Neu and Pho

# Set up things for each object
settings = {
        "fnames": {
                    "el": "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_4.slcio"}, 
        "labelname": {  "ph": "Photon",
                        "mu": "Muon",
                        "el": "Electron"},
        "plotdir":{ "ph": "photons",
                    "mu": "muons",
                    "el": "electrons"},
        "pdgid":  { "ph": 22,
                    "mu": 13,
                    "el": 11},
        "mass":   { "ph": 0,
                    "mu": 0.106,
                    "el": 0.000511}
}
print("Running on", settings["labelname"][obj_type])

# Gather input files
samples = glob.glob(settings["fnames"][obj_type])
fnames = []
for s in samples:
    if s.endswith(".slcio"):
        fnames.append(s)
    else:
        fnames += glob.glob(f"{s}/*.slcio")
print("Found %i files."%len(fnames))

# Define good particle
def isGood(tlv):
    if abs(tlv.Eta()) < 2:
        return True
    return False

# Perform matching between two TLVs
def isMatched(tlv1, tlv2, req_pt = True):
    if tlv1.DeltaR(tlv2) > 0.1: return False
    if req_pt:
        drelpt = abs(tlv1.Perp()-tlv2.Perp())/tlv2.Perp()
        if drelpt > 0.1*tlv2.Perp()/100: return False 
    return True

def getClusterEta(cluster):
    theta = cluster.getITheta()
    return -1*math.ln(math.tan(theta/2))

# ############## CREATE EMPTY HISTOGRAM OBJECTS #############################
variables = {}
variables["E"] = {"nbins": 50, "xmin": 0, "xmax": 50, "title": "E [GeV]"}
hists = {}

objects = {}
objects["mcp"] = f"True {settings['labelname'][obj_type]}"
objects["sim"] = "Sim Calorimeter"
objects["dig"] = "Digi Calorimeter"
objects["rec"] = "Reco Calorimeter"
objects["clu"] = "Matched Cluster"
objects["pfo"] = f"Reconstructed {settings['labelname'][obj_type]}"
objects["neu"] = "Neutron Fragments" 
objects["pho"] = "Photon Fragments" # Added for photon study

for obj in objects:
    for var in variables:
        hists[obj+"_"+var] = ROOT.TH1F(obj+"_"+var, objects[obj], variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"])

ranges = ["_0to1p1", "_1p1to1p2", "_1p2to2"]

hists2d = {}
for obj in objects:
    for var in variables:
        if obj == "mcp": continue
        for r in ranges:
            hists2d[obj+"_v_mcp_"+var+r] = ROOT.TH2F(obj+"_v_mcp_"+var+r, obj+"_v_mcp_"+var+r, variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"], variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"])

# ############## LOOP OVER EVENTS AND FILL HISTOGRAMS #############################
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "ECalBarrelCollection", "ECalEndcapCollection", "EcalBarrelCollectionDigi", "EcalEndcapCollectionDigi", "EcalBarrelCollectionRec", "EcalEndcapCollectionRec", "PandoraClusters"])
i = 0
for f in fnames:
    reader.open(f)
    if max_events > 0 and i >= max_events: break

    for event in reader:
        if max_events > 0 and i >= max_events: break
        if i%100 == 0: print("Processing event %i."%i)

        mcpCollection = event.getCollection("MCParticle")
        pfoCollection = event.getCollection("PandoraPFOs")
        simCollection_b = event.getCollection("ECalBarrelCollection")
        simCollection_e = event.getCollection("ECalEndcapCollection")
        try: digCollection_b = event.getCollection("EcalBarrelCollectionDigi")
        except: digCollection_b = None
        try: digCollection_e = event.getCollection("EcalEndcapCollectionDigi")
        except: digCollection_e = None
        try: recCollection_b = event.getCollection("EcalBarrelCollectionRec")
        except: recCollection_b = None
        try: recCollection_e = event.getCollection("EcalEndcapCollectionRec")
        except: recCollection_e = None
        cluCollection = event.getCollection("PandoraClusters")

        n_mcp_ob = 0
        n_pfo_ob = 0
        n_clu_ob = 0
        has_mcp_ob = False
        has_pfo_ob = False
        has_clu_ob = False
        my_pfo_ob = None
        my_clu_ob = None
        my_mcp_ob = None

        neu_E = 0
        pho_E = 0 # Added photon energy counter

        for mcp in mcpCollection:
            mcp_tlv = getTLV(mcp)
            if abs(mcp.getPDG())==settings['pdgid'][obj_type] and mcp.getGeneratorStatus()==1 and isGood(mcp_tlv):
                has_mcp_ob = True
                n_mcp_ob += 1
                my_mcp_ob = mcp_tlv

        for clu in cluCollection:
            if has_mcp_ob:
                n_clu_ob += 1
                has_clu_ob = True
                if n_clu_ob == 1:
                    my_clu_ob = clu
                elif n_clu_ob > 1 and clu.getEnergy() > my_clu_ob.getEnergy():
                    my_clu_ob = clu

        for pfo in pfoCollection:
            pfo_tlv = getTLV(pfo)

            # Match the main Electron PFO
            if abs(pfo.getType())==settings['pdgid'][obj_type]:
                if has_mcp_ob:
                    n_pfo_ob += 1
                    has_pfo_ob = True
                    if n_pfo_ob == 1:
                        my_pfo_ob = pfo_tlv
                    elif n_pfo_ob > 1 and pfo_tlv.E() > my_pfo_ob.E():
                        my_pfo_ob = pfo_tlv
            
            # Logic for Neutron Fragments
            elif abs(pfo.getType()) == 2112 and has_mcp_ob:
                if my_mcp_ob.DeltaR(pfo_tlv) < 0.2:
                    pfo_clusters = pfo.getClusters()
                    if pfo_clusters:
                        for cluster in pfo_clusters:
                            if cluster:
                                for hit in cluster.getCalorimeterHits():
                                    if hit: neu_E += hit.getEnergy()

            # Logic for Photon Fragments (Misidentified bits of the electron shower)
            elif abs(pfo.getType()) == 22 and has_mcp_ob:
                if my_mcp_ob.DeltaR(pfo_tlv) < 0.2:
                    pfo_clusters = pfo.getClusters()
                    if pfo_clusters:
                        for cluster in pfo_clusters:
                            if cluster:
                                for hit in cluster.getCalorimeterHits():
                                    if hit: pho_E += hit.getEnergy()

        sim_E = 0
        if simCollection_b:
            for sim in simCollection_b: sim_E += sim.getEnergy()*sampling_scaling
        if simCollection_e:
            for sim in simCollection_e: sim_E += sim.getEnergy()*sampling_scaling

        dig_E = 0
        if digCollection_b:
            for dig in digCollection_b: dig_E += dig.getEnergy()*calibration_mip_to_reco
        if digCollection_e:
            for dig in digCollection_e: dig_E += dig.getEnergy()*calibration_mip_to_reco

        rec_E = 0
        if recCollection_b:
            for rec in recCollection_b: rec_E += rec.getEnergy()
        if recCollection_e:
            for rec in recCollection_e: rec_E += rec.getEnergy()

        if n_pfo_ob > 1: print("Found multiple PFOs:", n_pfo_ob)

        if has_mcp_ob:
            hists["mcp_E"].Fill(my_mcp_ob.E())
            hists["neu_E"].Fill(neu_E) 
            hists["pho_E"].Fill(pho_E) # Fill Photon 1D
            if has_pfo_ob:
                hists["pfo_E"].Fill(my_pfo_ob.E())
            if has_clu_ob:
                hists["clu_E"].Fill(my_clu_ob.getEnergy())
            hists["sim_E"].Fill(sim_E)
            hists["dig_E"].Fill(dig_E)
            hists["rec_E"].Fill(rec_E)

            for r in ranges:
                r1 = r.replace("p", ".").strip("_")
                low_eta = r1.split("to")[0]
                high_eta = r1.split("to")[1]
                selection_string = f"my_mcp_ob.Eta()>={low_eta} and my_mcp_ob.Eta()<{high_eta}"
                if eval(selection_string):
                    hists2d["neu_v_mcp_E"+r].Fill(my_mcp_ob.E(), neu_E) 
                    hists2d["pho_v_mcp_E"+r].Fill(my_mcp_ob.E(), pho_E) # Fill Photon 2D
                    if has_pfo_ob: hists2d["pfo_v_mcp_E"+r].Fill(my_mcp_ob.E(), my_pfo_ob.E())
                    if has_clu_ob: hists2d["clu_v_mcp_E"+r].Fill(my_mcp_ob.E(), my_clu_ob.getEnergy())
                    hists2d["sim_v_mcp_E"+r].Fill(my_mcp_ob.E(), sim_E)
                    hists2d["dig_v_mcp_E"+r].Fill(my_mcp_ob.E(), dig_E)
                    hists2d["rec_v_mcp_E"+r].Fill(my_mcp_ob.E(), rec_E)

        i+=1
    reader.close()

# ############## SAVE HISTOGRAMS #############################
for var in variables:
    h_to_plot = {}
    for obj in objects:
        h_to_plot[obj] = hists[obj+"_"+var]
    plotHistograms(h_to_plot, f"plots/calo/comp_{var}_{append}.png", variables[var]["title"], "Count")

for hist in hists2d:
    c = ROOT.TCanvas("c_%s"%hist, "c")
    hists2d[hist].Draw("colz")
    var = hist.split("_")[-2]
    obj = hist.split("_")[0]
    hists2d[hist].GetXaxis().SetTitle("True "+settings['labelname'][obj_type]+" "+variables[var]["title"])
    hists2d[hist].GetYaxis().SetTitle(objects[obj]+" "+variables[var]["title"])
    c.SetRightMargin(0.18)
    c.SetLogz()
    c.SaveAs(f"plots/calo/{hist}_{append}.png")

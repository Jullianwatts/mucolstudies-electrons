import math
import glob
import ROOT
import pyLCIO

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1

# File paths - adjust as needed
samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v5/reco/electronGun*")

files = {}
slices = ["0_50", "50_250", "250_1000", "1000_5000"]
for s in slices: 
    files[f"electronGun_pT_{s}"] = []

for s in samples:
    sname = s.split("/")[-1]
    if sname not in files: continue
    files[sname] = glob.glob(f"{s}/*.slcio")

# Helper for display names
label_map = {
    "electronGun_pT_0_50": "0-50 GeV",
    "electronGun_pT_50_250": "50-250 GeV", 
    "electronGun_pT_250_1000": "250-1000 GeV",
    "electronGun_pT_1000_5000": "1000-5000 GeV"
}

# Set up histograms for efficiency calculation
hists = {}
for s in files:
    hists[s] = {}
    # MCP electrons (denominator)
    hists[s]["mcp_el_eta"] = ROOT.TH1F(f"{s}_mcp_el_eta", s, 20, -3, 3)
    hists[s]["mcp_el_E"] = ROOT.TH1F(f"{s}_mcp_el_E", s, 30, 0, 5000)
    hists[s]["mcp_el_pt"] = ROOT.TH1F(f"{s}_mcp_el_pt", s, 30, 0, 3000)
    
    # Successfully reconstructed as electrons (numerator)
    hists[s]["reco_el_match_eta"] = ROOT.TH1F(f"{s}_reco_el_match_eta", s, 20, -3, 3)
    hists[s]["reco_el_match_E"] = ROOT.TH1F(f"{s}_reco_el_match_E", s, 30, 0, 5000)
    hists[s]["reco_el_match_pt"] = ROOT.TH1F(f"{s}_reco_el_match_pt", s, 30, 0, 3000)
    
    # Missed electrons reconstructed as photons
    hists[s]["missed_as_photon_eta"] = ROOT.TH1F(f"{s}_missed_as_photon_eta", s, 20, -3, 3)
    hists[s]["missed_as_photon_E"] = ROOT.TH1F(f"{s}_missed_as_photon_E", s, 30, 0, 5000)
    hists[s]["missed_as_photon_pt"] = ROOT.TH1F(f"{s}_missed_as_photon_pt", s, 30, 0, 3000)
    
    # Missed electrons reconstructed as other
    hists[s]["missed_as_other_eta"] = ROOT.TH1F(f"{s}_missed_as_other_eta", s, 20, -3, 3)
    hists[s]["missed_as_other_E"] = ROOT.TH1F(f"{s}_missed_as_other_E", s, 30, 0, 5000)
    hists[s]["missed_as_other_pt"] = ROOT.TH1F(f"{s}_missed_as_other_pt", s, 30, 0, 3000)
    
    # Completely unrecognized
    hists[s]["missed_unrecognized_eta"] = ROOT.TH1F(f"{s}_missed_unrecognized_eta", s, 20, -3, 3)
    hists[s]["missed_unrecognized_E"] = ROOT.TH1F(f"{s}_missed_unrecognized_E", s, 30, 0, 5000)
    hists[s]["missed_unrecognized_pt"] = ROOT.TH1F(f"{s}_missed_unrecognized_pt", s, 30, 0, 3000)

# Perform matching between two TLVs
def isMatched(tlv1, tlv2):
    if tlv1.DeltaR(tlv2) < 0.1:
        return True
    return False

# Create a reader object
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "PandoraPFOs"])

# Loop over the different samples
for s in files:
    print("Working on sample", s)
    i = 0

    # Loop over the files in a sample
    for f in files[s]:
        if max_events > 0 and i >= max_events: break
        reader.open(f)

        # Loop over events in each file
        for event in reader:
            if max_events > 0 and i >= max_events: break
            if i % 100 == 0: print("\tProcessing event:", i)

            # Get the collections we care about
            try:
                mcps = event.getCollection("MCParticle")
            except:
                mcps = []
                print("No MCP")
                continue
                
            try:
                pfos = event.getCollection("PandoraPFOs")
            except:
                pfos = []
                continue

            # Collect MCP electrons
            mcp_electrons = []
            for mcp in mcps:
                if not mcp.getGeneratorStatus() == 1: continue
                if abs(mcp.getPDG()) != 11: continue  # Only electrons
                
                mcp_tlv = getTLV(mcp)
                if mcp_tlv.E() < 10: continue  # Energy cut
                if abs(mcp_tlv.Eta()) > 2.5: continue  # Eta cut
                
                mcp_electrons.append(mcp_tlv)
                # Fill denominator histograms
                hists[s]["mcp_el_eta"].Fill(mcp_tlv.Eta())
                hists[s]["mcp_el_E"].Fill(mcp_tlv.E())
                hists[s]["mcp_el_pt"].Fill(mcp_tlv.Perp())

            # Only process events that had at least one electron in fiducial region
            if len(mcp_electrons) == 0: continue

            # Collect reconstructed particles by type
            reco_electrons = []
            reco_photons = []
            reco_others = []
            
            for pfo in pfos:
                pfo_tlv = getTLV(pfo)
                if pfo_tlv.E() < 5: continue  # Low energy cut
                
                pfo_type = pfo.getType()
                if abs(pfo_type) == 11:  # Electron
                    reco_electrons.append(pfo_tlv)
                elif abs(pfo_type) == 22:  # Photon
                    reco_photons.append(pfo_tlv)
                else:
                    reco_others.append(pfo_tlv)

            # One-to-one matching for each MCP electron
            for mcp_el in mcp_electrons:
                # Check if matched to any reconstructed electron
                matched_as_electron = any(isMatched(mcp_el, reco_e) for reco_e in reco_electrons)
                
                if matched_as_electron:
                    # Successfully reconstructed as electron
                    hists[s]["reco_el_match_eta"].Fill(mcp_el.Eta())
                    hists[s]["reco_el_match_E"].Fill(mcp_el.E())
                    hists[s]["reco_el_match_pt"].Fill(mcp_el.Perp())
                else:
                    # Not reconstructed as electron - check what happened
                    matched_as_photon = any(isMatched(mcp_el, reco_p) for reco_p in reco_photons)
                    matched_as_other = any(isMatched(mcp_el, reco_o) for reco_o in reco_others)
                    
                    if matched_as_photon:
                        # Reconstructed as photon
                        hists[s]["missed_as_photon_eta"].Fill(mcp_el.Eta())
                        hists[s]["missed_as_photon_E"].Fill(mcp_el.E())
                        hists[s]["missed_as_photon_pt"].Fill(mcp_el.Perp())
                    elif matched_as_other:
                        # Reconstructed as something else
                        hists[s]["missed_as_other_eta"].Fill(mcp_el.Eta())
                        hists[s]["missed_as_other_E"].Fill(mcp_el.E())
                        hists[s]["missed_as_other_pt"].Fill(mcp_el.Perp())
                    else:
                        # Not reconstructed at all
                        hists[s]["missed_unrecognized_eta"].Fill(mcp_el.Eta())
                        hists[s]["missed_unrecognized_E"].Fill(mcp_el.E())
                        hists[s]["missed_unrecognized_pt"].Fill(mcp_el.Perp())

            i += 1

        reader.close()

# ELECTRON EFFICIENCY CALCULATION
print("\n=== CALCULATING ELECTRON EFFICIENCY ===")

# Efficiency vs eta
eff_eta = {}
for s in hists:
    if hists[s]["mcp_el_eta"].GetEntries() == 0:
        continue
    
    num = hists[s]["reco_el_match_eta"]
    denom = hists[s]["mcp_el_eta"]
    eff = ROOT.TGraphAsymmErrors()
    eff.BayesDivide(num, denom, "B")
    eff_eta[label_map[s]] = eff

plotEfficiencies(eff_eta, "plots/electron_efficiency_vs_eta.png",
                xlabel="#eta", ylabel="Electron Reconstruction Efficiency")

# Efficiency vs energy
eff_energy = {}
for s in hists:
    if hists[s]["mcp_el_E"].GetEntries() == 0:
        continue
    
    num = hists[s]["reco_el_match_E"]
    denom = hists[s]["mcp_el_E"]
    eff = ROOT.TGraphAsymmErrors()
    eff.BayesDivide(num, denom, "B")
    eff_energy[label_map[s]] = eff

plotEfficiencies(eff_energy, "plots/electron_efficiency_vs_energy.png",
                xlabel="Energy [GeV]", ylabel="Electron Reconstruction Efficiency")

# Efficiency vs pt
eff_pt = {}
for s in hists:
    if hists[s]["mcp_el_pt"].GetEntries() == 0:
        continue
    
    num = hists[s]["reco_el_match_pt"]
    denom = hists[s]["mcp_el_pt"]
    eff = ROOT.TGraphAsymmErrors()
    eff.BayesDivide(num, denom, "B")
    eff_pt[label_map[s]] = eff

plotEfficiencies(eff_pt, "plots/electron_efficiency_vs_pt.png",
                xlabel="p_{T} [GeV]", ylabel="Electron Reconstruction Efficiency")

# BREAKDOWN OF MISSED ELECTRONS
print("\n=== CALCULATING BREAKDOWN OF MISSED ELECTRONS ===")

# What happens to missed electrons vs eta
missed_breakdown_eta = {}
for s in hists:
    # Total missed electrons
    total_missed = (hists[s]["missed_as_photon_eta"].GetEntries() + 
                   hists[s]["missed_as_other_eta"].GetEntries() + 
                   hists[s]["missed_unrecognized_eta"].GetEntries())
    
    if total_missed == 0:
        continue
    
    # Create efficiency plots for each category
    for category in ["missed_as_photon", "missed_as_other", "missed_unrecognized"]:
        num = hists[s][f"{category}_eta"]
        # Use total MCP electrons as denominator to show fraction of all electrons
        denom = hists[s]["mcp_el_eta"]
        eff = ROOT.TGraphAsymmErrors()
        eff.BayesDivide(num, denom, "B")
        
        category_label = category.replace("missed_as_", "").replace("_", " ").title()
        key = f"{label_map[s]} - {category_label}"
        missed_breakdown_eta[key] = eff

plotEfficiencies(missed_breakdown_eta, "plots/electron_breakdown_vs_eta.png",
                xlabel="#eta", ylabel="Fraction of MC Electrons")

# Print summary
print("\nElectron Reconstruction Efficiency Summary:")
for s in hists:
    matched = hists[s]["reco_el_match_eta"].GetEntries()
    total = hists[s]["mcp_el_eta"].GetEntries()
    efficiency = matched/total if total > 0 else 0
    
    missed_photon = hists[s]["missed_as_photon_eta"].GetEntries()
    missed_other = hists[s]["missed_as_other_eta"].GetEntries() 
    missed_unrecognized = hists[s]["missed_unrecognized_eta"].GetEntries()
    
    print(f"\n{label_map.get(s, s)}:")
    print(f"  Total electrons: {int(total)}")
    print(f"  Reconstructed as electrons: {int(matched)} ({efficiency:.3f})")
    print(f"  Missed electrons:")
    print(f"    → Reconstructed as photons: {int(missed_photon)} ({missed_photon/total:.3f} of total)")
    print(f"    → Reconstructed as other: {int(missed_other)} ({missed_other/total:.3f} of total)")
    print(f"    → Not reconstructed: {int(missed_unrecognized)} ({missed_unrecognized/total:.3f} of total)")

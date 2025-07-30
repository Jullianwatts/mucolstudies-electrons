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
    
    # Missed electrons reconstructed as specific particles
    particle_types = ["photon", "pion_charged", "neutron", "other"]
    for ptype in particle_types:
        hists[s][f"missed_as_{ptype}_eta"] = ROOT.TH1F(f"{s}_missed_as_{ptype}_eta", s, 20, -3, 3)
        hists[s][f"missed_as_{ptype}_E"] = ROOT.TH1F(f"{s}_missed_as_{ptype}_E", s, 30, 0, 5000)
        hists[s][f"missed_as_{ptype}_pt"] = ROOT.TH1F(f"{s}_missed_as_{ptype}_pt", s, 30, 0, 3000)
    
    # Completely unrecognized
    hists[s]["missed_unrecognized_eta"] = ROOT.TH1F(f"{s}_missed_unrecognized_eta", s, 20, -3, 3)
    hists[s]["missed_unrecognized_E"] = ROOT.TH1F(f"{s}_missed_unrecognized_E", s, 30, 0, 5000)
    hists[s]["missed_unrecognized_pt"] = ROOT.TH1F(f"{s}_missed_unrecognized_pt", s, 30, 0, 3000)

# Perform matching between two TLVs
def isMatched(tlv1, tlv2):
    if tlv1.DeltaR(tlv2) < 0.1:
        return True
    return False

# Function to classify reconstructed particles by PDG ID
def classifyParticle(pfo_type):
    """
    Classify PFO particles based on PDG ID
    Returns: particle category string (only electron, photon, pion, neutron, other)
    """
    abs_type = abs(pfo_type)
    
    if abs_type == 11:  # Electron/positron
        return "electron"
    elif abs_type == 22:  # Photon
        return "photon"
    elif abs_type == 211:  # Charged pion
        return "pion_charged"
    elif abs_type == 111:  # Neutral pion
        return "pion_neutral"
    elif abs_type == 2112:  # Neutron
        return "neutron"
    else:
        return "other"

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
                if abs(mcp_tlv.Eta()) > 2.4: continue  # Eta cut
                
                mcp_electrons.append(mcp_tlv)
                # Fill denominator histograms
                hists[s]["mcp_el_eta"].Fill(mcp_tlv.Eta())
                hists[s]["mcp_el_E"].Fill(mcp_tlv.E())
                hists[s]["mcp_el_pt"].Fill(mcp_tlv.Perp())

            # Only process events that had at least one electron in fiducial region
            if len(mcp_electrons) == 0: continue

            # Collect reconstructed particles by type
            reco_particles = {
                "electron": [],
                "photon": [],
                "pion_charged": [],
                "neutron": [],
                "other": []
            }
            
            for pfo in pfos:
                pfo_tlv = getTLV(pfo)
                if pfo_tlv.E() < 5: continue  # Low energy cut
                
                pfo_type = pfo.getType()
                particle_class = classifyParticle(pfo_type)
                reco_particles[particle_class].append(pfo_tlv)

            # One-to-one matching for each MCP electron
            for mcp_el in mcp_electrons:
                # Check if matched to any reconstructed electron
                matched_as_electron = any(isMatched(mcp_el, reco_e) for reco_e in reco_particles["electron"])
                
                if matched_as_electron:
                    # Successfully reconstructed as electron
                    hists[s]["reco_el_match_eta"].Fill(mcp_el.Eta())
                    hists[s]["reco_el_match_E"].Fill(mcp_el.E())
                    hists[s]["reco_el_match_pt"].Fill(mcp_el.Perp())
                else:
                    # Not reconstructed as electron - check what happened
                    matched_category = None
                    
                    # Check each specific particle type in order of priority
                    for category in ["photon", "pion_charged", "neutron", "other"]:
                        if any(isMatched(mcp_el, reco_p) for reco_p in reco_particles[category]):
                            matched_category = category
                            break
                    
                    if matched_category:
                        # Fill the appropriate histogram
                        hists[s][f"missed_as_{matched_category}_eta"].Fill(mcp_el.Eta())
                        hists[s][f"missed_as_{matched_category}_E"].Fill(mcp_el.E())
                        hists[s][f"missed_as_{matched_category}_pt"].Fill(mcp_el.Perp())
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
                xlabel="#eta", ylabel="Electron Reconstruction Efficiency", with_bib=False)

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
                xlabel="Energy [GeV]", ylabel="Electron Reconstruction Efficiency", with_bib=False)

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
                xlabel="p_{T} [GeV]", ylabel="Electron Reconstruction Efficiency", with_bib=False)

# BREAKDOWN OF MISSED ELECTRONS - COMBINED ACROSS ALL pT SLICES
print("\n=== CALCULATING BREAKDOWN OF MISSED ELECTRONS (COMBINED) ===")

# Combine all pT slices for each particle type
categories = ["reco_el_match", "missed_as_photon", "missed_as_pion_charged", 
             "missed_as_neutron", "missed_as_other", "missed_unrecognized"]

# Create combined histograms for each category
combined_hists = {}
for category in categories:
    combined_hists[f"{category}_eta"] = ROOT.TH1F(f"combined_{category}_eta", "Combined", 20, -3, 3)

# Also create combined denominator (total MCP electrons)
combined_hists["mcp_el_eta"] = ROOT.TH1F("combined_mcp_el_eta", "Combined", 20, -3, 3)

# Add all pT slices together
for s in hists:
    for category in categories:
        combined_hists[f"{category}_eta"].Add(hists[s][f"{category}_eta"])
    combined_hists["mcp_el_eta"].Add(hists[s]["mcp_el_eta"])

# Create efficiency plots for each particle type (combined across all pT slices)
missed_breakdown_eta = {}
for category in categories:
    if combined_hists[f"{category}_eta"].GetEntries() == 0:
        continue
        
    num = combined_hists[f"{category}_eta"]
    denom = combined_hists["mcp_el_eta"]
    eff = ROOT.TGraphAsymmErrors()
    eff.BayesDivide(num, denom, "B")
    
    # Clean up category labels
    if category == "reco_el_match":
        category_label = "Electron"
    elif "pion_charged" in category:
        category_label = "Charged Pion"
    elif "pion_neutral" in category:
        category_label = "Neutral Pion"
    else:
        category_label = category.replace("missed_as_", "").replace("missed_", "").replace("_", " ").title()
    
    missed_breakdown_eta[category_label] = eff

plotEfficiencies(missed_breakdown_eta, "plots/electron_breakdown_vs_eta_combined.png",
                xlabel="#eta", ylabel="Fraction of MC Electrons", with_bib=False)

# Print summary
print("\nElectron Reconstruction Efficiency Summary:")
for s in hists:
    matched = hists[s]["reco_el_match_eta"].GetEntries()
    total = hists[s]["mcp_el_eta"].GetEntries()
    efficiency = matched/total if total > 0 else 0
    
    missed_photon = hists[s]["missed_as_photon_eta"].GetEntries()
    missed_pion_charged = hists[s]["missed_as_pion_charged_eta"].GetEntries()
    missed_neutron = hists[s]["missed_as_neutron_eta"].GetEntries()
    missed_other = hists[s]["missed_as_other_eta"].GetEntries()
    missed_unrecognized = hists[s]["missed_unrecognized_eta"].GetEntries()
    
    print(f"\n{label_map.get(s, s)}:")
    print(f"  Total electrons: {int(total)}")
    print(f"  Reconstructed as electrons: {int(matched)} ({efficiency:.3f})")
    print(f"  Missed electrons:")
    print(f"    → Reconstructed as photons: {int(missed_photon)} ({missed_photon/total:.3f} of total)")
    print(f"    → Reconstructed as charged pions: {int(missed_pion_charged)} ({missed_pion_charged/total:.3f} of total)")
    print(f"    → Reconstructed as neutrons: {int(missed_neutron)} ({missed_neutron/total:.3f} of total)")
    print(f"    → Reconstructed as other: {int(missed_other)} ({missed_other/total:.3f} of total)")
    print(f"    → Not reconstructed: {int(missed_unrecognized)} ({missed_unrecognized/total:.3f} of total)")

# Print combined summary
print("\n=== COMBINED SUMMARY (ALL pT SLICES) ===")
total_combined = combined_hists["mcp_el_eta"].GetEntries()
matched_combined = combined_hists["reco_el_match_eta"].GetEntries()
efficiency_combined = matched_combined/total_combined if total_combined > 0 else 0

missed_photon_combined = combined_hists["missed_as_photon_eta"].GetEntries()
missed_pion_charged_combined = combined_hists["missed_as_pion_charged_eta"].GetEntries()
missed_neutron_combined = combined_hists["missed_as_neutron_eta"].GetEntries()
missed_other_combined = combined_hists["missed_as_other_eta"].GetEntries()
missed_unrecognized_combined = combined_hists["missed_unrecognized_eta"].GetEntries()

print(f"Total electrons: {int(total_combined)}")
print(f"Reconstructed as electrons: {int(matched_combined)} ({efficiency_combined:.3f})")
print(f"Missed electrons:")
print(f"  → Reconstructed as photons: {int(missed_photon_combined)} ({missed_photon_combined/total_combined:.3f} of total)")
print(f"  → Reconstructed as charged pions: {int(missed_pion_charged_combined)} ({missed_pion_charged_combined/total_combined:.3f} of total)")
print(f"  → Reconstructed as neutrons: {int(missed_neutron_combined)} ({missed_neutron_combined/total_combined:.3f} of total)")
print(f"  → Reconstructed as other: {int(missed_other_combined)} ({missed_other_combined/total_combined:.3f} of total)")
print(f"  → Not reconstructed: {int(missed_unrecognized_combined)} ({missed_unrecognized_combined/total_combined:.3f} of total)")

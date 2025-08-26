import math
import glob
import ROOT
import pyLCIO

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

max_events = -1

# ---------- INPUT: only these two slices ----------
samples = glob.glob("/data/fmeloni/DataMuC_MAIA_v0/v5/reco/pionGun*")

files = {}
slices = ["0_50", "50_250"]  # ONLY these two slices
for s in slices:
    files[f"pionGun_pT_{s}"] = []

for s in samples:
    sname = s.split("/")[-1]
    if sname not in files: 
        # silently ignore other directories (e.g., 250_1000, 1000_5000)
        continue
    files[sname] = glob.glob(f"{s}/*.slcio")

label_map = {
    "pionGun_pT_0_50": "0-50 GeV",
    "pionGun_pT_50_250": "50-250 GeV",
}

def getSubdetectorEnergyFractions(pfo):
    try:
        clusters = pfo.getClusters()
        if not clusters or len(clusters) == 0:
            return None, None, None, None, {"total_clusters": 0, "ecal_energy": 0, "hcal_energy": 0}

        total_ecal_energy = 0.0
        total_hcal_energy = 0.0
        processed_clusters = 0
        
        for cluster in clusters:
            try:
                subdet_energies = cluster.getSubdetectorEnergies()
                if not subdet_energies or len(subdet_energies) < 2:
                    continue
                
                ecal_energy = subdet_energies[0]
                hcal_energy = subdet_energies[1]
                if ecal_energy is None or hcal_energy is None: continue
                if math.isnan(ecal_energy) or math.isnan(hcal_energy): continue
                if ecal_energy < 0 or hcal_energy < 0: continue
                
                total_ecal_energy += ecal_energy
                total_hcal_energy += hcal_energy
                processed_clusters += 1
            except Exception:
                continue

        total_subdet_energy = total_ecal_energy + total_hcal_energy
        if total_subdet_energy <= 0:
            return None, None, None, None, {
                "total_clusters": len(clusters),
                "ecal_energy": total_ecal_energy,
                "hcal_energy": total_hcal_energy
            }

        ecal_frac = total_ecal_energy / total_subdet_energy
        hcal_frac = total_hcal_energy / total_subdet_energy
        hcal_ecal_ratio = total_hcal_energy / total_ecal_energy if total_ecal_energy > 0 else None

        energy_stats = {
            "total_clusters": len(clusters),
            "processed_clusters": processed_clusters,
            "ecal_energy": total_ecal_energy,
            "hcal_energy": total_hcal_energy,
            "total_subdet_energy": total_subdet_energy,
            "pfo_energy": pfo.getEnergy(),
            "hcal_ecal_ratio": hcal_ecal_ratio
        }
        return ecal_frac, hcal_frac, hcal_ecal_ratio, total_subdet_energy, energy_stats

    except Exception:
        return None, None, None, None, {"total_clusters": 0, "ecal_energy": 0, "hcal_energy": 0}

# ---------- HISTOGRAMS ----------
hists = {}
for s in files:
    hists[s] = {}
    # True (MC) charged pions
    hists[s]["mcp_pi_eta"] = ROOT.TH1F(f"{s}_mcp_pi_eta", s, 20, -3, 3)
    hists[s]["mcp_pi_E"]   = ROOT.TH1F(f"{s}_mcp_pi_E",   s, 30, 0, 5000)
    hists[s]["mcp_pi_pt"]  = ROOT.TH1F(f"{s}_mcp_pi_pt",  s, 30, 0, 3000)

    # Reconstructed matched-as-pion
    hists[s]["reco_pi_match_eta"] = ROOT.TH1F(f"{s}_reco_pi_match_eta", s, 20, -3, 3)
    hists[s]["reco_pi_match_E"]   = ROOT.TH1F(f"{s}_reco_pi_match_E",   s, 30, 0, 5000)
    hists[s]["reco_pi_match_pt"]  = ROOT.TH1F(f"{s}_reco_pi_match_pt",  s, 30, 0, 3000)

    # HCAL/ECAL distributions
    hists[s]["reco_el_hcal_ecal_ratio"]     = ROOT.TH1F(f"{s}_reco_el_hcal_ecal_ratio",     s, 50, 0, 5)
    hists[s]["reco_pion_hcal_ecal_ratio"]   = ROOT.TH1F(f"{s}_reco_pion_hcal_ecal_ratio",   s, 50, 0, 5)
    hists[s]["reco_photon_hcal_ecal_ratio"] = ROOT.TH1F(f"{s}_reco_photon_hcal_ecal_ratio", s, 50, 0, 5)

    # Missed categories (explicitly include electron & neutral pion)
    particle_types = ["electron", "photon", "pion_neutral", "neutron", "other"]
    for ptype in particle_types:
        hists[s][f"missed_as_{ptype}_eta"] = ROOT.TH1F(f"{s}_missed_as_{ptype}_eta", s, 20, -3, 3)
        hists[s][f"missed_as_{ptype}_E"]   = ROOT.TH1F(f"{s}_missed_as_{ptype}_E",   s, 30, 0, 5000)
        hists[s][f"missed_as_{ptype}_pt"]  = ROOT.TH1F(f"{s}_missed_as_{ptype}_pt",  s, 30, 0, 3000)

    # Completely unrecognized
    hists[s]["missed_unrecognized_eta"] = ROOT.TH1F(f"{s}_missed_unrecognized_eta", s, 20, -3, 3)
    hists[s]["missed_unrecognized_E"]   = ROOT.TH1F(f"{s}_missed_unrecognized_E",   s, 30, 0, 5000)
    hists[s]["missed_unrecognized_pt"]  = ROOT.TH1F(f"{s}_missed_unrecognized_pt",  s, 30, 0, 3000)

def isMatched(tlv1, tlv2):
    return tlv1.DeltaR(tlv2) < 0.1

def classifyParticle(pfo_type):
    abs_type = abs(pfo_type)
    if abs_type == 11:   return "electron"
    if abs_type == 22:   return "photon"
    if abs_type == 211:  return "pion_charged"
    if abs_type == 111:  return "pion_neutral"
    if abs_type == 2112: return "neutron"
    return "other"

# ---------- READER ----------
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "PandoraClusters"])

global_subdet_stats = {
    "pfos_processed": 0,
    "pfos_with_clusters": 0,
    "pfos_with_valid_fractions": 0,
    "total_clusters_processed": 0
}

# For per-slice totals (to print a clean summary)
per_slice_counts = {s: {
    "true_pi": 0, "pi_to_pi": 0, "pi_to_e": 0, "pi_to_ph": 0,
    "pi_to_pi0": 0, "pi_to_n": 0, "pi_to_other": 0, "pi_unrec": 0
} for s in files}

# ---------- EVENT LOOP ----------
for s in files:
    print("Working on sample", s)
    i = 0

    for f in files[s]:
        if max_events > 0 and i >= max_events: break
        reader.open(f)

        for event in reader:
            if max_events > 0 and i >= max_events: break
            if i % 100 == 0: print("\tProcessing event:", i)

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

            try:
                clusters = event.getCollection("PandoraClusters")
            except:
                clusters = []
                print("No PandoraClusters")
                continue

            # ---- select true charged pions (PDG = ±211) ----
            mcp_pions = []
            for mcp in mcps:
                if not mcp.getGeneratorStatus() == 1: continue
                if abs(mcp.getPDG()) != 211: continue
                mcp_tlv = getTLV(mcp)
                if mcp_tlv.E() < 10: continue
                if abs(mcp_tlv.Eta()) > 2.4: continue
                mcp_pions.append(mcp_tlv)
                hists[s]["mcp_pi_eta"].Fill(mcp_tlv.Eta())
                hists[s]["mcp_pi_E"].Fill(mcp_tlv.E())
                hists[s]["mcp_pi_pt"].Fill(mcp_tlv.Perp())
            if len(mcp_pions) == 0: 
                continue

            # ---- bucket reconstructed PFOs by class ----
            reco_particles = {
                "electron": [], "photon": [], "pion_charged": [],
                "pion_neutral": [], "neutron": [], "other": []
            }
            reco_pfos = {
                "electron": [], "photon": [], "pion_charged": [],
                "pion_neutral": [], "neutron": [], "other": []
            }
            for pfo in pfos:
                pfo_tlv = getTLV(pfo)
                if pfo_tlv.E() < 5: continue
                pfo_type = pfo.getType()
                particle_class = classifyParticle(pfo_type)
                reco_particles[particle_class].append(pfo_tlv)
                reco_pfos[particle_class].append(pfo)

            # ---- subdetector energy summaries ----
            subdet_stats = {"total_pfos": 0, "pfos_with_clusters": 0, "valid_fractions": 0}
            for particle_type in ["electron", "photon", "pion_charged"]:
                for pfo in reco_pfos[particle_type]:
                    subdet_stats["total_pfos"] += 1
                    global_subdet_stats["pfos_processed"] += 1
                    ecal_frac, hcal_frac, hcal_ecal_ratio, total_energy, energy_stats = getSubdetectorEnergyFractions(pfo)
                    if energy_stats["total_clusters"] > 0:
                        global_subdet_stats["pfos_with_clusters"] += 1
                        global_subdet_stats["total_clusters_processed"] += energy_stats["total_clusters"]
                        subdet_stats["pfos_with_clusters"] += 1
                    if hcal_ecal_ratio is not None:
                        subdet_stats["valid_fractions"] += 1
                        global_subdet_stats["pfos_with_valid_fractions"] += 1
                        if particle_type == "electron":
                            hists[s]["reco_el_hcal_ecal_ratio"].Fill(hcal_ecal_ratio)
                        elif particle_type == "photon":
                            hists[s]["reco_photon_hcal_ecal_ratio"].Fill(hcal_ecal_ratio)
                        elif particle_type == "pion_charged":
                            hists[s]["reco_pion_hcal_ecal_ratio"].Fill(hcal_ecal_ratio)

            if i % 100 == 0 and subdet_stats["total_pfos"] > 0:
                print(f"\t  PFO subdetector analysis: {subdet_stats['pfos_with_clusters']}/{subdet_stats['total_pfos']} PFOs with clusters")
                print(f"\t  Valid energy fractions: {subdet_stats['valid_fractions']}")

            # ---- match each true pion to a reconstructed object ----
            for mcp_pi in mcp_pions:
                per_slice_counts[s]["true_pi"] += 1
                matched_as_pion = any(isMatched(mcp_pi, reco_pi) for reco_pi in reco_particles["pion_charged"])

                if matched_as_pion:
                    hists[s]["reco_pi_match_eta"].Fill(mcp_pi.Eta())
                    hists[s]["reco_pi_match_E"].Fill(mcp_pi.E())
                    hists[s]["reco_pi_match_pt"].Fill(mcp_pi.Perp())
                    per_slice_counts[s]["pi_to_pi"] += 1
                else:
                    matched_category = None
                    for category in ["electron", "photon", "pion_neutral", "neutron", "other"]:
                        if any(isMatched(mcp_pi, reco_p) for reco_p in reco_particles[category]):
                            matched_category = category
                            break

                    if matched_category:
                        hists[s][f"missed_as_{matched_category}_eta"].Fill(mcp_pi.Eta())
                        hists[s][f"missed_as_{matched_category}_E"].Fill(mcp_pi.E())
                        hists[s][f"missed_as_{matched_category}_pt"].Fill(mcp_pi.Perp())
                        if matched_category == "electron":       per_slice_counts[s]["pi_to_e"] += 1
                        elif matched_category == "photon":        per_slice_counts[s]["pi_to_ph"] += 1
                        elif matched_category == "pion_neutral":  per_slice_counts[s]["pi_to_pi0"] += 1
                        elif matched_category == "neutron":       per_slice_counts[s]["pi_to_n"] += 1
                        else:                                     per_slice_counts[s]["pi_to_other"] += 1
                    else:
                        hists[s]["missed_unrecognized_eta"].Fill(mcp_pi.Eta())
                        hists[s]["missed_unrecognized_E"].Fill(mcp_pi.E())
                        hists[s]["missed_unrecognized_pt"].Fill(mcp_pi.Perp())
                        per_slice_counts[s]["pi_unrec"] += 1

            i += 1
        reader.close()

# ---------- PRINT SUBDETECTOR SUMMARY ----------
print("\n=== SUBDETECTOR ENERGY ANALYSIS STATISTICS ===")
print(f"Total PFOs processed: {global_subdet_stats['pfos_processed']}")
print(f"PFOs with clusters: {global_subdet_stats['pfos_with_clusters']}")
print(f"PFOs with valid ECAL/HCAL energy fractions: {global_subdet_stats['pfos_with_valid_fractions']}")
print(f"Total clusters processed: {global_subdet_stats['total_clusters_processed']}")
if global_subdet_stats['pfos_processed'] > 0:
    cluster_availability = global_subdet_stats['pfos_with_clusters'] / global_subdet_stats['pfos_processed']
    print(f"Cluster availability rate: {cluster_availability:.3f}")
    if global_subdet_stats['pfos_with_clusters'] > 0:
        fraction_success = global_subdet_stats['pfos_with_valid_fractions'] / global_subdet_stats['pfos_with_clusters']
        print(f"Energy fraction calculation success rate: {fraction_success:.3f}")

# ---------- EFFICIENCIES (per-slice) ----------
eff_eta = {}
for s in hists:
    if hists[s]["mcp_pi_eta"].GetEntries() == 0: continue
    num = hists[s]["reco_pi_match_eta"]
    denom = hists[s]["mcp_pi_eta"]
    eff = ROOT.TGraphAsymmErrors(); eff.BayesDivide(num, denom, "B")
    eff_eta[label_map[s]] = eff
plotEfficiencies(eff_eta, "plots/pion_efficiency_vs_eta.png",
                 xlabel="#eta", ylabel="Pion Reconstruction Efficiency", with_bib=False)

eff_energy = {}
for s in hists:
    if hists[s]["mcp_pi_E"].GetEntries() == 0: continue
    num = hists[s]["reco_pi_match_E"]; denom = hists[s]["mcp_pi_E"]
    eff = ROOT.TGraphAsymmErrors(); eff.BayesDivide(num, denom, "B")
    eff_energy[label_map[s]] = eff
plotEfficiencies(eff_energy, "plots/pion_efficiency_vs_energy.png",
                 xlabel="Energy [GeV]", ylabel="Pion Reconstruction Efficiency", with_bib=False)

eff_pt = {}
for s in hists:
    if hists[s]["mcp_pi_pt"].GetEntries() == 0: continue
    num = hists[s]["reco_pi_match_pt"]; denom = hists[s]["mcp_pi_pt"]
    eff = ROOT.TGraphAsymmErrors(); eff.BayesDivide(num, denom, "B")
    eff_pt[label_map[s]] = eff
plotEfficiencies(eff_pt, "plots/pion_efficiency_vs_pt.png",
                 xlabel="p_{T} [GeV]", ylabel="Pion Reconstruction Efficiency", with_bib=False)

# ---------- HCAL/ECAL RATIO OVERLAYS ----------
particle_types_ratio = ["electron", "photon", "pion_charged"]
for ptype in particle_types_ratio:
    ratio_hists = {}
    for s in hists:
        hist_name = f"reco_{ptype.replace('_charged', '')}_hcal_ecal_ratio"
        if hist_name in hists[s] and hists[s][hist_name].GetEntries() > 0:
            ratio_hists[label_map[s]] = hists[s][hist_name]
    if len(ratio_hists) > 0:
        plotEfficiencies(ratio_hists, f"plots/{ptype}_hcal_ecal_ratio_vs_pt.png",
                         xlabel="HCAL/ECAL Energy Ratio", ylabel="Normalized Entries", with_bib=False)

# ---------- COMBINED + PER-SLICE SUMMARIES ----------
categories = [
    "reco_pi_match",
    "missed_as_electron", "missed_as_photon",
    "missed_as_pion_neutral",
    "missed_as_neutron", "missed_as_other",
    "missed_unrecognized"
]
combined_hists = {}
for category in categories:
    combined_hists[f"{category}_eta"] = ROOT.TH1F(f"combined_{category}_eta", "Combined", 20, -3, 3)
combined_hists["mcp_pi_eta"] = ROOT.TH1F("combined_mcp_pi_eta", "Combined", 20, -3, 3)

for s in hists:
    for category in categories:
        combined_hists[f"{category}_eta"].Add(hists[s][f"{category}_eta"])
    combined_hists["mcp_pi_eta"].Add(hists[s]["mcp_pi_eta"])

# Per-slice printout
print("\n=== SUMMARY (per slice) ===")
for s in files:
    lab = label_map[s]
    c = per_slice_counts[s]
    eff = c["pi_to_pi"]/c["true_pi"] if c["true_pi"]>0 else 0.0
    print(f"{lab}:")
    print(f"  True charged pions:        {c['true_pi']}")
    print(f"  Reco as charged pions:     {c['pi_to_pi']} (eff = {eff:.3f})")
    print(f"  Missed as electrons:       {c['pi_to_e']}")
    print(f"  Missed as photons:         {c['pi_to_ph']}")
    print(f"  Missed as neutral pions:   {c['pi_to_pi0']}")
    print(f"  Missed as neutrons:        {c['pi_to_n']}")
    print(f"  Missed as other:           {c['pi_to_other']}")
    print(f"  Not reconstructed:         {c['pi_unrec']}")

# Combined breakdown plot
missed_breakdown_eta = {}
for category in categories:
    if combined_hists[f"{category}_eta"].GetEntries() == 0:
        continue
    num = combined_hists[f"{category}_eta"]
    denom = combined_hists["mcp_pi_eta"]
    eff = ROOT.TGraphAsymmErrors()
    eff.BayesDivide(num, denom, "B")
    if category == "reco_pi_match":
        category_label = "Charged Pion"
    elif "pion_neutral" in category:
        category_label = "Neutral Pion"
    else:
        category_label = category.replace("missed_as_", "").replace("missed_", "").replace("_", " ").title()
    missed_breakdown_eta[category_label] = eff

plotEfficiencies(missed_breakdown_eta, "plots/pion_breakdown_vs_eta_combined.png",
                 xlabel="#eta", ylabel="Fraction of MC Charged Pions", with_bib=False)

# Combined counts
total_combined = combined_hists["mcp_pi_eta"].GetEntries()
matched_combined = combined_hists["reco_pi_match_eta"].GetEntries()
efficiency_combined = matched_combined/total_combined if total_combined > 0 else 0
missed_el_combined   = combined_hists["missed_as_electron_eta"].GetEntries()
missed_ph_combined   = combined_hists["missed_as_photon_eta"].GetEntries()
missed_pi0_combined  = combined_hists["missed_as_pion_neutral_eta"].GetEntries()
missed_n_combined    = combined_hists["missed_as_neutron_eta"].GetEntries()
missed_other_combined= combined_hists["missed_as_other_eta"].GetEntries()
missed_unrec_combined= combined_hists["missed_unrecognized_eta"].GetEntries()

print("\n=== SUMMARY (combined) ===")
print(f"Total true charged pions: {int(total_combined)}")
print(f"Reconstructed as charged pions: {int(matched_combined)}  (eff = {efficiency_combined:.3f})")
print(f"Missed as electrons:      {int(missed_el_combined)}")
print(f"Missed as photons:        {int(missed_ph_combined)}")
print(f"Missed as neutral pions:  {int(missed_pi0_combined)}")
print(f"Missed as neutrons:       {int(missed_n_combined)}")
print(f"Missed as other:          {int(missed_other_combined)}")
print(f"Not reconstructed:        {int(missed_unrec_combined)}")


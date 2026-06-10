import math
import os
import pyLCIO
import ROOT
from pyLCIO import EVENT, IOIMPL
from plotHelper import *


file_path     = "/scratch/jwatts/mucol/v2.11/reco/electronGun_pT_0_50/10kelectron0to50_recoCalib2.slcio"
output_dir    = "plots2026/diagnosis"
os.makedirs(output_dir, exist_ok=True)
output_prefix = f"{output_dir}/electronPFOdiag"
DR_MATCH      = 0.2
MAX_EVENTS    = 10000
MC_PT_CUT     = 10.0   # GeV — MC truth electron pT cut
TRK_PT_CUT    = 2.0    # GeV — min pT on tracks attached to PFO candidates
CLUS_E_CUT    = 2.0    # GeV — min energy on clusters attached to PFO candidates


def dR(tlv1, tlv2):
    deta = tlv1.Eta() - tlv2.Eta()
    dphi = tlv1.DeltaPhi(tlv2)
    return math.sqrt(deta**2 + dphi**2)

def getPfoMomentum(pfo):
    # For charged PFOs: use track momentum directly (avoids massless TLV issue)
    # For neutral PFOs: fall back to TLV (no track anyway)
    tracks = pfo.getTracks()
    if tracks.size() > 0:
        return getP(tracks[0])   # getP from plotHelper uses omega + tanLambda
    return getTLV(pfo).P()

def getPfoEnergy(pfo):
    # Sum cluster energies directly — more reliable than TLV energy for E/p
    clusters = pfo.getClusters()
    if clusters.size() > 0:
        return sum(c.getEnergy() for c in clusters)
    return getTLV(pfo).E()

def getEcalFrac(pfo):
    # Use getSubdetectorEnergies(): index 0 = ECAL, index 1 = HCAL
    ecal_e, hcal_e = 0.0, 0.0
    for c in pfo.getClusters():
        subdet = c.getSubdetectorEnergies()
        if subdet.size() > 1:
            ecal_e += subdet[0]
            hcal_e += subdet[1]
        else:
            # fallback: use cluster hit positions — hits closer to IP are ECAL
            for hit in c.getCalorimeterHits():
                pos = hit.getPosition()
                r   = math.sqrt(pos[0]**2 + pos[1]**2)
                if r < 1800:   # ~ECAL outer radius in mm for MAIA
                    ecal_e += hit.getEnergy()
                else:
                    hcal_e += hit.getEnergy()
    total = ecal_e + hcal_e
    return ecal_e / total if total > 0 else 0.0


categories   = ["correct_ID", "wrong_ID", "neutral_misID", "no_match"]
matched_cats = categories[:3]

h_pt       = {c: ROOT.TH1F(f"pt_{c}",    "", 25,   0,   55)  for c in categories}
h_eta      = {c: ROOT.TH1F(f"eta_{c}",   "", 20,  -3,    3)  for c in categories}
h_energy   = {c: ROOT.TH1F(f"energy_{c}","", 30,   0,  400)  for c in categories}
h_eop      = {c: ROOT.TH1F(f"eop_{c}",   "", 60,   0,    3)  for c in matched_cats}
h_ecalfrac = {c: ROOT.TH1F(f"ecalf_{c}", "", 25,   0, 1.05)  for c in matched_cats}
h_ntracks  = {c: ROOT.TH1F(f"ntrk_{c}",  "",  6, -0.5, 5.5) for c in matched_cats}
h_nclus    = {c: ROOT.TH1F(f"nclus_{c}", "",  6, -0.5, 5.5) for c in matched_cats}

h_pfo_type = ROOT.TH1F("pfo_type", "", 30, -0.5, 29.5)
h_e_vs_eta = ROOT.TH2F("e_vs_eta", "", 20,  -3,    3, 30, 0, 200)

match_strats     = ["closest_dR", "highest_pT"]
h_match_pfo_type = {s: ROOT.TH1F(f"pfo_type_{s}", "", 30, -0.5, 29.5) for s in match_strats}
h_match_eop      = {s: ROOT.TH1F(f"eop_{s}",      "", 60,   0,     3) for s in match_strats}
h_match_dr       = {s: ROOT.TH1F(f"dr_{s}",        "", 50,   0,  0.15) for s in match_strats}


reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(file_path)

event_count  = 0
mc_e_count   = 0
match_counts = {c: 0 for c in categories}
pfo_types_seen = set()   # debug: track all PFO types seen in file

for event in reader:

    event_count += 1

    try:
        mc_coll  = event.getCollection("MCParticle")
        pfo_coll = event.getCollection("PandoraPFOs")
    except:
        continue

    # debug: print subdetector energy indices on first event
    if event_count == 1:
        for pfo in pfo_coll:
            for c in pfo.getClusters():
                subdet = c.getSubdetectorEnergies()
                print(f"[DEBUG] cluster energy={c.getEnergy():.2f}  subdetectorEnergies size={subdet.size()}  values={[subdet[i] for i in range(subdet.size())]}")
                break
            break

    def pfo_passes(pfo):
        trk_ok  = any(getPt(t) > TRK_PT_CUT for t in pfo.getTracks())
        clus_ok = any(c.getEnergy() > CLUS_E_CUT for c in pfo.getClusters())
        return trk_ok or clus_ok

    pfo_list = [(pfo, getTLV(pfo)) for pfo in pfo_coll if pfo_passes(pfo)]

    # collect all PFO types for debug summary
    for pfo, _ in pfo_list:
        pfo_types_seen.add(abs(pfo.getType()))

    for mc in mc_coll:

        if abs(mc.getPDG()) != 11 or mc.getGeneratorStatus() != 1:
            continue

        mc_tlv = getTLV(mc)
        mc_pt  = mc_tlv.Perp()
        mc_eta = mc_tlv.Eta()
        mc_e   = mc_tlv.E()

        if mc_pt < MC_PT_CUT:
            continue

        mc_e_count += 1
        h_e_vs_eta.Fill(mc_eta, mc_e)

        candidates = [(pfo, tlv, dR(mc_tlv, tlv)) for pfo, tlv in pfo_list if dR(mc_tlv, tlv) < DR_MATCH]

        if len(candidates) == 0:
            cat = "no_match"
            h_pt[cat].Fill(mc_pt); h_eta[cat].Fill(mc_eta); h_energy[cat].Fill(mc_e)
            match_counts[cat] += 1
            continue

        pfo_closestDR, _, best_dr  = min(candidates, key=lambda x: x[2])
        pfo_highestPT, _, dr_highPT = max(candidates, key=lambda x: x[1].Perp())

        # match strategy comparison
        for strat, pfo, dr_val in [("closest_dR", pfo_closestDR, best_dr),
                                    ("highest_pT", pfo_highestPT, dr_highPT)]:
            h_match_pfo_type[strat].Fill(abs(pfo.getType()))
            h_match_dr[strat].Fill(dr_val)
            p = getPfoMomentum(pfo)
            e = getPfoEnergy(pfo)
            if p > 0:
                h_match_eop[strat].Fill(e / p)

        # primary match: closest dR
        best_pfo  = pfo_closestDR
        pfo_type  = abs(best_pfo.getType())
        h_pfo_type.Fill(pfo_type)

        if   pfo_type == 11:              cat = "correct_ID"
        elif pfo_type in (22, 2112, 130): cat = "neutral_misID"
        else:                             cat = "wrong_ID"

        h_pt[cat].Fill(mc_pt); h_eta[cat].Fill(mc_eta); h_energy[cat].Fill(mc_e)
        match_counts[cat] += 1

        # E/p: track momentum for charged, cluster energy for numerator
        p = getPfoMomentum(best_pfo)
        e = getPfoEnergy(best_pfo)
        if p > 0:
            h_eop[cat].Fill(e / p)

        h_ecalfrac[cat].Fill(getEcalFrac(best_pfo))
        h_ntracks[cat].Fill(best_pfo.getTracks().size())
        h_nclus[cat].Fill(best_pfo.getClusters().size())

    if event_count >= MAX_EVENTS:
        break

reader.close()


print(f"\n({event_count} events, {mc_e_count} MC electrons)")
print(f"PFO types seen in file: {sorted(pfo_types_seen)}")
print(f"\n{'Category':<20} | {'Count':<8} | {'Fraction'}")
for cat in categories:
    n = match_counts[cat]
    print(f"{cat:<20} | {n:<8} | {n/mc_e_count:.3f}")


atext  = ["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4"]
labels = {
    "correct_ID":    "Correct e^{-} ID",
    "wrong_ID":      "Wrong charged ID",
    "neutral_misID": "Neutral misID",
    "no_match":      "No PFO match"
}

def hmap(d, keys): return {labels[k]: d[k] for k in keys}

plotHistograms(hmap(h_pt,  categories),     f"{output_prefix}_pt_byCategory.png",
    xlabel="MC Electron p_{T} [GeV]",  ylabel="Electrons / bin", atltext=atext)

plotHistograms(hmap(h_eta, categories),     f"{output_prefix}_eta_byCategory.png",
    xlabel="MC Electron #eta",         ylabel="Electrons / bin", atltext=atext)

plotHistograms(hmap(h_energy, categories),  f"{output_prefix}_energy_byCategory.png",
    xlabel="MC Electron Energy [GeV]", ylabel="Electrons / bin", logy=True, atltext=atext)

plotHistograms(hmap(h_eop, matched_cats),   f"{output_prefix}_eop.png",
    xlabel="E / p", ylabel="Electrons / bin", atltext=atext)

plotHistograms(hmap(h_ecalfrac, matched_cats), f"{output_prefix}_ecalFrac.png",
    xlabel="ECAL Energy Fraction", ylabel="Electrons / bin", atltext=atext)

plotHistograms(hmap(h_ntracks, matched_cats),  f"{output_prefix}_nTracks.png",
    xlabel="N Tracks on PFO", ylabel="Electrons / bin", atltext=atext)

plotHistograms(hmap(h_nclus, matched_cats),    f"{output_prefix}_nClusters.png",
    xlabel="N Clusters on PFO", ylabel="Electrons / bin", atltext=atext)

plotHistograms({"Matched electrons": h_pfo_type}, f"{output_prefix}_pfoTypeAssigned.png",
    xlabel="PFO PDG type (abs)", ylabel="Count", atltext=atext)

# 2D energy vs eta
can2d = ROOT.TCanvas("can2d", "can2d", 900, 700)
can2d.SetLeftMargin(0.12); can2d.SetRightMargin(0.15); can2d.SetBottomMargin(0.12)
ROOT.gStyle.SetOptStat(0)
h_e_vs_eta.SetTitle("")
h_e_vs_eta.GetXaxis().SetTitle("MC Electron #eta")
h_e_vs_eta.GetYaxis().SetTitle("MC Electron Energy [GeV]")
h_e_vs_eta.Draw("COLZ")
can2d.SaveAs(f"{output_prefix}_energy_vs_eta_2D.png")
can2d.Close()

plotHistograms(
    {"Closest #DeltaR": h_match_pfo_type["closest_dR"], "Highest p_{T} in cone": h_match_pfo_type["highest_pT"]},
    f"{output_prefix}_matchStrat_pfoType.png",
    xlabel="PFO PDG type (abs)", ylabel="Count", atltext=atext)

plotHistograms(
    {"Closest #DeltaR": h_match_eop["closest_dR"], "Highest p_{T} in cone": h_match_eop["highest_pT"]},
    f"{output_prefix}_matchStrat_eop.png",
    xlabel="E / p", ylabel="Electrons / bin", atltext=atext)

plotHistograms(
    {"Closest #DeltaR": h_match_dr["closest_dR"], "Highest p_{T} in cone": h_match_dr["highest_pT"]},
    f"{output_prefix}_matchStrat_dR.png",
    xlabel="#DeltaR(MC e^{-}, matched PFO)", ylabel="Electrons / bin", atltext=atext)

print(f"\nPlots saved with prefix: {output_prefix}_*.png")

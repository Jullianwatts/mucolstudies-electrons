import os, math, ROOT, glob, pyLCIO
from pyLCIO import IOIMPL, EVENT, UTIL

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()
os.makedirs("plots", exist_ok=True)

# Config: Only 0-50 GeV slices, B=5T
B_FIELD, MAX_EVENTS, DR_CUT = 5, -1, 0.2
base_path = "/data/fmeloni/DataMuC_MAIA_v0/v5/reco/"
samples = ["electronGun_pT_0_50", "pionGun_pT_0_50"]

# --- Profile Discrepancy Functions ---
def generate_expected_em_profile(energy, num_layers=50, X0_per_layer=0.6286):
    a = 1.0 + 0.5 * math.log(max(energy, 0.01) / 0.01)
    b = 0.5 
    expected = {}
    total = 0.0
    for l in range(num_layers):
        t = l * X0_per_layer
        try: val = (b**a * t**(a-1) * math.exp(-b*t)) / math.gamma(a) if t > 0 else 0
        except: val = 0
        expected[l] = val
        total += val
    return {l: v/total for l, v in expected.items()} if total > 0 else {l: 0 for l in range(num_layers)}

def get_max_profile_discrepancy(hit_energies_by_layer, total_energy):
    if total_energy <= 0: return 0.0
    expected = generate_expected_em_profile(total_energy)
    max_disc = 0.0
    for l in range(50):
        obs_frac = hit_energies_by_layer.get(l, 0.0) / total_energy
        max_disc = max(max_disc, abs(obs_frac - expected.get(l, 0.0)))
    return max_disc

# --- Initialization ---
hists = {s: {
    "ep": ROOT.TH1F(f"ep_{s}", f"{s};E/p;Entries", 100, 0, 2.0),
    "disc": ROOT.TH1F(f"disc_{s}", f"{s};Max Profile Discrepancy;Entries", 100, 0, 0.5)
} for s in samples}

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

for s in samples:
    files = glob.glob(f"{base_path}{s}/*.slcio")
    count = 0
    for f in files:
        if MAX_EVENTS > 0 and count >= MAX_EVENTS: break
        reader.open(f)
        for event in reader:
            if MAX_EVENTS > 0 and count >= MAX_EVENTS: break
            try:
                mcps = event.getCollection("MCParticle")
                # Updated to SiTracks_Refitted
                tracks = event.getCollection("SiTracks_Refitted")
                clusters = event.getCollection("PandoraClusters")
                ecal_coll = event.getCollection("EcalBarrelCollectionRec")
                decoder = UTIL.BitField64(ecal_coll.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding))
            except: continue

            for mcp in mcps:
                pdg = abs(mcp.getPDG())
                if pdg not in [11, 211] or mcp.getGeneratorStatus() != 1: continue
                mcp_tlv = getTLV(mcp)
                if mcp_tlv.P() < 1.0 or abs(mcp_tlv.Eta()) > 2.4: continue
                
                best_p, min_dr_t = 0, DR_CUT
                for trk in tracks:
                    trk_tlv = getTrackTLV(trk, m=0.000511 if pdg==11 else 0.139, b_field=B_FIELD)
                    dr = trk_tlv.DeltaR(mcp_tlv)
                    if dr < min_dr_t: min_dr_t, best_p = dr, getP(trk, b_field=B_FIELD)
                
                best_clu, min_dr_c = None, DR_CUT
                for clu in clusters:
                    clu_tlv = getClusterTLV(clu)
                    dr = clu_tlv.DeltaR(mcp_tlv)
                    if dr < min_dr_c: min_dr_c, best_clu = dr, clu

                if best_p > 0 and best_clu:
                    e_total = best_clu.getEnergy()
                    hists[s]["ep"].Fill(e_total / best_p)
                    layers = {}
                    hits = best_clu.getCalorimeterHits()
                    for h in hits:
                        decoder.setValue(int(h.getCellID0()))
                        l = decoder["layer"].value()
                        layers[l] = layers.get(l, 0.0) + h.getEnergy()
                    hists[s]["disc"].Fill(get_max_profile_discrepancy(layers, e_total))
            count += 1
        reader.close()

for s in samples: print(f"{s} | E/p Mean: {hists[s]['ep'].GetMean():.4f} | Disc Mean: {hists[s]['disc'].GetMean():.4f}")

plotHistograms({s: hists[s]["ep"] for s in samples}, "plots/matched_ep.png", xlabel="E/p", ylabel="Entries", atltext=["Muon Collider", "SiTracks_Refitted", "0-50 GeV"])
plotHistograms({s: hists[s]["disc"] for s in samples}, "plots/profile_discrepancy.png", xlabel="Max Profile Discrepancy", ylabel="Entries", atltext=["Muon Collider", "EM Profile Match", "Electron vs Pion"])

import os, math, ROOT, glob, pyLCIO
from pyLCIO import IOIMPL, EVENT, UTIL; from math import gamma, exp

exec(open("./plotHelper.py").read()); ROOT.gROOT.SetBatch(); os.makedirs("plots", exist_ok=True)

# CONFIG & HISTOS
TARGET_PDG, B_FIELD, MAX_EVENTS, DR_CUT = 11, 5, -1, 0.2
file_path = "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_4.slcio"
h_ep_match = ROOT.TProfile("h_ep_match", "E/p vs Eta; #eta; E/p", 50, -2.5, 2.5)
h_ep_pfo = ROOT.TProfile("h_ep_pfo", "E/p vs Eta; #eta; E/p", 50, -2.5, 2.5)
h_disc = ROOT.TH1F("h_disc", "Profile Discrepancy; Profile Discrepancy; Entries", 100, 0, 1.0)

def generate_expected_em_profile(num_layers=60, energy=10.0, X0_scale=0.628):
    # Calculate a (shower max) and b (tail)
    a, b, expected_profile, norm = 1.0 + 0.5 * math.log(energy / 0.01), 0.5, {}, 0.0
    for i in range(num_layers): 
        t = i * X0_scale 
        f = (b * (b * t)**(a-1) * exp(-b * t)) / gamma(a)
        expected_profile[i] = f
        norm += f
    
    if norm > 0:
        res = {i: f/norm for i, f in expected_profile.items()}
        # PRINT TO TERMINAL: Find which layer has the highest energy fraction
        peak_layer = max(res, key=res.get)
        print(f"[DEBUG] E={energy:.1f} GeV | X0={X0_scale} | Peak Layer={peak_layer} | Peak Val={res[peak_layer]:.4f}")
        return res
    return {}

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader(); reader.open(file_path); count = 0
for event in reader:
    if MAX_EVENTS > 0 and count >= MAX_EVENTS: break
    try:
        mcps, tracks = event.getCollection("MCParticle"), event.getCollection("SiTracks_Refitted")
        clusters, ecal_coll = event.getCollection("PandoraClusters"), event.getCollection("EcalBarrelCollectionRec")
        pfos = event.getCollection("PandoraPFOs")
        decoder = UTIL.BitField64(ecal_coll.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding))
    except: continue

    # PFO Loop
    for pfo in pfos:
        if abs(pfo.getType()) == 11:
            e_calo = sum([c.getEnergy() for c in pfo.getClusters()])
            # Updated to use plotHelper getTLV and .P() for momentum
            pfo_tlv = getTLV(pfo)
            p_track = pfo_tlv.P()
            
            if p_track > 0:
                h_ep_pfo.Fill(pfo_tlv.Eta(), e_calo / p_track)

    # Truth-Matching Loop
    for mcp in [m for m in mcps if abs(m.getPDG()) == TARGET_PDG and m.getGeneratorStatus() == 1]:
        mcp_tlv = getTLV(mcp)
        if abs(mcp_tlv.Eta()) > 2.4: continue
        
        m_trk, min_dr_t = None, DR_CUT
        for trk in tracks:
            dr = getTrackTLV(trk, 0.000511, B_FIELD).DeltaR(mcp_tlv)
            if dr < min_dr_t: min_dr_t, m_trk = dr, trk

        m_clu, min_dr_c = None, DR_CUT
        for clu in clusters:
            dr = getClusterTLV(clu).DeltaR(mcp_tlv)
            if dr < min_dr_c: min_dr_c, m_clu = dr, clu

        if m_trk and m_clu:
            p, e, lyrs, total_e = getP(m_trk, B_FIELD), m_clu.getEnergy(), {}, 0.0
            if p > 0: h_ep_match.Fill(mcp_tlv.Eta(), e/p)
            for h in m_clu.getCalorimeterHits():
                decoder.setValue(int(h.getCellID0())); l = decoder["layer"].value()
                lyrs[l] = lyrs.get(l, 0.0) + h.getEnergy(); total_e += h.getEnergy()
            
            if total_e > 0 and len(lyrs) >= 3:
                # Using the scale we calculated
                expected = generate_expected_em_profile(energy=total_e, X0_scale=0.628)
                max_layer_disc = 0.0
                for lyr, obs_e in lyrs.items():
                    layer_disc = abs((obs_e / total_e) - expected.get(lyr, 0.0))
                    if layer_disc > max_layer_disc: max_layer_disc = layer_disc
                h_disc.Fill(max_layer_disc)
    count += 1

reader.close()
plotHistograms({"Truth-Matched": h_ep_match, "Pandora PFOs": h_ep_pfo}, "plots/analysis_ep.png", xlabel="#eta", ylabel="<E/p>", atltext=["Muon Collider", "Electron Gun", "E/p Comparison"])
plotHistograms({"0-50 GeV": h_disc}, "plots/analysis_disc.png", xlabel="Profile Discrepancy", ylabel="Entries", logy=True, atltext=["Muon Collider", "Simulation, no BIB", "|#eta| < 2.4", "MAIA Detector Concept"])

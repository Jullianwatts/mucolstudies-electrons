import pyLCIO, ROOT, os, plotHelper, math
ROOT.gROOT.SetBatch()

TARGET_EVENT, FILE_PATH = 29, "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_9.slcio"
OUT_DIR = "plots2026"
if not os.path.exists(OUT_DIR): os.makedirs(OUT_DIR)

def getDR(v1, v2): return math.sqrt((v1.Eta()-v2.Eta())**2 + v1.DeltaPhi(v2)**2)

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(FILE_PATH)

for i_evt, event in enumerate(reader):
    if i_evt != TARGET_EVENT: continue
    
    mc_coll = event.getCollection("MCParticle")
    electrons = [m for m in mc_coll if abs(m.getPDG()) == 11 and m.getGeneratorStatus() == 1]
    if not electrons: break
    
    mcp = electrons[0]; p = mcp.getMomentum()
    ref = ROOT.TLorentzVector(); ref.SetPxPyPzE(p[0], p[1], p[2], mcp.getEnergy())
    
    print(f"\n--- EVENT {i_evt} | TRUTH: E={ref.E():.2f} Eta={ref.Eta():.2f} ---")
    print(f"{'Type':<10} | {'E':<8} | {'Eta':<8} | {'Phi':<8} | {'DR_Truth':<10} | {'PDG':<5}\n" + "-"*70)

    pfo_tlvs = []
    pfo_coll = event.getCollection("PandoraPFOs")
    for pfo in pfo_coll:
        tlv = ROOT.TLorentzVector()
        p = pfo.getMomentum()
        tlv.SetPxPyPzE(p[0], p[1], p[2], pfo.getEnergy())
        dr_truth = getDR(ref, tlv)
        
        if dr_truth < 0.5:
            pfo_tlvs.append(tlv)
            print(f"{'PFO':<10} | {tlv.E():<8.2f} | {tlv.Eta():<8.2f} | {tlv.Phi():<8.2f} | {dr_truth:<10.3f} | {pfo.getType():<5}")
    min_pfo_dr = 999.0
    if len(pfo_tlvs) > 1:
        for i in range(len(pfo_tlvs)):
            for j in range(i + 1, len(pfo_tlvs)):
                dist = getDR(pfo_tlvs[i], pfo_tlvs[j])
                if dist < min_pfo_dr: min_pfo_dr = dist
        print(f"\nhow close are the pfos in this event? : {min_pfo_dr:.4f}")
    else:
        print(f"\nhow close are the pfos in this event? : N/A (Only {len(pfo_tlvs)} PFO near truth)")

    print(f"Summary: Tracks: {event.getCollection('SelectedTracks').getNumberOfElements()}, Total PFOs: {pfo_coll.getNumberOfElements()}")
    break

reader.close()

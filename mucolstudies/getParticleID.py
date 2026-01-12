import math
import pyLCIO
import ROOT
from pyLCIO import EVENT, IOIMPL

# --- Configuration ---
file_path = "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_3.slcio"
pdg_map = {11: "Electron", 13: "Muon", 22: "Photon", 211: "Pion", 2112: "Neutron"}

# Stats storage: stats[region][pdg] = [sum_ep, count]
stats = {
    "Barrel": {},
    "Endcap": {}
}

def get_tlv_from_mcp(mcp):
    v = ROOT.TLorentzVector()
    v.SetPxPyPzE(mcp.getMomentum()[0], mcp.getMomentum()[1], mcp.getMomentum()[2], mcp.getEnergy())
    return v

reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(file_path)

print(f"{'Evt':<4} | {'Region':<8} | {'Type':<8} | {'E_calo':<8} | {'P_track':<8} | {'E/p':<7} | {'dR'}")
print("-" * 85)

for i, event in enumerate(reader):
    if i >= 1000: break # Process up to 1000 events for statistics

    mcp_coll = event.getCollection("MCParticle")
    pfo_coll = event.getCollection("PandoraPFOs")

    # 1. Identify Truth Electron for Region and Matching
    truth_el_mcp = None
    for mcp in mcp_coll:
        if abs(mcp.getPDG()) == 11 and mcp.getGeneratorStatus() == 1:
            truth_el_mcp = mcp
            break
    if not truth_el_mcp: continue
    
    truth_tlv = get_tlv_from_mcp(truth_el_mcp)

    # Calculate Theta manually to avoid AttributeError
    p_truth = truth_el_mcp.getMomentum()
    px, py, pz = p_truth[0], p_truth[1], p_truth[2]
    pt = math.sqrt(px**2 + py**2)
    theta_rad = math.atan2(pt, pz)
    theta_deg = math.degrees(theta_rad)

    # MAIA Geometry Transition: ~38.8 degrees
    # Barrel is the central region; Endcaps are the 'plugs' at the ends
    if (38.8 < theta_deg < 141.2):
        region = "Barrel"
    else:
        region = "Endcap"

    # 2. Analyze PFOs matched to truth (within cone)
    for pfo in pfo_coll:
        pfo_mom = pfo.getMomentum()
        pfo_tlv = ROOT.TLorentzVector()
        pfo_tlv.SetPxPyPzE(pfo_mom[0], pfo_mom[1], pfo_mom[2], pfo.getEnergy())

        dr = truth_tlv.DeltaR(pfo_tlv)
        
        # Only look at PFOs in the electron's path (shower + fragments)
        if dr < 0.4:
            p_type = abs(pfo.getType())
            name = pdg_map.get(p_type, f"ID({p_type})")
            
            # RAW CALO ENERGY: Summing clusters associated with PFO
            e_calo = sum([c.getEnergy() for c in pfo.getClusters()])
            
            # TRACK MOMENTUM: Magnitude of the PFO momentum vector (linked to track)
            p_track = math.sqrt(pfo_mom[0]**2 + pfo_mom[1]**2 + pfo_mom[2]**2)
            
            # Print first 20 events line-by-line for debugging
            if i < 20:
                ep_val = f"{e_calo/p_track:.3f}" if p_track > 0 else "N/A"
                print(f"{i:<4} | {region:<8} | {name:<8} | {e_calo:<8.2f} | {p_track:<8.2f} | {ep_val:<7} | {dr:.3f}")

            # Collect stats for the Final Summary
            if p_track > 0:
                if p_type not in stats[region]:
                    stats[region][p_type] = [0.0, 0] # [sum_ep, count]
                stats[region][p_type][0] += (e_calo / p_track)
                stats[region][p_type][1] += 1

reader.close()

# --- FINAL SUMMARY ---
print("\n" + "="*60)
print("MAIA DETECTOR CALIBRATION SUMMARY (E/p)")
print("="*60)
for reg in ["Barrel", "Endcap"]:
    print(f"\nREGION: {reg}")
    print(f"{'Particle Type':<15} | {'Average E/p':<12} | {'Count'}")
    print("-" * 50)
    for p_type in sorted(stats[reg].keys()):
        data = stats[reg][p_type]
        avg_ep = data[0] / data[1]
        name = pdg_map.get(p_type, f"ID({p_type})")
        print(f"{name:<15} | {avg_ep:<12.4f} | {data[1]}")

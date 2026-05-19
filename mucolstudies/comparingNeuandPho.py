import math, ROOT, pyLCIO, os
exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

INPUT_FILE = "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_7.slcio"
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(INPUT_FILE)

print(f"{'EVT':<5} | {'E-RECO?':<8} | {'MIS-RECO':<10} | {'TRUE E':<8} | {'MIS-E':<8} | {'ECAL(0)':<8} | {'HCAL(1)':<8}")
print("-" * 80)

for ievt, event in enumerate(reader):
    mcps = [m for m in event.getCollection("MCParticle") if m.getGeneratorStatus() == 1 and abs(m.getPDG()) == 11 and m.getEnergy() > 10]
    if not mcps: continue
    true_e = mcps[0].getEnergy()

    pfos = event.getCollection("PandoraPFOs")
    has_electron = any(abs(p.getType()) == 11 for p in pfos)
    e_status = "YES" if has_electron else "NO"

    for pfo in pfos:
        pdg = abs(pfo.getType())
        # Skip the electron; show photons, neutrons, and pions > 2 GeV
        if pdg == 11 or pfo.getEnergy() < 2: continue
        
        clusters = pfo.getClusters()
        if clusters.size() == 0: continue
        
        p_label = "PHOTON" if pdg == 22 else ("NEUTRON" if pdg == 2112 else ("PION" if pdg == 211 else str(pdg)))
        cl = clusters[0]
        sub_e = cl.getSubdetectorEnergies()
        e0 = sub_e[0] if sub_e.size() > 0 else 0.0
        e1 = sub_e[1] if sub_e.size() > 1 else 0.0

        print(f"{ievt:<5} | {e_status:<8} | {p_label:<10} | {true_e:<8.1f} | {pfo.getEnergy():<8.1f} | {e0:<8.2f} | {e1:<8.2f}")

reader.close()

import math
import ROOT
import pyLCIO
import os

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

INPUT_FILE = "/scratch/jwatts/mucol/v11Container/reco/1kelectron25_reco.slcio"
#INPUT_FILE = "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_9.slcio"
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(INPUT_FILE)

categories = [
    "Nothing Reconstructed", "Just Electrons", "Just Pions", "Just Neutrons",
    "Just Photons", "Just Photons + Neutrons", "Ele + Photon", "Ele + Neutron", "Ele + Photon + Neutron"
]
counts = {cat: 0 for cat in categories}
total_events = 0

for event in reader:
    mcps = [m for m in event.getCollection("MCParticle") if m.getGeneratorStatus() == 1 and abs(m.getPDG()) == 11 and m.getEnergy() > 0]
    if not mcps: continue

    total_events += 1
    pfos = event.getCollection("PandoraPFOs")

    e = any(abs(p.getType()) == 11 and p.getEnergy() > 0 for p in pfos)
    pi = any(abs(p.getType()) == 211 and p.getEnergy() > 0 for p in pfos)
    n = any(abs(p.getType()) == 2112 and p.getEnergy() > 0 for p in pfos)
    ph = any(abs(p.getType()) == 22 and p.getEnergy() > 0 for p in pfos)

    if not e and not pi and not n and not ph: counts["Nothing Reconstructed"] += 1
    elif e and not pi and not n and not ph:     counts["Just Electrons"] += 1
    elif pi and not e and not n and not ph:    counts["Just Pions"] += 1
    elif n and not e and not pi and not ph:     counts["Just Neutrons"] += 1
    elif ph and not e and not pi and not n:    counts["Just Photons"] += 1
    elif ph and n and not e and not pi:        counts["Just Photons + Neutrons"] += 1
    elif e and ph and not n and not pi:        counts["Ele + Photon"] += 1
    elif e and n and not ph and not pi:        counts["Ele + Neutron"] += 1
    elif e and ph and n and not pi:            counts["Ele + Photon + Neutron"] += 1

reader.close()

if total_events > 0:
    print(f"\n===== EVENT ANALYSIS (Total Events: {total_events}) =====")
    print(f"{'Category':<30} | {'Freq (%)':<10}")
    print("-" * 45)
    for cat, count in counts.items():
        print(f"{cat:<30} | {(count / total_events) * 100:>8.1f}%")

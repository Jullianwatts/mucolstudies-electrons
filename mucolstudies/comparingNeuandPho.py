import math, ROOT, pyLCIO, os
exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)
INPUT_FILE = "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_7.slcio"

hists = {}
for p in ["photon", "neutron"]:
    hists[f"{p}_nhits"] = ROOT.TH1F(f"{p}_nhits", ";Number of Hits;Counts", 60, 0, 300)
    hists[f"{p}_dr"]    = ROOT.TH1F(f"{p}_dr", ";#DeltaR(PFO, MC);Counts", 50, 0, 0.2)
    hists[f"{p}_E"]     = ROOT.TH1F(f"{p}_E", ";Energy [GeV];Counts", 50, 0, 60)

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(INPUT_FILE)
first_event = True

for event in reader:
    if first_event:
        try:
            names = event.getCollection("PandoraClusters").getParameters().getStringVals("ClusterSubdetectorNames")
            for i, name in enumerate(names): print(f"Index {i}: {name}")
        except: print("Mapping not found")
        first_event = False

    mcps = [(abs(m.getPDG()), getTLV(m)) for m in event.getCollection("MCParticle") if m.getGeneratorStatus() == 1 and abs(m.getPDG()) in [22, 2112]]
    pfos = event.getCollection("PandoraPFOs")

    for pfo in pfos:
        pdg = abs(pfo.getType())
        if pdg not in [22, 2112] or pfo.getEnergy() < 2: continue
        
        p_label = "photon" if pdg == 22 else "neutron"
        pfo_tlv = getTLV(pfo)
        cl = pfo.getClusters()
        if cl.size() == 0: continue
        
        hists[f"{p_label}_nhits"].Fill(len(cl[0].getCalorimeterHits()))
        hists[f"{p_label}_E"].Fill(pfo_tlv.E())

        drs = [pfo_tlv.DeltaR(m[1]) for m in mcps if m[0] == pdg]
        if drs: hists[f"{p_label}_dr"].Fill(min(drs))

reader.close()

for m, xl in [("nhits", "Hits"), ("dr", "#DeltaR"), ("E", "Energy [GeV]")]:
    plotHistograms({"Photons": hists[f"photon_{m}"], "Neutrons": hists[f"neutron_{m}"]}, 
                   f"{PLOT_DIR}/CMP_{m}.png", xlabel=xl, ylabel="Entries", atltext="MAIA PFO Study")

for p in ["photon", "neutron"]:
    h = hists[f"{p}_nhits"]
    print(f"{p.upper()}: {int(h.GetEntries())} entries, {h.GetMean():.1f} avg hits")

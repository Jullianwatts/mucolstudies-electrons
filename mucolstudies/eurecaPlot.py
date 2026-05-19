import math, pyLCIO, ROOT, os, numpy as np
from pyLCIO import EVENT
pdg_map = {11:"Electrons", 13:"Muons", 22:"Photons", 211:"Charged Pions", 130:"Neutral Kaons", 2112:"Neutrons", 2212:"Protons"}
f_in, p_out = "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_7.slcio", "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(p_out, exist_ok=True)
ROOT.gROOT.SetBatch()
hists, ev_count, n_bins, xmin, xmax = {}, 0, 50, 2.0, 200.0
log_bins = np.logspace(np.log10(xmin), np.log10(xmax), n_bins + 1)
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
try:
    reader.open(f_in)
    for ev in reader:
        if any(abs(m.getPDG())==11 and m.getGeneratorStatus()==1 and m.getEnergy()>10 for m in ev.getCollection("MCParticle")):
            ev_count += 1
            if "PandoraPFOs" in ev.getCollectionNames():
                for p in ev.getCollection("PandoraPFOs"):
                    p_type, en = abs(p.getType()), p.getEnergy()
                    if p_type != 11 and en >= 2.0:
                        nm = pdg_map.get(p_type, f"Other_{p_type}")
                        if nm not in hists: hists[nm] = ROOT.TH1F(f"h_{nm}", ";Energy [GeV];Counts", n_bins, log_bins)
                        hists[nm].Fill(en)
    reader.close()
except Exception as e: print(e)
c = ROOT.TCanvas("c", "c", 800, 600)
c.SetLogx(); c.SetLogy()
s = ROOT.THStack("s", "")
l = ROOT.TLegend(0.55, 0.65, 0.85, 0.85); l.SetBorderSize(0); l.SetFillStyle(0)
cols = [ROOT.kRed+1, ROOT.kAzure+1, ROOT.kSpring-5, ROOT.kOrange+7, ROOT.kMagenta-4]
for i, (nm, h) in enumerate(hists.items()):
    h.SetLineColor(cols[i%5]); h.SetLineWidth(2); s.Add(h); l.AddEntry(h, nm, "l")
s.Draw("nostack"); s.SetMinimum(0.8); s.SetMaximum(s.GetMaximum("nostack")*20)
ax, ay = s.GetXaxis(), s.GetYaxis()
ax.SetTitle("Energy [GeV]"); ay.SetTitle("Counts")
ax.SetMoreLogLabels(); ax.SetNoExponent()
l.Draw(); c.SaveAs(os.path.join(p_out, "pfo_energy_compact_reco7.svg"))

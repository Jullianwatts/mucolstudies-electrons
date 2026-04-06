import math, ROOT, pyLCIO, os
exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR, exist_ok=True)
INPUT_FILE = "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_7.slcio"

DR_CONE = 0.1
ECAL_IDX = 0
HCAL_IDX = 1
TOLERANCE = 0.05 

n_events = 0
n_real_em = 0
n_with_gam = 0
n_with_neu = 0

h_Ratio = ROOT.TH1F("h_Ratio", "Conservation; (E_{e} + E_{#gamma}) / E_{truth}; Events", 100, 0, 2.0)

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(INPUT_FILE)

print(f"{'Event':<6} | {'TruthE':<7} | {'ECAL':<7} | {'HCAL':<7} | {'NeuEM/Had':<12} | {'Ratio':<5}")
print("-" * 65)

for event in reader:
    n_events += 1
    mcps = [m for m in event.getCollection("MCParticle") if m.getGeneratorStatus() == 1 and abs(m.getPDG()) == 11 and m.getEnergy() > 10]
    pfos = event.getCollection("PandoraPFOs")

    for m in mcps:
        e_true = m.getEnergy()
        mc_tlv = getTLV(m)
        e_ele, e_gam, e_neu_em, e_neu_had, s_ecal, s_hcal = 0, 0, 0, 0, 0, 0
        has_gam, has_neu = False, False

        for pfo in pfos:
            if getTLV(pfo).DeltaR(mc_tlv) < DR_CONE:
                pdg = abs(pfo.getType())
                p_e = pfo.getEnergy()
                p_ecal, p_hcal = 0, 0
                if pfo.getClusters().size() > 0:
                    cl = pfo.getClusters()[0]
                    p_ecal = cl.getSubdetectorEnergies()[ECAL_IDX]
                    p_hcal = cl.getSubdetectorEnergies()[HCAL_IDX]
                s_ecal += p_ecal
                s_hcal += p_hcal
                if pdg == 11: e_ele += p_e
                elif pdg == 22: 
                    e_gam += p_e
                    has_gam = True
                elif pdg == 2112:
                    e_neu_em += p_ecal
                    e_neu_had += p_hcal
                    has_neu = True

        ratio = (e_ele + e_gam) / e_true if e_true > 0 else 0
        h_Ratio.Fill(ratio)
        if 1.0-TOLERANCE <= ratio <= 1.0+TOLERANCE: n_real_em += 1
        if has_gam: n_with_gam += 1
        if has_neu: n_with_neu += 1

        if n_events <= 20:
            neu_s = f"{e_neu_em:.1f}/{e_neu_had:.1f}"
            print(f"{n_events:<6} | {e_true:<7.1f} | {s_ecal:<7.1f} | {s_hcal:<7.1f} | {neu_s:<12} | {ratio:<5.2f}")

reader.close()

real_freq = (n_real_em / n_events * 1000) if n_events > 0 else 0
print(f"\nSummary per 1000 events:")
print(f"- Conserved EM (Real): {real_freq:.1f}")
print(f"- Events w/ Photons:  {n_with_gam}")
print(f"- Events w/ Neutrons: {n_with_neu}")

plotHistograms({"EM_Balance": h_Ratio}, f"{PLOT_DIR}/balance.png", xlabel="E(reco)/E(truth)")

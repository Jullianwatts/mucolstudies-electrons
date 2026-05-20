import math,glob,ROOT,pyLCIO,os

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

PLOT_DIR="/scratch/jwatts/mucol/mucolstudies/plots2026"
os.makedirs(PLOT_DIR,exist_ok=True)

samples=glob.glob("/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_7.slcio")
files={"electronGun_pT_0_50":samples}

DR_MATCH=0.2
REC_MIN=0.8
REC_MAX=1.2

BARREL_MIN=0.698
BARREL_MAX=2.443
ENDCAP1_MIN=0.139
ENDCAP1_MAX=0.698
ENDCAP2_MIN=2.443
ENDCAP2_MAX=3.002

hists={}

for s in files:

    hists[s]={}

    for obj in ["mcp_el","baseline","recovered"]:

        hists[s][obj+"_barrel"]=ROOT.TH1F(f"{s}_{obj}_barrel",";#theta [rad];Counts",50,BARREL_MIN,BARREL_MAX)
        hists[s][obj+"_endcap"]=ROOT.TH1F(f"{s}_{obj}_endcap",";#theta [rad];Counts",50,ENDCAP1_MIN,ENDCAP2_MAX)

reader=pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

for s in files:

    for f in files[s]:

        print(f"Opening {f}")

        reader.open(f)

        for event in reader:

            mcps=event.getCollection("MCParticle")
            trks=event.getCollection("SelectedTracks")
            clusters=event.getCollection("PandoraClusters")

            mcp_electrons=[]

            for mcp in mcps:

                if mcp.getGeneratorStatus()!=1 or abs(mcp.getPDG())!=11: continue

                tlv=getTLV(mcp)

                if tlv.E()<10: continue

                theta=tlv.Theta()

                if BARREL_MIN<theta<BARREL_MAX:

                    hists[s]["mcp_el_barrel"].Fill(theta)
                    suffix="_barrel"

                elif (ENDCAP1_MIN<theta<ENDCAP1_MAX) or (ENDCAP2_MIN<theta<ENDCAP2_MAX):

                    hists[s]["mcp_el_endcap"].Fill(theta)
                    suffix="_endcap"

                else:
                    continue

                mcp_electrons.append((tlv,tlv.E(),suffix))

            trk_tlvs=[]

            for t in trks:

                tlv=getTrackTLV(t)

                if tlv.Pt()>2: trk_tlvs.append(tlv)

            cluster_tlvs=[]

            for c in clusters:

                if c.getEnergy()>2: cluster_tlvs.append(getClusterTLV(c))

            for mcp_el,trueE,suffix in mcp_electrons:

                has_track=False
                has_cluster=False

                for t in trk_tlvs:

                    if mcp_el.DeltaR(t)<DR_MATCH:
                        has_track=True
                        break

                for c in cluster_tlvs:

                    if mcp_el.DeltaR(c)<DR_MATCH:
                        has_cluster=True
                        break

                if has_track and has_cluster:

                    hists[s]["baseline"+suffix].Fill(mcp_el.Theta())
                    hists[s]["recovered"+suffix].Fill(mcp_el.Theta())

                    continue

                sumE=0

                for c in cluster_tlvs:
                    sumE+=c.E()

                ratio=sumE/trueE if trueE>0 else 0

                if REC_MIN<ratio<REC_MAX:
                    hists[s]["recovered"+suffix].Fill(mcp_el.Theta())

        reader.close()

print("\n--- Efficiency Results ---")

for s in hists:

    for region in ["barrel","endcap"]:

        denom=hists[s]["mcp_el_"+region]
        num_base=hists[s]["baseline_"+region]
        num_rec=hists[s]["recovered_"+region]

        n_den=denom.GetEntries()

        base_eff=num_base.GetEntries()/n_den if n_den>0 else 0
        rec_eff=num_rec.GetEntries()/n_den if n_den>0 else 0

        region_name="0.698<θ<2.443" if region=="barrel" else "0.139<θ<0.698 and 2.443<θ<3.002"

        print(f"Region: {region_name} | Baseline: {base_eff:.4f} | Recovered: {rec_eff:.4f}")

        teff1=ROOT.TEfficiency(num_base,denom)
        teff2=ROOT.TEfficiency(num_rec,denom)

        plotEfficiencies(
            {"Baseline":teff1.CreateGraph(),"Recovered":teff2.CreateGraph()},
            os.path.join(PLOT_DIR,f"ELECTRON_RECOVERY_{region}.png"),
            xlabel="#theta [rad]",
            ylabel="Efficiency"
        )

print(f"\nPlots saved to {PLOT_DIR}")

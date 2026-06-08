import math
import glob
import ROOT
import pyLCIO
import os

exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()

samples = glob.glob(
    "/scratch/jwatts/mucol/data/reco/electronGun_pT_0_50/electronGun_pT_0_50_reco_7.slcio"
)

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

total_mcp = 0
matched_mcp = 0
unmatched_mcp = 0

unmatched_with_tracks = 0
unmatched_without_tracks = 0

for f in samples:

    print(f"\nOpening {f}")

    reader.open(f)

    for ievt, event in enumerate(reader):

        mcps = event.getCollection("MCParticle")
        tracks = event.getCollection("SelectedTracks")

        # Build SelectedTrack TLVs
        track_tlvs = []

        for trk in tracks:

            tlv = getTrackTLV(trk)

            if tlv.E() < 10:
                continue

            if abs(tlv.Eta()) > 2.4:
                continue

            track_tlvs.append(tlv)

        # Loop over MCP electrons
        for mcp in mcps:

            if mcp.getGeneratorStatus() != 1:
                continue

            if abs(mcp.getPDG()) != 11:
                continue

            mcp_tlv = getTLV(mcp)

            if mcp_tlv.E() < 10:
                continue

            if abs(mcp_tlv.Eta()) > 2.4:
                continue

            total_mcp += 1

            matched = False

            for trk_tlv in track_tlvs:

                if mcp_tlv.DeltaR(trk_tlv) < 0.2:

                    matched = True
                    matched_mcp += 1
                    break

            # If NO matched SelectedTrack
            if not matched:

                unmatched_mcp += 1

                if len(track_tlvs) > 0:
                    unmatched_with_tracks += 1

                    print(f"Event {ievt}")
                    print("Unmatched MCP electron")
                    print(f"BUT event HAS {len(track_tlvs)} SelectedTracks")

                else:
                    unmatched_without_tracks += 1

    reader.close()

print(f"Total MCP electrons           : {total_mcp}")
print(f"Matched MCP electrons         : {matched_mcp}")
print(f"Unmatched MCP electrons       : {unmatched_mcp}")
print("")
print(f"Unmatched + event has tracks  : {unmatched_with_tracks}")
print(f"Unmatched + no tracks in evt  : {unmatched_without_tracks}")

if unmatched_mcp > 0:

    frac_with = unmatched_with_tracks / unmatched_mcp
    frac_without = unmatched_without_tracks / unmatched_mcp

    print("")
    print(f"Fraction with tracks    : {frac_with:.4f}")
    print(f"Fraction without tracks : {frac_without:.4f}")


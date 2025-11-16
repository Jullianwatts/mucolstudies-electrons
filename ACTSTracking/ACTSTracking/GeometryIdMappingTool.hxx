#pragma once

#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>

#include <UTIL/CellIDDecoder.h>

#include <unordered_map>

namespace ACTSTracking {

//! \brief Maps DD4hep cell ID's to ACTS geometry ID's
/**
 * Lots of hardcoded values that match up to mod5 geometry
 * dumped to TGeo.
 *
 * @author Karol Krizka
 * @version $Id$
 */
class GeometryIdMappingTool {
 public:
  enum class DetSchema : char { MuColl_v1, MAIA_v0, MuSIC_v1, MuSIC_v2 };

  using modules_map = std::unordered_map<uint32_t, uint32_t>;
  using det_mod_map = std::unordered_map<DetSchema, modules_map>;

  /** Create a mapping tool using the provided encoderString to
   * interpret cell ID's.
   */
  GeometryIdMappingTool(const std::string& encoderString,
                        DetSchema dType = DetSchema::MuColl_v1);

  /** Decode hit
   */
  uint64_t getGeometryID(const lcio::SimTrackerHit* hit);
  uint64_t getGeometryID(const lcio::TrackerHit* hit);

  uint64_t getGeometryID(uint32_t systemID, uint32_t layerID, int32_t sideID,
                         uint32_t ladderID, uint32_t moduleID);

 private:
  std::string _encoderString;
  const DetSchema det_type;

  static const std::unordered_map<int32_t, uint32_t> VolumeMap;

  static const int32_t VertexEndCapNegative;
  static const int32_t VertexBarrel;
  static const int32_t VertexEndCapPositive;
  static const int32_t InnerTrackerEndCapNegative;
  static const int32_t InnerTrackerBarrel;
  static const int32_t InnerTrackerEndCapPositive;
  static const int32_t OuterInnerTrackerEndCapNegative;
  static const int32_t OuterInnerTrackerBarrel;
  static const int32_t OuterInnerTrackerEndCapPositive;
  static const int32_t OuterTrackerEndCapNegative;
  static const int32_t OuterTrackerBarrel;
  static const int32_t OuterTrackerEndCapPositive;

  // Modules in phi ladder per layer
  static const det_mod_map NLad_VertexBarrel;
  static const det_mod_map NLad_InnerTrackerBarrel;
  static const det_mod_map NLad_OuterInnerTrackerBarrel;
  static const det_mod_map NLad_OuterTrackerBarrel;

  // Modules in ring per layer
  static const det_mod_map NRng_VertexEndCap;
  static const det_mod_map NRng_InnerTrackerEndCap;
  static const det_mod_map NRng_OuterInnerTrackerEndCap;
  static const det_mod_map NRng_OuterTrackerEndCap;
};

}  // namespace ACTSTracking

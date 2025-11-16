#include "ACTSProcBase.hxx"

#include <TGeoManager.h>

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/Detector.h>

#include <UTIL/LCTrackerConf.h>

#include <Acts/Definitions/Units.hpp>
#include <Acts/Geometry/CylinderVolumeBuilder.hpp>
#include <Acts/Geometry/CylinderVolumeHelper.hpp>
#include <Acts/Geometry/ITrackingVolumeBuilder.hpp>
#include <Acts/Geometry/LayerArrayCreator.hpp>
#include <Acts/Geometry/LayerCreator.hpp>
#include <Acts/Geometry/ProtoLayerHelper.hpp>
#include <Acts/Geometry/SurfaceArrayCreator.hpp>
#include <Acts/Geometry/TrackingGeometryBuilder.hpp>
#include <Acts/Geometry/TrackingVolumeArrayCreator.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/Plugins/Json/JsonMaterialDecorator.hpp>
#include <Acts/Plugins/TGeo/TGeoDetectorElement.hpp>
#include <Acts/Plugins/TGeo/TGeoLayerBuilder.hpp>

#include "Helpers.hxx"

using namespace ACTSTracking;
using DetSchema = GeometryIdMappingTool::DetSchema;

ACTSProcBase::ACTSProcBase(const std::string& procname) : Processor(procname) {
  // configuration
  registerProcessorParameter(
      "MatFile", "Path to the material description JSON file. Can be empty.",
      _matFile, _matFile);

  registerProcessorParameter("TGeoFile", "Path to the tracker geometry file.",
                             _tgeoFile, _tgeoFile);

  registerProcessorParameter("TGeoDescFile", "Path to the JSON file describing the subdetectors.",
                             _tgeodescFile, _tgeodescFile);

  registerProcessorParameter("DetectorSchema", "Detector schema name (MuColl_v1, MAIA_v0, MuSIC_v1, MuSIC_v2).",
                             _detSchema, _detSchema);
}

std::shared_ptr<GeometryIdMappingTool> ACTSProcBase::geoIDMappingTool() const {
  return _geoIDMappingTool;
}

const Acts::MagneticFieldContext& ACTSProcBase::magneticFieldContext() const {
  return _magneticFieldContext;
}

const Acts::GeometryContext& ACTSProcBase::geometryContext() const {
  return _geometryContext;
}

const Acts::CalibrationContext& ACTSProcBase::calibrationContext() const {
  return _calibrationContext;
}

std::shared_ptr<Acts::MagneticFieldProvider> ACTSProcBase::magneticField()
    const {
  return _magneticField;
}

std::shared_ptr<const Acts::TrackingGeometry> ACTSProcBase::trackingGeometry()
    const {
  return _trackingGeometry;
}

const Acts::Surface* ACTSProcBase::findSurface(
    const EVENT::TrackerHit* hit) const {
  uint64_t moduleGeoId = _geoIDMappingTool->getGeometryID(hit);
  return _trackingGeometry->findSurface(moduleGeoId);
}

void ACTSProcBase::init() {
  // Parse parameters
  _matFile = findFile(_matFile);
  _tgeoFile = findFile(_tgeoFile);
  _tgeodescFile = findFile(_tgeodescFile);

  // Print the initial parameters
  printParameters();

  // Load geometry
  streamlog_out(MESSAGE) << " -------------------------------------"
                         << std::endl;

  streamlog_out(MESSAGE) << " -- Building magnetic field" << std::endl;
  buildBfield();
  streamlog_out(MESSAGE) << " -- Building tracking detector" << std::endl;
  buildDetector();

  streamlog_out(MESSAGE)  // << " ---- instantiated  geometry for detector " <<
                          // theDetector.header().name()  << std::endl
      << " -------------------------------------" << std::endl;

  // Initialize mapping tool
  DetSchema dSchema;
  if (_detSchema == "MuSIC_v1") {
    dSchema = DetSchema::MuSIC_v1;
  } else if (_detSchema == "MuSIC_v2") {
    dSchema = DetSchema::MuSIC_v2;
  } else if (_detSchema == "MAIA_v0") {
    dSchema = DetSchema::MAIA_v0;
  } else if (_detSchema == "MuColl_v1") {
    dSchema = DetSchema::MuColl_v1;
  } else {
    dSchema = DetSchema::MuColl_v1;
  }

  _geoIDMappingTool = std::make_shared<GeometryIdMappingTool>(
      lcio::LCTrackerCellID::encoding_string(), dSchema);
}

void ACTSProcBase::processRunHeader(LCRunHeader* run) {}

void ACTSProcBase::processEvent(LCEvent* evt) {}

void ACTSProcBase::check(LCEvent* evt) {}

void ACTSProcBase::end() {}

void ACTSProcBase::buildDetector() {
  // Logging
  Acts::Logging::Level surfaceLogLevel = Acts::Logging::INFO;
  Acts::Logging::Level layerLogLevel = Acts::Logging::INFO;
  Acts::Logging::Level volumeLogLevel = Acts::Logging::INFO;

  // Material description
  std::shared_ptr<const Acts::IMaterialDecorator> matDeco = nullptr;
  if (!_matFile.empty()) {
    // Set up the converter first
    Acts::MaterialMapJsonConverter::Config jsonGeoConvConfig;
    // Set up the json-based decorator
    matDeco = std::make_shared<const Acts::JsonMaterialDecorator>(
        jsonGeoConvConfig, _matFile, Acts::Logging::INFO);
  }

  // Geometry
  TGeoManager* gGeoManagerOld = nullptr;
  if (!_tgeoFile.empty()) {
    // Save current geometry. This is needed by all the other Processors
    gGeoManagerOld = gGeoManager;
    gGeoManager = nullptr;  // prevents it from being deleted

    // Load new geometry
    TGeoManager::Import(_tgeoFile.c_str());
  }

  // configure surface array creator
  Acts::SurfaceArrayCreator::Config sacConfig;
  auto surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>(
      sacConfig,
      Acts::getDefaultLogger("SurfaceArrayCreator", surfaceLogLevel));

  // configure the proto layer helper
  Acts::ProtoLayerHelper::Config plhConfig;
  auto protoLayerHelper = std::make_shared<const Acts::ProtoLayerHelper>(
      plhConfig, Acts::getDefaultLogger("ProtoLayerHelper", layerLogLevel));

  // configure the layer creator that uses the surface array creator
  Acts::LayerCreator::Config lcConfig;
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator = std::make_shared<const Acts::LayerCreator>(
      lcConfig, Acts::getDefaultLogger("LayerCreator", layerLogLevel));

  // configure the layer array creator
  Acts::LayerArrayCreator::Config lacConfig;
  auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
      lacConfig, Acts::getDefaultLogger("LayerArrayCreator", layerLogLevel));

  // tracking volume array creator
  Acts::TrackingVolumeArrayCreator::Config tvacConfig;
  auto tVolumeArrayCreator =
      std::make_shared<const Acts::TrackingVolumeArrayCreator>(
          tvacConfig,
          Acts::getDefaultLogger("TrackingVolumeArrayCreator", volumeLogLevel));

  // configure the cylinder volume helper
  Acts::CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper =
      std::make_shared<const Acts::CylinderVolumeHelper>(
          cvhConfig,
          Acts::getDefaultLogger("CylinderVolumeHelper", volumeLogLevel));

  //-------------------------------------------------------------------------------------
  // list the volume builders
  std::list<std::shared_ptr<const Acts::ITrackingVolumeBuilder>> volumeBuilders;

  /**
   * Disable for now
  // Create a beam pipe if configured to do so
  auto beamPipeParameters =
      vm["geo-tgeo-bp-parameters"].template as<read_range>();
  if (beamPipeParameters.size() > 2) {
    /// configure the beam pipe layer builder
    Acts::PassiveLayerBuilder::Config bplConfig;
    bplConfig.layerIdentification = "BeamPipe";
    bplConfig.centralLayerRadii = std::vector<double>(1, beamPipeParameters[0]);
    bplConfig.centralLayerHalflengthZ =
        std::vector<double>(1, beamPipeParameters[1]);
    bplConfig.centralLayerThickness =
        std::vector<double>(1, beamPipeParameters[2]);
    auto beamPipeBuilder = std::make_shared<const Acts::PassiveLayerBuilder>(
        bplConfig,
        Acts::getDefaultLogger("BeamPipeLayerBuilder", layerLogLevel));
    // create the volume for the beam pipe
    Acts::CylinderVolumeBuilder::Config bpvConfig;
    bpvConfig.trackingVolumeHelper = cylinderVolumeHelper;
    bpvConfig.volumeName = "BeamPipe";
    bpvConfig.layerBuilder = beamPipeBuilder;
    bpvConfig.layerEnvelopeR = {1. * Acts::UnitConstants::mm,
                                1. * Acts::UnitConstants::mm};
    bpvConfig.buildToRadiusZero = true;
    auto beamPipeVolumeBuilder =
        std::make_shared<const Acts::CylinderVolumeBuilder>(
            bpvConfig,
            Acts::getDefaultLogger("BeamPipeVolumeBuilder", volumeLogLevel));
    // add to the list of builders
    volumeBuilders.push_back(beamPipeVolumeBuilder);
  }
  */

  //
  // Detector definition
  std::vector<Acts::TGeoLayerBuilder::Config> layerBuilderConfigs;

  // Check if the geometry has been defined
  if (_tgeodescFile.empty()) {
    throw std::runtime_error("Required geometry description file (TGeoDescFile) missing.");
  }

  // Open the description
  nlohmann::json tgeodesc;
  std::ifstream tgeodescFile(_tgeodescFile, std::ifstream::in | std::ifstream::binary);
  if(!tgeodescFile.is_open()) {
    throw std::runtime_error("Unable to open TGeo description file: "+_tgeodescFile);
  }
  tgeodescFile >> tgeodesc;

  // Helper parsing functions
  auto range_from_json = [](const nlohmann::json& jsonpair) -> std::pair<double, double> {
    return {
      jsonpair["lower"], 
      jsonpair["upper"]
    };
  };

  // Loop over volumes to define sub-detectors
  for(const auto& volume : tgeodesc["Volumes"]) {
    // Volume information
    Acts::TGeoLayerBuilder::Config layerBuilderConfig;
    layerBuilderConfig.configurationName = volume["geo-tgeo-volume-name"];
    layerBuilderConfig.unit = 1 * Acts::UnitConstants::cm;
    layerBuilderConfig.autoSurfaceBinning = true;

    // AutoBinning
    std::vector<std::pair<double, double>> binTolerances{(int)Acts::binValues,
                                                         {0., 0.}};
    binTolerances[Acts::binR] = range_from_json(volume["geo-tgeo-sfbin-r-tolerance"]);
    binTolerances[Acts::binZ] = range_from_json(volume["geo-tgeo-sfbin-z-tolerance"]);
    binTolerances[Acts::binPhi] = range_from_json(volume["geo-tgeo-sfbin-phi-tolerance"]);
    layerBuilderConfig.surfaceBinMatcher =
        Acts::SurfaceBinningMatcher(binTolerances);

    // Loop over subvolumes (two endcaps and one barrel)
    std::array<std::string, 3> subvolumeNames = {"negative","central","positive"}; // List of possible subvolume names. Order corresponds to layerConfigurations.
    for(std::size_t idx = 0; idx < 3; idx++) {
      const std::string& subvolumeName = subvolumeNames[idx];
      if(!volume["geo-tgeo-volume-layers"][subvolumeName]) {
        // Skip disabled volume
        continue;
      }

      // Create the layer config object and fill it
      Acts::TGeoLayerBuilder::LayerConfig lConfig;
      lConfig.volumeName = volume["geo-tgeo-subvolume-names"][subvolumeName];
      lConfig.sensorNames = volume["geo-tgeo-sensitive-names"][subvolumeName];
      lConfig.localAxes = volume["geo-tgeo-sensitive-axes"][subvolumeName];
      lConfig.envelope = std::pair<double, double>(
          0.1 * Acts::UnitConstants::mm, 0.1 * Acts::UnitConstants::mm);

      // Fill the parsing restrictions in r
      lConfig.parseRanges.push_back({Acts::binR, range_from_json(volume["geo-tgeo-layer-r-ranges"][subvolumeName])});

      // Fill the parsing restrictions in z
      lConfig.parseRanges.push_back({Acts::binZ, range_from_json(volume["geo-tgeo-layer-z-ranges"][subvolumeName])});

      // Fill the layer splitting parameters in z
      float rsplit = volume["geo-tgeo-layer-r-split"][subvolumeName];
      if(rsplit > 0) {
        lConfig.splitConfigs.push_back({Acts::binR, rsplit});
      }

      // Fill the layer splitting parameters in z
      float zsplit = volume["geo-tgeo-layer-z-split"][subvolumeName];
      if(zsplit > 0) {
        lConfig.splitConfigs.push_back({Acts::binZ, zsplit});
      }

      // Save
      layerBuilderConfig.layerConfigurations[idx].push_back(lConfig);
    }

    // Save
    layerBuilderConfigs.push_back(layerBuilderConfig);
  }

  // remember the layer builders to collect the detector elements
  std::vector<std::shared_ptr<const Acts::TGeoLayerBuilder>> tgLayerBuilders;

  for (auto& lbc : layerBuilderConfigs) {
    std::shared_ptr<const Acts::LayerCreator> layerCreatorLB = nullptr;

    if (lbc.autoSurfaceBinning) {
      // Configure surface array creator (optionally) per layer builder
      // (in order to configure them to work appropriately)
      Acts::SurfaceArrayCreator::Config sacConfigLB;
      sacConfigLB.surfaceMatcher = lbc.surfaceBinMatcher;
      auto surfaceArrayCreatorLB =
          std::make_shared<const Acts::SurfaceArrayCreator>(
              sacConfigLB, Acts::getDefaultLogger(
                               lbc.configurationName + "SurfaceArrayCreator",
                               surfaceLogLevel));

      // configure the layer creator that uses the surface array creator
      Acts::LayerCreator::Config lcConfigLB;
      lcConfigLB.surfaceArrayCreator = surfaceArrayCreatorLB;
      layerCreatorLB = std::make_shared<const Acts::LayerCreator>(
          lcConfigLB,
          Acts::getDefaultLogger(lbc.configurationName + "LayerCreator",
                                 layerLogLevel));
    }

    // Configure the proto layer helper
    Acts::ProtoLayerHelper::Config plhConfigLB;
    auto protoLayerHelperLB = std::make_shared<const Acts::ProtoLayerHelper>(
        plhConfigLB,
        Acts::getDefaultLogger(lbc.configurationName + "ProtoLayerHelper",
                               layerLogLevel));

    //-------------------------------------------------------------------------------------
    lbc.layerCreator =
        (layerCreatorLB != nullptr) ? layerCreatorLB : layerCreator;
    lbc.protoLayerHelper =
        (protoLayerHelperLB != nullptr) ? protoLayerHelperLB : protoLayerHelper;

    auto layerBuilder = std::make_shared<const Acts::TGeoLayerBuilder>(
        lbc, Acts::getDefaultLogger(lbc.configurationName + "LayerBuilder",
                                    layerLogLevel));
    // remember the layer builder
    tgLayerBuilders.push_back(layerBuilder);

    // build the pixel volume
    Acts::CylinderVolumeBuilder::Config volumeConfig;
    volumeConfig.trackingVolumeHelper = cylinderVolumeHelper;
    volumeConfig.volumeName = lbc.configurationName;
    volumeConfig.buildToRadiusZero = (volumeBuilders.size() == 0);
    volumeConfig.layerEnvelopeR = {1. * Acts::UnitConstants::mm,
                                   5. * Acts::UnitConstants::mm};
    auto ringLayoutConfiguration =
        [&](const std::vector<Acts::TGeoLayerBuilder::LayerConfig>& lConfigs)
        -> void {
      for (const auto& lcfg : lConfigs) {
        for (const auto& scfg : lcfg.splitConfigs) {
          if (scfg.first == Acts::binR and scfg.second > 0.) {
            volumeConfig.ringTolerance =
                std::max(volumeConfig.ringTolerance, scfg.second);
            volumeConfig.checkRingLayout = true;
          }
        }
      }
    };
    ringLayoutConfiguration(lbc.layerConfigurations[0]);
    ringLayoutConfiguration(lbc.layerConfigurations[2]);
    volumeConfig.layerBuilder = layerBuilder;
    auto volumeBuilder = std::make_shared<const Acts::CylinderVolumeBuilder>(
        volumeConfig,
        Acts::getDefaultLogger(lbc.configurationName + "VolumeBuilder",
                               volumeLogLevel));
    // add to the list of builders
    volumeBuilders.push_back(volumeBuilder);
  }

  //-------------------------------------------------------------------------------------
  // create the tracking geometry
  Acts::TrackingGeometryBuilder::Config tgConfig;
  // Add the builders
  tgConfig.materialDecorator = matDeco;

  for (auto& vb : volumeBuilders) {
    tgConfig.trackingVolumeBuilders.push_back(
        [=](const auto& gcontext, const auto& inner, const auto&) {
          return vb->trackingVolume(gcontext, inner);
        });
  }
  // Add the helper
  tgConfig.trackingVolumeHelper = cylinderVolumeHelper;
  auto cylinderGeometryBuilder =
      std::make_shared<const Acts::TrackingGeometryBuilder>(
          tgConfig,
          Acts::getDefaultLogger("TrackerGeometryBuilder", volumeLogLevel));
  // get the geometry
  _trackingGeometry =
      cylinderGeometryBuilder->trackingGeometry(_geometryContext);
  // collect the detector element store
  for (auto& lBuilder : tgLayerBuilders) {
    auto detElements = lBuilder->detectorElements();
    _detectorStore.insert(_detectorStore.begin(), detElements.begin(),
                          detElements.end());
  }

  //
  // Restore old gGeoManager
  if (gGeoManagerOld != nullptr) {
    gGeoManager = gGeoManagerOld;
  }
}

void ACTSProcBase::buildBfield() {
  // Get the magnetic field
  dd4hep::Detector& lcdd = dd4hep::Detector::getInstance();
  const double position[3] = {
      0, 0,
      0};  // position to calculate magnetic field at (the origin in this case)
  double magneticFieldVector[3] = {
      0, 0, 0};  // initialise object to hold magnetic field
  lcdd.field().magneticField(
      position,
      magneticFieldVector);  // get the magnetic field vector from DD4hep

  // Build ACTS representation of field
  // Note:
  //  magneticFieldVector[2] = 3.57e-13
  //  dd4hep::tesla = 1e-13
  //  Acts::UnitConstants::T = 0.000299792
  _magneticField = std::make_shared<Acts::ConstantBField>(Acts::Vector3(
      magneticFieldVector[0] / dd4hep::tesla * Acts::UnitConstants::T,
      magneticFieldVector[1] / dd4hep::tesla * Acts::UnitConstants::T,
      magneticFieldVector[2] / dd4hep::tesla * Acts::UnitConstants::T));
}

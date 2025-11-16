########################################################
# cmake file for building Marlin ACTSTracking Package
# @author Karol Krizka
CMAKE_MINIMUM_REQUIRED(VERSION 3.15 FATAL_ERROR)
########################################################

# project name
PROJECT( ACTSTracking )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 3 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

### BUILD OPTIONS ###########################################################

set(PACKAGE_INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(PACKAGE_INSTALL_INC_DIR include CACHE PATH "Installation directory for headers")

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
INCLUDE_DIRECTORIES( ${Marlin_INCLUDE_DIRS} )
LINK_LIBRARIES( ${Marlin_LIBRARIES} )
ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

FIND_PACKAGE( DD4hep REQUIRED COMPONENTS DDCore )
INCLUDE_DIRECTORIES( SYSTEM ${DD4hep_INCLUDE_DIRS} )
LINK_LIBRARIES( ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES}  )

FIND_PACKAGE( LCIO REQUIRED )
INCLUDE_DIRECTORIES( ${LCIO_INCLUDE_DIRS} )
LINK_LIBRARIES( ${LCIO_LIBRARIES} )

FIND_PACKAGE( ROOT REQUIRED)
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} )

FIND_PACKAGE( Acts REQUIRED COMPONENTS Core PluginJson PluginTGeo )

INCLUDE(GNUInstallDirs)

### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()

### DATA ###########################################################
ADD_SUBDIRECTORY(data)
ADD_SUBDIRECTORY(example)

### LIBRARY #################################################################
ADD_DEFINITIONS( "-Wno-long-long -Wno-effc++ -Wno-unused-parameter -Wno-narrowing" )

# sources
CONFIGURE_FILE(src/config.h.in config.h)
ADD_LIBRARY(ACTSTracking SHARED
  ./src/ACTSProcBase.cxx
  ./src/ACTSTruthTrackingProc.cxx
  ./src/ACTSTruthCKFTrackingProc.cxx
  ./src/ACTSSeededCKFTrackingProc.cxx
  ./src/ACTSDuplicateRemoval.cxx
  ./src/TrackTruthProc.cxx
  ./src/GeometryIdSelector.cxx
  ./src/GeometryIdMappingTool.cxx
  ./src/Helpers.cxx
  )
TARGET_INCLUDE_DIRECTORIES(ACTSTracking
  PRIVATE
    ./ACTSTracking
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Marlin_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS}
    ${DD4hep_INCLUDE_DIRS}
    ${LCIO_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(ACTSTracking
  ${ROOT_LIBRARIES} DD4hep::DDCore
  ActsCore ActsPluginTGeo ActsPluginJson)

INSTALL( TARGETS ACTSTracking DESTINATION lib )
INSTALL( DIRECTORY ACTSTracking DESTINATION include FILES_MATCHING PATTERN "*.hxx" )



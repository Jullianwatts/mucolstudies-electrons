# Install script for directory: /scratch/jwatts/mucol/MyBIBUtils

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/scratch/jwatts/mucol/MyBIBUtils")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/scratch/jwatts/mucol/MyBIBUtils/build/data/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMyBIBUtils.so.0.1.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMyBIBUtils.so.0.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "$ORIGIN/../lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/marlin-1.19.5-la7poeudkvkc2tgvw4o2jfw5djzm6dh3/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/lcio-2.22.5-v3x5uiwfim3nefmeizxn5xooqa4ptgki/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/gear-1.9.5-3bjxgebzqs4bw37h7rncz3nzfmlrfogb/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/clhep-2.4.7.1-pi7m3sbrmz5ex4audlodc5vaqc2o4zg7/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/ilcutil-1.7.3-rbqyt64imuq24jaou6uehkomwkrxombo/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/marlinutil-1.18.2-vegnxpkvhgsqx6ti6jcdcvy3t5hldwmy/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/ced-1.10-7hwb3ptti7nxcks4uyuvx6alnnr4il43/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/dd4hep-master-fzr6nyeeos2hpkhcx4w3gsptbvlsbhkt/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/root-6.32.08-ywhwc7irxiytd63xmlgxbunobvc5yjoe/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/sio-0.2-ypovfe4xkzponiqfb627l6ny5wbgnelp/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/xerces-c-3.3.0-7udplkbcss57oci6gkpochnq4yatplgw/lib")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE FILES
    "/scratch/jwatts/mucol/MyBIBUtils/build/lib/libMyBIBUtils.so.0.1.0"
    "/scratch/jwatts/mucol/MyBIBUtils/build/lib/libMyBIBUtils.so.0.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMyBIBUtils.so.0.1.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMyBIBUtils.so.0.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/marlin-1.19.5-la7poeudkvkc2tgvw4o2jfw5djzm6dh3/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/lcio-2.22.5-v3x5uiwfim3nefmeizxn5xooqa4ptgki/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/gear-1.9.5-3bjxgebzqs4bw37h7rncz3nzfmlrfogb/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/clhep-2.4.7.1-pi7m3sbrmz5ex4audlodc5vaqc2o4zg7/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/ilcutil-1.7.3-rbqyt64imuq24jaou6uehkomwkrxombo/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/marlinutil-1.18.2-vegnxpkvhgsqx6ti6jcdcvy3t5hldwmy/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/ced-1.10-7hwb3ptti7nxcks4uyuvx6alnnr4il43/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/dd4hep-master-fzr6nyeeos2hpkhcx4w3gsptbvlsbhkt/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/root-6.32.08-ywhwc7irxiytd63xmlgxbunobvc5yjoe/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/sio-0.2-ypovfe4xkzponiqfb627l6ny5wbgnelp/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/xerces-c-3.3.0-7udplkbcss57oci6gkpochnq4yatplgw/lib:::::::::::::::"
           NEW_RPATH "$ORIGIN/../lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/marlin-1.19.5-la7poeudkvkc2tgvw4o2jfw5djzm6dh3/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/lcio-2.22.5-v3x5uiwfim3nefmeizxn5xooqa4ptgki/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/gear-1.9.5-3bjxgebzqs4bw37h7rncz3nzfmlrfogb/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/clhep-2.4.7.1-pi7m3sbrmz5ex4audlodc5vaqc2o4zg7/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/ilcutil-1.7.3-rbqyt64imuq24jaou6uehkomwkrxombo/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/marlinutil-1.18.2-vegnxpkvhgsqx6ti6jcdcvy3t5hldwmy/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/ced-1.10-7hwb3ptti7nxcks4uyuvx6alnnr4il43/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/dd4hep-master-fzr6nyeeos2hpkhcx4w3gsptbvlsbhkt/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/root-6.32.08-ywhwc7irxiytd63xmlgxbunobvc5yjoe/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/sio-0.2-ypovfe4xkzponiqfb627l6ny5wbgnelp/lib:/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/xerces-c-3.3.0-7udplkbcss57oci6gkpochnq4yatplgw/lib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE FILES "/scratch/jwatts/mucol/MyBIBUtils/build/lib/libMyBIBUtils.so")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/scratch/jwatts/mucol/MyBIBUtils/build/CMakeFiles/MyBIBUtils.dir/install-cxx-module-bmi-RelWithDebInfo.cmake" OPTIONAL)
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/scratch/jwatts/mucol/MyBIBUtils/build/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/scratch/jwatts/mucol/MyBIBUtils/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()

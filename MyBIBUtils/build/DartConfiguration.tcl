# This file is configured by CMake automatically as DartConfiguration.tcl
# If you choose not to use CMake, this file may be hand configured, by
# filling in the required variables.


# Configuration directories and files
SourceDirectory: /scratch/jwatts/mucol/MyBIBUtils
BuildDirectory: /scratch/jwatts/mucol/MyBIBUtils/build

# Where to place the cost data store
CostDataFile: 

# Site is something like machine.domain, i.e. pragmatic.crd
Site: ap23.uc.osg-htc.org

# Build name is osname-revision-compiler, i.e. Linux-2.4.2-2smp-c++
BuildName: Linux-g++

# Subprojects
LabelsForSubprojects: 

# Submission information
SubmitURL: http://
SubmitInactivityTimeout: 

# Dashboard start time
NightlyStartTime: 00:00:00 EDT

# Commands for the build/test/submit cycle
ConfigureCommand: "/opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/cmake-3.31.4-u7bitlxdtcr7tlb64k22k5x5dbbyt5od/bin/cmake" "/scratch/jwatts/mucol/MyBIBUtils"
MakeCommand: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/cmake-3.31.4-u7bitlxdtcr7tlb64k22k5x5dbbyt5od/bin/cmake --build . --config "${CTEST_CONFIGURATION_TYPE}" -- -i
DefaultCTestConfigurationType: Release

# version control
UpdateVersionOnly: 

# CVS options
# Default is "-d -P -A"
CVSCommand: 
CVSUpdateOptions: 

# Subversion options
SVNCommand: 
SVNOptions: 
SVNUpdateOptions: 

# Git options
GITCommand: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/git-2.47.0-v5r5iyrovrcxzip5xp4ugwtkzzmqu2uc/bin/git
GITInitSubmodules: 
GITUpdateOptions: 
GITUpdateCustom: 

# Perforce options
P4Command: 
P4Client: 
P4Options: 
P4UpdateOptions: 
P4UpdateCustom: 

# Generic update command
UpdateCommand: /opt/spack/opt/spack/linux-almalinux9-x86_64/gcc-11.5.0/git-2.47.0-v5r5iyrovrcxzip5xp4ugwtkzzmqu2uc/bin/git
UpdateOptions: 
UpdateType: git

# Compiler info
Compiler: /usr/bin/g++
CompilerVersion: 11.5.0

# Dynamic analysis (MemCheck)
PurifyCommand: 
ValgrindCommand: 
ValgrindCommandOptions: 
DrMemoryCommand: 
DrMemoryCommandOptions: 
CudaSanitizerCommand: 
CudaSanitizerCommandOptions: 
MemoryCheckType: 
MemoryCheckSanitizerOptions: 
MemoryCheckCommand: MEMORYCHECK_COMMAND-NOTFOUND
MemoryCheckCommandOptions: 
MemoryCheckSuppressionFile: 

# Coverage
CoverageCommand: /usr/bin/gcov
CoverageExtraFlags: -l

# Testing options
# TimeOut is the amount of time in seconds to wait for processes
# to complete during testing.  After TimeOut seconds, the
# process will be summarily terminated.
# Currently set to 25 minutes
TimeOut: 1500

# During parallel testing CTest will not start a new test if doing
# so would cause the system load to exceed this value.
TestLoad: 

TLSVerify: 
TLSVersion: 

UseLaunchers: 
CurlOptions: 
# warning, if you add new options here that have to do with submit,
# you have to update cmCTestSubmitCommand.cxx

# For CTest submissions that timeout, these options
# specify behavior for retrying the submission
CTestSubmitRetryDelay: 5
CTestSubmitRetryCount: 3

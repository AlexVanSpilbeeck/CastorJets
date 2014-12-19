#!/bin/bash          

export ROOTSYS=/user/cmssoft/root_5.32.00/root
export PATH=$PATH:$ROOTSYS/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/lib


source $VO_CMS_SW_DIR/cmsset_default.sh 
cd /localgrid/avanspil/Castor_Analysis/CMSSW_4_2_10_patch2/src/UACastor/CastorTree/Analysis/.
eval `scram runtime -sh`

cd /scratch
if [ ! -d CMSSW_4_2_10_patch2 ]
then
	cp -r /user/avanspil/Castor_Analysis/CMSSW_4_2_10_patch2 .
fi
cd CMSSW_4_2_10_patch2/src/UACastor/CastorTree/Analysis

./Run Pythia6Z2star_new 7000 JetAnalyzer 

cp *root /user/avanspil/Castor_Analysis/CMSSW_4_2_10_patch2/src/UACastor/CastorTree/Analysis/.

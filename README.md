## 2018 HLT rate, timing and efficiency and integration studies. 

Commands for each step ----

## 1. HLT rate:

hltGetConfiguration /users/oozcelik/Displaced-Jpsi/HLT --globaltag 100X_dataRun2_relval_ForTSG_v1 --input root://eoscms.cern.ch//eos/cms/tier0/store/data/Run2017E/EphemeralHLTPhysics1/RAW/v1/000/304/777/00000/00175E91-D0AD-E711-A24F-02163E01451E.root - -data --process MYHLT --offline --customise HLTrigger/Configuration/customizeHLTforCMSSW.customiseFor2017DtUnpacking --setup /dev/CMSSW_10_0_0/HLT --prescale none --max-events 10 --output none > hlt.py

edmConfigDump hlt.py > hlt_config.py

-- Test HLT config file on 5 events
cmsRun hlt_config.py &> log_hlt_config

-- Edit run_steamflow_cfg.py and set this variable to False:
-- switchL1PS=False

-- create the list of input files with all the EphemeralHLTPhysics datasets for run 304777 (Run2017E)
-- => edit run_steamflow_cfg.py and set nEvents=-1

./cmsBatch.py 10 run_steamflow_cfg.py -o OUT -q 8nh -r /eos/cms/store/user/oozcelik/HLTrigger/OUT/hlt.root

cd $CMSSW_BASE/src/SteamRatesEdmWorkflow/Rates/

ls /eos/cms/store/user/oozcelik/HLTrigger/OUT/hlt_*.root >> filesInput.py

cp /afs/cern.ch/user/n/ndaci/public/STEAM/Production/Xudong_HLTv4/json_HLTPhysicsL1v4_2p0e34.txt .

-- "Edit config_makeBatchJobs.py and set:\nEDM file folder, CMSSW base directory, json file and LSF settings"

python config_makeBatchJobs.py

./sub_total.jobb

python config_mergeOutputs.py

## 2. HLT Timing :
(using Ephemeral 2017E dataset)

hltGetConfiguration /users/oozcelik/Displaced-Jpsi/HLT/V4 --globaltag 100X_dataRun2_relval_ForTSG_v1 --input root://eoscms.cern.ch//eos/cms/tier0/store/data/Run2017E/EphemeralHLTPhysics1/RAW/v1/000/304/777/00000/00175E91-D0AD-E711-A24F-02163E01451E.root --data --process TIMING --offline --l1 L1Menu_Collisions2017_v4_m6 --unprescale --max-events 100 --setup /dev/CMSSW_10_0_0/GRun --no-output --timing > hlt.py

cmsRun hlt.py >& full.log
cmsRun harvesting.py >& harvesting.log
root -l DQM_V0001_R000304777__HLT__FastTimerService__All.root

and go to QMData/Run 304777/HLT/Run summary/TimerService/process TIMING paths/path HLT_Mu3_PFJet200CSV_1p5_v15/ directory, you can explore, e.g, the path time real plot and the module timing plot of that path

## 3. HLT integration test

-- to get the updated L1T menu circulated July 24th (L1Menu_Collisions2017_dev_r9) --
(it might be run w/o L1 menu, not so crucial.)
git clone https://github.com/cms-l1-dpg/2017-pp-menu-dev -b 2017-07-24 ../2017-pp-menu-dev
mkdir -p L1Trigger/L1TGlobal/data/Luminosity/startup
cp ../2017-pp-menu-dev/Apr12/*.xml L1Trigger/L1TGlobal/data/Luminosity/startup/

hltIntegrationTests -s /dev/CMSSW_10_0_0/GRun -i root://xrootd-cms.infn.it//store/mc/RunIISummer17DRStdmix/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/00916118-3AA9-E711-9619-008CFAC93DC0.root --mc -n 1000 -x "--l1Xml L1Menu_Collisions2017_dev_r9.xml --globaltag 94X_mc2017_realistic_TSG_2017_12_19_13_49_40 --unprescale --customise HLTrigger/Configuration/customizeHLTforCMSSW.customiseFor2017DtUnpacking"  



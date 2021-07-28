import ROOT
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
from DataFormats.FWLite import Events, Handle, Runs
from math import *
import pandas as pd
import numpy as np
import re
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

#file_list = 'ken_ttz_v1_files.txt'
file_list  = 'ttbb_eft.txt'
file_list = '/cms/data/store/user/bcaraway/miniAOD_eft/2018/TTJets_13TeV_TopEFT_MINIAOD_v01_2018-v1.txt'
#file_list = '/cms/data/store/user/bcaraway/miniAOD_eft/2018/TTBBJet_13TeV_TopEFT_MINIAOD_2018-v1.txt '
#rfile     = 'TOP-RunIIFall18wmLHEG_10.root'
#rfile      = 'test_4_ken_ttbb_1_27.root'
rfile      = 'ttzjet_cpq3_m7p0.root'
rfile      = 'ttz_eft_sm.root'
#rfile = '/cms/data/store/user/hatake/TopEFT/TTBB_13TeV_TopEFT_v12_2018/TTBB_13TeV_TopEFT_MINIAOD_v12_2018-v1/step_miniaodsim_1.root'

#rfile  = '/cms/data/store/user/hatake/TopEFT/TTHJet_13TeV_TopEFT_v11hs/TTHJet_13TeV_TopEFT_MINIAOD_v11hs-v3/step_miniaodsim_5321.root'
#rfile  = '/cms/data/store/user/hatake/TopEFT/TTBB_13TeV_TopEFT_v12_2018/TTBB_13TeV_TopEFT_MINIAOD_v12_2018-v1/step_miniaodsim_1622.root'

#rfile = '/cms/data/store/user/hatake/TopEFT/TTHJet_13TeV_TopEFT_v11_2017/TTHJet_13TeV_TopEFT_2017_MINIAOD_v11_2017-v1/step_miniaodsim_5610.root' # 2017 ttH
#rfile = '/cms/data/store/user/hatake/TopEFT/TTHJet_13TeV_TopEFT_v11hs/TTHJet_13TeV_TopEFT_MINIAOD_v11hs-v3/step_miniaodsim_9958.root' # 2018 ttH
rfile = '/cms/data/store/user/hatake/TopEFT/TTHJet_13TeV_TopEFT_v11_2016/TTHJet_13TeV_TopEFT_MINIAOD_v11_2016-v2/step_miniaodsim_1000.root'

#rfile = '/cms/data/store/user/hatake/TopEFT/TTZJet_13TeV_TopEFT_v20_2017/TTZJet_13TeV_TopEFT_MINIAOD_v20_2017-v1/step_miniaodsim_8844.root' # 2017 ttZ
#rfile = '/cms/data/store/user/hatake/TopEFT/TTZJet_13TeV_TopEFT_v20/TTZJet_13TeV_TopEFT_MINIAOD_v20-v1/step_miniaodsim_407.root' # 2018 ttZ


f = open(file_list)
#tth_list = [line.strip('\n') for line in f.readlines()]
tth_list = f.readlines()[0].strip('\n')
#events = Events (["/home/bcaraway/CombineTool/CMSSW_10_2_9/src/EFT_Analyzer/HIG-RunIIFall17MiniAOD-00821ND_298299.root"])
#rfile = "/home/bcaraway/CombineTool/CMSSW_10_2_9/src/EFT_Analyzer/step_miniaodsim_TTHJet_v11_test.root" 
#rfile = "/cms/data/store/user/hatake/TopEFT/TTHJet_13TeV_TopEFT_v11_2017/TTHJet_13TeV_TopEFT_2017_MINIAOD_v11_2017-v1/step_miniaodsim_11.root" 
#rfile = "/cms/data/store/user/hatake/TopEFT/TTHJet_13TeV_TopEFT_v11hs/TTHJet_13TeV_TopEFT_MINIAOD_v11hs-v3/step_miniaodsim_111.root"
#rfile = "/cms/data/store/user/hatake/TopEFT/TTZJet_13TeV_TopEFT_v18/TTZJet_13TeV_TopEFT_MINIAOD_v18-v1/step_miniaodsim_1111.root"
#rfile = "/cms/data/store/user/hatake/TopEFT/TTZJet_13TeV_TopEFT_v18_2016/TTZJet_13TeV_TopEFT_MINIAOD_v18_2016-v1/step_miniaodsim_1111.root"
#rfile  = "/cms/data/store/user/hatake/TopEFT/TTHJet_13TeV_TopEFT_v11_2016/TTHJet_13TeV_TopEFT_MINIAOD_v11_2016-v1/step_miniaodsim_111.root"
#rfile  = "/cms/data/store/user/hatake/TopEFT/TTHJet_13TeV_TopEFT_v11_2016/TTHJet_13TeV_TopEFT_MINIAOD_v11_2016-v2/step_miniaodsim_10001.root"
#rfile = "/cms/data/hatake/ana/TTX/Prod_sl6/2017/step_miniaodsim.test.root"
#rfile = "/home/bcaraway/CombineTool/CMSSW_10_2_9/src/EFT_Analyzer/step_aodsim.root" 
events = Events ([rfile])
#runs = Runs(     [rfile])
#events= Events(tth_list)
'''
Notes on file from Andrew W:

Each event is a ttH sample, with multi-lepton final state (verify final state) H->WW->stuff
Additionally, each event has 184 EFT weights corresponding to varying 16 Wilson Coefficients
'''

handlePruned  = Handle ("std::vector<reco::GenParticle>")
handlePacked  = Handle ("std::vector<pat::PackedGenParticle>")
labelPruned = ("prunedGenParticles")
#labelPacked = ("packedGenParticles")
labelGen   = ("genParticles")

handleLHE  = Handle("LHEEventProduct")
handleRLHE = Handle("LHERunInfoProduct")
#labelLHE   = ("externalLHEProducer") # GEN
labelLHE = ("LHEFile") # miniAOD

# loop over events

#runs.getByLabel("source","",labelLHE, handleRLHE)
#runs.getByLabel(labelLHE, handleRLHE)
#lher = handleRLHE.product()
#for i in lher:
#    print(i)


aux_df = pd.DataFrame()
daughter_list = []

for event in events:
    #event.getByLabel( labelLHE,"","",handleLHE)       # GEN
    event.getByLabel ("source","",labelLHE, handleLHE) # miniAOD
    #event.getByLabel ("source","",labelLHE, handleRLHE)
    lhe  = handleLHE.product()
    #help(lhe)
    #lher = handleRLHE.product()
    for w in lhe.weights():
        pass
        #if 'rwgt_' in w.id:
        print(w.id, w.wgt)

    aux_df = pd.DataFrame(
        data=   [[float(wc.split('_')[1]) for wc in re.findall(r'c\w+_-?\d+\.\d+', w.id)] for w in lhe.weights() if 'EFT' in w.id],
        columns=[str(wc.split('_')[0])    for wc in re.findall(r'c\w+_-?\d+\.\d+', lhe.weights()[180].id) if 'EFT' in lhe.weights()[180].id],  # have to hard-code to read a EFT weight (180)
        index=[re.findall(r'EFTrwgt\d+',w.id)[0] for w in lhe.weights() if 'EFT' in w.id]
    )
    #sm_df = pd.DataFrame(
    #    data = [np.zeros(16)],
    #    columns=[str(wc.split('_')[0])    for wc in re.findall(r'c\w+_-?\d+\.\d+', lhe.weights()[180].id) if 'EFT' in lhe.weights()[180].id],
    #    index=['EFTrwgt183']
    #)
    #aux_df = aux_df.append(sm_df)
    break
    print('\n')

#print(file_list)
print(aux_df)


data = []
for count, event in enumerate(events):
    if count % 1000 == 0:
        print('Event #:{0:6} of {1:7}'.format(count,events.size()))
    #event.getByLabel (labelPacked, handlePacked)

    event.getByLabel (labelPruned, handlePruned)
    event.getByLabel ("source","",labelLHE, handleLHE) # miniAOD
    #event.getByLabel(labelGen, handlePruned)
    #event.getByLabel (labelLHE, handleLHE)
    # get the product
    #packed = handlePacked.product()
    pruned = handlePruned.product()
    lhe = handleLHE.product()
    
    
    hasHiggs = False
    n_daughters = {'n':0, 'n_d':0}
    def getAllDaughters(h, iters):
        #if iters == 20: return
        for i in xrange(0, h.numberOfDaughters() ):
            #if h.daughter(i).pdgId() not in daughter_list:
            if ((h.pdgId() == 25) and (h.status() == 62) and (n_daughters['n'] != 2)):
                #if h.daughter(i).pdgId() == 22 or h.daughter(i).pdgId() == 23: print(h.daughter(i).pdgId())
                daughter_list.append(h.daughter(i).pdgId())
                n_daughters['n'] += 1
            if h.daughter(i).pdgId() == 25:
                getAllDaughters(h.daughter(i), iters+1)
    for p in pruned : # pruned
        break
        if p.pdgId() == 25:
            hasHiggs = True
            data.append(
                [w.wgt for w in lhe.weights() if 'EFT' in w.id or 'SM' in w.id]+[p.pt(),p.eta(),p.phi(),p.mass()]
            )
            #print "PdgId : %s   pt : %s  eta : %s   phi : %s" %(p.pdgId(),p.pt(),p.eta(),p.phi())    
            getAllDaughters(p,0)
            break    
            #df = pd.concat([df,temp_df], ignore_index=True)
    #if n_daughters['n'] != 2 : print(n_daughters['n'])
    data.append([w.wgt for w in lhe.weights() if 'EFT' in w.id or 'SM' in w.id])

#daughter_list , d_count =np.unique(daughter_list, return_counts=True)
#daughter_df = pd.Series(d_count, index=daughter_list)
#print(daughter_df.sort_values(ascending=False))
df = pd.DataFrame(
    data = data,
    #columns = np.append(aux_df.index.values,['H_pt','H_eta','H_phi','H_mass'])
    columns = aux_df.index.values
)
print(df)
#import matplotlib.pyplot as plt
#plt.hist(df['EFTrwgt183'].clip(0,50), bins = 20)

#aux_df.to_pickle('aux_EFT_ken_tthv2.pkl')    
#df.to_pickle('eventInfo_EFT_ken_tthv2.pkl')

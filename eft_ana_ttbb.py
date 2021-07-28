import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
from DataFormats.FWLite import Events, Handle, Runs
from math import *
import pandas as pd
import numpy as np
import re
import ROOT
from ROOT import TLorentzVector as tlv
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

file_list  = 'ttbb_2018_test.txt'
file_list  = '/cms/data/store/user/bcaraway/miniAOD_eft/2018/TTBBJet_13TeV_TopEFT_MINIAOD_2018-v12.txt'
#rfile = '/cms/data/store/user/hatake/TopEFT/TTBB/GenOnly/TTBB_4f_cbW/TOP-RunIIFall18wmLHEG.root'
rfile = 'ttbb_2018_test.root'

#tth_list = [line.strip('\n') for line in f.readlines()]
#rfile_list = open(file_list).readlines()[0].strip('\n')
#events = Events ([rfile])
events= Events([line.strip('\n') for line in open(file_list).readlines()])

#runs = Runs(     [rfile])
'''
Notes on file from Andrew W:

Each event is a ttH sample, with multi-lepton final state (verify final state) H->WW->stuff
Additionally, each event has 184 EFT weights corresponding to varying 16 Wilson Coefficients
'''

handlePruned  = Handle ("std::vector<reco::GenParticle>")
handlePacked  = Handle ("std::vector<pat::PackedGenParticle>")
labelGen = ("prunedGenParticles") # miniAOD
#labelPacked = ("packedGenParticles")
#labelGen   = ("genParticles") # GEN

handleLHE  = Handle("LHEEventProduct")
###handleRLHE = Handle("LHERunInfoProduct")
#labelLHE   = ("externalLHEProducer")
labelLHE = ("LHEFile")

# -- miniAOD
#Type                                  Module                      Label             Process     
#------------------------------------------------------------------------------------------------
#LHEEventProduct                       "source"                    ""                "LHEFile"

# loop over events


aux_df = pd.DataFrame()
daughter_list = []

for event in events:
    #event.getByLabel( labelLHE,"","",handleLHE)       # GEN
    event.getByLabel ("source","",labelLHE, handleLHE) # miniAOD
    lhe  = handleLHE.product()
    
    isold = False
    for w in lhe.weights():
        if 'rwgt' in w.id:
            print(w.id, w.wgt)
            if 'EFT' in w.id: 
                isold=True
                break

    #exit()
    if isold:
        aux_df = pd.DataFrame(
            data=   [[float(wc.split('_')[1]) for wc in re.findall(r'c\w+_-?\d+\.\d+', w.id)] for w in lhe.weights() if 'EFT' in w.id],
            columns=[str(wc.split('_')[0])    for wc in re.findall(r'c\w+_-?\d+\.\d+', lhe.weights()[180].id) if 'EFT' in lhe.weights()[180].id],  # have to hard-code to read a EFT weight (180)
            index=[re.findall(r'EFTrwgt\d+',w.id)[0] for w in lhe.weights() if 'EFT' in w.id]
        )
    else:
        aux_df = pd.DataFrame(
            data=   [[float(wc.split('_')[1]) for wc in re.findall(r'_-?\d+\.?\d*', w.id.replace("min",'-').replace('p','.'))] for w in lhe.weights() if 'rwgt' in w.id],
            columns=[str(wc.split('_')[0])    for wc in re.findall(r'(c\w+|SM)_', lhe.weights()[150].id) if 'rwgt' in lhe.weights()[150].id],  # have to hard-code to read a EFT weight (180)
            index=[re.findall(r'rwgt\w+',w.id)[0] for w in lhe.weights() if 'rwgt' in w.id]
    )
        aux_df.loc['rwgt_SM','cbW'] = 0
    print('\n')
    break


#print(file_list)
print(aux_df)
data = []
for count, event in enumerate(events):
    if count % 1000 == 0:
        print('Event #:{0:6} of {1:7}'.format(count,events.size()))

    event.getByLabel(labelGen, handlePruned)
    #event.getByLabel (labelLHE, handleLHE) # GEN
    event.getByLabel ("source","",labelLHE, handleLHE) # miniAOD
    # get the product
    pruned = handlePruned.product()
    lhe = handleLHE.product()
    
    _bb = []
    bb_count = 0
    for p in pruned:
        if abs(p.pdgId()) == 5 and p.isHardProcess():
            if abs(p.mother().pdgId()) != 6 and abs(p.mother().pdgId()) != 24:
                b =tlv()
                b.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy())
                _bb.append(b)
                bb_count += 1

    if bb_count == 3:
        for p in pruned:
            if p.isHardProcess():
                try:
                    print(p.pdgId(), p.mother().pdgId())
                except:
                    pass
        exit()
        
    extra_bb = _bb[0]+_bb[1]
    data.append([w.wgt for w in lhe.weights() if 'rwgt' in w.id]+[bb_count,extra_bb.Pt(),extra_bb.M(),_bb[0].DeltaR(_bb[1])])
    if count == 200000:
        break

# end of event loop
df = pd.DataFrame(data = data, columns = np.append(aux_df.index.values, ['bb_count','bb_pt', 'bb_m','bb_dr']))
print(df)
print(np.unique(df['bb_count'],return_counts=True))


#daughter_list , d_count =np.unique(daughter_list, return_counts=True)
#daughter_df = pd.Series(d_count, index=daughter_list)
#print(daughter_df.sort_values(ascending=False))
#df = pd.DataFrame(
#    data = data,
#    columns = np.append(aux_df.index.values,['H_pt','H_eta','H_phi','H_mass'])
#)
#print(df)

#aux_df.to_pickle('aux_EFT_ken_ttbb.pkl')    
df.to_pickle('eventInfo_EFT_ken_ttbb_2018_v12.pkl')

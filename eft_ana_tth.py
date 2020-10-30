import ROOT
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
from DataFormats.FWLite import Events, Handle
from math import *
import pandas as pd
import numpy as np
import re
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()


file_list = 'tth_file_list.txt'
f = open(file_list)
tth_list = [line.strip('\n') for line in f.readlines()]
#events = Events (["/home/bcaraway/CombineTool/CMSSW_10_2_9/src/EFT_Analyzer/HIG-RunIIFall17MiniAOD-00821ND_298299.root"])
events= Events(tth_list)
'''
Notes on file from Andrew W:

Each event is a ttH sample, with multi-lepton final state (verify final state) H->WW->stuff
Additionally, each event has 184 EFT weights corresponding to varying 16 Wilson Coefficients
'''

handlePruned  = Handle ("std::vector<reco::GenParticle>")
handlePacked  = Handle ("std::vector<pat::PackedGenParticle>")
labelPruned = ("prunedGenParticles")
labelPacked = ("packedGenParticles")
handleLHE  = Handle("LHEEventProduct")
labelLHE   = ("externalLHEProducer")
# loop over events

aux_df = pd.DataFrame()
daughter_list = []

for event in events:
    event.getByLabel (labelLHE, handleLHE)
    lhe = handleLHE.product()
    for w in lhe.weights():
        print(w.id, w.wgt)
    #exit()
    aux_df = pd.DataFrame(
        data=   [[float(wc.split('_')[1]) for wc in re.findall(r'c\w+_-?\d+\.\d+', w.id)] for w in lhe.weights() if 'EFT' in w.id],
        columns=[str(wc.split('_')[0])    for wc in re.findall(r'c\w+_-?\d+\.\d+', lhe.weights()[1].id) if 'EFT' in lhe.weights()[1].id], 
        index=[re.findall(r'EFTrwgt\d+',w.id)[0] for w in lhe.weights() if 'EFT' in w.id]
    )
    break


data = []
for count, event in enumerate(events):
    if count % 1000 == 0:
        print('Event #:{0:6} of {1:7}'.format(count,events.size()))
    event.getByLabel (labelPacked, handlePacked)
    event.getByLabel (labelPruned, handlePruned)
    event.getByLabel (labelLHE, handleLHE)
    # get the product
    packed = handlePacked.product()
    pruned = handlePruned.product()
    
    lhe = handleLHE.product()
    
    
    hasHiggs = False
    #n_daughters = {'n':0, 'n_d':0}
    #def getAllDaughters(h, iters):
    #    #if iters == 20: return
    #    for i in xrange(0, h.numberOfDaughters() ):
    #        #if h.daughter(i).pdgId() not in daughter_list:
    #        if ((h.pdgId() == 25) and (h.status() == 62) and (n_daughters['n'] != 2)):
    #            if h.daughter(i).pdgId() == 22 or h.daughter(i).pdgId() == 23: print(h.daughter(i).pdgId())
    #            daughter_list.append(h.daughter(i).pdgId())
    #            n_daughters['n'] += 1
    #        if h.daughter(i).pdgId() == 25:
    #            getAllDaughters(h.daughter(i), iters+1)

    for p in pruned : # pruned
        if p.pdgId() == 25:
            hasHiggs = True
            data.append(
                [w.wgt for w in lhe.weights() if 'EFT' in w.id]+[p.pt(),p.eta(),p.phi(),p.mass()]
            )
            #print "PdgId : %s   pt : %s  eta : %s   phi : %s" %(p.pdgId(),p.pt(),p.eta(),p.phi())    
            #getAllDaughters(p,0)
            break    
            #df = pd.concat([df,temp_df], ignore_index=True)
    #if n_daughters['n'] != 2 : print(n_daughters['n'])

daughter_list , d_count =np.unique(daughter_list, return_counts=True)
daughter_df = pd.Series(d_count, index=daughter_list)
print(daughter_df.sort_values(ascending=False))
df = pd.DataFrame(
    data = data,
    columns = np.append(aux_df.index.values,['H_pt','H_eta','H_phi','H_mass'])
)
print(df)

#aux_df.to_pickle('aux_EFT.pkl')    
#df.to_pickle('eventInfo_EFT_tth.pkl')


    

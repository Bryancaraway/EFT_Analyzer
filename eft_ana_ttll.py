import ROOT
from ROOT import TLorentzVector
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

file_list = 'ttll_file_list.txt'
f = open(file_list)
ttll_list = [line.strip('\n') for line in f.readlines()]
#events = Events (["/home/bcaraway/CombineTool/CMSSW_10_2_9/src/EFT_Analyzer/ttllND_67924.root"])
events = Events(ttll_list)

'''
Notes on file from Andrew W:

weird tt+ll sample 21->ll what?! thats weird
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
mother_list = []

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
print(aux_df)
exit()
#
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
    #temp_df = pd.DataFrame(
    #    data=   [[w.wgt for w in lhe.weights() if 'EFT' in w.id]],
    #    columns=[re.findall(r'EFTrwgt\d+',w.id)[0] for w in lhe.weights() if 'EFT' in w.id]
    #) 
    
    hasHiggs = False
    n_daughters = {'n':0, 'n_d':0}
    def getAllDaughters(h, iters):
        #if iters == 20: return
        for i in xrange(0, h.numberOfDaughters() ):
            #if h.daughter(i).pdgId() not in daughter_list:
            if ((h.pdgId() == 25) and (h.status() == 62) and (n_daughters['n'] != 2)):
                if h.daughter(i).pdgId() == 22 or h.daughter(i).pdgId() == 23: print(h.daughter(i).pdgId())
                daughter_list.append(h.daughter(i).pdgId())
                n_daughters['n'] += 1
            if h.daughter(i).pdgId() == 25:
                getAllDaughters(h.daughter(i), iters+1)

    for i,p1 in enumerate(pruned) : # pruned
        if (abs(p1.pdgId()) <= 16 and (abs(p1.pdgId()) >= 11)):
            foundZ = False
            for j, p2 in enumerate(pruned) :
                if (p1.pdgId() == -1*p2.pdgId()) and p1.mother() and p2.mother(): 
                    if (p1.mother() == p2.mother()) and p1.mother().status() == 62:
                        tlv_p1 = TLorentzVector(p1.px(), p1.py(), p1.pz(), p1.energy() )
                        tlv_p2 = TLorentzVector(p2.px(), p2.py(), p2.pz(), p2.energy() )
                        tlv_mom = tlv_p1+tlv_p2
                        if abs(p1.mother().pdgId()) == 21: 
                            #print(p1.mother() , p2.mother())
                            pass
                            #temp = pd.DataFrame(np.array([[p.pdgId() for p in pruned], [p.status() for p in pruned]]).T, columns=['pdgId','status'],
                            #                    index=[int(p.mother().pdgId())  if p.mother() else np.nan for p in pruned])
                            #print(temp.head(50))
                            #print(frmt.format(*[p.pdgId() for p in pruned]))
                            #print(frmt.format(*[p.mother().pdgId() for p in pruned]))
                            
                            #exit()
                        if abs(tlv_mom.M() - 91.2) <= 10.0:
                            data.append([w.wgt for w in lhe.weights() if 'EFT' in w.id]+
                                        [tlv_mom.Pt(),tlv_mom.M()])
                                        #[tlv_mom.Pt(),tlv_mom.M(), abs(p1.pdgId()), p1.mother().pdgId(), p1.mother().status()])
                            mother_list.append(p1.mother().pdgId())
                            foundZ = True
                            break
            if foundZ: break
        
    #if n_daughters['n'] != 2 : print(n_daughters['n'])

mother_list , m_count = np.unique(mother_list, return_counts=True)
mother_df = pd.Series(m_count, index=mother_list)
print(mother_df.sort_values(ascending=False))
df = pd.DataFrame(
    data = data,
    #columns = np.append(aux_df.index.values,['ll_pt','ll_mass', 'll_id','ll_mom', 'll_mom_status'])
    columns = np.append(aux_df.index.values,['Z_pt','Z_mass'])
)
df.to_pickle('eventInfo_EFT.pkl')
import matplotlib.pyplot as plt
#plt.hist(np.clip(df['ll_pt'].values, 0, 1000), bins = 50, range=(0,1000))
#plt.show()
plt.hist(np.clip(df['Z_mass'].values, 0, 500), bins = 50, range=(75,125))
plt.show()

#aux_df.to_pickle('aux_EFT.pkl')    



    

import sys
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt

def main():
    aux_df = pd.read_pickle('aux_EFT.pkl')    
    df     = pd.read_pickle('eventInfo_EFT.pkl')
    
    bins = np.arange(0,1050,50)
    h_pt = np.clip(df['H_pt'].to_numpy(), bins[0], bins[-1])
    h_pt = np.repeat(h_pt[:,np.newaxis],3,axis=1)
    
    #print(df)
    #print(aux_df)
    #print(df)
    
    sm     = aux_df.index.values[-1]
    
    def plotWCvsWgt(wc_str):
        f_x = aux_df.sort_values(by=wc_str)[wc_str]
        #
        row1 = df.filter(regex=r'EFTrwgt\d+').iloc[0]
        x    = row1/row1['EFTrwgt183']
        #
        fx_df = pd.concat([f_x,x],axis='columns')
        plt.scatter(y=fx_df.iloc[:,1], x=fx_df.loc[:,wc_str])
        plt.axhline(1, color='r')
        plt.ylim(0)
        plt.grid(True)
        plt.title(wc_str)
        plt.xlabel(f'WC {wc_str} ')
        plt.ylabel('EFT/SM (weight)' )
        plt.show()
    #
    
    #plotWCvsWgt('ctZ')    
    
    #
    
    def giveHiLoWC(wc_str):
        hi = aux_df.sort_values(by=wc_str , ascending=False).index.values[0]
        lo = aux_df.sort_values(by=wc_str , ascending=True) .index.values[0]
        return hi, lo
    
    def plot_EFT_scale(wc_str):
        
        hi, lo = giveHiLoWC(wc_str)
        
        #
        w_hi, w_lo, w_sm = df[hi], df[lo], df[sm]
        
        plt.hist(h_pt,
                 histtype='step',
                 bins=bins,
                 weights=[w_hi,w_lo,w_sm], 
                 label=[(lambda x : '{0} = {1:5.3f}'.format(wc_str,aux_df.loc[x,wc_str]))(y) for y in [hi,lo,sm]])#str(aux_df.loc[hi,wc_str]),str(aux_df.loc[lo,wc_str]),str(aux_df.loc[sm,wc_str])])
        plt.title(wc_str)
        plt.xlabel('GEN H pt (GeV)')
    
        plt.grid(True)
        plt.xlim(bins[0],bins[-1])
    
        plt.legend()
    
        plt.savefig('pdf/{}.pdf'.format(wc_str))
        #plt.show()
        plt.clf()
    
    #
    eft_wc    = ( # all 16 considered WC 
        'ctW', 'ctp', 'cpQM', 'ctli', 'cQei', 'ctZ', 'cQlMi', 'cQl3i', 'ctG', 'ctlTi', 'cbW', 'cpQ3', 'ctei', 'cpt', 'ctlSi', 'cptb'
    )
    eft_hb_wc = ( # two heavy + boson Wilson Coefficients
        'ctp', 'cpQM', 'cpQ3', 'cpt', 'cptb', 'ctW', 'ctZ', 'cbW', 'ctG'
    )
    eft_2l_wc = ( # two heavy + two lepton Wilson Coefficients
        'cQl3i', 'ctli', 'cQei', 'cQlMi', 'ctlTi', 'ctei', 'ctlSi'
    )
    
    for wc in eft_hb_wc:
        plot_EFT_scale(wc)
    
    
    def computeBeta():
        # taken from Jon's code
        # Build the experiment matrix 
        x =[np.ones(184)]
        count = 0
        for i in range(16):
            for j in range(i+1):
                count+=1
                print(count, i , j)
                x.append(aux_df.iloc[:,i].values * aux_df.iloc[:,j].values)
            count+=1
            print(count,i)
            #if (i == 0 ): print( count, aux_df.iloc[:,i])
            x.append(aux_df.iloc[:,i].values)
        x = np.matrix(x).T
        # Build the result matrix y
        y = np.asmatrix(df.iloc[:,:-3].to_numpy()).T
        # Compute beta matrix
        beta = ((x.T * x).I * x.T * y).A
        return beta, x, y
    
    beta, x, y = computeBeta()
    print(list(map(lambda _: _.shape, [beta, x, y])))
    
    print(aux_df.columns[0])
    
    w_c0 = lambda c : beta[1]*c*c + beta[2]*c + beta[0]
    
    plt.hist(h_pt,
             histtype='step',
             bins=bins,
             weights=[w_c0(1)/np.sum(w_c0(1)),w_c0(-1)/np.sum(w_c0(-1)),beta[0]/np.sum(beta[0])],
             density=True,
             label=[f'{aux_df.columns.values[0]} = 1',f'{aux_df.columns.values[0]} = -1', 'SM'])
    plt.yscale('log')
    plt.xlim(0)
    plt.legend()
    plt.show()

#def get_root_hist(roofile):
#    from ROOT import TFile, TH1F
#    f = TFile(roofile,'READ')
#    hist_content = []
#    for i in f.GetListOfKeys():
#        hist = i.ReadObj()
#        hist_content.append([hist.GetBinContent(j+1) for j in range(hist.GetNbinsX()+1)])
#    hist_content = np.array(hist_content)
#    hist_content = np.sum(hist_content, axis=0)
#    hist_binned = [np.sum(hist_content[0:4]),np.sum(hist_content[4:6]), np.sum(hist_content[6:8]), np.sum(hist_content[8:])]
#    return hist_binned/np.sum(hist_binned)

def get_root_hist(roofile):
    import uproot
    roo = uproot.open(roofile)
    hist_content = []
    for hist in roo:
        hist_content.append(roo[hist].values)
    hist_content = np.array(hist_content)
    hist_content = np.sum(hist_content, axis=0)
    hist_binned = [np.sum(hist_content[0:4]),np.sum(hist_content[4:6]), np.sum(hist_content[6:8]), np.sum(hist_content[8:])]
    return hist_binned/np.sum(hist_binned)

class EFT_DC_Prep:
    
    gen_bins = {
        'pt': np.arange(0,1050,50),
        '400inc': [0,200,300,400,500],
        '300inc': [0,200,300,400]
    } # this will be clipped, i.e 400 - inf
    gen_bins_labels = {
        '400inc': ['tt_HZ_bin1','tt_HZ_bin2','tt_HZ_bin3','tt_HZ_bin4'],
        '300inc': ['tt_HZ_bin1','tt_HZ_bin2','tt_HZ_bin3']
    }
    
    # heavy flavour + boson: 'ctp', 'cpQM', 'cpQ3', 'cpt', 'cptb', 'ctW', 'ctZ', 'cbW', 'ctG'
    # 2 sigma interval (sm) for AN2019_011
    an2019_11 = {
        'ctp' : [-14.12,-1.48,32.30,44.48],
        'cpQM': [-3.45,3.33],
        'cpQ3': [-7.21,2.25],
        'cpt' : [-20.91,-14.10,-6.52,4.24],
        'cptb': [-9.87, 9.67],
        'ctW' : [-2.15,-0.29,0.21,1.96],
        'ctZ' : [-2.14,2.19],
        'cbW' : [-4.12,4.09],
        'ctG' : [-1.26,-0.69,0.08,0.79]
     }
    
    def __init__(self, part, aux_file, wgt_file, eft_params):
        self.part   = part
        self.aux_df = pd.read_pickle(aux_file)
        self.wgt_df = pd.read_pickle(wgt_file).loc[:,self.aux_df.index.values]
        self.gen_df = pd.read_pickle(wgt_file).filter(regex=r'[HZ]_\w+') # will need to generalize
        self.eft    = eft_params
        self.beta   = self.compute_beta()
        print(sum(self.beta['ctZ_ctZ']*4)/sum(self.beta['SM'])+sum(self.beta['ctZ']*2)/sum(self.beta['SM'])+1)

    def compute_beta(self):
        # taken from Jon's code
        # Build the experiment matrix 
        x =[np.ones(len(self.aux_df.index.values))]
        beta_cols = ['SM']
        for i in range(len(self.aux_df.columns.values)):
            for j in range(i+1):
                x.append(self.aux_df.iloc[:,i].values * self.aux_df.iloc[:,j].values)
                beta_cols.append(f'{self.aux_df.columns.values[i]}_{self.aux_df.columns.values[j]}')
            x.append(self.aux_df.iloc[:,i].values)
            beta_cols.append(f'{self.aux_df.columns.values[i]}')
        x = np.matrix(x).T
        # Build the result matrix y
        y = np.asmatrix(self.wgt_df.to_numpy()).T
        # Compute beta matrix
        beta = ((x.T * x).I * x.T * y).A
        return pd.DataFrame(data = beta.T, columns=beta_cols)

    def plot_all_eft_scale(self, k_type='pt', opt='all', tag=None):
        kinem = np.clip(self.gen_df[self.part+'_'+k_type], self.gen_bins[k_type][0], self.gen_bins[k_type][-1])
        #h_SM, bins = np.histogram(kinem, bins=self.gen_bins[k_type], weights=self.beta['SM'])
        #for i, row in 
        eft_calc_opt = {'all':self.calc_eft_weight,
                        'q': (lambda x,y,z: self.calc_eft_weight_Q(x,y,z) + self.beta['SM']),
                        'p': (lambda x,y,z: self.calc_eft_weight_P(x,y,z) + self.beta['SM'])}
        w_SM  = self.beta['SM']
        for param in self.eft:
            #w_eft = [eft_calc_opt[opt]({param:i}, kinem, k_type) for i in [-1,1]]
            w_eft  = [eft_calc_opt[opt]({param:i}, kinem, k_type) for i in self.an2019_11[param]]
            fig, (ax,ax2) = plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios':[3,1]})
            fig.subplots_adjust(hspace=0.0)
            n, bins,_ = ax.hist([kinem for i in range(len(w_eft)+1)],
                                #[kinem,kinem,kinem],
                                bins=self.gen_bins[k_type],                 
                                histtype='step',
                                weights = w_eft+[w_SM],#[w_eft[0],w_eft[1],w_SM],
                                #label=np.append([param +'='+str(i) for i in [-1,1]],'SM'))
                                label=np.append([param +'='+str(i) for i in self.an2019_11[param]],'SM'))
            ax.set_yscale('log')
            ax.set_ylabel('Events')
            bin_c = (self.gen_bins[k_type][1:]+self.gen_bins[k_type][:-1])/2
            bin_w = (self.gen_bins[k_type][1:]-self.gen_bins[k_type][:-1])

            for i in range(len(n)-1):
                ax2.errorbar(x=bin_c, xerr=bin_w/2,
                             y=n[i]/n[-1],
                             fmt='.')
                #ax2.errorbar(x=bin_c, xerr=bin_w/2,
                #             y=n[1]/n[2],
                #             fmt='.')
            ax2.axhline(1, color='k', linewidth='1', linestyle='--', dashes=(4,8), snap=True)
            ax2.set_ylim(.5,1.5)
            ax2.set_ylabel('EFT/SM')
            plt.xlim(self.gen_bins[k_type][0],self.gen_bins[k_type][-1])
            plt.xlabel(f'GEN {self.part} pt [GeV]')
            ax.legend()
            #plt.show()
            #plt.close(fig)
            #
        import matplotlib.backends.backend_pdf as matpdf
        pdf = matpdf.PdfPages(f"pdf/{self.part}_WC_calc{opt}{'' if tag is None else '_'+tag}.pdf")
        for fig_ in range(1, plt.gcf().number+1):
            pdf.savefig( fig_ )
        pdf.close()
        plt.close('all')

    def plot_xsec_eft_scale(self, tag=None):
        w_SM = self.beta['SM']
        n_bins = 50
        for param in self.eft:
            low_edge, hi_edge, step = self.an2019_11[param][0], self.an2019_11[param][-1], (self.an2019_11[param][-1]-self.an2019_11[param][0])/n_bins
            x = np.arange(low_edge, hi_edge+step,step)
            y = [np.sum(self.calc_eft_weight({param: i}))/np.sum(w_SM) for i in x]

            plt.figure()
            plt.plot(x, y, label=f'{param}')
            plt.xlabel(f'{param} value')
            plt.ylabel('xsec ratio (EFT/SM)')
            plt.ylim(0.5, 1.5)
            plt.title(f'xsec ratio vs WC {param}')
            #plt.legend()
            #plt.show()
            #plt.clf()
        import matplotlib.backends.backend_pdf as matpdf
        pdf = matpdf.PdfPages(f"pdf/{self.part}_WC_xsec_ratio{'' if tag is None else '_'+tag}.pdf")
        for fig_ in range(1, plt.gcf().number+1):
            pdf.savefig( fig_ )
        pdf.close()
        plt.close('all')            



    def plot_eft_scale(self, eft_dict, k_type='pt'):
        kinem = np.clip(self.gen_df[self.part+'_'+k_type], self.gen_bins[k_type][0], self.gen_bins[k_type][-1])
        #h_SM, bins = np.histogram(kinem, bins=self.gen_bins[k_type], weights=self.beta['SM'])
        #for i, row in 
        w_SM  = self.beta['SM']
        w_eft = self.calc_eft_weight(eft_dict, kinem, k_type)
        #print(w_SM,w_eft)
        print(w_eft)
        plt.hist([kinem,kinem],
                 bins=self.gen_bins[k_type],                 
                 histtype='step',
                 weights = [w_eft,w_SM],
                 label=np.append([k +'='+str(eft_dict[k]) for k in eft_dict],'SM'))
        #plt.yscale('log')
        plt.xlim(self.gen_bins[k_type][0],self.gen_bins[k_type][-1])
        plt.legend()
        plt.show()
        
    def calc_eft_weight_P(self, eft_dict, kinem=None, k_type=None):
        p   = np.sum([self.beta[k1+'_'+k2]*v1*v2\
                      for k1,v1 in eft_dict.items() for k2,v2 in eft_dict.items()], axis=0)
        return p

    def calc_eft_weight_Q(self, eft_dict, kinem=None, k_type=None):
        q   = np.sum([self.beta[k]*v\
                      for k,v in eft_dict.items()], axis=0)
        return q

    def calc_eft_weight(self, eft_dict, kinem=None, k_type=None):
        p   = self.calc_eft_weight_P(eft_dict,kinem,k_type)
              
        q   = self.calc_eft_weight_Q(eft_dict,kinem,k_type) 
             
        r   = self.beta['SM']
        #print(p,q,r)
        return p + q + r

    def create_dc_Df(self, inc='400inc'):
        pt_dist = np.clip(self.gen_df[self.part+'_pt'], self.gen_bins['400inc'][0],self.gen_bins['400inc'][-1])
        pqr_df   = self.beta.groupby(
            pd.cut(pt_dist,self.gen_bins['400inc'], labels=list(map((lambda x : x.replace('_HZ_', self.part)), self.gen_bins_labels[inc])))
        ).sum()
        #
        pqr_list = ['SM'] + list(self.eft) + \
                   list(set([self.eft[i]+'_'+self.eft[j] for i in range(len(self.eft)) for j in range(i,len(self.eft))] +\
                   [self.eft[j]+'_'+self.eft[i] for i in range(len(self.eft)) for j in range(i,len(self.eft))]))
        pqr_df = pqr_df.filter(items=pqr_list, axis='columns')
        pqr_df = pqr_df.divide(pqr_df['SM'], axis='rows')
        #
        wc_min_max = pd.DataFrame({key:[self.aux_df[key].min(), self.aux_df[key].max()] for key in self.aux_df.filter(items=self.eft).columns}, index=pd.CategoricalIndex(['min','max']))
        #
        pqr_df.to_pickle('pqr_df_'+self.part+'.pkl')
        wc_min_max.to_pickle('wc_minmax_df_'+self.part+'.pkl')
    
    def getMCshape(self):
        pt_dist = np.clip(self.gen_df[self.part+'_pt'], self.gen_bins['400inc'][0],self.gen_bins['400inc'][-1])
        h, bins = np.histogram(pt_dist, bins=self.gen_bins['400inc'], weights=self.wgt_df['EFTrwgt183'])
        h_unw, _ = np.histogram(pt_dist, bins=self.gen_bins['400inc'])
        return h,bins, h_unw

    def __getitem__(self,key):
        return getattr(self,key)

def plot4Andrew(mc_hist,nd_hist,nd_bins,raw_cm, c_name='KH_files'):
    fig, (ax0,ax1) = plt.subplots(nrows=2)
    
    nd_err = (np.sqrt(raw_c)/raw_c)*nd_hist/np.sum(nd_hist)
    nd_hist = nd_hist/np.sum(nd_hist)
    bins  = [100,250,350,450]
    bin_w = [200,100,100,100]
    #
    ax0.bar(bins,mc_hist, width=bin_w, label='MC_central', fill=False, edgecolor='tab:red')
    ax0.bar(bins,nd_hist, width=bin_w, yerr=nd_err, ecolor='tab:blue', label=f'{c_name}', fill=False, edgecolor='tab:blue')
    ax1.bar(bins,nd_hist/mc_hist, width=bin_w, yerr=nd_err/mc_hist, ecolor='tab:blue', label=f'{c_name}/central_sample', fill=False, edgecolor='tab:blue')
    ax0.legend()
    ax0.set_yscale('log')
    ax0.set_ylabel('% of Total')
    ax1.set_ylabel(f'{c_name}/central_sample')
    ax1.set_ylim(.5,1.5)
    plt.grid(True)
    plt.xlim(0,500)
    plt.show()

if __name__ == '__main__':
    #main()
    sfile = sys.argv[1]
    #
    pklDir   = 'pkl_files/'
    #aux_file = pklDir+'aux_EFT.pkl'
    aux_file = pklDir+'aux_EFT_ken_tthv2.pkl'
    #aux_file = pklDir+'aux_EFT_ken_tth.pkl'
    file_dict = {'H': pklDir+'eventInfo_EFT_tth.pkl',
                 'Z': pklDir+'eventInfo_EFT_ttZ.pkl',
                 'K': pklDir+'eventInfo_EFT_ken_tthv2.pkl',
                 'ND': pklDir+'eventInfo_EFT_ND_tthv2.pkl'}
    #wgt_file = file_dict['K']
    part = {'K':'H',
            'ND':'H'}
    wgt_file = file_dict[sfile]

    eft_wc    = ( # all 16 considered WC 
        'ctW', 'ctp', 'cpQM', 'ctli', 'cQei', 'ctZ', 'cQlMi', 'cQl3i', 'ctG', 'ctlTi', 'cbW', 'cpQ3', 'ctei', 'cpt', 'ctlSi', 'cptb'
    )
    eft_hb_wc = ( # two heavy + boson Wilson Coefficients
        'ctp', 'cpQM', 'cpQ3', 'cpt', 'cptb', 'ctW', 'ctZ', 'cbW', 'ctG'
    )
    eft_2l_wc = ( # two heavy + two lepton Wilson Coefficients
        'cQl3i', 'ctli', 'cQei', 'cQlMi', 'ctlTi', 'ctei', 'ctlSi'
    )
    

    eft = EFT_DC_Prep(part.get(sfile,sfile), aux_file, wgt_file, eft_hb_wc)
    #print(eft.wgt_df)
    #print(eft['beta']['ctp'])
    #print(eft['beta']['ctp_ctp'])
    #print(eft['beta']['SM'])
    #print(pd.DataFrame.from_dict({'ctZ':[1,-1], 'ctW': 0}))
    
    eft.plot_all_eft_scale( tag='ND_H')
    eft.plot_xsec_eft_scale(tag='ND_H')

    #eft.plot_all_eft_scale(opt='q')
    #eft.plot_all_eft_scale(opt='p')
    #eft.create_dc_Df()
    ##
    #mc_hist = get_root_hist('h_ttz.root')
    #mc_hist = get_root_hist('h_tth.root')
    #c_hist, c_bins, raw_c = eft.getMCshape()
    #print(c_hist)
    ###
    #plot4Andrew(mc_hist,c_hist,c_bins,raw_c, c_name='KH_ND files')
    

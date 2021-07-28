import pandas as pd
import numpy as np
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
from matplotlib import rc
rc("figure", figsize=(6, 6*(6./8.)), dpi=200)

def main():
    ttbb_v12 = compute_beta()
    ttbb_cbw = pd.read_pickle("eventInfo_EFT_ken_ttbb_2018_cbWonly.pkl")
    fig, ax =plt.subplots()
    _makeplots(ttbb_v12,'v12',ax)
    _makeplots(ttbb_cbw,'cbw_only',ax)
    endplt(fig,ax)
    plt.show()

def _makeplots(ttbb,_info,ax):
    #ttbb = pd.read_pickle("eventInfo_EFT_ken_ttbb_2018_cbW.pkl")
    color_dict = {
        'v12': plt.cm.Oranges(np.linspace(0,1,10))[::-1][:6],
        'cbw_only': plt.cm.Blues(np.linspace(0,1,10))[::-1][:6],
    }
    c = color_dict[_info]
    edges = plot_impact(ttbb,f'no_dr_cut {_info}',ax, c[0:2])
    ttbb = ttbb.loc[ttbb['bb_dr']<1.2]
    _= plot_impact(ttbb,f'dr(bb)<1.2 {_info}',ax, c[2:4])
    ttbb = ttbb.loc[ttbb['bb_dr']<0.8]
    _= plot_impact(ttbb,f'dr(bb)<0.8 {_info}',ax, c[4:6])
    ax.set_xlim(edges[0],edges[-1])
    #
    
def plot_impact(ttbb,add_info,ax, c):
    cbw7, edges = np.histogram(ttbb['bb_pt'].clip(0,650), bins = [0,200,300,450,600], weights=ttbb['rwgt_cbW_7'])
    cbwmin7, _ = np.histogram(ttbb['bb_pt'].clip(0,650), bins = [0,200,300,450,600],  weights=ttbb['rwgt_cbW_min7p0'])
    sm, _ = np.histogram(ttbb['bb_pt'].clip(0,650), bins = [0,200,300,450,600], weights=ttbb['rwgt_SM'])
    #
    ax.step(x=edges,y=np.append(cbw7/sm,0), where='post', color=c[0], alpha=.75, label='cbW = 7 {}'.format(add_info))
    ax.step(x=edges,y=np.append(cbwmin7/sm,0), where='post', color=c[1], alpha=.75, label='cbW = -7 {}'.format(add_info))
    return edges
    
def endplt(fig,ax):
    ax.axhline(1,c='k',ls='--')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.legend(framealpha = 0, ncol=2, fontsize='xx-small')

    ax.set_ylim(.75,1.75)
    ax.set_xlabel('bb pT (GeV)')
    fig.suptitle('cbW impact on ttbb')
    ax.grid()

def compute_beta():
    # taken from Jon's code
    # Build the experiment matrix 
    df     = pd.read_pickle('eventInfo_EFT_ken_ttbb_2018_v12.pkl')
    aux_df = pd.read_pickle('pkl_files/aux_EFT.pkl')
    wgt_df = df.filter(items=aux_df.index.values)
    x =[np.ones(len(aux_df.index.values))]
    beta_cols = ['SM']
    for i in range(len(aux_df.columns.values)):
        for j in range(i+1):
            x.append(aux_df.iloc[:,i].values * aux_df.iloc[:,j].values)
            beta_cols.append(f'{aux_df.columns.values[i]}_{aux_df.columns.values[j]}')
        x.append(aux_df.iloc[:,i].values)
        beta_cols.append(f'{aux_df.columns.values[i]}')
    x = np.matrix(x).T
    # Build the result matrix y
    y = np.asmatrix(wgt_df.to_numpy()).T
    # Compute beta matrix
    beta = ((x.T * x).I * x.T * y).A
    beta_df = pd.DataFrame(data = beta.T, columns=beta_cols)
    df['rwgt_cbW_min7p0'] = calc_rwgt(-7,beta_df,'cbW')
    df['rwgt_cbW_7']      = calc_rwgt(7,beta_df,'cbW')
    df['rwgt_SM']         = beta_df['SM']
    return df

def calc_rwgt(v,df=None,wc=None):
    p = df[f'{wc}_{wc}']*v*v
    q = df[f'{wc}']*v
    r = df['SM']
    return p + q + r


if __name__ == '__main__':
    main()

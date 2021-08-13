##
# Code used to make histogram of lambda 2 values and plot of lam2 v splitting intesity difference
import discre


Deng = discre.Tester('/Users/ja17375/DiscrePy/Sheba/Results/Deng/Deng_SKS_SKKS.stk','/Users/ja17375/DiscrePy/Sheba/Runs/Deng_events_2')
Jacks_Split = discre.Tester('/Users/ja17375/DiscrePy/Sheba/Results/Jacks_Split_2/Accepted_SKS_SKKS_all.stk','/Users/ja17375/DiscrePy/Sheba/Runs/Jacks_Split_2')

Dd = pd.read_csv('Deng_diffs.stk',delim_whitespace=True)
Dm = pd.read_csv('Deng_matches.stk',delim_whitespace=True)
Jm = pd.read_csv('../Jacks_Split_2/Accepted_SKS_SKKS_d_SI_matches.stk',delim_whitespace=True)
Jd = pd.read_csv('../Jacks_Split_2/Accepted_SKS_SKKS_d_SI_diffs.stk',delim_whitespace=True)

diff = Jacks_Split.stk_diffs.LAM2.append(Deng.stk_diffs.LAM2)
match = Jacks_Split.stk_matches.LAM2.append(Deng.stk_matches.LAM2)

m_lam2 = Jm.LAM2.append(Dm.LAM2)
m_si = Jm.D_SI.append(Dm.D_SI)
d_lam2 = Jd.LAM2.append(Dd.LAM2)
d_si = Jd.D_SI.append(Dd.D_S)

fig,(ax1,ax2) = plt.subplots(1,2,figsize=(14,6))
    ...:
    ...: ax1.hist([match,diff],bins, histtype='bar', stacked=True,label=["'Matching'","'Discrepent'"])
    ...: ax1.set_xlabel(r'$\lambda _2$ values')
    ...: ax1.set_ylabel('Frequency')
    ...: ax1.legend()
    ...: ax1.set_xlim([0,1])
    ...: ax2.plot(m_lam2,m_si,color='blue',marker='.',ls='None',label="'Matching'")
    ...: ax2.plot(d_lam2,d_si,color='darkorange',marker='.',ls='None',label="'Discrepant'")
    ...: ax2.plot([0, 1.0],[0.4,0.4],ls='dashed',color='black')
    ...: ax2.set_xlabel(r'$\lambda _2$ values')
    ...: ax2.set_ylabel(r'$\Delta$ SI')
    ...: ax2.legend()
    ...: ax2.set_xlim([0,1])
    ...: ax2.set_ylim([0,4.0])
    ...: plt.savefig('/Users/ja17375/DiscrePy/Figures/hist_and_dVs.eps',format='eps',dpi=1000)
    ...: plt.show()

from pathlib import Path
import pandas as pd, numpy as np, json, hashlib
from scipy.stats import norm,t
BASE=Path(__file__).resolve().parent
p=pd.read_csv(BASE/'repo/classifications/country_year_regimes.csv').sort_values(['country','year'])
panel=pd.read_csv(BASE/'repo/data/derived/analytical_panel_2000_2024.csv')
D='Fossil displacement'
assert len(p)==475 and p.country.nunique()==19 and not p.duplicated(['country','year']).any()
assert not p[['renewables_change_twh','fossil_change_twh','fossil_electricity']].isna().any().any()
assert ((p.renewables_change_twh>0)&(p.fossil_change_twh<0)).equals(p.state_primary.eq(D))
assert p.loc[p.year==2000,'state_primary'].ne(D).all()
sp=[]
for c,g in p.groupby('country'):
 g=g.set_index('year')
 for y in g.index:
  if g.loc[y,'state_primary']==D and (y==2000 or g.loc[y-1,'state_primary']!=D):
   end=y
   while end<2024 and g.loc[end+1,'state_primary']==D: end+=1
   sp.append(dict(country=c,start=y,end=end,duration=end-y+1,censored=end==2024))
ep=pd.DataFrame(sp); ep.to_csv(BASE/'reconstructed_displacement_episodes.csv',index=False)
rng=np.random.default_rng(20260901); countries=sorted(p.country.unique()); draws=rng.integers(0,19,(20000,19))
def rates(df,cols):
 out={}
 for col in cols:
  stat=df.groupby('country')[col].agg(['sum','count']).reindex(countries,fill_value=0).to_numpy()
  num,den=stat[:,0],stat[:,1]
  boot=num[draws].sum(1)/den[draws].sum(1)
  loco=(num.sum()-num)/(den.sum()-den)
  out[col]={'events':int(num.sum()),'n':int(den.sum()),'pct':100*num.sum()/den.sum(),'bootstrap_95_pct':(np.quantile(boot,[.025,.975])*100).tolist(),'loco_pct':[float(loco.min()*100),float(loco.max()*100)]}
 return out
res={'commit':'e518a95082ba6f2b4df3b3eea433fb6752f35aac','n_panel':len(p),'n_countries':p.country.nunique(),'state_counts':p.state_primary.value_counts().to_dict(),'episodes':len(ep),'episode_length_counts':ep.duration.value_counts().sort_index().to_dict(),'censored':int(ep.censored.sum()),'analyses':{}}
rows=[]
for H in [2,3]:
 for e in ep[ep.start<=2025-H].itertuples():
  g=p[p.country==e.country].set_index('year'); s=e.start; end=s+H-1
  w=g.loc[s:end]; f0=g.loc[s-1,'fossil_electricity']; f1=g.loc[end,'fossil_electricity']
  rows.append(dict(country=e.country,start=s,horizon=H,endpoint=end,joint_continuity=bool(w.state_primary.eq(D).all()),fossil_continuity=bool(w.fossil_change_twh.lt(0).all()),endpoint_retention=bool(f1<f0),fossil_baseline_twh=f0,fossil_endpoint_twh=f1,change_twh=f1-f0,change_pct=100*(f1-f0)/f0,prior_total_generation_twh=g.loc[s-1,'electricity_generation']))
a=pd.DataFrame(rows);a.to_csv(BASE/'matched_horizon_episode_results.csv',index=False)
for H,g in a.groupby('horizon'):
 assert (g.joint_continuity<=g.fossil_continuity).all() and (g.fossil_continuity<=g.endpoint_retention).all()
 fail=g[~g.joint_continuity]
 g=g.copy(); g['fossil_minus_joint']=g.fossil_continuity.astype(int)-g.joint_continuity.astype(int); g['endpoint_minus_joint']=g.endpoint_retention.astype(int)-g.joint_continuity.astype(int); g['endpoint_minus_fossil']=g.endpoint_retention.astype(int)-g.fossil_continuity.astype(int); out=rates(g,['joint_continuity','fossil_continuity','endpoint_retention','fossil_minus_joint','endpoint_minus_joint','endpoint_minus_fossil']); out['retention_given_joint_failure']=rates(fail,['endpoint_retention'])['endpoint_retention']
 out['magnitude']={'change_pct_q25_median_q75':g.change_pct.quantile([.25,.5,.75]).tolist(),'sum_episode_change_twh':g.change_twh.sum(),'failure_change_pct_q25_median_q75':fail.change_pct.quantile([.25,.5,.75]).tolist()}
 # OLS FE, CR1 covariance normal reference for manuscript replication and t18 sensitivity
 trend=[]
 for omit2020 in [False,True]:
  z=g[g.start.ne(2020)].copy() if omit2020 else g.copy()
  X=np.column_stack([np.ones(len(z)),z.start.to_numpy()-2000,pd.get_dummies(z.country,drop_first=True,dtype=float).to_numpy()]); y=z.joint_continuity.astype(float).to_numpy(); inv=np.linalg.pinv(X.T@X); b=inv@X.T@y; u=y-X@b; meat=np.zeros((X.shape[1],X.shape[1])); G=z.country.nunique();n=len(z);k=np.linalg.matrix_rank(X)
  for c in z.country.unique():
   mask=z.country.eq(c).to_numpy();v=X[mask].T@u[mask];meat+=np.outer(v,v)
  cov=G/(G-1)*(n-1)/(n-k)*inv@meat@inv;se=np.sqrt(cov[1,1]);stat=b[1]/se
  trend.append(dict(exclude_2020=omit2020,n=n,country_clusters=G,estimate_pp=100*b[1],se_pp=100*se,p_normal=2*norm.sf(abs(stat)),p_t_clusters_minus_1=2*t.sf(abs(stat),G-1)))
 out['trend']=trend
 out['retention_sensitivity']={'exclude_2020':rates(g[g.start.ne(2020)],['endpoint_retention'])['endpoint_retention']}
 for threshold in [.001,.0025]:
  material=g.assign(material_retention=(-g.change_twh > threshold*g.prior_total_generation_twh))
  out['retention_sensitivity'][str(threshold)]=rates(material,['material_retention'])['material_retention']
 res['analyses'][str(H)]=out
# Exact directional transitions
nxt=p.groupby('country').state_primary.shift(-1)
origin=p.state_primary.eq(D)&p.year.lt(2024)&nxt.ne('Boundary')
res['D_next_directional']=nxt[origin].value_counts().to_dict();res['D_origin_directional_n']=int(origin.sum())
res['demand']={'displacement_n':int(p.state_primary.eq(D).sum()),'nondecreasing_n':int((p.state_primary.eq(D)&p.electricity_demand_growth_pct.ge(0)).sum()),'missing_demand_growth_in_D':int(p.loc[p.state_primary.eq(D),'electricity_demand_growth_pct'].isna().sum())}
res['period_displacement_pct']=p.assign(period=((p.year-2000)//5)*5+2000,isD=p.state_primary.eq(D)).groupby('period').isD.mean().mul(100).to_dict()
# independent aggregate level context, not ratio from unweighted country shares
agg=p.groupby('year')[['renewables_electricity','fossil_electricity','electricity_generation']].sum()
res['aggregate_2000_2024']={'2000':agg.loc[2000].to_dict(),'2024':agg.loc[2024].to_dict(),'change_pct':((agg.loc[2024]/agg.loc[2000]-1)*100).to_dict()}
# Kaplan-Meier reaches H after surviving completed failures at lengths below H.
def km(e,H):
 survival=1.
 for duration in range(1,H):
  at=int(e.duration.ge(duration).sum()); fail=int((e.duration.eq(duration)&~e.censored).sum())
  if at: survival*=1-fail/at
 return survival*100
res['kaplan_meier_primary']={str(H):km(ep,H) for H in [2,3]}
extra=[]
for column in ['state_threshold_010','state_threshold_025','demand_conditioned']:
 z=p.copy(); isD=(z.state_primary.eq(D)&z.electricity_demand_growth_pct.ge(0)) if column=='demand_conditioned' else z[column].eq(D)
 z['isD']=isD; episodes=[]
 for country,g in z.groupby('country'):
  g=g.set_index('year')
  for year in g.index:
   if g.loc[year,'isD'] and (year==2000 or not g.loc[year-1,'isD']):
    end=year
    while end<2024 and g.loc[end+1,'isD']: end+=1
    episodes.append(dict(country=country,start=year,end=end,duration=end-year+1,censored=end==2024))
 e=pd.DataFrame(episodes);extra.append(dict(specification=column,n_years=int(isD.sum()),n_episodes=len(e),km_SD2=km(e,2),km_SD3=km(e,3)))
res['duration_sensitivities']=extra
# Export manuscript-ready numeric tables (percent or percentage-point units).
summary=[]
for horizon,values in res['analyses'].items():
 for measure,stats in values.items():
  if isinstance(stats,dict) and 'events' in stats:
   summary.append(dict(horizon=int(horizon),measure=measure,events=stats['events'],n=stats['n'],estimate_pct=stats['pct'],bootstrap_low=stats['bootstrap_95_pct'][0],bootstrap_high=stats['bootstrap_95_pct'][1],loco_low=stats['loco_pct'][0],loco_high=stats['loco_pct'][1]))
pd.DataFrame(summary).to_csv(BASE/'nested_horizon_summary.csv',index=False)
(BASE/'audit_results.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))

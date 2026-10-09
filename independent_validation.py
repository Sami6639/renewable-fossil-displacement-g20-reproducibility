import pandas as pd,numpy as np,json
from pathlib import Path
from collections import Counter

from scipy.stats import norm,t
P=Path(__file__).resolve().parent
p=pd.read_csv(P/'repo/data/derived/analytical_panel_2000_2024.csv').sort_values(['country','year'])
r={}
r['quality']={'rows':len(p),'countries':p.country.nunique(),'duplicates':int(p.duplicated(['country','year']).sum()),'years':sorted(p.year.unique().tolist()),'required_missing':p[['renewables_electricity','fossil_electricity','renewables_change_twh','fossil_change_twh','prior_generation_twh']].isna().sum().to_dict()}
rows=[];episodes=[];trans=[]
for country,g in p.groupby('country'):
 g=g.set_index('year')
 assert list(g.index)==list(range(2000,2025))
 def state(y,c=0):
  x=g.loc[y];a=x.renewables_change_twh;b=x.fossil_change_twh;thr=c*x.prior_generation_twh
  if abs(a)<=thr or abs(b)<=thr:return 'B'
  return ('D' if b<0 else 'A') if a>0 else ('C' if b<0 else 'R')
 for y in range(2001,2025):
  for level,change in [('renewables_electricity','renewables_change_twh'),('fossil_electricity','fossil_change_twh')]:
   assert np.isclose(g.loc[y,level]-g.loc[y-1,level],g.loc[y,change],atol=.00001)
 for y in range(2000,2025):
  if state(y)!='D':continue
  if y<2024 and state(y+1)!='B':trans.append(state(y+1))
  if y>2000 and state(y-1)=='D':continue
  assert y>2000
  n=1
  while y+n<=2024 and state(y+n)=='D':n+=1
  episodes.append({'country':country,'start':y,'duration':n,'censored':y+n-1==2024})
  for H in [2,3]:
   endpoint=y+H-1
   if endpoint>2024:continue
   levels=g.loc[y-1:endpoint,'fossil_electricity'].values
   rows.append({'country':country,'start':y,'H':H,'joint':all(state(k)=='D' for k in range(y,endpoint+1)),'fossil':all(np.diff(levels)<0),'retained':levels[-1]<levels[0],'baseline':levels[0],'change_pct':100*(levels[-1]/levels[0]-1)})
e=pd.DataFrame(episodes);a=pd.DataFrame(rows)
r['episodes']={'n':len(e),'durations':e.duration.value_counts().sort_index().to_dict(),'censored':int(e.censored.sum()),'earliest_start':int(e.start.min())};r['transitions']=dict(Counter(trans));r['horizons']={}
# Separate multinomial country multiplicity bootstrap; same inferential design, independent random stream.
rng=np.random.default_rng(931106);countries=sorted(p.country.unique());weights=rng.multinomial(19,[1/19]*19,size=100000)
for H,g in a.groupby('H'):
 vals={}
 g=g.copy()
 for name in ['joint','fossil','retained']:g[name]=g[name].astype(int)
 g['fossil_minus_joint']=g.fossil-g.joint;g['retained_minus_joint']=g.retained-g.joint;g['retained_minus_fossil']=g.retained-g.fossil
 for name in ['joint','fossil','retained','fossil_minus_joint','retained_minus_joint','retained_minus_fossil']:
  den=g.groupby('country').size().reindex(countries,fill_value=0).values
  num=g.groupby('country')[name].sum().reindex(countries,fill_value=0).values
  boot=(weights@num)/(weights@den)*100
  loco=100*(num.sum()-num)/(den.sum()-den)
  vals[name]={'events':int(num.sum()),'n':len(g),'pct':100*g[name].mean(),'bootstrap_independent_100k':np.quantile(boot,[.025,.975]).tolist(),'loco':[min(loco),max(loco)]}
 vals['conditional_retention']={'events':int(g.loc[g.joint==0,'retained'].sum()),'n':int((g.joint==0).sum())}
 vals['median_change_pct']=g.change_pct.median();vals['trend']=[]
 for no2020 in [False,True]:
  z=g[g.start!=2020] if no2020 else g
  # Frisch-Waugh-Lovell within-country slope and cluster sandwich, without explicit dummy inversion.
  x=(z.start-z.groupby('country').start.transform('mean')).to_numpy()
  y=(z.joint-z.groupby('country').joint.transform('mean')).to_numpy()
  b=np.dot(x,y)/np.dot(x,x); u=y-b*x; G=z.country.nunique(); N=len(z); K=G+1
  score=pd.Series(x*u,index=z.index).groupby(z.country).sum()
  variance=(score**2).sum()/(np.dot(x,x)**2)*G/(G-1)*(N-1)/(N-K)
  se=np.sqrt(variance);crit=t.ppf(.975,G-1)
  vals['trend'].append({'exclude_2020':no2020,'b_pp':100*b,'se_pp':100*se,'normal_p':2*norm.sf(abs(b/se)),'t_p':2*t.sf(abs(b/se),G-1),'t_df':G-1,'t_ci_pp':[100*(b-crit*se),100*(b+crit*se)]})
 vals['nested_violations']=int(((g.joint>g.fossil)|(g.fossil>g.retained)).sum())
 r['horizons'][int(H)]=vals
(P/'independent_validation.json').write_text(json.dumps(r,indent=2));a.to_csv(P/'independent_episode_endpoints.csv',index=False)
print(json.dumps(r,indent=2))
# Exact-seed reconciliation against worker's reported endpoints, from independently built outcomes.
exact=np.random.default_rng(20260901).integers(0,19,size=(20000,19))
W=np.array([np.bincount(x,minlength=19) for x in exact])
reported=pd.read_csv(P/'nested_horizon_summary.csv');checks=[]
for row in reported.itertuples():
 z=a[a.H==row.horizon].copy(); mapping={'joint_continuity':'joint','fossil_continuity':'fossil','endpoint_retention':'retained'}
 if row.measure=='retention_given_joint_failure':z=z[~z.joint]; v=z.retained.astype(int)
 elif row.measure in mapping:v=z[mapping[row.measure]].astype(int)
 else:
  left,right={'fossil_minus_joint':('fossil','joint'),'endpoint_minus_joint':('retained','joint'),'endpoint_minus_fossil':('retained','fossil')}[row.measure]
  v=z[left].astype(int)-z[right].astype(int)
 sums=v.groupby(z.country).sum().reindex(countries,fill_value=0).values
 counts=z.groupby('country').size().reindex(countries,fill_value=0).values
 ci=np.quantile(100*(W@sums)/(W@counts),[.025,.975]);assert np.allclose(ci,[row.bootstrap_low,row.bootstrap_high])
 checks.append((row.horizon,row.measure,'PASS'))
print('EXACT BOOTSTRAP',checks)
# Classifications recreated from signs and threshold boundary rules.
cl=pd.read_csv(P/'repo/classifications/country_year_regimes.csv')
for threshold,col in [(0,'state_primary'),(.001,'state_threshold_010'),(.0025,'state_threshold_025')]:
 dr=cl.renewables_change_twh;df=cl.fossil_change_twh;boundary=(dr.abs()<=threshold*cl.prior_generation_twh)|(df.abs()<=threshold*cl.prior_generation_twh)
 computed=np.where(boundary,'Boundary',np.where(dr>0,np.where(df<0,'Fossil displacement','Additive expansion'),np.where(df<0,'Joint contraction','Fossil resurgence')))
 print(col,'mismatches',int((computed!=cl[col]).sum()),'boundary',int(boundary.sum()))
 assert (computed==cl[col]).all()
print('No fossil zero baseline:',bool((a.baseline>0).all()))

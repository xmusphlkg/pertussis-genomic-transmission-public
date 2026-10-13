#!/usr/bin/env python3
"""Recalculate the Letter's numerical evidence from frozen metadata; no network calls."""
from pathlib import Path
import json
import pandas as pd
R=Path(__file__).resolve().parents[1]
O=R/'results/recomputed'; O.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(R/'data/accession_metadata.tsv',sep='\t')
assert not d.tree_sample_id.duplicated().any()
f=d[(d.epidemic_period=='resurgence') & d.country_iso3.isin(['AUS','CHN','JPN'])].copy()
t=f[f.primary_model_lineage_id=='L1_02.07'].copy()
assert t.groupby('country_iso3').size().to_dict()=={'AUS':3,'CHN':148,'JPN':54}
for col in ['collection_lower','collection_upper_effective','public_date']:
 f[col]=pd.to_datetime(f[col]);t[col]=pd.to_datetime(t[col])
def counts(g,date):
 return dict(definitely_collected=int((g.collection_upper_effective<=date).sum()),possibly_collected=int((g.collection_lower<=date).sum()),public=int((g.public_date<=date).sum()))
case=pd.read_csv(R/'data/country_month_cases.tsv',sep='\t');case['model_month']=pd.to_datetime(case.model_month)
ms=[]; project_ms=[]; lag=[]; curves=[]; ks=[]; projects=[]; matches=[]
for c,g in t.groupby('country_iso3'):
 cc=case[(case.country_iso3==c)&(case.case_data_available==True)].sort_values('model_month')
 base=cc[cc.model_month.dt.year==2019];post=cc[cc.model_month.dt.year>=2023]
 mx=base.cases.max();med=base.cases.median()
 above=post.cases.gt(mx); consecutive=above & above.shift(-1,fill_value=False)
 dates={'maximum_2019':post.loc[above,'model_month'].iloc[0],'median_2019':post.loc[post.cases.gt(med),'model_month'].iloc[0],'two_consecutive':post.loc[consecutive,'model_month'].iloc[0],'peak':post.loc[post.cases.idxmax(),'model_month']}
 for name,date in dates.items():
  ms.append(dict(country_iso3=c,definition=name,date=str(date.date()),**counts(g,date)))
  for p,pg in g.groupby('project_id'):project_ms.append(dict(country_iso3=c,project_id=p,definition=name,date=str(date.date()),**counts(pg,date)))
 lag.append(dict(country_iso3=c,n=len(g),median_min_days=g.lag_min_days.median(),median_max_days=g.lag_max_days.median()))
 for date in sorted(set(g.collection_lower)|set(g.collection_upper_effective)|set(g.public_date)|set(dates.values())):
  curves.append(dict(country_iso3=c,date=str(date.date()),**counts(g,date)))
 for k in range(1,len(g)+1):
  lo=g.collection_lower.sort_values().iloc[k-1];hi=g.collection_upper_effective.sort_values().iloc[k-1];pub=g.public_date.sort_values().iloc[k-1]
  ks.append(dict(country_iso3=c,k=k,collection_lower=str(lo.date()),collection_upper=str(hi.date()),public_date=str(pub.date()),displacement_min_days=(pub-hi).days,displacement_max_days=(pub-lo).days))
 for p,pg in g.groupby('project_id'):
  vc=pg.public_date.value_counts();n=int(vc.max());dt=min(vc[vc==n].index)
  projects.append(dict(country_iso3=c,project_id=p,n=len(pg),modal_date=str(dt.date()),modal_n=n,modal_pct=round(n/len(pg)*100,1),median_min_days=pg.lag_min_days.median(),median_max_days=pg.lag_max_days.median()))
for (c,p,y),g in f.groupby(['country_iso3','project_id','year']):
 a=g[g.primary_model_lineage_id=='L1_02.07'];b=g[g.primary_model_lineage_id!='L1_02.07']
 if len(a)==0 or len(b)==0:continue
 amin,amax=a.lag_min_days.median(),a.lag_max_days.median();bmin,bmax=b.lag_min_days.median(),b.lag_max_days.median()
 def modal(g):
  v=g.public_date.value_counts();return min(v[v==v.max()].index)
 matches.append(dict(country_iso3=c,project_id=p,year=y,n_target=len(a),n_comparator=len(b),target_min=amin,target_max=amax,comparator_min=bmin,comparator_max=bmax,shared_modal=modal(a)==modal(b),direction='target_shorter' if amax<bmin else 'target_longer' if amin>bmax else 'intervals_overlap'))
for name,rows in [('milestones',ms),('project_milestones',project_ms),('lag_summary',lag),('cumulative_counts',curves),('count_displacement',ks),('project_release',projects),('matched_comparators',matches)]:pd.DataFrame(rows).to_csv(O/f'{name}.tsv',sep='\t',index=False)
x=pd.DataFrame(projects);m=pd.DataFrame(matches)
assert len(x)==13 and (x.modal_pct==100).sum()==11
assert len(m)==10 and m.n_target.sum()==81 and m.n_comparator.sum()==44 and m.shared_modal.all()
assert m.direction.value_counts().to_dict()=={'intervals_overlap':5,'target_shorter':4,'target_longer':1}
assert x.loc[x.project_id=='PRJNA1182239','modal_n'].iloc[0]==68
assert x.loc[x.project_id=='PRJNA1182239','modal_pct'].iloc[0]==82.9
expected={'CHN':[(15,29,0),(89,106,0)],'JPN':[(38,39,0),(48,49,16)]}
for c,vals in expected.items():
 for definition,counts_expected in zip(['maximum_2019','peak'],vals):
  row=next(v for v in ms if v['country_iso3']==c and v['definition']==definition)
  assert tuple(row[k] for k in ['definitely_collected','possibly_collected','public'])==counts_expected
summary={'target_records':len(t),'projects':len(x),'single_date_projects':int((x.modal_pct==100).sum()),'larger_projects':int((x.n>=5).sum()),'larger_single_date_projects':int(((x.n>=5)&(x.modal_pct==100)).sum()),'matched_strata':len(m),'all_numerical_checks_passed':True}
(O/'audit.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

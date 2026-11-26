import glob, json, re
import numpy as np
import pandas as pd
import rdata
from scipy.stats import spearmanr, mannwhitneyu, chi2
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, FixedLocator
genes=['CCNE1','CCNE2']
o='out/'
w='work/'
raw=[pd.read_csv(f'upload/BRCA_{i}_50_50.csv').set_index('Patient') for i in (898,9134)]
clin=pd.read_csv(glob.glob('upload/TCGA*')[0],low_memory=False).set_index('bcr_patient_barcode')
assert set(raw[0].index)==set(raw[1].index)
t=raw[0][['Days','Status']].join(raw[1][['Days','Status']],rsuffix='2',validate='one_to_one').join(clin,rsuffix='_clinical',validate='one_to_one')
print('endpoint checks',len(t),[(a,b,int((t[a]!=t[b]).sum())) for a,b in [('Days','Days2'),('Status','Status2'),('Days','Days_clinical'),('Status','Status_clinical')]])
assert all((t[a]==t[b]).all() for a,b in [('Days','Days2'),('Status','Status2'),('Days','Days_clinical'),('Status','Status_clinical')])
for g,x in zip(genes,raw):
 t[g]=np.log2(x.Expression+1)
 t[g+'high']=x.Group.eq('High').astype(int)
 print(g,x.groupby('Group').Expression.agg(['min','max','size']))
t['time']=t.Days/365.25
t['event']=t.Status.map({'Alive':0,'Dead':1})
t['age']=pd.to_numeric(t.age_at_diagnosis,errors='coerce')/10
t['er']=t.er_status_by_ihc.map({'Negative':0,'Positive':1})
t['stage']=t.ajcc_pathologic_tumor_stage.str.extract(r'^Stage (IV|III|II|I)(?:[ABC])?$')[0].map({'I':1,'II':2,'III':3,'IV':4})
t['female']=t.gender.eq('FEMALE')
p=pd.read_csv(w+'met_patient.tsv',sep='\t',comment='#')
s=pd.read_csv(w+'met_sample.tsv',sep='\t',comment='#')
m=s.merge(p,on='PATIENT_ID',validate='many_to_one').set_index('SAMPLE_ID')
for gene,i in zip(genes,(898,9134)):
 e=pd.read_json(w+f'met{i}.json').set_index('sampleId')
 assert e.index.is_unique and e.entrezGeneId.eq(i).all()
 m[gene]=e.value
m=m[m.SAMPLE_TYPE.eq('Primary')]
m['time']=m.OS_MONTHS/12
m['event']=m.OS_STATUS.map({'0:LIVING':0,'1:DECEASED':1})
m['er']=m.ER_STATUS.map({'Negative':0,'Positive':1})
m['age']=m.AGE_AT_DIAGNOSIS/10
m['stage']=m.TUMOR_STAGE.where(m.TUMOR_STAGE.isin([1,2,3,4]))
m['grade']=m.GRADE.where(m.GRADE.isin([1,2,3]))
x=pd.read_csv(w+'scanexpr.csv',index_col=0).T
assert set(x.columns)==set(genes) and x.index.is_unique
b=next(iter(rdata.read_rda(w+'scan.rda').values())).set_index('title')
assert b.index.is_unique and set(x.index)==set(b.index)
b=b.loc[~b.isRepl.astype(bool)].join(x,validate='one_to_one')
assert len(b)==3273
b['time']=pd.to_numeric(b.ovrallSurvDays,errors='coerce')/365.25
b['event']=pd.to_numeric(b.ovrallSurvEvent,errors='coerce').where(lambda x:x.isin([0,1]))
b['age']=pd.to_numeric(b.ageAtDx,errors='coerce')/10
b['er']=b.er_Status.astype(str).map({'0':0,'1':1})
b['size']=pd.to_numeric(b.tumor_size,errors='coerce')/10
b['node']=b.lymphNode_Status.map({'NodeNegative':0,'NodePositive':1})
b['grade']=b.nhg.astype(str).map({'G1':1,'G2':2,'G3':3})
sets={'TCGA-BRCA':t,'METABRIC':m,'SCAN-B':b}
base={'TCGA-BRCA':'age + er + C(stage)','METABRIC':'age + er + C(stage)','SCAN-B':'age + er + size + node'}
adjust={'TCGA-BRCA':['age','er','stage'],'METABRIC':['age','er','stage'],'SCAN-B':['age','er','size','node']}
records=[]; desc=[]; klin=[]; summary=[]; shape=[]; background=[]
def fit(d,form,cohort,model,terms):
 f=CoxPHFitter().fit(d,'time','event',formula=form)
 ph=proportional_hazard_test(f,d,time_transform='rank').summary
 for g in terms:
  r=f.summary.loc[g]
  records.append(dict(cohort=cohort,model=model,gene=g,n=len(d),events=int(d.event.sum()),HR=r['exp(coef)'],lower=r['exp(coef) lower 95%'],upper=r['exp(coef) upper 95%'],p=r.p,ph=ph.loc[g,'p'],c=f.concordance_index_))
 return f
plt.rcParams.update({'font.family':'serif','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(8.1,3.2))
for name,d in sets.items():
 d=d.dropna(subset=['time','event']+genes).query('time >= 0').copy()
 assert d.index.is_unique
 for g in genes:
  if name!='TCGA-BRCA': d[g+'high']=(d[g]>d[g].median()).astype(int)
  d[g]=(d[g]-d[g].mean())/d[g].std()
 cc=d.dropna(subset=adjust[name]).copy()
 rev=KaplanMeierFitter().fit(d.time,1-d.event).median_survival_time_
 summary.append(dict(cohort=name,n=len(d),events=int(d.event.sum()),complete=len(cc),complete_events=int(cc.event.sum()),reverseKM=rev,rho=spearmanr(d.CCNE1,d.CCNE2).statistic))
 for g in genes:
  a=d[d[g+'high']==0]; h=d[d[g+'high']==1]
  lr=logrank_test(h.time,a.time,h.event,a.event)
  desc.append(dict(cohort=name,gene=g,nlow=len(a),dlow=int(a.event.sum()),nhigh=len(h),dhigh=int(h.event.sum()),p=lr.p_value))
  q=d.dropna(subset=['er'])
  klin.append(dict(cohort=name,gene=g,n=len(q),meanERneg=q.loc[q.er.eq(0),g].mean(),meanERpos=q.loc[q.er.eq(1),g].mean(),p=mannwhitneyu(q.loc[q.er.eq(0),g],q.loc[q.er.eq(1),g]).pvalue))
  fit(d,g,name,'Unadjusted',[g])
  fit(cc,g,name,'Complete-case unadjusted',[g])
  main=fit(cc,g+' + '+base[name],name,'Adjusted',[g])
  spline=CoxPHFitter().fit(cc,'time','event',formula=f'bs({g}, df=3) + '+base[name])
  clinical=CoxPHFitter().fit(cc,'time','event',formula=base[name])
  shape.append(dict(cohort=name,gene=g,nonlinear_p=chi2.sf(2*(spline.log_likelihood_-main.log_likelihood_),2),delta_c=main.concordance_index_-clinical.concordance_index_))
  late=cc[cc.time.gt(5)].assign(entry=5.)
  if name=='METABRIC':
   f=CoxPHFitter().fit(late,'time','event',entry_col='entry',formula=g+' + '+base[name])
   r=f.summary.loc[g]
   records.append(dict(cohort=name,model='Adjusted after 5 years',gene=g,n=len(late),events=int(late.event.sum()),HR=r['exp(coef)'],lower=r['exp(coef) lower 95%'],upper=r['exp(coef) upper 95%'],p=r.p))
  if name!='TCGA-BRCA':
   col='CLAUDIN_SUBTYPE' if name=='METABRIC' else 'pam50_subtype'
   sub=cc[cc[col].astype(str).isin(['LumA','LumB','Basal','Her2','Normal','claudin-low','luminal A','luminal B','HER2-enriched','Luminal A','Luminal B','Normal-like','Basal-like'])].copy()
   if len(sub)>100:
    sub['subtype']=sub[col].astype(str)
    fit(sub,g+' + '+base[name]+' + C(subtype)',name,'Adjusted with subtype',[g])
  fit(d.dropna(subset=['age','er']),g+' + age + er',name,'Age and ER',[g])
  fit(cc,g+'high + '+base[name],name,'Adjusted median',[g+'high'])
  if name!='TCGA-BRCA':
   extra=cc.dropna(subset=['grade'])
   fit(extra,g+' + '+base[name]+' + C(grade)',name,'Adjusted with grade',[g])
  else:
   fit(cc[cc.female],g+' + '+base[name],name,'Female only',[g])
  short=cc.assign(event=cc.event*(cc.time<=5),time=cc.time.clip(upper=5))
  fit(short,g+' + '+base[name],name,'Adjusted through 5 years',[g])
  inter=fit(cc,g+' * er + '+base[name],name,'ER interaction',[g,g+':er'])
  if name=='TCGA-BRCA':
   ax=axes[genes.index(g)]
   for label,k,color in [('Low',a,'#376894'),('High',h,'#b94b38')]:
    km=KaplanMeierFitter().fit(k.time,k.event,label=label)
    km.plot_survival_function(ax=ax,color=color,show_censors=True,censor_styles={'ms':2,'marker':'+'},ci_alpha=.12)
    ax.text(-.105,-(.34 if label=='Low' else .44),label,transform=ax.transAxes,ha='right',color=color)
    for year in [0,2,4,6,8,10]:
     ax.text(year,-(.34 if label=='Low' else .44),str((k.time>=year).sum()),transform=ax.get_xaxis_transform(),ha='center',color=color)
   ax.set(xlim=(0,10.15),ylim=(0,1.03),xlabel='Years since diagnosis',ylabel='Overall survival probability',title=f'{g}   Log-rank P = {lr.p_value:.3g}',xticks=[0,2,4,6,8,10])
   ax.legend(loc='lower left',frameon=False)
   ax.xaxis.set_label_coords(.5,-.17)
 fit(cc,'CCNE1 + CCNE2 + '+base[name],name,'Joint adjusted',genes)
fig.subplots_adjust(left=.10,right=.98,bottom=.40,top=.88,wspace=.35)
fig.savefig(o+'Figure1.png',dpi=330); plt.close(fig)
res=pd.DataFrame(records)
for _,inds in res.groupby(['cohort','model']).groups.items():
 p=res.loc[inds,'p'].to_numpy()
 order=np.argsort(p); q=np.empty(len(p))
 q[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(max=1)
 res.loc[inds,'q']=q
ext=res[(res.model=='Adjusted')&(res.cohort!='TCGA-BRCA')].p
order=ext.sort_values().index
res.loc[order,'external_q']=np.minimum.accumulate((res.loc[order,'p'].to_numpy()*4/np.arange(1,5))[::-1])[::-1].clip(max=1)
for name,tab in [('cox',res),('cohorts',pd.DataFrame(summary)),('logrank',pd.DataFrame(desc)),('er_association',pd.DataFrame(klin)),('model_shape',pd.DataFrame(shape))]: tab.to_csv(o+name+'.csv',index=False)
fig,ax=plt.subplots(figsize=(7.7,3.0))
for i,(_,r) in enumerate(res[res.model.eq('Adjusted')].sort_values(['gene','cohort'],key=lambda x:x.map({'TCGA-BRCA':'1','METABRIC':'2','SCAN-B':'3'}).fillna(x)).iterrows()):
 y=5-i+(.45 if i<3 else 0)
 ax.errorbar(r.HR,y,xerr=[[r.HR-r.lower],[r.upper-r.HR]],fmt='o',color=['#376894','#b94b38'][i//3],capsize=3)
 ax.text(.69,y,r.gene+'  '+r.cohort,ha='right',va='center')
 ax.text(2.04,y,f'{r.HR:.2f} ({r.lower:.2f}–{r.upper:.2f})',ha='left',va='center')
ax.set(xscale='log',xlim=(.70,1.99),ylim=(-.6,6),yticks=[],xlabel='Adjusted hazard ratio per cohort SD')
ax.set_xticks([.75,1,1.25,1.5,1.75],labels=['0.75','1','1.25','1.50','1.75'])
ax.xaxis.set_minor_formatter(NullFormatter())
ax.axvline(1,color='.5',ls='--',lw=1)
ax.spines['left'].set_visible(False)
fig.subplots_adjust(left=.28,right=.77,bottom=.20)
fig.savefig(o+'Figure2.png',dpi=330);plt.close(fig)
gs=[]; ranks=[]
for g,i in zip(genes,(898,9134)):
 folder=w+'CCNE/'+g+'/'
 r=pd.read_csv(glob.glob(folder+'*.rnk')[0],sep='\t',header=None,index_col=0)[1].rename(g)
 assert list(r[r.eq(1)].index)==[i]
 ranks.append(r)
 table=pd.read_csv(glob.glob(folder+'enrichment*.txt')[0],sep='\t')
 table.insert(0,'gene',g); gs.append(table)
 html=open(glob.glob(folder+'Report*.html')[0]).read()
 m=re.search(r'var background = (\[)',html)
 bg=pd.DataFrame(json.JSONDecoder().raw_decode(html[m.start(1):])[0])
 bg.insert(0,'gene',g); background.append(pd.concat([table[bg.columns.drop('_row')],bg.drop(columns='_row')],ignore_index=True).drop_duplicates('geneSet'))
pd.concat(background).to_csv(o+'kegg_all_pathways.csv',index=False)
bg1,bg2=background
bgs=bg1.merge(bg2,on='geneSet',suffixes=('_CCNE1','_CCNE2'))
print('KEGG FDR counts',[(int(x.FDR.lt(.05).sum()),int((x.FDR.lt(.05)&x.normalizedEnrichmentScore.gt(0)).sum())) for x in background])
print('Shared FDR<.05',bgs.query('FDR_CCNE1 < .05 and FDR_CCNE2 < .05')[['geneSet','normalizedEnrichmentScore_CCNE1','normalizedEnrichmentScore_CCNE2']].to_string(index=False))
print('Cell cycle',bgs[bgs.geneSet.eq('hsa04110')].to_string(index=False))
full=pd.concat(gs); full.to_csv(o+'kegg_exports.csv',index=False)
common=gs[0].merge(gs[1],on=['geneSet','description'],suffixes=('_CCNE1','_CCNE2'))
def overlap(row):
 a,b=[set(row['leadingEdgeId_'+g].split(';')) for g in genes]
 return pd.Series(dict(shared=len(a&b),union=len(a|b),jaccard=len(a&b)/len(a|b)))
common=common.join(common.apply(overlap,axis=1)); common.to_csv(o+'kegg_overlap.csv',index=False)
rp=pd.concat(ranks,axis=1).dropna().drop(index=[898,9134],errors='ignore')
print('rank concordance',len(rp),spearmanr(rp.CCNE1,rp.CCNE2).statistic)
print('cross rank',ranks[0].get(9134),ranks[1].get(898))
selected={
 'hsa04110':'Cell cycle','hsa03030':'DNA replication','hsa03013':'Nucleocytoplasmic transport',
 'hsa03440':'Homologous recombination','hsa03460':'Fanconi anemia pathway',
 'hsa03430':'Mismatch repair','hsa03050':'Proteasome','hsa04340':'Hedgehog signaling'
}
fig,ax=plt.subplots(figsize=(7.7,3.5))
for j,(g,bg) in enumerate(zip(genes,background)):
 z=bg.set_index('geneSet').loc[list(selected)]
 y=np.arange(len(z))+(j-.5)*.31
 ax.barh(y,z.normalizedEnrichmentScore,height=.27,label=g,color=['#376894','#b94b38'][j])
 for yy,(_,r) in zip(y,z.iterrows()):
  if r.FDR<.05:
   ax.text(r.normalizedEnrichmentScore+(.055 if r.normalizedEnrichmentScore>0 else -.055),yy,'*',
           va='center',ha='left' if r.normalizedEnrichmentScore>0 else 'right')
ax.set(yticks=range(len(selected)),yticklabels=list(selected.values()),xlim=(-2.5,3.6),xlabel='Normalized enrichment score')
ax.invert_yaxis()
ax.axvline(0,color='black',lw=.6)
ax.legend(loc='lower right',frameon=False)
fig.tight_layout();fig.savefig(o+'Figure3.png',dpi=330);plt.close(fig)

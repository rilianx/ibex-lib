import json,collections,sys
R=collections.defaultdict(dict); setof={}
for f,s in [('results/local/runs-closed.jsonl','223'),('results/local/runs-open.jsonl','75')]:
    for l in open(f):
        d=json.loads(l)
        if d['arm'] in ('astatsc','acomb'): R[d['instance']][d['arm']]=d['result']; setof[d['instance']]=s
st=collections.Counter(); same=0; tot=0; diff=[]
for i,v in R.items():
    a=v.get('astatsc'); 
    if not a: continue
    st[(setof[i],a.get('status'))]+=1
    if 'acomb' in v and a.get('status')==v['acomb'].get('status')=='complete':
        tot+=1; same+= a['nodes']==v['acomb']['nodes']
        if a['nodes']!=v['acomb']['nodes']: diff.append((i,v['acomb']['nodes'],a['nodes']))
print('statuses',dict(st)); print('astatsc vs acomb nodes, both closed: %d/%d same'%(same,tot), diff[:10])
K=['cvx_nodes','cvx_with_dual','cvx_certified','cvx_close','cvx_better_lb','cvx_center_indef','cvx_raised','cvx_pruned']
for s in ('223','75','all'):
    T=collections.Counter(); n=0; inst=collections.Counter()
    for i,v in R.items():
        a=v.get('astatsc')
        if not a or (s!='all' and setof[i]!=s) or 'cvx_nodes' not in a: continue
        n+=1
        for k in K: T[k]+=a.get(k,0)
        if a.get('cvx_certified',0)>0: inst['certified']+=1
        if a.get('cvx_better_lb',0)>0: inst['better_lb']+=1
        if a.get('cvx_close',0)>0: inst['close']+=1
        if a.get('cvx_pruned',0)>0: inst['pruned']+=1
    N=T['cvx_nodes'] or 1
    print('\n## set %s: %d instances, %d nodes'%(s,n,T['cvx_nodes']))
    print('| counter | total | % of nodes | instances with > 0 |'); print('|---|---|---|---|')
    for k in K: print('| %s | %d | %.2f | %s |'%(k,T[k],100*T[k]/N, inst.get({'cvx_certified':'certified','cvx_better_lb':'better_lb','cvx_close':'close','cvx_pruned':'pruned'}.get(k,''),'')))
# per instance: share of nodes certified, top instances
rows=[]
for i,v in R.items():
    a=v.get('astatsc')
    if a and a.get('cvx_nodes'): rows.append((a.get('cvx_certified',0)/a['cvx_nodes'],i,setof[i],a['cvx_nodes'],a.get('cvx_certified',0),a.get('cvx_close',0),a.get('cvx_pruned',0),a.get('status'),round(a.get('time',0),1)))
rows.sort(reverse=True)
print('\n| instance | set | nodes | certified | % | would close | pruned | status | time |'); print('|---|---|---|---|---|---|---|---|---|')
for r in rows:
    if r[4]>0: print('| %s | %s | %d | %d | %.0f | %d | %d | %s | %s |'%(r[1],r[2],r[3],r[4],100*r[0],r[5],r[6],r[7],r[8]))

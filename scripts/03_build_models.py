import re,pickle,os,textwrap
from collections import defaultdict

# ---------- genome ----------
G={}; ACC={}; order=[]
name=None
for line in open('/home/claude/U25.fa'):
    if line.startswith('>'):
        acc=line[1:].split()[0]; c=int(re.search(r'chromosome (\d+)',line).group(1))
        name=acc; G[acc]=[]; ACC[acc]=c; order.append((c,acc))
    else: G[name].append(line.strip())
G={k:''.join(v).upper() for k,v in G.items()}
order.sort()

# ---------- source annotations (symbols) ----------
recs=pickle.load(open('/home/claude/recs.pkl','rb'))
meta={}
for r in recs:
    k=r['lt']
    if k not in meta or len(r['tr'])>len(meta[k]['tr']): meta[k]=r

# ---------- parse miniprot gff ----------
def parse(g):
    out={}; cur=None
    for line in open(g):
        if line.startswith('#'): continue
        f=line.rstrip('\n').split('\t')
        if len(f)<9: continue
        at=dict(kv.split('=',1) for kv in f[8].split(';') if '=' in kv)
        if f[2]=='mRNA':
            if at.get('Rank')!='1': cur=None; continue
            t=at['Target'].split()[0]
            cur=dict(tgt=t,chrom=f[0],s=int(f[3]),e=int(f[4]),score=int(f[5]),strand=f[6],
                     ident=float(at['Identity']),cds=[])
            out[t]=cur
        elif f[2]=='CDS' and cur is not None:
            cur['cds'].append((int(f[3]),int(f[4]),int(f[7])))
    return out
A=parse('/home/claude/hap1.gff'); B=parse('/home/claude/hap2.gff')

# ---------- pair & choose ----------
chosen=[]
usedB=set()
for a,ma in A.items():
    b='ACKKBF_B'+a.split('_A')[1]
    mb=B.get(b)
    if mb is None:
        chosen.append((ma,a,None,ma['ident'],None,'A-only')); continue
    usedB.add(b)
    if mb['score']>ma['score']: chosen.append((mb,a,b,ma['ident'],mb['ident'],'B'))
    elif ma['score']>mb['score']: chosen.append((ma,a,b,ma['ident'],mb['ident'],'A'))
    else: chosen.append((ma,a,b,ma['ident'],mb['ident'],'tie'))
for b,mb in B.items():
    if b in usedB: continue
    a='ACKKBG_A'+b.split('_B')[1]
    if a in A: continue
    chosen.append((mb,None,b,None,mb['ident'],'B-only'))

# ---------- resolve genomic overlaps ----------
chosen.sort(key=lambda x:-x[0]['score'])
occ=defaultdict(list); keep=[]
def ov(s1,e1,s2,e2): return max(0,min(e1,e2)-max(s1,s2)+1)
for item in chosen:
    m=item[0]; L=m['e']-m['s']+1; clash=False
    for (s,e,st) in occ[m['chrom']]:
        if st==m['strand'] and ov(m['s'],m['e'],s,e) > 0.5*min(L,e-s+1): clash=True; break
    if clash: continue
    occ[m['chrom']].append((m['s'],m['e'],m['strand'])); keep.append(item)
keep.sort(key=lambda x:(ACC[x[0]['chrom']],x[0]['s']))
print("採用モデル:",len(keep))

# ---------- translate ----------
CODON={}
bases='TCAG'; aas='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
i=0
for b1 in bases:
    for b2 in bases:
        for b3 in bases:
            CODON[b1+b2+b3]=aas[i]; i+=1
def rc(s): return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def translate(nt):
    return ''.join(CODON.get(nt[i:i+3],'X') for i in range(0,len(nt)-2,3))

final=[]; n=0
stats=dict(ok=0,internal_stop=0,no_met=0,frameshift=0)
for m,a,b,ia,ib,src in keep:
    blocks=sorted(m['cds'])
    seq=''.join(G[m['chrom']][s-1:e] for s,e,p in blocks)
    if m['strand']=='-':
        seq=rc(seq); ph=blocks[-1][2]
    else: ph=blocks[0][2]
    seq=seq[ph:]
    if len(seq)%3: stats['frameshift']+=1
    prot=translate(seq)
    if prot.endswith('*'): prot=prot[:-1]
    if '*' in prot: stats['internal_stop']+=1
    if not prot.startswith('M'): stats['no_met']+=1
    if '*' not in prot and prot.startswith('M'): stats['ok']+=1
    n+=5
    src_lt=b if src in ('B','B-only') else a
    md=meta.get(src_lt,{})
    sym=(md.get('gene') or '') if md else ''
    prod=(md.get('product') or 'hypothetical protein') if md else 'hypothetical protein'
    final.append(dict(lt='APU25_%05d'%n,chrom=m['chrom'],chrnum=ACC[m['chrom']],s=m['s'],e=m['e'],
        strand=m['strand'],blocks=blocks,phase=ph,prot=prot,nt=seq,sym=sym,prod=prod,
        src=src_lt,srchap=src,identA=ia,identB=ib,ident=m['ident'],score=m['score']))
print(stats)
pickle.dump(final,open('/home/claude/final.pkl','wb'))

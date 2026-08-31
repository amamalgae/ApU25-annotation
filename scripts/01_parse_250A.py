import re,glob,os,pickle
recs=[]
for hap in ['haplotype1','haplotype2']:
    for f in glob.glob(f'/home/claude/ann/UTEX_250-A_{hap}_*/*.gb'):
        chrn=re.search(r'chromosome (\d+),',os.path.basename(f)).group(1)
        txt=open(f,encoding='utf-8',errors='replace').read()
        feat=txt.split('FEATURES')[1]
        for m in re.finditer(r'^     CDS +(\S[^\n]*(?:\n {21}[^/][^\n]*)*)\n((?:^ {21}/.*\n(?: {21}[^/][^\n]*\n)*)+)',feat,re.M):
            loc=re.sub(r'\s','',m.group(1)); q=m.group(2)
            def get(k):
                mm=re.search(r'/%s="((?:[^"]|\n)*)"'%k,q)
                return re.sub(r'\s+',' ',mm.group(1)).strip() if mm else ''
            def gett(k):
                mm=re.search(r'/%s="((?:[^"]|\n)*)"'%k,q)
                return re.sub(r'\s','',mm.group(1)) if mm else ''
            recs.append(dict(hap=hap,chr=chrn,loc=loc,lt=get('locus_tag'),gene=get('gene'),
                             prod=get('product'),pid=get('protein_id'),tr=gett('translation')))
pickle.dump(recs,open('/home/claude/recs.pkl','wb'))
print("CDS total:",len(recs))
from collections import Counter
c=Counter(r['prod'] for r in recs)
print("distinct products:",len(c))
print("hypothetical:",c['hypothetical protein'], "=%.1f%%"%(100*c['hypothetical protein']/len(recs)))
print("\n--- top 30 named products ---")
for p,n in c.most_common(35):
    print("%5d  %s"%(n,p[:100]))

import pickle,re,textwrap,os
from collections import defaultdict
final=pickle.load(open('/home/claude/final.pkl','rb'))
G={}; name=None
hdr={}
for line in open('/home/claude/U25.fa'):
    if line.startswith('>'):
        acc=line[1:].split()[0]; c=int(re.search(r'chromosome (\d+)',line).group(1))
        name=acc; G[acc]=[]; hdr[acc]=c
    else: G[name].append(line.strip())
G={k:''.join(v).upper() for k,v in G.items()}
bychr=defaultdict(list)
for f in final: bychr[f['chrom']].append(f)

OUT='/mnt/user-data/outputs/UTEX25_annotation'
os.makedirs(OUT,exist_ok=True)

def wrap_loc(loc):
    out=[];cur=''
    for tok in loc.replace(',',', ').split(' '):
        if not tok: continue
        if len(cur)+len(tok)>58: out.append(cur); cur=tok
        else: cur=cur+tok
    out.append(cur)
    return ('\n'+' '*21).join(x.rstrip() for x in out)

def qual(k,v,quote=True):
    s='/%s=%s'%(k,'"%s"'%v if quote else v)
    return '\n'.join(textwrap.wrap(s,58,initial_indent=' '*21,subsequent_indent=' '*21,
                                   break_long_words=True,break_on_hyphens=False))

def locstr(f):
    parts=['%d..%d'%(s,e) for s,e,p in sorted(f['blocks'])]
    l=parts[0] if len(parts)==1 else 'join(%s)'%','.join(parts)
    if f['strand']=='-': l='complement(%s)'%l
    return l

DATE='31-AUG-2026'
tot=0
for acc in sorted(G,key=lambda a:hdr[a]):
    c=hdr[acc]; seq=G[acc]; genes=sorted(bychr[acc],key=lambda x:x['s'])
    p=os.path.join(OUT,'Auxenochlorella protothecoides isolate UTEX 25 chromosome %d, whole genome shotgun sequence.gb'%c)
    with open(p,'w') as o:
        o.write('LOCUS       %-16s%9d bp    DNA     linear   CON %s\n'%(acc.split('.')[0],len(seq),DATE))
        o.write('DEFINITION  Auxenochlorella protothecoides isolate UTEX 25 chromosome %d, whole\n            genome shotgun sequence.\n'%c)
        o.write('ACCESSION   %s\n'%acc.split('.')[0])
        o.write('VERSION     %s\n'%acc)
        o.write('DBLINK      BioProject: PRJNA1328465\n')
        o.write('KEYWORDS    WGS.\n')
        o.write('SOURCE      Auxenochlorella protothecoides\n')
        o.write('  ORGANISM  Auxenochlorella protothecoides\n')
        o.write('            Eukaryota; Viridiplantae; Chlorophyta; core chlorophytes;\n')
        o.write('            Trebouxiophyceae; Chlorellales; Chlorellaceae; Auxenochlorella.\n')
        o.write('COMMENT     ##Genome-Annotation-Data-START##\n')
        o.write('            Annotation Provider   :: homology projection (in-house)\n')
        o.write('            Annotation Method     :: miniprot 0.18-r281, protein-to-genome\n')
        o.write('            Annotation Source     :: UTEX 250-A haplotype1 (PRJNA1195245) and\n')
        o.write('                                     haplotype2 (PRJNA1195244) proteins\n')
        o.write('            Assembly              :: UTEX 25 (PRJNA1328465)\n')
        o.write('            Note                  :: models are projected, NOT experimentally\n')
        o.write('                                     validated; no RNA-seq/IsoSeq evidence\n')
        o.write('            ##Genome-Annotation-Data-END##\n')
        o.write('FEATURES             Location/Qualifiers\n')
        o.write('     source          1..%d\n'%len(seq))
        o.write(qual('organism','Auxenochlorella protothecoides')+'\n')
        o.write(qual('mol_type','genomic DNA')+'\n')
        o.write(qual('isolate','UTEX 25')+'\n')
        o.write(qual('chromosome',str(c))+'\n')
        o.write(qual('db_xref','taxon:3075')+'\n')
        for f in genes:
            gl='%d..%d'%(f['s'],f['e'])
            if f['strand']=='-': gl='complement(%s)'%gl
            o.write('     gene            %s\n'%wrap_loc(gl))
            o.write(qual('locus_tag',f['lt'])+'\n')
            if f['sym']: o.write(qual('gene',f['sym'])+'\n')
            L=locstr(f)
            o.write('     mRNA            %s\n'%wrap_loc(L))
            o.write(qual('locus_tag',f['lt'])+'\n')
            o.write(qual('product',f['prod'])+'\n')
            o.write('     CDS             %s\n'%wrap_loc(L))
            o.write(qual('locus_tag',f['lt'])+'\n')
            if f['sym']: o.write(qual('gene',f['sym'])+'\n')
            o.write(qual('codon_start',str(f['phase']+1),quote=False)+'\n')
            o.write(qual('product',f['prod'])+'\n')
            o.write(qual('inference','similar to AA sequence:INSDC:%s'%f['src'])+'\n')
            flags=[]
            if '*' in f['prot']: flags.append('internal stop')
            if not f['prot'].startswith('M'): flags.append('no start Met')
            if len(f['nt'])%3: flags.append('length not multiple of 3')
            note='projected from UTEX 250-A %s (subgenome %s); miniprot identity %.4f'%(
                f['src'],f['srchap'],f['ident'])
            if flags: note+='; QC flags: '+', '.join(flags)
            o.write(qual('note',note)+'\n')
            o.write(qual('translation',f['prot'])+'\n')
        o.write('ORIGIN      \n')
        for i in range(0,len(seq),60):
            chunk=seq[i:i+60].lower()
            o.write('%9d %s\n'%(i+1,' '.join(chunk[j:j+10] for j in range(0,len(chunk),10))))
        o.write('//\n')
    tot+=os.path.getsize(p)
    print('chr%-3d %-14s %8d bp  genes=%-5d  %6.1f MB'%(c,acc,len(seq),len(genes),os.path.getsize(p)/1e6))
print('total %.1f MB'%(tot/1e6))

# Builds LaTeX tables for Section S9 from the *_summary_table.txt files written by report_v2_E*.R.
# Usage: python3 gen_tables.py <results dir> <out.tex>
import sys, re
D = sys.argv[1]; out = open(sys.argv[2], 'w')
MOD = {'A':'(M1)','B':'(M2)',"C'":'(M3)','C':'(M3)','C(N(0,I))':'(M3)','C2':'(M4)','D':'(M5)','E':'(M6)'}
def lines(f): return [l.rstrip('\n') for l in open(f'{D}/{f}', encoding='utf-8')]
def blocks(L, first):
    """split a wrapped R table into blocks: each block = header + rows; first = first header token"""
    B=[]; cur=None
    for l in L:
        t=l.split()
        if not t: 
            if cur: B.append(cur); cur=None
            continue
        if t[0]==first or (cur and len(t) and not re.match(r"^[-0-9.!]+$", t[0]) and t[0] not in MOD and not l.startswith(' ') ):
            pass
        if cur is None or (t[0] in ('rf','MST','MCA') and l.startswith(' ')) and cur and all(re.match(r"^[-0-9.!]+$",x) for x in cur[-1].split()[:1]):
            pass
    return B
def cell(x):
    x=x.replace('!','*')
    if x=='-': return '--'
    return x.replace('*','$^*$')
def tab(caption,label,head,rows,align=None,size='\\footnotesize'):
    n=len(head)
    align=align or 'l'+'c'*(n-1)
    out.write('\\begin{table}[p]\n\\renewcommand{\\baselinestretch}{1}'+size+'\n\\caption{'+caption+'\\label{'+label+'}}\n\\begin{center}\n\\setlength{\\tabcolsep}{3pt}\n'+('\\resizebox{\\textwidth}{!}{%\n' if n>8 else ''))
    out.write('\\begin{tabular}{'+align+'}\n\\toprule\n'+' & '.join(head)+' \\\\\n\\midrule\n')
    for r in rows: out.write(' & '.join(r)+' \\\\\n')
    out.write('\\bottomrule\n\\end{tabular}'+('}' if n>8 else '')+'\n\\end{center}\n\\end{table}\n\n')

# ---------- E1: three wrapped blocks of 24 rows
L=[l for l in lines('E1_summary_table.txt')[1:] if l.strip()]
hdr_idx=[i for i,l in enumerate(L) if l.split()[0] in ('model','rf','MCA')]
b=[L[hdr_idx[j]+1:(hdr_idx[j+1] if j+1<len(hdr_idx) else len(L))] for j in range(len(hdr_idx))]
r1=[];r2=[]
for x,y,z in zip(b[0],b[1],b[2]):
    t=x.split()
    # model may contain spaces? 'C(N(0,I))' no spaces
    model,k,n,d=t[0],t[1],t[2],t[3]; vals=t[6:]  # skip reps and range
    y=y.split(); z=z.split()
    lab=MOD[model]; nn='equal' if n=='300' else '150, 300, 600'
    r1.append([lab,k,nn,d]+[cell(v) for v in [vals[0]]+[z[0],z[1]]+vals[1:4]+[y[0],y[1],y[2],y[3],y[4]]])
    r2.append([lab,k,nn,d]+[cell(v) for v in [y[5],y[7],y[8],y[9],y[10],y[11]]])
tab('Empirical sizes of MCA-kDET and its variants under $H_0$ (all settings). Columns: MCA-kDET at $\\alpha = 0.05$, 0.01 and 0.10; single learners; two-sided (2s), sum-type, combined and pairwise variants, at $\\alpha = 0.05$. \\Qnote{200 replications per setting; asterisks mark values outside $[0.020, 0.080]$}',
    'tab:S-size1',['Null','$k$','$n$','$d$','MCA','$\\alpha{=}.01$','$\\alpha{=}.10$','l1','quad','knn','rf','2s','Sum','Comb','Pair'],r1)
tab('Empirical sizes at $\\alpha = 0.05$ of the shared network, the multiclass accuracy test and the competitors (all settings). \\Qnote{200 replications per setting}',
    'tab:S-size2',['Null','$k$','$n$','$d$','Shared','Accuracy','DISCO','MMD','MST','Schott'],r2)

# ---------- E2: power, two blocks
L=[l for l in lines('E2_summary_table.txt')[1:]]
L=L[:L.index('')] if '' in L else L
L=[l for l in L if l.strip()]
h=[i for i,l in enumerate(L) if l.split()[0] in ('model','MST')]
b1=L[h[0]+1:h[1]]; b2=L[h[1]+1:]
rows=[]
for x,y in zip(b1,b2):
    t=x.split(); y=y.split()
    lab=MOD[t[0]]; beta=t[1]
    par=f'$s = {round(float(beta)*100)}$' if t[0]=='D' else f'$\\beta = {beta}$'
    v=t[3:]  # l1 quad knn rf MCA sum comb pair shared macc DISCO MMD
    rows.append([lab,par]+[v[4],v[0],v[1],v[2],v[3],v[5],v[6],v[7],v[8],v[9],v[10],v[11],y[0],y[1]])
tab('Power at $\\alpha = 0.05$ in models (M1)--(M5), $k = 3$ (the data of Figure~\\ref{M-fig:power} and Table~\\ref{tab:variants}). \\Qnote{100 replications per setting}',
    'tab:S-power',['Model','Setting','MCA','l1','quad','knn','rf','Sum','Comb','Pair','Shared','Accuracy','DISCO','MMD','MST','Schott'],rows)

# ---------- E4: model E (M6)
L=[l for l in lines('E4_summary_table.txt')[1:] if l.strip()]
h=[i for i,l in enumerate(L) if l.split()[0] in ('beta','MMD')]
b1=L[h[0]+1:h[1]]; b2=L[h[1]+1:]; b2=[l for l in b2 if re.match(r'^\s*[0-9.]+\s+[0-9.]+\s*$',l)]
rows=[]
for x,y in zip(b1,b2):
    t=x.split(); y=y.split()
    v=t[2:]  # l1 quad knn rf MCA 2s sum comb shared macc DISCO
    rows.append([t[0],v[4],v[0],v[1],v[2],v[3],v[5],v[6],v[7],v[8],v[9],v[10],y[0],y[1]])
tab('Rejection rates at $\\alpha = 0.05$ in model (M6), $k = 3$, $d = 100$; $\\beta = 0$ is $H_0$. \\Qnote{200 replications at $\\beta = 0$, 100 otherwise}',
    'tab:S-M6',['$\\beta$','MCA','l1','quad','knn','rf','2s','Sum','Comb','Shared','Accuracy','DISCO','MMD','MST'],rows)

# ---------- E6: dimension
L=[l for l in lines('E6_summary_table.txt') if l.strip()]
h=[i for i,l in enumerate(L) if l.split()[0] in ('model','MMD')]
b1=L[h[0]+1:h[0]+9]; b2=L[h[1]+1:h[1]+9]
rows=[]
for x,y in zip(b1,b2):
    t=x.split(); y=y.split()
    v=t[5:]  # l1 quad knn rf MCA sum comb shared macc DISCO
    rows.append([MOD[t[0]],t[1],v[4],v[0],v[1],v[2],v[3],v[5],v[6],v[7],v[8],v[9],y[0],y[1],y[2]])
tab('Power at $\\alpha = 0.05$ against the dimension, $k = 3$, $\\beta = 0.8$; model (M3) has ten dependent pairs (five when $d = 10$). \\Qnote{100 replications per setting}',
    'tab:S-dim',['Model','$d$','MCA','l1','quad','knn','rf','Sum','Comb','Shared','Accuracy','DISCO','MMD','MST','Schott'],rows)

# ---------- E5: tau and multiple splits
L=[l for l in lines('E5_summary_table.txt')[1:] if l.strip()]
h=[i for i,l in enumerate(L) if l.split()[0] in ('model','2x')]
b1=L[h[0]+1:h[1]]; b2=L[h[1]+1:h[1]+13]
rows=[]
for x,y in zip(b1,b2):
    t=x.split(); y=y.split()
    rows.append([MOD[t[0]],t[1],t[2],t[4],t[5],t[6],t[7],t[8],t[9],t[10],y[0],y[1]])
tab('Rejection rates at $\\alpha = 0.05$ for split fractions $\\tau = 0.3$, 0.5, 0.7, $k = 3$; (M1) with $d = 100$ and (M3) with $d = 20$; $\\beta = 0$ is $H_0$. The last two columns combine five splits by twice the median and twice the mean of the p-values (Remark~\\ref{M-rem:split}). \\Qnote{200 replications under $H_0$, 100 otherwise}',
    'tab:S-tau',['Model','$\\beta$','$\\tau$','MCA','l1','quad','knn','rf','Sum','Comb','$2\\times$median','$2\\times$mean'],rows)

# ---------- E7: time
L=lines('E7_summary_table.txt')
rows=[]
for l in L:
    m=re.match(r'^(\s*)(.*?)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*$',l)
    if m and 'k=3' not in l:
        name=m.group(2).strip()
        if m.group(1): name='\\quad '+name
        rows.append([name,m.group(3),m.group(4),m.group(5)])
tab('Mean computing time in seconds per replication, model (M1), $d = 100$, $n = 300$ per population, $B = 200$, one core per replication. \\Qnote{5 replications per $k$; to be re-measured on an idle machine}',
    'tab:S-time',['Component','$k = 3$','$k = 5$','$k = 10$'],rows,align='lccc',size='\\small')
out.close()

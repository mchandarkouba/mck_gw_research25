inpath = "/Users/mck/Desktop/palmese_research/corr_scripts/corr_output/gw_dataFrame_clean.csv"
outpath = "/Users/mck/Desktop/palmese_research/corr_scripts/gw_dataFrame_overleaf.txt"

###############################################################################

def formatted(s0:str):
    s1 = ''
    for i,c in enumerate(s0):
        if c in ('_',): 
            s1 = s1 + f"$\{c}$"
        else:
            s1 = s1 + c
    return s1

###############################################################################

with open(inpath, 'r') as csv:
    lines = csv.readlines()
        
cells = [line.split(',') for line in lines]
header = cells[0]

cinclude = ['grace_id', 
            'mass_1_source', 
            'mass_2_source', 
            'final_mass_source_upper',
            ]
rexclude = {"grace_id": [#"GW231123_135430", 
                         #"GW200308_173609", 
                         #"GW200322_091133", 
                         #"GW190929_012149",
                         ]
           }

inds_ex = []
for k in rexclude:
    c = header.index(k)
    col = [row[c] for row in cells]
    inds_ex += [col.index(x) for x in rexclude[k]]

inds_in = {header.index(x) : n
           for n,x in enumerate(cinclude)
           }

finlist = []
for r,row in enumerate(cells):
    if r in inds_ex: continue

    rowlist = ['',]*len(cinclude)
    for c,cell in enumerate(row):
        if c in inds_in: rowlist[inds_in[c]] = formatted(cell)
    
    rowstr = (" & ").join(rowlist)
    if r==0:
        rowstr = "\hline \\\\" + rowstr + "\\\\ \hline"
    finlist += [rowstr]

finstr = (" \\\\ ").join(finlist)

with open(outpath, 'w') as txt:
    txt.write(finstr)
    txt.close()
            
            
        
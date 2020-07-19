# %% 
import numpy as np
import glob 

prefix = glob.glob("*.export")[0].split('.')[0]

# read k list
klist = []
with open("klist.dat") as f:
    for line in f.readlines():
        if line.strip() == "":
            continue
        if len(line.strip().split()) == 3:
            klist.append(np.array(line.strip().split(),dtype=float))

# %%
for ki, k in enumerate(klist):
    # read g list
    glist = []
    with open(f"{prefix}.export/grid.{ki+1}") as f:
        for line in f.readlines():
            if line.strip() == "":
                continue
            if not line.strip().startswith("<") and len(line.strip().split()) == 3:
                glist.append(np.array(line.strip().split(),dtype=float))
    # %%
    nbnd = 0
    ng=0
    u_k = np.zeros([len(glist),144],dtype=complex)
    with open(f"{prefix}.export/wfc.{ki+1}") as f:
        for line in f.readlines():
            if line.strip() == "":
                continue
            if line.strip().startswith("<Wfc"):
                nbnd += 1
                ng=0
            elif nbnd > 0 and len(line.strip().split(","))==2:
                r,i = np.array(line.strip().split(","),dtype=float)
                u_k[ng,nbnd-1]=r+1j*i
                ng+=1
                
    print("# kpoint",ki)
    print("#\t",k)
    for gi,g in enumerate(glist):
        print(*g, *u_k[gi])


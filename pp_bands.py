#%%
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from glob import glob

# %%
# IO functions
def read_bands(fname):
    E_ki = []
    E_i =[]
    with open(fname) as f:
        k_ind=-1
        for line in f:
            line_s = line.split()
            if line.strip().startswith('k = '):
                k_ind +=1
                if len(E_i) > 0:
                    E_ki.append(np.array(E_i))
                E_i =[]                
            elif len(line_s)==8:
                E_i=np.concatenate((E_i,np.asarray(line_s,dtype=float)))
                
    if len(E_i) >0:
        E_ki.append(E_i)
    return np.stack(E_ki)


def read_klist(fname):
    def extract_k(s):
        k1,k2,k3 = s[s.find('=')+1:s.find('(')].split()
        return np.array([float(k1),float(k2),float(k3)])
    k_list = []
    
    with open(fname) as f:
        for line in tqdm(f):
            if line.strip().startswith('k = '):
                k_list.append(extract_k(line))
            
    return k_list

out_file = 'bands.out'
E_ki=read_bands(out_file)
klist=read_klist(out_file)

# write klist
f_klist = open("klist.dat",'w')

for ki, k in enumerate(klist):
    print(*k,2/27)
    f_klist.write(f"{k[0]} {k[1]} {k[2]}\n")

    for i,E in enumerate(E_ki[ki]):
        print(i,E)

f_klist.close()

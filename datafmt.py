# %% 
import numpy as np
from tqdm import tqdm

# %%

def read_bands(fname):
    E_ki = []
    E_i =[]
    with open(fname) as f:
        k_ind=-1
        for line in f:
            line_s = line.split()
            if len(line_s) == 4:
                k_ind +=1
                if len(E_i) > 0:
                    E_ki.append(np.array(E_i))
                E_i =[]                
            elif len(line_s)==2:
                E_i.append(float(line_s[-1]))
    if len(E_i) >0:
        E_ki.append(E_i)
    return np.stack(E_ki)

def read_u(fname):
    def extract_k(s):
        k1,k2,k3 = s[s.find('[')+1:s.find(']')].split()
        return np.array([float(k1),float(k2),float(k3)])
    k_list = []
    g_list = []
    u_kgi = []
    g_k = []
    u_k = []
    
    with open(fname) as f:
        for line in tqdm(f):
            if line.startswith('#') and 'kpoint' not in line:
                k_list.append(extract_k(line))
                
                if len(g_k) > 0:
                    g_list.append(g_k)
                    u_kgi.append(np.stack(u_k))
                    g_k=[]
                    u_k=[]
            elif '#' not in line:
                line_s=line.split()
                g_k.append([int(float(x)) for x in line_s[:3]])
                u_k.append([complex(x) for x in line_s[3:]])
        if len(g_k) > 0:
            g_list.append(g_k)
            u_kgi.append(np.stack(u_k))
    return k_list,g_list,u_kgi

# %%

k_list,g_list,u_kgi= read_u('u_i.dat')
E_ki = read_bands('bands.dat') 

# %%

# %%
with open('Ek.dat','w') as f:
    lines = []
    for E_i in E_ki:
        for E in E_i:
            lines.append(f"{E}\n")
    lines[-1]=lines[-1].strip()
    f.writelines(lines)

with open('klist.dat','w') as f:
    lines = []
    for k in k_list:
        lines.append(f"{k[0]} {k[1]} {k[2]}\n")
    lines[-1]=lines[-1].strip()
    f.writelines(lines)
# %%
with open('glist.dat','w') as f:
    lines=[]
    for g_k in g_list:
        lines.append(f"{len(g_k)}\n")
        for g in g_k:
            lines.append(f"{g[0]} {g[1]} {g[2]}\n")
    lines[-1]=lines[-1].strip()
    f.writelines(lines)


# %%
with open('ui.dat','w') as f:
    lines = []
    for u_gi in u_kgi:
        for u_i in u_gi:
            for u in u_i:
                lines.append(f"({u.real},{u.imag})\n")
    lines[-1]=lines[-1].strip()
    f.writelines(lines)


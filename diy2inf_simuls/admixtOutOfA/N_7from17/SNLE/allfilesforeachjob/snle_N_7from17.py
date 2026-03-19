#!/usr/bin/env python
# coding: utf-8

# In[12]:


import os
import torch, numpy as np
from sbi.inference import NLE
from sbi.utils import BoxUniform
import time
import pickle
import sys
if hasattr(sys,"ps1"): # interactive run
    fromrep = 0
    torep = 1
else:
    fromrep = int(sys.argv[1])-1   
    torep = int(sys.argv[2])   
    


# In[13]:


from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    workdir = "".join([os.getcwd(),"/"]) # subdir of "/home/frousset/Infusionplus/sbi/paper-infusion-SNLE/admixture/"
    diyabccall = './diyabc-RF-linux-v1.1.36-o -p ./ -R "ALL" -t 10 -m -o params.txt'
elif platform == "darwin":
    # OS X 
    workdir = "/Users/marin/Desktop/paper-infusion-SNLE/admixture/"
    diyabccall = "diyabc-RF-mac -p /Users/marin/Desktop/paper-infusion-SNLE/admixture/ -o params.txt"
elif platform == "win32":
    # Windows... '32' is confusing, but anyway
    workdir = "d:/home/francois/travail/stats/Infusionplus/sbi/paper-infusion-SNLE/admixture/"
    # cf https://github.com/diyabc/diyabc for the optional arguments; -t <thread nbr> -m for parallelization
    diyabccall = 'diyabc-RF-windows-v1.1.36-o.exe -p ./ -R "ALL" -t 10 -m -o params.txt'

os.chdir(workdir)


# In[14]:


torch.manual_seed(123)
device = torch.device("cpu")

log10 = lambda x: torch.log10(torch.tensor(x, device=device))

low = torch.tensor([
    3.0,            # log10_N2
    1.0,            # t1
    log10(50.0),    # log10_t12
    0.0,            # log10_(1+t23)
    5.0,            # Nbn34
    0.0,            # log10_(1+t34)
    2.0             # log10_Na
], device=device)

high = torch.tensor([
    5.0,            # log10_N2
    30.0,           # t1
    log10(5000.0),  # log10_t12
    log10(5001.0),  # log10_(1+t23)
    500.0,          # Nbn34
    log10(5001.0),  # log10_(1+t34)
    4.0             # log10_Na
], device=device)


prior = BoxUniform(low=low, high=high, device=device)


# In[15]:


def simulator(theta):
    batch = theta.shape[0]
    device = theta.device    
    N2    = 10.0 ** theta[:, 0]
    t1    = theta[:, 1]
    t12   = 10.0 ** theta[:, 2]
    t23   = 10.0 ** theta[:, 3] - 1.0
    Nbn34 = theta[:, 4]
    t34   = 10.0 ** theta[:, 5] - 1.0
    Na    = 10.0 ** theta[:, 6]
    
    # N2, t1, t12, t23, Nbn34, t34, Na = theta.T
    Scenario = torch.full((batch,), 1, device=device)
    N1   = torch.full((batch,), 12589.0, device=device)
    N3   = torch.full((batch,), 19953.0, device=device)
    N4   = torch.full((batch,), 3162.0, device=device)
    ra   = torch.full((batch,), 0.2, device=device)
    d3   = torch.full((batch,), 42.0, device=device)
    Nbn3 = torch.full((batch,), 160.0, device=device)
    d4   = torch.full((batch,), 9.0, device=device)
    Nbn4 = torch.full((batch,), 98.0, device=device)
    N34  = torch.full((batch,), 1259.0, device=device)
    d34  = torch.full((batch,), 24.0, device=device)
    theta_full = torch.stack([Scenario, N1, N2, N3, N4, t1, ra, t12,
    d3, Nbn3, d4, Nbn4, N34, t23, d34, Nbn34, t34, Na],dim=1)
    theta_np = theta_full.detach().cpu().numpy()
    theta_np[:, 0] = theta_np[:, 0].astype(int)
    # print(theta_full.size())
    np.savetxt("".join([workdir,"params.txt"]),
               theta_np, fmt="%d " + " ".join(["%.6f"]* (theta_np.shape[1]-1)))
    if os.path.isfile("statfile.txt"):
        os.remove("statfile.txt")
    chk = os.system(diyabccall)
    if chk != 0:
        print("diyabc simulator returned an error")
    x_np = np.loadtxt("".join([workdir,"statfile.txt"]))
    x_np = x_np[:, 1:]
    x = torch.from_numpy(x_np).to(theta)
    # print(x.size())
    return x

def estimate(x0, prior, rounds, num_sims):
    inference = NLE(prior=prior, density_estimator="maf", device=device)
    proposal = prior
    for r in range(rounds):
        print(f"  Round {r}")
        if r == 0:
            theta = prior.sample((num_sims,))
        else:
            theta = proposal.sample((num_sims,), x=x0, show_progress_bars=False)
        # print(theta.size())
        x = simulator(theta)
        # print(x.size())
        lf = inference.append_simulations(theta, x).train(
            training_batch_size=128,
            show_train_summary=False
        )
        posterior = inference.build_posterior(lf)
        proposal  = posterior.set_default_x(x0)
    samples = posterior.sample((10000,), x=x0, show_progress_bars=False)
    samples[:, 0] = 10.0 ** samples[:, 0]
    samples[:, 2] = 10.0 ** samples[:, 2]
    samples[:, 3] = 10.0 ** samples[:, 3] - 1.0
    samples[:, 5] = 10.0 ** samples[:, 5] - 1.0
    samples[:, 6] = 10.0 ** samples[:, 6]    
    posterior_mean = samples[1:1000,:].mean(dim=0)
    posterior_q1 = samples[1:1000,:].quantile(0.025, dim=0)
    posterior_q2 = samples[1:1000,:].quantile(0.975, dim=0)   
    posterior_mean_10K = samples.mean(dim=0)
    posterior_q1_10K = samples.quantile(0.025, dim=0)
    posterior_q2_10K = samples.quantile(0.975, dim=0)   
    return posterior_mean, posterior_q1, posterior_q2, posterior_mean_10K, posterior_q1_10K, posterior_q2_10K, inference


# In[16]:


x0_np = np.loadtxt("".join([workdir,"obs.txt"]))
# n_rep = x0_np.shape[0]

n_rep = 1
theta_dim = 7

posterior_means = torch.empty((n_rep, theta_dim), device=device)
posterior_q1s = torch.empty((n_rep, theta_dim), device=device)
posterior_q2s = torch.empty((n_rep, theta_dim), device=device)

posterior_means_10K = torch.empty((n_rep, theta_dim), device=device)
posterior_q1s_10K = torch.empty((n_rep, theta_dim), device=device)
posterior_q2s_10K = torch.empty((n_rep, theta_dim), device=device)


# In[17]:


for rep in range(fromrep,torep):
    print(f"\n=== Replicate {rep+1} ===")
    start = time.time()
    x0 = torch.tensor(x0_np[rep], dtype=torch.float32, device=device).unsqueeze(0)
    mean_rep, q1_rep, q2_rep, mean_rep_10k, q1_rep_10k, q2_rep_10k, inference = estimate(x0, prior, 10, 2000)
    posterior_means[0, :] = mean_rep
    posterior_q1s[0, :] = q1_rep
    posterior_q2s[0, :] = q2_rep
    posterior_means_10K[0, :] = mean_rep_10k
    posterior_q1s_10K[0, :] = q1_rep_10k
    posterior_q2s_10K[0, :] = q2_rep_10k
    end = time.time()
    print("Time spent:", end - start, "seconds")
    with open(".".join(["inference",str(rep),"pkl"]), "wb") as handle:
        pickle.dump(inference, handle)
    np.savetxt(".".join(["posterior_means",str(rep),"csv"]), posterior_means.cpu().numpy(), delimiter=";")
    np.savetxt(".".join(["posterior_q1s",str(rep),"csv"]), posterior_q1s.cpu().numpy(), delimiter=";")
    np.savetxt(".".join(["posterior_q2s",str(rep),"csv"]), posterior_q2s.cpu().numpy(), delimiter=";")
    np.savetxt(".".join(["posterior_means_10K",str(rep),"csv"]), posterior_means_10K.cpu().numpy(), delimiter=";")
    np.savetxt(".".join(["posterior_q1s_10K",str(rep),"csv"]), posterior_q1s_10K.cpu().numpy(), delimiter=";")
    np.savetxt(".".join(["posterior_q2s_10K",str(rep),"csv"]), posterior_q2s_10K.cpu().numpy(), delimiter=";")


# In[23]:


true_params = np.array([3162, 6, 158, 31, 65, 999, 501], dtype=np.float32)
arr = np.round(posterior_means.cpu().numpy(), 2).ravel()
true = true_params.ravel()

print("=== Comparaison posterior_mean vs true_params ===\n")
for a, t in zip(arr, true):
    print(f"{a:>8.2f}   {t:>8.2f}")


# In[ ]:





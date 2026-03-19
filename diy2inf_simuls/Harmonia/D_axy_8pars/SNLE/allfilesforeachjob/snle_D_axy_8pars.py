#!/usr/bin/env python
# coding: utf-8

# In[2]:


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
    workdir = "".join([os.getcwd(),"/"]) # subdir of "/home/frousset/Infusionplus/sbi/paper-infusion-SNLE/D_axy_8pars/"
    diyabccall = './diyabc-RF-linux-v1.1.36-o -p ./ -R "ALL" -t 10 -m -o params.txt'
elif platform == "darwin":
    # OS X 
    workdir = "/Users/marin/Desktop/paper-infusion-SNLE/D_axy_8pars/"
    diyabccall = "diyabc-RF-mac -p /Users/marin/Desktop/paper-infusion-SNLE/D_axy_8pars/ -o params.txt"
elif platform == "win32":
    # Windows... '32' is confusing, but anyway
    workdir = "d:/home/francois/travail/stats/Infusionplus/sbi/paper-infusion-SNLE/D_axy_8pars/"
    # cf https://github.com/diyabc/diyabc for the optional arguments; -t <thread nbr> -m for parallelization
    diyabccall = 'diyabc-RF-windows-v1.1.36-o.exe -p ./ -R "ALL" -t 10 -m -o params.txt'

os.chdir(workdir)


# In[14]:


torch.manual_seed(123)
device = torch.device("cpu")

log10 = lambda x: torch.log10(torch.tensor(x, device=device))

low = torch.tensor([
    -2.0,            # logTh1
    -2.0,            # logTh2
    -2.0,            # logTh3
    -2.0,            # logTh4
    3.0,            # T1
    0.05,            # ar
    -5.0,            # logMu
    0.01,            # MEANP
], device=device)

high = torch.tensor([
    1.0,            # logTh1
    1.0,            # logTh2
    1.0,            # logTh3
    1.0,            # logTh4
    45.0,            # T1
    0.95,            # ar
    -2.0,            # logMu
    0.5,            # MEANP
], device=device)


prior = BoxUniform(low=low, high=high, device=device)


# In[21]:


def simulator(theta):
    batch = theta.shape[0]
    device = theta.device    
    T1    = theta[:, 4]
    ar   = theta[:, 5]
    MEANMU   = 10.0 ** theta[:, 6]
    MEANP = theta[:, 7]
    Th1    = 10.0 ** theta[:, 0] / MEANMU
    Th2    = 10.0 ** theta[:, 1] / MEANMU
    Th3    = 10.0 ** theta[:, 2] / MEANMU
    Th4    = 10.0 ** theta[:, 3] / MEANMU
    
    Scenario = torch.full((batch,), 1, device=device)
    theta_full = torch.stack([Scenario, Th1, Th2, Th3, Th4, T1, ar, MEANMU, MEANP],dim=1)
    theta_np = theta_full.detach().cpu().numpy()
    theta_np[:, 0] = theta_np[:, 0].astype(int)
    # print(theta_np[0,:])
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


# In[23]:


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
    posterior_mean = samples[1:1000,:].mean(dim=0)
    posterior_q1 = samples[1:1000,:].quantile(0.025, dim=0)
    posterior_q2 = samples[1:1000,:].quantile(0.975, dim=0)   
    posterior_mean_10K = samples.mean(dim=0)
    posterior_q1_10K = samples.quantile(0.025, dim=0)
    posterior_q2_10K = samples.quantile(0.975, dim=0)   
    return posterior_mean, posterior_q1, posterior_q2, posterior_mean_10K, posterior_q1_10K, posterior_q2_10K, inference


# In[16]:


x0_np = np.loadtxt("".join([workdir,"S_obs_table.txt"]))
# n_rep = x0_np.shape[0]

n_rep = 1
theta_dim = 8

posterior_means = torch.empty((n_rep, theta_dim), device=device)
posterior_q1s = torch.empty((n_rep, theta_dim), device=device)
posterior_q2s = torch.empty((n_rep, theta_dim), device=device)

posterior_means_10K = torch.empty((n_rep, theta_dim), device=device)
posterior_q1s_10K = torch.empty((n_rep, theta_dim), device=device)
posterior_q2s_10K = torch.empty((n_rep, theta_dim), device=device)


# In[17]:


# In[ ]:


for rep in range(fromrep,torep):
    print(f"\n=== Replicate {rep+1} ===")
    start = time.time()
    x0 = torch.tensor(x0_np[rep], dtype=torch.float32, device=device).unsqueeze(0)
    mean_rep, q1_rep, q2_rep, mean_rep_10k, q1_rep_10k, q2_rep_10k, inference = estimate(x0, prior, 10, 2300)
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


# In[25]:


true_params = np.array([0.00, 0.00, 0.00, 0.00, 13.00, 0.50, -3.00, 0.25 ], dtype=np.float32)
arr = np.round(posterior_means.cpu().numpy(), 2).ravel()
true = true_params.ravel()

print("=== Comparaison posterior_mean vs true_params ===\n")
for a, t in zip(arr, true):
    print(f"{a:>8.2f}   {t:>8.2f}")


# In[ ]:


# In[ ]:





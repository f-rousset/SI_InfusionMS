#!/usr/bin/env python
# coding: utf-8

# In[3]:


# --------------------------------------------------------
# SNLE for 15-parameters Cholesky example 
# --------------------------------------------------------

import os
import time
import pickle
import torch, numpy as np
from sbi.inference import NLE
from sbi.utils import BoxUniform

torch.manual_seed(123)
device = torch.device("cpu")


# In[4]:


import sys
if not hasattr(sys,"ps1"): # interactive run
    fromrep = int(sys.argv[1])-1   
    torep = int(sys.argv[2])   
else:
    fromrep = 0
    torep = 1


# In[96]:


from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    workdir = "".join([os.getcwd(),"/"]) # subdir of "/home/frousset/Infusionplus/sbi/paper-infusion-SNLE/covariance/"
elif platform == "darwin":
    # OS X 
    workdir = "/Users/marin/Desktop/paper-infusion-SNLE/covariance/"
elif platform == "win32":
    # Windows... '32' is confusing, but anyway
    workdir = "d:/home/francois/travail/stats/Infusionplus/sbi/paper-infusion-SNLE/covariance/"
    # cf https://github.com/diyabc/diyabc for the optional arguments; -t <thread nbr> -m for parallelization

os.chdir(workdir)


# In[5]:


dim = 5
theta_dim = dim * (dim+1) // 2
n_obs = 50

def upper_tri_indices(d):
    idx = []
    for i in range(d):
        for j in range(i, d):
            idx.append((i, j))
    return idx

UT_IDX = upper_tri_indices(dim)

low  = torch.empty(theta_dim, device=device)
high = torch.empty(theta_dim, device=device)

k = 0
for i in range(dim):
    for j in range(i, dim):
        if i == j:
            low[k]  = 0.5
            high[k] = 6.0
        else:
            low[k]  = -4.0
            high[k] = 4.0
        k += 1

prior = BoxUniform(low=low, high=high, device=device)

def theta_to_U(theta):
    batch = theta.shape[0]
    U = torch.zeros(batch, dim, dim, device=theta.device)
    for k,(i,j) in enumerate(UT_IDX):
        U[:,i,j] = theta[:,k]
    return U

def cov_to_summaries(C):
    stats = []
    for (i,j) in UT_IDX:
        stats.append(C[:,i,j])
    return torch.stack(stats, dim=1)

def simulator(theta):
    batch = theta.shape[0]
    U = theta_to_U(theta)
    V = U.transpose(1,2) @ U
    loc = torch.zeros(batch, dim, device=theta.device)
    mvn = torch.distributions.MultivariateNormal(loc=loc, covariance_matrix=V)
    X = mvn.rsample((n_obs,))
    X = X.permute(1,0,2)
    Xc = X - X.mean(dim=1,keepdim=True)
    C = (Xc.transpose(1,2) @ Xc)/(n_obs-1)
    return cov_to_summaries(C)


def estimate(x0, prior, rounds, num_sims):
# def estimate(theta_star, prior, rounds, num_sims):
#   x0 = simulator(theta_star).squeeze(0)
    inference = NLE(prior=prior, density_estimator="maf", device=device)
    proposal = prior
    for r in range(rounds):
        print(f"  Round {r}")
        if r == 0:
            theta = prior.sample((num_sims,))
        else:
            theta = proposal.sample((num_sims,), x=x0, show_progress_bars=False)
        x = simulator(theta)
        lf = inference.append_simulations(theta, x).train(
            training_batch_size=128,
            show_train_summary=False)
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


# In[ ]:


x0_np = np.loadtxt("".join([workdir,"obs.txt"]))
n_rep = x0_np.shape[0]

n_rep = 1

posterior_means = torch.empty((n_rep, theta_dim), device=device)
posterior_q1s = torch.empty((n_rep, theta_dim), device=device)
posterior_q2s = torch.empty((n_rep, theta_dim), device=device)

posterior_means_10K = torch.empty((n_rep, theta_dim), device=device)
posterior_q1s_10K = torch.empty((n_rep, theta_dim), device=device)
posterior_q2s_10K = torch.empty((n_rep, theta_dim), device=device)

for rep in range(fromrep,torep):
    print(f"\n=== Replicate {rep+1} ===")
    start = time.time()
    x0 = torch.tensor(x0_np[rep], dtype=torch.float32, device=device).unsqueeze(0)
    mean_rep, q1_rep, q2_rep, mean_rep_10k, q1_rep_10k, q2_rep_10k, inference = estimate(x0, prior, 10, 4400)
    end = time.time()
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


# In[10]:


true_params =  torch.tensor([2.603, -1.689, -0.1089, 1.280, 2.528, 
                          3.670, -0.1172, -1.727, 0.5490, 3.522, 1.690, 
                          0.2382, 2.578, -1.048, 2.592], 
                          device=device).unsqueeze(0)

arr = np.round(posterior_means.cpu().numpy(), 2).ravel()
true = true_params.ravel()

print("=== Comparaison posterior_mean vs true_params ===\n")
for a, t in zip(arr, true):
    print(f"{a:>8.2f}   {t:>8.2f}")


# In[ ]:





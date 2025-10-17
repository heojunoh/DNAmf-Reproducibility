import GPy
import numpy as np
np.bool = np.bool_ # use this command if your numpy >1.23.1
from matplotlib import pyplot as plt
import matplotlib.mlab as ml
import matplotlib.patches as mpatches
import scipy.stats as stats
import rpy2.robjects as robjects

import time

''' load data '''
r = robjects.r
loaded_data = r.readRDS("tmp_data.rds")
X1=np.array(loaded_data[0])
X2=np.array(loaded_data[1])
X3=np.array(loaded_data[2])
Y1=np.array(loaded_data[3])
Y2=np.array(loaded_data[4])
Y3=np.array(loaded_data[5])
Xtest=np.array(loaded_data[6])
Ytest=np.array(loaded_data[7])

dim = Xtest.shape[1]
ntest = Xtest.shape[0]
active_dimensions = np.arange(0,dim)
max_iter = 10 # for optimization
MCsample = 100 # number of MC samples

''' Train level 1 '''
start = time.time()
k1 = GPy.kern.RBF(dim, ARD = True)
m1 = GPy.models.GPRegression(X=X1, Y=Y1, kernel=k1)

m1[".*Gaussian_noise"] = m1.Y.var()*0.01
m1[".*Gaussian_noise"].fix()

m1.optimize(max_iters = max_iter)

m1[".*Gaussian_noise"].unfix()
m1[".*Gaussian_noise"].constrain_positive()

m1.optimize_restarts(10, optimizer = "bfgs",  max_iters = max_iter)

mu1, v1 = m1.predict(X2)


''' Train level 2 '''
XX = np.hstack((X2, mu1))

k2 = GPy.kern.RBF(1, active_dims = [dim])*GPy.kern.RBF(dim, active_dims = active_dimensions, ARD = True) \
+ GPy.kern.RBF(dim, active_dims = active_dimensions, ARD = True)

m2 = GPy.models.GPRegression(X=XX, Y=Y2, kernel=k2)

m2[".*Gaussian_noise"] = m2.Y.var()*0.01
m2[".*Gaussian_noise"].fix()

m2.optimize(max_iters = max_iter)

m2[".*Gaussian_noise"].unfix()
m2[".*Gaussian_noise"].constrain_positive()

m2.optimize_restarts(10, optimizer = "bfgs",  max_iters = max_iter)


# Prepare for level 3: sample f_1 at X3
mu0, C0 = m1.predict(X3, full_cov=True)
Z = np.random.multivariate_normal(mu0.flatten(),C0,MCsample)
tmp_m = np.zeros((MCsample,X3.shape[0]))
tmp_v = np.zeros((MCsample,X3.shape[0]))

# push samples through f_2
for i in range(0,MCsample):
 mu, v = m2.predict(np.hstack((X3, Z[i,:][:,None])))
 tmp_m[i,:] = mu.flatten()
 tmp_v[i,:] = v.flatten()

# get mean and variance at X3
mu2 = np.mean(tmp_m, axis = 0)
v2 = np.mean(tmp_v, axis = 0) + np.var(tmp_m, axis = 0)
mu2 = mu2[:,None]
v3 = np.abs(v2[:,None])


''' Train level 3 '''
XX = np.hstack((X3, mu2))

k3 = GPy.kern.RBF(1, active_dims = [dim])*GPy.kern.RBF(dim, active_dims = active_dimensions, ARD = True) \
+ GPy.kern.RBF(dim, active_dims = active_dimensions, ARD = True)

m3 = GPy.models.GPRegression(X=XX, Y=Y3, kernel=k3)

m3[".*Gaussian_noise"] = m3.Y.var()*0.01
m3[".*Gaussian_noise"].fix()

m3.optimize(max_iters = max_iter)

m3[".*Gaussian_noise"].unfix()
m3[".*Gaussian_noise"].constrain_positive()

m3.optimize_restarts(10, optimizer = "bfgs",  max_iters = max_iter)

# Compute posterior mean and variance for level 3 evaluated at the test points

# sample f_1 at Xtest
mu0, C0 = m1.predict(Xtest, full_cov=True)
Z = np.random.multivariate_normal(mu0.flatten(),C0,MCsample)

# push samples through f_2 and f_3
tmp_m = np.zeros((MCsample**2,ntest))
tmp_v = np.zeros((MCsample**2,ntest))
cnt = 0
for i in range(0,MCsample):
 mu, C = m2.predict(np.hstack((Xtest, Z[i,:][:,None])), full_cov=True)
 Q = np.random.multivariate_normal(mu.flatten(),C,MCsample)
 for j in range(0,MCsample):
  mu, v = m3.predict(np.hstack((Xtest, Q[j,:][:,None])))
  tmp_m[cnt,:] = mu.flatten()
  tmp_v[cnt,:] = v.flatten()
  cnt = cnt + 1


# get f_2 posterior mean and variance at Xtest
mu3 = np.mean(tmp_m, axis = 0)
v3 = np.mean(tmp_v, axis = 0) + np.var(tmp_m, axis = 0)
mu3 = mu3[:,None]
v3 = np.abs(v3[:,None])
end = time.time()

error = np.sqrt(np.mean((mu3-Ytest)**2))
crps = np.mean(-np.sqrt(v3)*(1/np.sqrt(np.pi)-2*stats.norm.pdf((Ytest-mu3)/np.sqrt(v3))-(Ytest-mu3)/np.sqrt(v3)*(2*stats.norm.cdf((Ytest-mu3)/np.sqrt(v3))-1)))
ctime = (end - start)

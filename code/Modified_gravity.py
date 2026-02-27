!pip install emcee corner

import numpy as np
import pandas as pd
from scipy.integrate import quad, solve_ivp
from scipy.optimize import approx_fprime
import emcee
import corner
import matplotlib.pyplot as plt

# Данные
filename = 'Pantheon+SH0ES (1).dat'
data = pd.read_csv(filename, sep='\s+', comment='#')
z_sn = data['zHD'].values[::50]
mu_obs = data['MU_SH0ES'].values[::50]
mu_err = data['MU_SH0ES_ERR_DIAG'].values[::50]
mask = (mu_err > 0) & (mu_err < 10) & (z_sn > 0.01)
z_sn, mu_obs, mu_err = z_sn[mask], mu_obs[mask], mu_err[mask]
print(f"Сверхновых: {len(z_sn)}")

# DESI BAO
z_bao = np.array([0.30, 0.51, 0.71, 0.93, 1.32, 1.49, 1.85])
dv_obs = np.array([8.467, 12.642, 16.648, 20.625, 26.923, 28.821, 33.476])
dv_err = np.array([0.167, 0.235, 0.330, 0.414, 0.614, 0.648, 0.892])

# fσ8
z_fs8 = np.array([0.02, 0.10, 0.15, 0.22, 0.30, 0.38, 0.51, 0.60, 0.70, 0.77, 0.85, 1.00, 1.20, 1.40, 1.60])
fs8_obs = np.array([0.428, 0.470, 0.490, 0.420, 0.440, 0.440, 0.455, 0.430, 0.430, 0.450, 0.460, 0.440, 0.430, 0.410, 0.400])
fs8_err = np.array([0.047, 0.060, 0.055, 0.070, 0.050, 0.040, 0.040, 0.040, 0.045, 0.040, 0.055, 0.050, 0.060, 0.070, 0.080])

class ModifiedGravity:
    def __init__(self, H0, Omega_m, z_tr, beta, alpha, Delta_z=1.5):
        self.H0 = H0
        self.Omega_m = Omega_m
        self.z_tr = z_tr
        self.beta = beta
        self.alpha = alpha
        self.Delta_z = Delta_z
        self.c = 299792.458
        self.Omega_sub = 1.0 - Omega_m
        self.eps = 1e-4
        self.sigma8_norm = 0.8
        
    def _derivative(self, func, x):
        def f(x0): return func(x0)
        return approx_fprime([x], f, self.eps)[0]
    
    def Phi(self, z):
        return 0.5 * (1 + np.tanh((self.z_tr - z) / self.Delta_z))
    
    def dPhi_dz(self, z):
        return -0.5 / self.Delta_z * (1 / np.cosh((self.z_tr - z) / self.Delta_z))**2
    
    def Geff(self, z):
        return 1.0 + self.beta * self.Phi(z)  # Усиление гравитации
    
    def E2(self, z):
        return self.Omega_m * (1+z)**3 + self.Omega_sub * (1+z)**self.alpha
    
    def H(self, z):
        return self.H0 * np.sqrt(self.E2(z))
    
    def growth_equation(self, z, y):
        delta, ddelta = y
        E = np.sqrt(self.E2(z))
        E_prime = self._derivative(lambda zz: np.sqrt(self.E2(zz)), z)
        Omega_m_z = self.Omega_m * (1+z)**3 / self.E2(z)
        
        A = 3/(z+1) + E_prime/E
        B = 1.5 * Omega_m_z * self.Geff(z)  # Модифицированная гравитация
        
        return [ddelta, -A*ddelta + B*delta]
    
    def growth_factor(self, z, z_init=10.0):
        if z >= z_init:
            return 1.0/(1+z)
        sol = solve_ivp(self.growth_equation, [z_init, z], [1.0, -1.0/(1+z_init)], 
                       method='RK45', rtol=1e-4)
        return sol.y[0,-1] if sol.success else 1.0/(1+z)
    
    def f(self, z):
        dz = 0.05
        z1, z2 = max(0.01, z-dz), z+dz
        D1, D2 = self.growth_factor(z1), self.growth_factor(z2)
        return -(1+z)*(np.log(D2)-np.log(D1))/(z2-z1)
    
    def fsigma8(self, z):
        D0 = self.growth_factor(0)
        Dz = self.growth_factor(z)
        fz = self.f(z)
        return fz * self.sigma8_norm * Dz / D0
    
    def distance_modulus(self, z):
        integral, _ = quad(lambda zp: 1/self.H(zp), 0, z, limit=50)
        return 5*np.log10(self.c*integral*(1+z)) + 25
    
    def D_V(self, z):
        dm = self.c * quad(lambda zp: 1/self.H(zp), 0, z, limit=50)[0]
        return (z * dm**2 / self.H(z))**(1/3)

def log_likelihood(params):
    H0, Omega_m, z_tr, beta, alpha = params
    
    if not (50 < H0 < 100 and 0.1 < Omega_m < 0.5 and 10 < z_tr < 100 
            and 0 < beta < 2 and 1 < alpha < 4):
        return -np.inf
    
    try:
        model = ModifiedGravity(H0, Omega_m, z_tr, beta, alpha)
        
        # BAO
        dv_model = np.array([model.D_V(z) for z in z_bao])
        chi2_bao = np.sum(((dv_obs - dv_model) / dv_err)**2)
        
        # Сверхновые
        chi2_sn = 0
        for i in range(len(z_sn)):
            chi2_sn += ((mu_obs[i] - model.distance_modulus(z_sn[i])) / mu_err[i])**2
        
        # fσ8
        fs8_model = np.array([model.fsigma8(z) for z in z_fs8])
        chi2_fs8 = np.sum(((fs8_obs - fs8_model) / fs8_err)**2)
        
        return -0.5 * (chi2_bao + chi2_sn/10 + chi2_fs8)
    except:
        return -np.inf

# MCMC
ndim, nwalkers = 5, 32
pos = np.array([70.0, 0.3, 50.0, 0.5, 2.0]) + 1e-3*np.random.randn(nwalkers, ndim)

print("Запуск MCMC для модифицированной гравитации...")
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood)
sampler.run_mcmc(pos, 500, progress=True)

samples = sampler.get_chain(discard=100, flat=True)
labels = ["H0", "Ω_m", "z_tr", "β", "α"]

print("\nРезультаты:")
for i, label in enumerate(labels):
    m, s = np.mean(samples[:, i]), np.std(samples[:, i])
    print(f"{label} = {m:.3f} ± {s:.3f}")

# График fσ8
model = ModifiedGravity(*np.mean(samples, axis=0))
z_plot = np.linspace(0, 1.8, 50)
fs8_model = [model.fsigma8(z) for z in z_plot]

plt.figure(figsize=(12, 8))
plt.errorbar(z_fs8, fs8_obs, yerr=fs8_err, fmt='o', color='red', label='Данные')
plt.plot(z_plot, fs8_model, 'b-', linewidth=2, label='Модифицированная гравитация')
plt.axhline(y=0.45, color='gray', linestyle='--', alpha=0.5)
plt.xlabel('z')
plt.ylabel('fσ₈(z)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

corner.corner(samples, labels=labels, quantiles=[0.16,0.5,0.84], show_titles=True)
plt.show()
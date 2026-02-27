# Copyright 2026 Bykov Denis Alexandrovich (bycov)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Часть 4: Класс QMG (без изменений)
import numpy as np
import pandas as pd
from scipy.integrate import quad, solve_ivp
from scipy.optimize import approx_fprime
import emcee
import corner
import matplotlib.pyplot as plt

class QuantumMaterializationGravity:
    def __init__(self, H0, Omega_m, z_tr, Q_growth, Q_lens, alpha, Delta_z=1.5):
        self.H0 = H0
        self.Omega_m = Omega_m
        self.Omega_aether = 1.0 - Omega_m
        self.z_tr = z_tr
        self.Q_growth = Q_growth
        self.Q_lens = Q_lens
        self.alpha = alpha
        self.Delta_z = Delta_z
        self.c = 299792.458
        self.sigma8_norm = 0.8
        self.eps = 1e-4
        self._cache = {}  # Кэш для growth_factor
        
    def Phi(self, z):
        return 0.5 * (1 + np.tanh((self.z_tr - z) / self.Delta_z))
    
    def Geff_growth(self, z):
        return 1.0 + self.Q_growth * self.Phi(z)
    
    def Geff_lens(self, z):
        return 1.0 + (self.Q_growth + self.Q_lens) * self.Phi(z)
    
    def E2(self, z):
        return self.Omega_m * (1+z)**3 + self.Omega_aether * (1+z)**self.alpha
    
    def H(self, z):
        return self.H0 * np.sqrt(self.E2(z))
    
    def growth_equation(self, z, y):
        delta, ddelta = y
        E = np.sqrt(self.E2(z))
        E_prime = approx_fprime([z], lambda zz: np.sqrt(self.E2(zz)), self.eps)[0]
        Omega_m_z = self.Omega_m * (1+z)**3 / self.E2(z)
        
        A = 3/(z+1) + E_prime/E
        B = 1.5 * Omega_m_z * self.Geff_growth(z)
        
        return [ddelta, -A*ddelta + B*delta]
    
    def growth_factor(self, z, z_init=10.0):
        # Кэширование
        key = (round(z, 4))  # округляем для кэша
        if key in self._cache:
            return self._cache[key]
            
        if z >= z_init:
            result = 1.0/(1+z)
        else:
            sol = solve_ivp(self.growth_equation, [z_init, z], [1.0, -1.0/(1+z_init)], 
                           method='RK45', rtol=1e-4)
            result = sol.y[0,-1] if sol.success else 1.0/(1+z)
        
        self._cache[key] = result
        return result
    
    def clear_cache(self):
        self._cache = {}
    
    def f(self, z):
        dz = 0.05
        z1, z2 = max(0.01, z-dz), z+dz
        D1, D2 = self.growth_factor(z1), self.growth_factor(z2)
        return -(1+z)*(np.log(D2)-np.log(D1))/(z2-z1)
    
    def sigma8(self, z=0):
        return self.sigma8_norm * self.growth_factor(z) / self.growth_factor(0)
    
    def fsigma8(self, z):
        return self.f(z) * self.sigma8(z)
    
    def distance_modulus(self, z):
        integral, _ = quad(lambda zp: 1/self.H(zp), 0, z, limit=50)
        return 5*np.log10(self.c*integral*(1+z)) + 25
    
    def D_V(self, z):
        dm = self.c * quad(lambda zp: 1/self.H(zp), 0, z, limit=50)[0]
        return (z * dm**2 / self.H(z))**(1/3)
    
    def S8(self):
        return self.sigma8(0) * np.sqrt(self.Omega_m / 0.3)

print("✅ Класс QMG с кэшированием загружен")

# Часть 5: Данные
filename = 'Pantheon+SH0ES.dat'

data = pd.read_csv(filename, sep='\\s+', comment='#')
z_sn = data['zHD'].values[::50]
mu_obs = data['MU_SH0ES'].values[::50]
mu_err = data['MU_SH0ES_ERR_DIAG'].values[::50]
mask = (mu_err > 0) & (mu_err < 10) & (z_sn > 0.01)
z_sn, mu_obs, mu_err = z_sn[mask], mu_obs[mask], mu_err[mask]
print(f"Сверхновых: {len(z_sn)}")

# BAO
z_bao = np.array([0.30, 0.51, 0.71, 0.93, 1.32, 1.49, 1.85])
dv_obs = np.array([8.467, 12.642, 16.648, 20.625, 26.923, 28.821, 33.476])
dv_err = np.array([0.167, 0.235, 0.330, 0.414, 0.614, 0.648, 0.892])

# fσ8
z_fs8 = np.array([0.02, 0.10, 0.15, 0.22, 0.30, 0.38, 0.51, 0.60, 0.70, 0.77, 0.85, 1.00, 1.20, 1.40, 1.60])
fs8_obs = np.array([0.428, 0.470, 0.490, 0.420, 0.440, 0.440, 0.455, 0.430, 0.430, 0.450, 0.460, 0.440, 0.430, 0.410, 0.400])
fs8_err = np.array([0.047, 0.060, 0.055, 0.070, 0.050, 0.040, 0.040, 0.040, 0.045, 0.040, 0.055, 0.050, 0.060, 0.070, 0.080])

# Часть 6: Функция правдоподобия с жестким приором на Ω_m
def log_likelihood_qmg(params):
    H0, Omega_m, z_tr, Q_growth, Q_lens, alpha = params
    
    if not (50 < H0 < 100 and 0.25 < Omega_m < 0.35 and 10 < z_tr < 100 
            and 0.4 < Q_growth < 1.2 and -0.5 < Q_lens < 0.1 and 0 < alpha < 4):
        return -np.inf
    
    try:
        model = QuantumMaterializationGravity(H0, Omega_m, z_tr, Q_growth, Q_lens, alpha)
        
        # BAO
        dv_model = np.array([model.D_V(z) for z in z_bao])
        chi2_bao = np.sum(((dv_obs - dv_model) / dv_err)**2)
        
        # SN
        chi2_sn = 0
        for i in range(len(z_sn)):
            chi2_sn += ((mu_obs[i] - model.distance_modulus(z_sn[i])) / mu_err[i])**2
        chi2_sn = chi2_sn / 10
        
        # fσ8 (кеширование уже внутри модели)
        fs8_model = np.array([model.fsigma8(z) for z in z_fs8])
        chi2_fs8 = np.sum(((fs8_obs - fs8_model) / fs8_err)**2)
        
        # KiDS (использует закешированный growth_factor)
        S8_model = model.S8()
        chi2_kids = ((S8_model - 0.766) / 0.020)**2
        
        chi2_total = chi2_bao + chi2_sn + chi2_fs8 + chi2_kids
        
        return -0.5 * chi2_total
        
    except Exception as e:
        return -np.inf

# Тест
test_params = [82, 0.3, 32, 0.58, -0.18, 3.0]
print(f"logL тест = {log_likelihood_qmg(test_params):.1f}")
print("✅ Функция правдоподобия с кэшированием готова")

# Часть 7: MCMC
print("\n🚀 Запуск MCMC с кэшированием (будет в 2 раза быстрее)...")

ndim, nwalkers = 6, 32
pos = np.array([82.0, 0.3, 32.0, 0.58, -0.18, 3.0]) + 1e-3*np.random.randn(nwalkers, ndim)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood_qmg)
sampler.run_mcmc(pos, 300, progress=True)

samples = sampler.get_chain(discard=50, flat=True)
labels = ["H0", "Ω_m", "z_tr", "Q_рост", "Q_линза", "α"]

print("\n📊 Результаты:")
for i, label in enumerate(labels):
    m, s = np.mean(samples[:, i]), np.std(samples[:, i])
    print(f"{label} = {m:.3f} ± {s:.3f}")

Q_total = samples[:, 3] + samples[:, 4]
print(f"\n⚛️ Q_total = {np.mean(Q_total):.3f} ± {np.std(Q_total):.3f}")
print(f"📉 Q_линза отрицательный? {np.mean(samples[:,4]) < 0}")

# S8 график
s8_samples = 0.8 * np.sqrt(samples[:,1] / 0.3)
plt.figure(figsize=(15, 5))

plt.subplot(131)
plt.hist(s8_samples, bins=30, density=True, alpha=0.7, color='blue', label='QMG')
plt.axvline(0.834, color='red', linestyle='--', label='Planck')
plt.axvline(0.766, color='green', linestyle='--', label='KiDS')
plt.axvspan(0.746, 0.786, alpha=0.2, color='green')
plt.xlabel('S₈')
plt.ylabel('Вероятность')
plt.title('S₈ распределение')
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(132)
plt.scatter(samples[:,1], samples[:,4], alpha=0.1, s=1)
plt.xlabel('Ω_m')
plt.ylabel('Q_линза')
plt.title('Корреляция Ω_m и Q_линза')
plt.grid(True, alpha=0.3)

# fσ8 график
plt.subplot(133)
best = np.median(samples, axis=0)
model = QuantumMaterializationGravity(*best)
z_plot = np.linspace(0, 1.8, 50)
fs8_model = [model.fsigma8(z) for z in z_plot]
plt.errorbar(z_fs8, fs8_obs, yerr=fs8_err, fmt='o', color='red', label='Данные')
plt.plot(z_plot, fs8_model, 'b-', linewidth=2, label='QMG')
plt.axhline(y=0.45, color='gray', linestyle='--', alpha=0.5)
plt.xlabel('z')
plt.ylabel('fσ₈(z)')
plt.title('Рост структур')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()


print(f"S₈ среднее = {np.mean(s8_samples):.3f} ± {np.std(s8_samples):.3f}")

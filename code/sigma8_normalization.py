import numpy as np
import camb

H0 = 82.9
omega_m = 0.30
target_sigma8 = 0.8

# Данные из вашего скана
As_values = np.array([1.50e-9, 1.61e-9, 1.72e-9, 1.83e-9, 1.94e-9, 2.06e-9, 2.17e-9, 2.28e-9, 2.39e-9, 2.50e-9])
sigma8_values = np.array([0.935, 0.969, 1.002, 1.033, 1.064, 1.094, 1.124, 1.152, 1.180, 1.207])

# Линейная интерполяция
from scipy.interpolate import interp1d
As_func = interp1d(sigma8_values, As_values, kind='linear', fill_value='extrapolate')
As_needed = As_func(target_sigma8)

print(f"As для σ₈=0.8: {As_needed:.2e}")

# Проверка
pars = camb.CAMBparams()
pars.set_cosmology(H0=H0, ombh2=0.0224, omch2=omega_m*H0**2/100**2 - 0.0224, mnu=0.06, omk=0)
pars.InitPower.set_params(As=As_needed, ns=0.965, r=0)
pars.set_matter_power(redshifts=[0.0], kmax=2.0)
results = camb.get_results(pars)
sigma8_check = results.get_sigma8()[0]
print(f"Проверка σ₈ = {sigma8_check:.3f}")
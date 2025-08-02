import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# PHYSICAL CONSTANTS AND CROSS-SECTIONS 
tho2_density = 10
uo2_density = 10.5
zr_density = 6.4
h2o_density = 0.7

tho2_Mr = 264
uo2_Mr = 265
zr_Mr = 91.22
h2o_Mr = 18

sga_Th = 7.56
sga_U233 = 573
sga_H2O = 0.66
sga_Zr = 0.185

sgs_Th = 12.6
sgs_U233 = 6
sgs_H2O = 103
sgs_Zr = 8

sgf_u233 = 524

# HELPER FUNCTIONS 
def Number(vol_frac, density, mr):
    return vol_frac * (density / mr) * 6.02e23

def SGA_Macroscopic(N, sigma):
    return N * (sigma / 1e24)

def SGS_Macroscopic(N, sigma):
    return N * (sigma / 1e24)

def SGF_Macroscopic(N, sigma):
    return N * (sigma / 1e24)

def SGTr_Macroscopic(SGA, SGS, one_minus_mu0):
    return SGA + one_minus_mu0 * SGS

def Dr(SGTr):
    return 1 / (3 * SGTr)

def KEff(SGF_Total, SGA_Total, v, Dr):
    return (v * SGF_Total) / (SGA_Total + 0.1)  # Empirical leakage term

def calculateCurrentKEff(N_th, N_u, vol_h2o, vol_zr, vol_uo2, vol_tho2):
    n_h2o = Number(vol_h2o, h2o_density, h2o_Mr)
    n_zr = Number(vol_zr, zr_density, zr_Mr)

    SGA_Th = SGA_Macroscopic(N_th, sga_Th)
    SGA_U233 = SGA_Macroscopic(N_u, sga_U233)
    SGA_H2O = SGA_Macroscopic(n_h2o, sga_H2O)
    SGA_Zr = SGA_Macroscopic(n_zr, sga_Zr)
    SGA_Total = SGA_Th + SGA_U233 + SGA_H2O + SGA_Zr

    SGS_Th = SGS_Macroscopic(N_th, sgs_Th)
    SGS_U233 = SGS_Macroscopic(N_u, sgs_U233)
    SGS_H2O = SGS_Macroscopic(n_h2o, sgs_H2O)
    SGS_Zr = SGS_Macroscopic(n_zr, sgs_Zr)

    SGTr_Th = SGTr_Macroscopic(SGA_Th, SGS_Th, 0.9971)
    SGTr_U233 = SGTr_Macroscopic(SGA_U233, SGS_U233, 0.9972)
    SGTr_H2O = SGTr_Macroscopic(SGA_H2O, SGS_H2O, 0.676)
    SGTr_Zr = SGTr_Macroscopic(SGA_Zr, SGS_Zr, 0.9927)
    SGTr_Total = SGTr_Th + SGTr_U233 + SGTr_H2O + SGTr_Zr

    SGF_Total = SGF_Macroscopic(N_u, sgf_u233)
    D = Dr(SGTr_Total)

    return KEff(SGF_Total, SGA_Total, 2.49, D)

# INITIAL VOLUME FRACTIONS 
vol_uo2 = 0.04
vol_tho2 = 0.26
vol_h2o = 0.55
vol_zr = 0.15
assert abs(vol_uo2 + vol_tho2 + vol_h2o + vol_zr - 1) < 1e-6, "Volume fractions must sum to 1"

# INITIAL NUMBER DENSITIES 
n_th = Number(vol_tho2, tho2_density, tho2_Mr)
n_u = Number(vol_uo2, uo2_density, uo2_Mr)
n_pa = 0

# THORIUM FUEL CYCLE DYNAMICS
sigma_th = 7.4e-24
sigma_pa = 200e-24
sigma_u = 530e-24
phi = 1e14
lambda_pa = 1.78e-5

days = np.linspace(0, 3000, 500)
time_sec = days * 24 * 3600
y0 = [n_th, n_pa, n_u]

def thorium_cycle(y, t, phi, sigma_th, lambda_pa, sigma_u):
    N_th, N_pa, N_u = y
    dN_th_dt = -sigma_th * phi * N_th
    dN_pa_dt = sigma_th * phi * N_th - lambda_pa * N_pa
    dN_u_dt = lambda_pa * N_pa - sigma_u * phi * N_u
    return [dN_th_dt, dN_pa_dt, dN_u_dt]

solution = odeint(thorium_cycle, y0, time_sec, args=(phi, sigma_th, lambda_pa, sigma_u))
N_th, N_pa, N_u = solution.T

# CALCULATE CRITICALITY OVER TIME 
target_days = [0] + [(x + 1) * 100 for x in range(30)]
keff_over_time = []

for target_day in target_days:
    idx = (np.abs(days - target_day)).argmin()
    thorium = N_th[idx]
    uranium = N_u[idx]
    keff = calculateCurrentKEff(thorium, uranium, vol_h2o, vol_zr, vol_uo2, vol_tho2)
    print(f"Criticality on day {target_day}: {keff}")
    keff_over_time.append((target_day, keff))

# PLOTTING 
plt.figure(figsize=(12, 6))
plt.plot(days, N_th, label='Th-232', color='red')
plt.plot(days, N_u, label='U-233', color='blue')
plt.xlabel('Time (days)')
plt.ylabel('Number Density (atoms/cm³)')
plt.title('Thorium Fuel Cycle: Number Density vs Time')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

days_k, keffs = zip(*keff_over_time)
plt.figure(figsize=(12, 6))
plt.plot(days_k, keffs, marker='o', color='green')
plt.xlabel('Time (days)')
plt.ylabel('k_eff')
plt.title('Reactor k_eff Over Time')
plt.grid(True)
plt.tight_layout()
plt.show()


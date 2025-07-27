import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# DIFFUSION CALCULATION
def Number(Vol_Fraction_Fuel, Density, Mr):
  avogadro_constant = (6.02)*(10**23)
  return Vol_Fraction_Fuel * (Density / Mr) * avogadro_constant
def SGA_Macroscopic(N, SGA_Microscopic):
  return N * (SGA_Microscopic / (10**24))
def SGS_Macroscopic(N, SGS_Microscopic):
  return N * (SGS_Microscopic / (10**24))
def SGF_Macroscopic(N, SGF_Microscopic):
  return N * (SGF_Microscopic / (10**24))
def SGTr_Macroscopic(SGA_Macroscopic, SGS_Macroscopic, one_minus_mu0):
  return SGA_Macroscopic + one_minus_mu0 * SGS_Macroscopic
def Dr(SGTr):
  return (1/(3*SGTr))
def KEff(SGF_Total, SGA_Total, v, Dr):
  numerator = v * SGF_Total
  denominator = SGA_Total + Dr * 3.14159 * 3.14159
  return numerator / denominator
def calculateCurrentKEff(currentThoriumNumber, currentUraniumNumber):
  current_SGA_Th = SGA_Macroscopic(currentThoriumNumber, 7.56)
  current_SGA_U233 = SGA_Macroscopic(currentUraniumNumber, 573)
  SGA_H2O = SGA_Macroscopic(Number(0.55, 0.7, 18), 0.66)
  SGA_Zr = SGA_Macroscopic(Number(0.15, 6.4, 91.22), 0.185)
  current_SGA_Total = current_SGA_Th + current_SGA_U233 + SGA_H2O + SGA_Zr

  current_SGS_Th = SGS_Macroscopic(currentThoriumNumber, 12.6)
  current_SGS_U233 = SGS_Macroscopic(currentUraniumNumber, 6)
  SGS_H2O = SGS_Macroscopic(Number(0.55, 0.7, 18), 103)
  SGS_Zr = SGS_Macroscopic(Number(0.15, 6.4, 91.22), 8)
  current_SGS_Total = current_SGS_Th + current_SGS_U233 + SGS_H2O + SGS_Zr


  current_SGTr_Th = SGTr_Macroscopic(current_SGA_Th, current_SGS_Th, 0.9971)
  current_SGTr_U233 = SGTr_Macroscopic(current_SGA_U233, current_SGS_U233, 0.9972)
  SGTr_H2O = SGTr_Macroscopic(SGA_H2O, SGS_H2O, 0.676)
  SGTr_Zr = SGTr_Macroscopic(SGA_Zr, SGS_Zr, 0.9927)
  current_SGTr_Total = current_SGTr_Th + current_SGTr_U233 + SGTr_H2O + SGTr_Zr

  SGF_Total = SGF_Macroscopic(currentUraniumNumber, 524)
  current_reactor_Dr = Dr(current_SGTr_Total)
  current_keff = KEff(SGF_Total, current_SGA_Total, 2.49, current_reactor_Dr)
  
  return current_keff

# Constants
avogadro_constant = (6.02)*(10**23)

volume_fraction_cladding = 0.15
volume_fraction_h2o = 0.55
volume_fraction_tho2 = 0.21
volume_fraction_uo2 = 0.09

tho2_density = 10
uo2_density = 10.5
zr_density = 6.4
h2o_density = 0.7

tho2_Mr = 264
Th_Mr = 232
Zr_Mr = 91.22
h2o_Mr = 18
uo2_Mr = 265

sga_Th = 7.56
sga_U233 = 573
sga_H2O = 0.66
sga_Zr = 0.185

sgs_Th = 12.6
sgs_U233 = 6
sgs_H2O = 103
sgs_Zr = 8

sgf_u233 = 524

# Numbers
n_tho2 = Number(volume_fraction_tho2, tho2_density, tho2_Mr)
n_uo2 = Number(volume_fraction_uo2, uo2_density, uo2_Mr)
n_zr = Number(volume_fraction_cladding, zr_density, Zr_Mr) # Number(0.15, 6.4, 91.22)
n_h2o = Number(volume_fraction_h2o, h2o_density, h2o_Mr) # Number(0.55, 0.7, 18)
n_th = n_tho2
n_u233 = n_uo2
n_h = n_h2o * 2
n_o = n_h2o + 2 * n_tho2

# Macroscopics
SGA_Th = SGA_Macroscopic(n_th, sga_Th)
SGA_U233 = SGA_Macroscopic(n_u233, sga_U233)
SGA_H2O = SGA_Macroscopic(n_h2o, sga_H2O) # SGA_Macroscopic(Number(0.55, 0.7, 18), 0.66)
SGA_Zr = SGA_Macroscopic(n_zr, sga_Zr) # SGA_Macroscopic(Number(0.15, 6.4, 91.22), 0.185)
SGA_Total = SGA_Th + SGA_U233 + SGA_H2O + SGA_Zr

SGS_Th = SGS_Macroscopic(n_th, sgs_Th)
SGS_U233 = SGS_Macroscopic(n_u233, sgs_U233)
SGS_H2O = SGS_Macroscopic(n_h2o, sgs_H2O)
SGS_Zr = SGS_Macroscopic(n_zr, sgs_Zr)

SGTr_Th = SGTr_Macroscopic(SGA_Th, SGS_Th, 0.9971)
SGTr_U233 = SGTr_Macroscopic(SGA_U233, SGS_U233, 0.9972)
SGTr_H2O = SGTr_Macroscopic(SGA_H2O, SGS_H2O, 0.676)
SGTr_Zr = SGTr_Macroscopic(SGA_Zr, SGS_Zr, 0.9927)
SGTr_Total = SGTr_Th + SGTr_U233 + SGTr_H2O + SGTr_Zr

SGF_Total = SGF_Macroscopic(n_u233, sgf_u233)

# KEff Calculation
Dr_Reactor = Dr(SGTr_Total)
k_eff = KEff(SGF_Total, SGA_Total, 2.49, Dr_Reactor)

# Printing
print("KEff: " + str(round(k_eff, 5)))

if k_eff > 1:
  print("Supercritical")
elif k_eff < 1:
  print("Subcritical")
elif k_eff == 1:
  print("Critical")



# GRAPHING PROCESS
sigma_th = 7.4e-24   # Microscopic capture cross-section of Th-232 in cm^2 (7.4 barns)
sigma_pa = 200e-24
sigma_u = 530e-24
phi = 1e14 # Neutron flux in n/cm^2/s
lambda_pa = 1.78e-5  # Decay constant of Pa-233 in 1/s (half-life ~27 days)

# Initial number densities in atoms/cm^3
N_th0 = n_th  # Initial Th-232
N_pa0 = 0     # Initial Pa-233
N_u0 = n_u233    # Initial U-233

# Time array (0 to 2000 days, converted to seconds)
days = np.linspace(0, 3000, 500)
time_sec = days * 24 * 3600

# Define system of Differential Equations
def thorium_cycle(y, t, phi, sigma_th, lambda_pa, sigma_u):
    N_th, N_pa, N_u = y
    dN_th_dt = -sigma_th * phi * N_th
    dN_pa_dt = sigma_th * phi * N_th - lambda_pa * N_pa
    dN_u_dt = lambda_pa * N_pa - sigma_u * phi * N_u
    return [dN_th_dt, dN_pa_dt, dN_u_dt]

# Initial conditions vector
y0 = [N_th0, N_pa0, N_u0]

# Solve ODEs
solution = odeint(thorium_cycle, y0, time_sec, args=(phi, sigma_th, lambda_pa, sigma_u))
N_th, N_pa, N_u = solution.T

# Get KEff at specific time points
target_days = [0]
for x in range(30):
  target_days.append((x+1)*100)
print(target_days)
print("\n--- k_eff at Specific Days ---")
for target_day in target_days:
    # Find the index closest to the desired day
    idx = (np.abs(days - target_day)).argmin()
    thorium = N_th[idx]
    uranium = N_u[idx]
    keff = calculateCurrentKEff(thorium, uranium)
    print(f"Day {target_day}: k_eff = {keff:.5f}")


# Plotting
plt.figure(figsize=(12, 6))
plt.plot(days, N_th, label='Th-232',color='red')
plt.plot(days, N_u, label='U-233',color='Blue')
plt.xlabel('Time (days)')
plt.ylabel('Number Density Of Th-232 / U-233')
plt.title('Thorium Fuel Cycle: Number Density vs Time')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

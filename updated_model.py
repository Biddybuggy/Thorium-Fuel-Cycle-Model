import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

NA = 6.022e23  
barn = 1e-24   

FISSILE_FRACTION_U233 = 0.25 

sga_Th_b   = 7.4
sga_Pa_b   = 200.0
sga_U233_b = 573.0
sga_H2O_b  = 0.6652
sga_Zr_b   = 0.185

sga_B10_b  = 3835.0   
sga_Gd_b   = 5000.0   

sgs_Th_b   = 3.2
sgs_U233_b = 7.2
sgs_H2O_b  = 48.0
sgs_Zr_b   = 6.0

sgf_u233_b = 524.0

nu = 2.49

# Converting barns to cm^2
def XS(x_barns):
    return x_barns * barn

sga_Th   = XS(sga_Th_b)
sga_Pa   = XS(sga_Pa_b)
sga_U233 = XS(sga_U233_b)
sga_H2O  = XS(sga_H2O_b)
sga_Zr   = XS(sga_Zr_b)
sga_B10  = XS(sga_B10_b)
sga_Gd   = XS(sga_Gd_b)

sgs_Th   = XS(sgs_Th_b)
sgs_U233 = XS(sgs_U233_b)
sgs_H2O  = XS(sgs_H2O_b)
sgs_Zr   = XS(sgs_Zr_b)

sgf_u233 = XS(sgf_u233_b)

tho2_density = 10.0
uo2_density  = 10.97
h2o_density  = 1.0
zr_density   = 6.52

tho2_Mr = 264.0
uo2_Mr  = 270.0
h2o_Mr  = 18.0
zr_Mr   = 91.224

vol_tho2 = 0.38
vol_uo2  = 0.07
vol_h2o  = 0.45
vol_zr   = 0.10

def number_density(volume_fraction, density, Mr):
    """Return atomic number density (#/cm^3)."""
    return (volume_fraction * density / Mr) * NA

def SGA_Macroscopic(N, sga): return N * sga
def SGS_Macroscopic(N, sgs): return N * sgs
def SGF_Macroscopic(N, sgf): return N * sgf

def SGTr_Macroscopic(SGA, SGS, mu):
    return SGS - mu * SGA

def Dr(SGTr):
    return 1.0 / (3.0 * SGTr)

def KEff(SGF_Total, SGA_Total, nu, D):
    return (nu * SGF_Total) / (SGA_Total + 0.1)

N_th   = number_density(vol_tho2, tho2_density, tho2_Mr)
N_uTot = number_density(vol_uo2,  uo2_density,  uo2_Mr)
N_u0   = FISSILE_FRACTION_U233 * N_uTot  

N_h2o  = number_density(vol_h2o, h2o_density, h2o_Mr)
N_zr   = number_density(vol_zr,  zr_density,  zr_Mr)

def boron_number_density_from_ppm(ppm):
    mass_fraction = ppm * 1e-6                
    mass_B_per_cm3 = mass_fraction * h2o_density  
    Mr_B = 10.81                             
    return (mass_B_per_cm3 / Mr_B) * NA

def gadolinium_number_density_from_wt_fraction(wt_frac):
    fuel_density = 10.5  
    mass_Gd_per_cm3 = wt_frac * fuel_density
    Mr_Gd = 157.25       
    return (mass_Gd_per_cm3 / Mr_Gd) * NA

B_ppm      = 200.0     
Gd_wt_frac = 2.0e-4     

N_B0  = boron_number_density_from_ppm(B_ppm)
N_Gd0 = gadolinium_number_density_from_wt_fraction(Gd_wt_frac)

phi = 3.0e13          
lambda_pa = 1.59e-6    

def derivatives(y, t):
    N_th, N_pa, N_u, N_B, N_Gd = y
    dN_th_dt = -sga_Th * phi * N_th
    dN_pa_dt = sga_Th * phi * N_th - lambda_pa * N_pa
    dN_u_dt  = lambda_pa * N_pa - sga_U233 * phi * N_u
    dN_B_dt  = -sga_B10 * phi * N_B
    dN_Gd_dt = -sga_Gd  * phi * N_Gd
    return [dN_th_dt, dN_pa_dt, dN_u_dt, dN_B_dt, dN_Gd_dt]


N_pa0 = 0.0
y0 = [N_th, N_pa0, N_u0, N_B0, N_Gd0]
def calculateCurrentKEff(N_th, N_pa, N_u, N_B, N_Gd):
    SGA_Th   = SGA_Macroscopic(N_th,   sga_Th)
    SGA_U    = SGA_Macroscopic(N_u,    sga_U233)
    SGA_Pa   = SGA_Macroscopic(N_pa,   sga_Pa)
    SGA_H2O  = SGA_Macroscopic(N_h2o,  sga_H2O)
    SGA_Zr   = SGA_Macroscopic(N_zr,   sga_Zr)
    SGA_B    = SGA_Macroscopic(N_B,    sga_B10)
    SGA_Gd   = SGA_Macroscopic(N_Gd,   sga_Gd)

    SGA_Total = SGA_Th + SGA_U + SGA_Pa + SGA_H2O + SGA_Zr + SGA_B + SGA_Gd

    SGS_Th   = SGS_Macroscopic(N_th,   sgs_Th)
    SGS_U    = SGS_Macroscopic(N_u,    sgs_U233)
    SGS_H2O  = SGS_Macroscopic(N_h2o,  sgs_H2O)
    SGS_Zr   = SGS_Macroscopic(N_zr,   sgs_Zr)

    SGTr_Total = SGTr_Macroscopic(SGA_Th,  SGS_Th,  0.9971) \
               + SGTr_Macroscopic(SGA_U,   SGS_U,   0.9972) \
               + SGTr_Macroscopic(SGA_H2O, SGS_H2O, 0.676)  \
               + SGTr_Macroscopic(SGA_Zr,  SGS_Zr,  0.9927)

    SGF_Total = SGF_Macroscopic(N_u, sgf_u233)

    D = Dr(SGTr_Total)
    return KEff(SGF_Total, SGA_Total, nu, D)

years = 5             
days  = years * 365
time  = np.linspace(0, days * 86400.0, 4000)  

solution = odeint(derivatives, y0, time)
N_th_t, N_pa_t, N_u_t, N_B_t, N_Gd_t = solution.T

k_eff_values = np.array([
    calculateCurrentKEff(N_th_t[i], N_pa_t[i], N_u_t[i], N_B_t[i], N_Gd_t[i])
    for i in range(len(time))
])

print(f"Initial k_eff = {k_eff_values[0]:.3f}")
idx_crit = np.where(k_eff_values <= 1.0)[0][0]
t_crit_days = time[idx_crit] / 86400.0
print(f"k_eff drops below 1 at ~{t_crit_days:.1f} days (~{t_crit_days/30.0:.1f} months)")

plt.figure(figsize=(8, 5))
plt.plot(time / 86400.0, k_eff_values, label=r"$k_{\mathrm{eff}}$")
plt.axhline(1.0, color="r", linestyle="--", label="critical (k=1)")
plt.xlabel("time (days)")
plt.ylabel(r"$k_{\mathrm{eff}}$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

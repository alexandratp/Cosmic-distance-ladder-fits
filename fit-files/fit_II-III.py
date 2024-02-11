import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#Lendo arquivos
Y_fits_path = 'Y_r22-pantheon+_unc.fits' #'Y_syn_w22.fits'
L_fits_path = 'L_II.fits'                #'L_III.fits'          
C_fits_path = 'C_r22-pantheon+_unc.fits'

Y = fits.open(Y_fits_path)[0].data
L = fits.open(L_fits_path)[0].data
C = fits.open(C_fits_path)[0].data

#São montados vetores de períodos e metalicidades de todas as 
#cefeidas a partir da matriz de quações
logPm1 = []
OH = []

for i in range(3130):
  logPm1.append(L[41, i])

for i in range(3130):
  OH.append(L[43, i])

#Aqui são carregados os valores do parâmetros ajustados pelo
#Riess et al. 2022
params = np.loadtxt( 'lstsq_results.txt', comments='#' )
Q_best = params[:,0]
i_Q_best = params[:,1]

#Carregando parâmetros de ajuste
M_B = Q_best[42]
sigma_M_B = i_Q_best[42]
M_HW = Q_best[38]
sigma_M_HW = i_Q_best[38]
Z_W = Q_best[43]
sigma_Z_W = i_Q_best[43]
b_W = -3.299
sigma_b_W = 0.015
mu_m31 = Q_best[40]
sigma_mu_m31 = i_Q_best[40]
delta_zp = Q_best[45]
sigma_delta_zp = i_Q_best[45]
delta_mu_NGC4258 = Q_best[37]
sigma_delta_mu_NGC4258 = i_Q_best[37]
delta_mu_LMC = Q_best[39]
sigma_delta_mu_LMC = i_Q_best[39]
fivelogH0 = Q_best[46]
sigma_fivelogH0 = i_Q_best[46]
c=299792
q0=-0.55
j0=1
mu_NGC4258 = 29.398
sigma_mu_NGC4258 = 0.032 
mu_LMC = 18.477
sigma_mu_LMC = 0.0263
beta_cal = 4.565
sigma_beta_cal = 0.384
beta = 3.09
sigma_beta = 0.04
alpha = 0.148
sigma_alpha = 0.004

#########################################################################################
#Tomando a matriz inteira C

C_ndiag = np.zeros( (len(C), len(C)) )

for i in range(len(C)):
  for j in range(len(C)):
      C_ndiag[i, j] = C[i, j]

#########################################################################################
#Reescrevendo a matriz L

L_r = np.zeros( (len(C), len(L)) )

for i in range(len(L)):
  for j in range(3492):
    L_r[j, i] = L[i, j]

#########################################################################################
#Inserindo os períodos das cefeidas no vetor Y

for i in range(3130):
  Y[i] = Y[i]+np.array(b_W)*logPm1[i]

#########################################################################################
#Fazendo o ajuste

C_inv = np.linalg.inv(C_ndiag)
cq_inv = L_r.T @ C_inv @ L_r
cq = np.linalg.inv(L_r.T @ C_inv @ L_r)
Q_best_fit = cq @ L_r.T @ C_inv @ Y

print('Delta_mu_NGC4258:', Q_best_fit[37], '+/-', np.sqrt(cq[37, 37]), '(', delta_mu_NGC4258, '+/-', sigma_delta_mu_NGC4258, ')')
print('M_HW:', Q_best_fit[38], '+/-', np.sqrt(cq[38, 38]), '(', M_HW, '+/-', sigma_M_HW, ')')
print('Delta_mu_LMC:', Q_best_fit[39], '+/-', np.sqrt(cq[39, 39]), '(', delta_mu_LMC, '+/-', sigma_delta_mu_LMC, ')')
print('mu_m31:', Q_best_fit[40], '+/-', np.sqrt(cq[40, 40]), '(', mu_m31, '+/-', sigma_mu_m31, ')')
print('b_W:', Q_best_fit[41], '+/-', np.sqrt(cq[41, 41]), '(', b_W, '+/-', sigma_b_W, ')')
print('M_B:', Q_best_fit[42], '+/-', np.sqrt(cq[42, 42]), '(', M_B, '+/-', sigma_M_B, ')')
print('Z_W:', Q_best_fit[43], '+/-', np.sqrt(cq[43, 43]), '(', Z_W, '+/-', sigma_Z_W, ')')
print('delta_zp:', Q_best_fit[45], '+/-', np.sqrt(cq[45, 45]), '(', delta_zp, ')')
print('5*log(H_0):', Q_best_fit[46], '+/-', np.sqrt(cq[46, 46]), '(', fivelogH0, '+/-', sigma_fivelogH0, ')')
print('beta_cal:', Q_best_fit[47], '+/-', np.sqrt(cq[47, 47]), '(', beta_cal, '+/-', sigma_beta_cal, ')')
print('beta:', Q_best_fit[48], '+/-', np.sqrt(cq[48, 48]), '(', beta, '+/-', sigma_beta, ')')
print('alpha:', Q_best_fit[49], '+/-', np.sqrt(cq[49, 49]), '(', alpha, '+/-', sigma_alpha, ')')
print('-------------------------------------------------------------------------------------------------------')
H0 = 10**(Q_best_fit[46]/5)
H0_error = (10**((Q_best_fit[46]/5))*np.log(10))*(np.sqrt(cq[46,46])/5)
print('A constante de Hubble é dada por:', H0, '+/-', H0_error)
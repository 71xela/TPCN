import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

n = np.array([10,20,50,70,100,150,200,500,700,1000])
ldlt = [0.0003, 0.0008, 0.0029, 0.0066, 0.0138, 0.0306, 0.061, 0.304, 0.837, 2.172]
lu = [0.0002, 0.0005, 0.0013, 0.0029, 0.0065, 0.0167, 0.039, 0.354, 1.242, 5.172]

param = linregress(n, np.sqrt(ldlt))
p = param.slope
b = param.intercept
p_err = param.stderr
b_err = param.intercept_stderr
cov = param.rvalue
print(f"slope = {p} +/- {p_err}")
print(f"intercept= {b} +/- {p_err}") 
print(f"correlation = {cov}")

N = (p*n + b)
fig = plt.figure(figsize=(10,8))
plt.plot(n, ldlt, marker='+', label='factorisation LDL^T')
# plt.plot(n, lu, marker='+', label='factorisation LU')
plt.plot(n, N**2, linestyle='--', label='Régression linéaire de racine cubique')
# plt.title("Temps d'execution en fonction de la taille de la matrice")
plt.xlabel("taille n")
plt.ylabel("Temps en secondes")
# plt.xscale('log')
# plt.yscale('log')
plt.legend(loc='best')
#plt.savefig("ldlt.png", dpi=300, format='png')
plt.show()

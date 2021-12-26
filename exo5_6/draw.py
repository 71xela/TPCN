import numpy as np
import matplotlib.pyplot as plt

l = []
with open("conv.txt", "r") as fic:
	rows = fic.readlines()[0].split()

rows = np.array(rows).astype(float)
k = np.arange(0,rows.size,1)

l2 = []
with open("comp.txt", "r") as fic2:
	rows2 = fic2.readlines()[0].split()

rows2 = np.array(rows2).astype(float)
k2 = np.arange(3,103,1)

fig = plt.figure(figsize=(10,8))
plt.plot(k, rows, label="||r||")
plt.xlabel("itération k")
plt.ylabel("Norme du résidu")
#plt.title("Historique de convergence de Jacobi pour Poisson 1D")
plt.yscale('log')
plt.legend(loc='best')
#plt.savefig("jaconv.png", dpi=300, format='png')
plt.show()

p = 15000
print(rows[p], rows[p+4759])

fig = plt.figure(figsize=(10,8))
plt.plot(k2, rows2, label="nombre d'itération")
plt.plot(k2, -8/np.log(np.cos(np.pi/(k2+1))), label="estimation théorique")
plt.xlabel("taille n")
plt.ylabel("Nombre d'itération")
# plt.title("Nombre d'itération nécessaire à Jacobi pour une taille n ")
# plt.yscale('log')
plt.legend(loc='best')
# plt.savefig("jaciter.png", dpi=300, format='png')
plt.show()

import matplotlib.pyplot as plt 
import numpy as np 
dados1 = np.loadtxt(f'./pos1000.dat',dtype=np.float64).T
L=20
angle_cell=np.zeros((L+2,L+2))
for j in range (0,int(L**2)):
	angle_cell[int(dados1[0,j]),int(dados1[1,j])] = dados1[2,j]
	
fig=plt.figure()
plt.imshow(angle_cell,origin='lower')
fig.savefig(r'./foto.png', format='png')
plt.clf()
plt.close()


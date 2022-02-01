import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

mfinpz = np.load("MFI_data.npz")
mfidata = mfinpz['MFI']

semiAng = np.array([15,20,25,30,35,40,45])
vf = np.array([0.1,0.2,0.3,0.5,0.7,1])

sloverdel = 1/vf-1

mfilist1 = np.zeros([semiAng.size,vf.size])
mfilist2 = np.zeros([semiAng.size,vf.size])
for ii in range(semiAng.size):
    for jj in range(vf.size):
        mfilist1[ii,jj] = semiAng[ii]
        mfilist2[ii,jj] = vf[jj]

print(mfilist1)
print(mfilist2)
print(mfidata)
        
fig1 = plt.figure()

f = interpolate.interp2d(semiAng, sloverdel, np.fliplr(np.flipud(np.transpose(mfidata))), kind='linear')

x = np.linspace(5,45,100)
y = np.linspace(min(sloverdel),max(sloverdel),100)
mfinew = f(x,y)

#h = plt.contourf(x, y, mfinew)
h = plt.contourf(semiAng, sloverdel, np.transpose(mfidata))
#plt.axis('scaled')
plt.legend()
plt.xlabel('$Hopper inclination angle, \mu[^o]$')
plt.ylabel('$\lambda/\Delta$')
plt.colorbar()
plt.xlim([30,45])
bbox_inches="tight"
plt.show()
fig1.savefig('designChart.png', dpi=700)

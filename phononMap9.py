import numpy as np
from numpy import polyfit
import matplotlib.pyplot as plt
import re
import time


is_sorted = lambda a: np.all(a[:-1] <= a[1:])

def contri():
	return 0

name="eqlbm.log"

xlo=-6.031575907312943 
xhi=246.03157590731337
ylo=-1.6290947109836174 
yhi=81.62909471098328
zlo=0
zhi=100

nulo=-0.5
nuhi=60
ratiolo=0
ratiohi=0.8

xparts=(xhi-xlo)*2//3
xslab=(xhi-xlo)/xparts
yparts=(yhi-ylo)*2//3
yslab=(yhi-ylo)/yparts
zparts=(zhi-zlo)*2//3
zslab=(zhi-zlo)/zparts

freqlo=0.32
freqhi=0.33

Ntotal=7997

modes = np.array([])
#eigvec = np.array([])
datapos = np.array([])
parti = np.array([])	
parti2 = np.array([])


######################## Reading Coordinates and Atomic ID

with open("str.data","r",encoding='utf-8') as fb:
	coordinates=np.array([])
	j=0
	for line in fb:
		if j==Ntotal:
			break
#		print(j)
		c0 = re.split(r'\s+|\s|: |, |:|,',line)
		c = [ele for ele in c0 if ele.strip()]
		c1 = np.array(c[:])
		coordinates = np.append(coordinates,np.array([c1[0].astype(int), c1[2].astype(np.float64), c1[3].astype(np.float64), c1[4].astype(np.float64)]))
		j += 1
		
coordinates=coordinates.reshape(-1,4)	

print(coordinates)
print('shape of coordinates array=',np.shape(coordinates))

coordinates=coordinates[coordinates[:, 0].argsort()]
print('New coordinates array=\n',coordinates)

######################## Reading Basis and Atomic ID mapping from phonon log file

with open(name,"r",encoding='utf-8') as fc:
	basis=np.array([])
	j=0
	for i in range(15):
#		print(i)
		next(fc)
	flag3=True
	for line in fc:
		if flag3:
			print('first line of basis=',line)
			flag3=False
		if j==Ntotal:
			break
#		print(j)
		c0 = re.split(r'\s+|\s|: |, |:|,',line)
		c = [ele for ele in c0 if ele.strip()]
		c1 = np.array(c[:])
#		print(c1)
		basis = np.append(basis,np.array([c1[3].astype(int), c1[4].astype(int)]))
		j += 1

basis=basis.reshape(-1,2)
print('shape of basis array=',np.shape(basis))

basis=basis[basis[:, 1].argsort()]
print('New basis array=\n',basis)

######################################################################## Map arrays

idsorted = np.append(coordinates,basis,axis=1)
idsorted = idsorted[:,:-1]
print('shape of combined mapping array=',np.shape(idsorted))
print('combined mapping array=',idsorted)

basissorted = idsorted[idsorted[:, 4].argsort()]
print('Basis sorted combined mapping array=',basissorted)


#####################################################################################
	
"""
using map.in and data.pos files in lammps Fixphonon simulation ensures that the atom basis number k and atomic IDs share the same sequence and there is no ambiguity - hence directly the atomic ID can be used to analyze the phonon output without confusion regarding whether to use the k sequence or the atomic ID sequence for it.
"""
startread = time.time()
print('started reading at time (s)=',startread)
with open("eigvec.dat","r",encoding='utf-8') as fa:
	flag = True
	lmda = 0
	f=0
	for line in fa:
		c0 = re.split(r'\s+|\s|: |, |:|,',line)
		c = [ele for ele in c0 if ele.strip()]
		c1 = np.array(c[:])
#		print(c1)
		if flag:
#			print(line)
#			print(c)
#			print(c1)#[-7])#.split(","))
			Ntotal = c1[-1].astype(int)
			sysdim = c1[-7].astype(int)
#			parti=np.zeros([sysdim*Ntotal])
#			modes=np.zeros([sysdim*Ntotal])
			k = 0
#			print(Ntotal,sysdim)
			flag = False
			flag2 = False
			continue
		if k==0:
			if lmda == sysdim*Ntotal:
				break			
			temp = c1[-1].astype(np.float64)
			temp2 = c1
			
			if (temp>freqlo) and (temp<freqhi):
				flag2=True
				contriroot=np.array([])
				f += 1
				selectedmode = temp
				
			partinv = 0
			partinum = 0
#			print(temp)
			for _ in range(1):
				next(fa)
			k = 1
			
#			print(k,lmda)
			continue

#		print(c1)
		n1 = c1.astype(np.float64)
		partinv = partinv + n1[-1]**4
#		partinum = partinum + n1[-1]**2
		if flag2:
			contriroot = np.append(contriroot,n1[-1])
			
#		eigvec = np.append(eigvec,n1)#,axis=0)
		k += 1
		if k==Ntotal+1:
			if flag2:
				flag2=False
			modes = np.append(modes,temp)#,axis=0)
			parti = np.append(parti,1.0/(Ntotal*partinv))
#			parti2 = np.append(parti2,(partinum**2)*1.0/(Ntotal*partinv))
			lmda += 1
			print(lmda)
			k = 0
			for _ in range(1):
				next(fa)
#		print(k)

endread = time.time()
print('time endread=',endread)
print('last frequency line=',temp2)

n1 = np.sum(np.array([1 for i in parti if i <= 500/Ntotal]))
#n2 = np.sum(np.array([1 for i in parti2 if i <= 0.05]))
		
print('are the modes sorted?',is_sorted(modes))
is_sorted(modes)
print('number of modes between freqlo and freqhi=',f)

contrisum=np.sum(contriroot**2)
print('Sum of all atomic contributions to the selected mode=',contrisum)

contri=contriroot**2

##################################################### Two variations

idmap=np.append(idsorted,contri.reshape(-1,1),axis=1)
print('shape of ID mapped array=',np.shape(idmap))
print('ID mapped array=\n',idmap)
basismap=np.append(basissorted,contri.reshape(-1,1),axis=1)
print('shape of basis mapped array=',np.shape(basismap))
print('basis mapped array=\n',basismap)

###################################################### Grouping along increasing x
idmap=idmap[idmap[:, 1].argsort()]
basismap=basismap[basismap[:, 1].argsort()]
slab=xslab
idmaptemp=np.array([])
basismaptemp=np.array([])

poscompare=xlo+slab
temparr=np.array([])
for i in idmap:
#	print(i)
	temp=i[1]
	if temp > poscompare:
#		print(temparr)
		temparr=temparr.reshape(-1,6)
		totalcontri=np.sum(temparr[:,-1])
		averagepos=np.average(temparr,axis=0)
		idmaptemp=np.append(idmaptemp,np.array([len(temparr[:,0])]))
		idmaptemp=np.append(idmaptemp,averagepos[1:4])
		idmaptemp=np.append(idmaptemp,np.array([totalcontri]))
		idmaptemp=np.append(idmaptemp,np.array([totalcontri/len(temparr[:,0])]))
		temparr=np.array([])
		poscompare=poscompare+slab
	temparr=np.append(temparr,i)
#print(idmaptemp)
#idmaptemp=idmaptemp.flatten
idmapx=idmaptemp.reshape(-1,6)

idmaptemp=np.array([])
basismaptemp=np.array([])

poscompare=xlo+slab
temparr=np.array([])
for i in basismap:
#	print(np.shape(i))
	temp=i[1]
	if temp > poscompare:
		temparr=temparr.reshape(-1,6)
		totalcontri=np.sum(temparr[:,-1])
		averagepos=np.average(temparr,axis=0)
		basismaptemp=np.append(basismaptemp,np.array([len(temparr[:,0])]))
		basismaptemp=np.append(basismaptemp,averagepos[1:4])
		basismaptemp=np.append(basismaptemp,np.array([totalcontri]))
		basismaptemp=np.append(basismaptemp,np.array([totalcontri/len(temparr[:,0])]))
		temparr=np.array([])
		poscompare=poscompare+slab
	temparr=np.append(temparr,i)
basismapx=basismaptemp.reshape(-1,6)

print(f'number of x slabs = {len(basismapx[:,0])}')

######################################################## along y

idmap=idmap[idmap[:, 2].argsort()]
basismap=basismap[basismap[:, 2].argsort()]
slab=yslab

idmaptemp=np.array([])
basismaptemp=np.array([])

poscompare=ylo+slab
temparr=np.array([])
for i in idmap:
#	print(i)
	temp=i[2]
	if temp > poscompare:
#		print(temparr)
		temparr=temparr.reshape(-1,6)
		totalcontri=np.sum(temparr[:,-1])
		averagepos=np.average(temparr,axis=0)
		idmaptemp=np.append(idmaptemp,np.array([len(temparr[:,0])]))
		idmaptemp=np.append(idmaptemp,averagepos[1:4])
		idmaptemp=np.append(idmaptemp,np.array([totalcontri]))
		idmaptemp=np.append(idmaptemp,np.array([totalcontri/len(temparr[:,0])]))
		temparr=np.array([])
		poscompare=poscompare+slab
	temparr=np.append(temparr,i)
#print(idmaptemp)
#idmaptemp=idmaptemp.flatten
idmapy=idmaptemp.reshape(-1,6)

idmaptemp=np.array([])
basismaptemp=np.array([])

poscompare=ylo+slab
temparr=np.array([])
for i in basismap:
#	print(np.shape(i))
	temp=i[2]
	if temp > poscompare:
		temparr=temparr.reshape(-1,6)
		totalcontri=np.sum(temparr[:,-1])
		averagepos=np.average(temparr,axis=0)
		basismaptemp=np.append(basismaptemp,np.array([len(temparr[:,0])]))
		basismaptemp=np.append(basismaptemp,averagepos[1:4])
		basismaptemp=np.append(basismaptemp,np.array([totalcontri]))
		basismaptemp=np.append(basismaptemp,np.array([totalcontri/len(temparr[:,0])]))
		temparr=np.array([])
		poscompare=poscompare+slab
	temparr=np.append(temparr,i)
basismapy=basismaptemp.reshape(-1,6)

print(f'number of y slabs = {len(basismapy[:,0])}')

######################################################################

idmap=idmap[idmap[:, 3].argsort()]
basismap=basismap[basismap[:, 3].argsort()]
slab=zslab
idmaptemp=np.array([])
basismaptemp=np.array([])

poscompare=zlo+slab
temparr=np.array([])
for i in idmap:
#	print(i)
	temp=i[3]
	if temp > poscompare:
#		print(temparr)
		temparr=temparr.reshape(-1,6)
		totalcontri=np.sum(temparr[:,-1])
		averagepos=np.average(temparr,axis=0)
		idmaptemp=np.append(idmaptemp,np.array([len(temparr[:,0])]))
		idmaptemp=np.append(idmaptemp,averagepos[1:4])
		idmaptemp=np.append(idmaptemp,np.array([totalcontri]))
		idmaptemp=np.append(idmaptemp,np.array([totalcontri/len(temparr[:,0])]))
		temparr=np.array([])
		poscompare=poscompare+slab
	temparr=np.append(temparr,i)
#print(idmaptemp)
#idmaptemp=idmaptemp.flatten
idmapz=idmaptemp.reshape(-1,6)

idmaptemp=np.array([])
basismaptemp=np.array([])

poscompare=zlo+slab
temparr=np.array([])
for i in basismap:
#	print(np.shape(i))
	temp=i[3]
	if temp > poscompare:
		temparr=temparr.reshape(-1,6)
		totalcontri=np.sum(temparr[:,-1])
		averagepos=np.average(temparr,axis=0)
		basismaptemp=np.append(basismaptemp,np.array([len(temparr[:,0])]))
		basismaptemp=np.append(basismaptemp,averagepos[1:4])
		basismaptemp=np.append(basismaptemp,np.array([totalcontri]))
		basismaptemp=np.append(basismaptemp,np.array([totalcontri/len(temparr[:,0])]))
		temparr=np.array([])
		poscompare=poscompare+slab
	temparr=np.append(temparr,i)
basismapz=basismaptemp.reshape(-1,6)

print(f'number of z slabs = {len(basismapz[:,0])}')

################## End of reading



#print('shape of participation ratio array=',np.shape(parti))
#print(parti)
#print('shape of participation ratio array #2=',np.shape(parti2))
#print(parti2)
#print('shape of frequency mode array=',np.shape(modes))

##################################### Plotting 1
nrows = 2
ncolumns = 1
fig, axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(12,8))
#axes = axes.flatten()

j=0
axes[j].scatter(modes, parti, label=f'Number of localized modes (P$_\lambda$ <= {500/Ntotal}) = {n1}', color='r', s=10)
axes[j].set_xlabel('Phonon mode Frequency, THz')
axes[j].set_ylabel(f'Participation Ratio (P$_\lambda$)')
axes[j].legend(loc='upper right')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')

j=1

axes[j].scatter(modes, parti, label=f'Plot limited to frequency range of {nulo} to {nuhi} THz \nand {ratiolo} to {ratiohi} participation ratios', color='r', s=10)
axes[j].set_xlabel('Phonon mode Frequency, THz')
axes[j].set_ylabel(f'Participation Ratio (P$_\lambda$)')
axes[j].legend(loc='upper right')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')
axes[j].set_xlim(nulo,nuhi)
axes[j].set_ylim(ratiolo,ratiohi)

plt.suptitle(f"Participation ratio")

plt.show()
fig.savefig('Participation ratio9.jpg',dpi=300, bbox_inches='tight')
plt.close()
##################################### Plotting 2

nrows = 3
ncolumns = 2

fig, axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(12,20))
axes = axes.flatten()


j=0
axes[j].plot(basismapx[:,1], basismapx[:,-1],'--bo', label=f'Slabwise phonon contribution \nper atom (based on basis) for \n$\lambda$ = {selectedmode: 0.4f} THz, \nalong x axis')
axes[j].set_xlabel('Position, A')
axes[j].set_ylabel(f'Local contribution')
axes[j].legend(loc='upper left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')

j=1
axes[j].plot(idmapx[:,1], idmapx[:,0], label=f'Number of atoms in each \n{xslab: 0.4f} angstroms slab \nalong x axis', color='r')
axes[j].set_xlabel('Position, A')
axes[j].set_ylabel(f'# of atoms')
axes[j].legend(loc='upper left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')







j=2
axes[j].plot(basismapy[:,2], basismapy[:,-1], '--bo', label=f'Slabwise phonon contribution \nper atom (based on basis) for \n$\lambda$ = {selectedmode: 0.4f} THz, \nalong y axis')
axes[j].set_xlabel('Position, A')
axes[j].set_ylabel(f'Local contribution')
axes[j].legend(loc='upper left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')

j=3

axes[j].plot(idmapy[:,2], idmapy[:,0], label=f'Number of atoms in each \n{yslab: 0.4f} angstroms slab \nalong y axis', color='r')
axes[j].set_xlabel('Position, A')
axes[j].set_ylabel(f'# of atoms')
axes[j].legend(loc='upper left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')




j=4
axes[j].plot(basismapz[:,3], basismapz[:,-1],'--bo', label=f'Slabwise phonon contribution \nper atom (based on basis) for \n$\lambda$ = {selectedmode: 0.4f} THz, \nalong z axis')
axes[j].set_xlabel('Position, A')
axes[j].set_ylabel(f'Local contribution')
axes[j].legend(loc='upper left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')

j=5
axes[j].plot(idmapz[:,3], idmapz[:,0], label=f'Number of atoms in each \n{yslab: 0.4f} angstroms slab \nalong z axis', color='r')
axes[j].set_xlabel('Position, A')
axes[j].set_ylabel(f'# of atoms')
axes[j].legend(loc='upper left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')





plt.suptitle(f"Local phonon contributions")

plt.show()
fig.savefig('Local phonon contributions9.jpg',dpi=300, bbox_inches='tight')
plt.close()


		

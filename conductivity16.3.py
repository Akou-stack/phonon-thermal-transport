import numpy as np
from numpy import polyfit
import matplotlib.pyplot as plt


ts=0.0002
dsteps=1000
slabsExpected=120
dstepsthermo=1000
runsteps=11000000

avpoints=100
av2points=500
avpointsthermo1=100
avpointsthermo2=500

timestart=1000
timeend=2100

distanceFromSource=6
di=int(6/0.5 + 1)

timeStartForAverage=timestart
timeEndForAverage=timeend

#xsize=395.91
xsize=1

################################## Read thermo data file
timethermo=np.array([])
temperatures1d=np.array([])
energy1d=np.array([])

with open("thermo3.dat","r",encoding='utf-8') as fa:
	for _ in range(1):
		next(fa)
	j=0
	totsteps=0
	temp=np.array([])	
	for line in fa:
		c = line.split(" ")
		c1 = np.array(c[:])
		if c1[0]=='#':
			continue
#		print(c1[:-1])
		n1 = c1[:-1].astype(np.float64)
#		print(np.shape(n1))
		temp=np.append(temp,n1,axis=0)
#		print(np.shape(temp))
		timethermo=np.append(timethermo,totsteps*ts)		
		totsteps=totsteps+dstepsthermo

		j +=1
#		print(temp)
#		print('total steps=',totsteps)
thermo2d=temp.reshape(j,11)
temperatures1d=np.copy(thermo2d[:,-1])
energy1d=np.copy(thermo2d[:,1])
ke1d=np.copy(thermo2d[:,3])
pe1d=np.copy(thermo2d[:,2])
xlength=np.average(thermo2d[int(timestart/(ts*dstepsthermo)+1),4]) 	#thermo entries start from 0 timestep
ylength=np.copy(thermo2d[int(timestart/(ts*dstepsthermo)+1),5])
zlength=3.4

stepsTruncated=int((len(energy1d[1:])//avpointsthermo1)*avpointsthermo1)

if stepsTruncated == len(energy1d[1:]):
	energy2d100=energy1d[1:].reshape(-1,avpointsthermo1)
elif stepsTruncated < len(energy1d[1:]):
	energy2d100=energy1d[1:stepsTruncated+1].reshape(-1,avpointsthermo1)
	
stepsTruncated=int((len(energy1d[1:])//avpointsthermo2)*avpointsthermo2)

if stepsTruncated == len(energy1d[1:]):
	energy2d500=energy1d[1:].reshape(-1,avpointsthermo2)
elif stepsTruncated < len(energy1d[1:]):
	energy2d500=energy1d[1:stepsTruncated+1].reshape(-1,avpointsthermo2)


timethermo100=timethermo[avpointsthermo1::avpointsthermo1]
timethermo500=timethermo[avpointsthermo2::avpointsthermo2]
energy100=np.average(energy2d100,axis=1)
energy100std=np.std(energy2d100,axis=1)
energy500=np.average(energy2d500,axis=1)
energy500std=np.std(energy2d500,axis=1)

print(f'\nx length is {xlength: 0.2f} A, y length (width) is {ylength: 0.2f} A, and thickness taken as {zlength: 0.2f} A\n')

################################## Read slab temperature data file

time=np.array([])
temperatures=np.array([])
coordinates=np.array([])
#time=np.append(time,0)

with open("tmp.dat","r",encoding='utf-8') as fa:
	for _ in range(3):
		next(fa)
	i=0
	totsteps=0
	for line in fa:
		c = line.split(" ")
		c1 = np.array(c[:])
		slabs=c1[1].astype(int)
		break
	if slabs==slabsExpected+1:
		slabsFlag = True
	elif slabs==slabsExpected:
		slabsFlag = False
	while True:
		j=0	
		temp=np.array([])
		loopFlag = True	
		for line in fa:
			if j==slabs:
#				print("totsteps=",totsteps)
#				next(fa)
				
				loopFlag = False
				break
			c = line.split(" ")
			c1 = np.array(c[:])
			n1 = c1[2:].astype(np.float64)
#			print(n1)
			temp=np.append(temp,n1,axis=0)
#			print(np.shape(temp))
			j +=1
			loopFlag = False
		if loopFlag:
#			print(line, j)
			break
#		print(temp)
#		print('total steps=',totsteps)
#		print(i)
		temp2d=temp.reshape(slabs,4)
#		print(temp2d)
		if slabsFlag:
			temperatures=np.append(temperatures,temp2d[:-1,3])
			coordinates=np.append(coordinates,temp2d[:-1,1])
			i += 1
		elif slabsFlag==False:
			temperatures=np.append(temperatures,temp2d[:,3])
			coordinates=np.append(coordinates,temp2d[:,1])
			i += 1
		totsteps=totsteps+dsteps
		time=np.append(time,totsteps*ts)
#		if totsteps==runsteps+dsteps:
#			break


#	time=np.append(time,totsteps*ts)
#	print(temperatures)
	stepsReal=totsteps
	print(f'total steps={stepsReal}')
	temperatures2d=temperatures.reshape(int(stepsReal/dsteps),slabsExpected)
	coordinates2d=coordinates.reshape(int(stepsReal/dsteps),slabsExpected)*xsize
	
	xforward=coordinates2d[-1,int(0+di):int(slabsExpected/2+1-di)]
	xreverse=coordinates2d[-1,int(slabsExpected/2+di):int(slabsExpected-di+1)]
	yforward=np.transpose(temperatures2d[:,int(0+di):int(slabsExpected/2+1-di)])
	yreverse=np.transpose(temperatures2d[:,int(slabsExpected/2+di):int(slabsExpected-di+1)])
	
	print(f'shape of final forward step coordinate array={np.shape(xforward)}\nshape of final reverse step coordinate array={np.shape(xreverse)}\nshape of forward range temperature profile={np.shape(yforward)}\nshape of reverse range temperature profile={np.shape(yreverse)}')
	
	slopeforward , dsforward = np.polyfit(xforward,yforward,1,cov=1)
	slopereverse , dsreverse = np.polyfit(xreverse,yreverse,1,cov=1)
	
	print(f'shape of forward slope array={np.shape(slopeforward)}\n\
shape of reverse slope array={np.shape(slopereverse)}\n\
shape of forward cov array={np.shape(dsforward)}\n\
shape of reverse cov array={np.shape(dsreverse)}')

	fa.close()
	
#	print('coordinates=',coordinates2d)



stepsTruncated=int(((stepsReal/dsteps)//avpoints)*avpoints)
print(stepsTruncated,(stepsReal/dsteps))

if stepsTruncated == int(stepsReal/dsteps):
	temp3d=temperatures2d[0:].reshape(-1,avpoints,slabsExpected)
	coord3d=coordinates2d[0:].reshape(-1,avpoints,slabsExpected)
	slopeforward2d=slopeforward[0,:].reshape(-1,avpoints)
	slopereverse2d=slopereverse[0,:].reshape(-1,avpoints)
elif stepsTruncated < int(stepsReal/dsteps):
	temp3d=temperatures2d[0:stepsTruncated+1].reshape(-1,avpoints,slabsExpected)
	coord3d=coordinates2d[0:stepsTruncated+1].reshape(-1,avpoints,slabsExpected)
	slopeforward2d=slopeforward[0,:stepsTruncated+1].reshape(-1,avpoints)
	slopereverse2d=slopereverse[0,:stepsTruncated+1].reshape(-1,avpoints)
	
stepsTruncated2=int(((stepsReal/dsteps)//av2points)*av2points)

if stepsTruncated2 == int(stepsReal/dsteps):
	temp3d2=temperatures2d[0:].reshape(-1,av2points,slabsExpected)
	coord3d2=coordinates2d[0:].reshape(-1,av2points,slabsExpected)
elif stepsTruncated2 < int(stepsReal/dsteps):
	temp3d2=temperatures2d[0:stepsTruncated2+1].reshape(-1,av2points,slabsExpected)
	coord3d2=coordinates2d[0:stepsTruncated2+1].reshape(-1,av2points,slabsExpected)
	
netTemp2d=temperatures2d[int(timestart/(ts*dsteps)):int(timeend/(ts*dsteps))+1,:]
netCoord2d=coordinates2d[int(timestart/(ts*dsteps)):int(timeend/(ts*dsteps))+1,:]

netslopeforwardarray=slopeforward[0,int(timestart/(ts*dsteps)):int(timeend/(ts*dsteps))+1]
netslopereversearray=slopeforward[0,int(timestart/(ts*dsteps)):int(timeend/(ts*dsteps))+1]

stdnetslopeforward=np.std(netslopeforwardarray)
stdnetslopereverse=np.std(netslopereversearray)

netTempProfile=np.average(netTemp2d,axis=0)
netCoordProfile=np.average(netCoord2d,axis=0)
stdNetTempProfile=np.std(netTemp2d,axis=0)
stdNetCoordProfile=np.std(netCoord2d,axis=0)

meanStdNetTempProfile=np.average(stdNetTempProfile)
netTdiff=np.max(netTempProfile)-np.min(netTempProfile)

#print(f'shape of temp array for averaging = {np.shape(temp3d)}')
#print('temperatures second to 12th entry (10000 TS to 110000 TS)=',temperatures2d[1:12])
#print('temperatures 3d first 2 entries=',temp3d[0:2])
avtemp2d=np.average(temp3d, axis=1)
print(f'shape of {dsteps*ts*avpoints} ps averaged=',np.shape(avtemp2d))
#print(f'shape of av temp array = {np.shape(avtemp2d)}')
#print('average temperatures 2d first entry=',avtemp2d[0])

stdtemp2d=np.std(temp3d, axis=1)
#print(f'shape of std array corresponding to av temperatures array = {np.shape(stdtemp2d)}')

avslopeforward=np.average(slopeforward2d, axis=1)
avslopereverse=np.average(slopereverse2d, axis=1)

stdslopeforward=np.std(slopeforward2d, axis=1)
stdslopereverse=np.std(slopereverse2d, axis=1)


ynetforward=np.transpose(netTempProfile[int(0+di):int(slabsExpected/2+1-di)])
ynetreverse=np.transpose(netTempProfile[int(slabsExpected/2+di):int(slabsExpected-di+1)])
	
slopenetforward , dsnetforward = np.polyfit(xforward,ynetforward,1,cov=1)
slopenetreverse , dsnetreverse = np.polyfit(xreverse,ynetreverse,1,cov=1)

#print(f'shape of coord array for averaging = {np.shape(coord3d)}')
#print('coords second to 12th entry (10000 TS to 110000 TS)=',coordinates2d[1:12])
#print('coords 3d first 2 entries=',coord3d[0:2])
avcoord2d=np.average(coord3d, axis=1)
#print(f'shape of av coord array = {np.shape(avcoord2d)}')
#print('average coord 2d first entry=',avcoord2d[0])

stdcoord2d=np.std(coord3d, axis=1)
#print(f'shape of std array corresponding to av coord array = {np.shape(stdcoord2d)}')

av2temp2d=np.average(temp3d2, axis=1)

std2temp2d=np.std(temp3d2, axis=1)

av2coord2d=np.average(coord3d2, axis=1)

std2coord2d=np.std(coord3d2, axis=1)



##################################### (instantaneous) Delta T calculations
#######################			Only the 'time' array is initialized at 0 ps (start of simulation

timeconvergeplot=np.copy(time[avpoints-1::avpoints])

##########
Tmiddle2d=np.copy(temp3d[:,:,int(slabsExpected/2)])
#print(f'half of slabs = {slabs/2} and Tmiddle array shape=',np.shape(Tmiddle2d))
avTmiddle1d=np.average(Tmiddle2d,axis=1)
stdTmiddle1d=np.std(Tmiddle2d,axis=1)
#print('averaged Tmiddle array shape=',np.shape(avTmiddle1d))
#print('averaged Tmiddle array=', avTmiddle1d)

Tfirst2d=np.copy(temp3d[:,:,0])
avTfirst1d=np.average(Tfirst2d,axis=1)
stdTfirst1d=np.std(Tfirst2d,axis=1)

Tinstantdiff2d=np.subtract(Tmiddle2d,Tfirst2d)
avTinstantdiff1d=abs(np.average(Tinstantdiff2d,axis=1))
stdTinstantdiff1d=np.std(Tinstantdiff2d,axis=1)
#############


avTdiff=np.max(avtemp2d,axis=1)-np.min(avtemp2d,axis=1)
spatialMeanStdAvT=np.average(stdtemp2d,axis=1)

	
##################  Plotting

nrows=4
ncolumns=4

fig, axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(18,25))
axes=axes.flatten()
#print(np.shape(axes))

	
j=0	

axes[j].plot(timethermo, temperatures1d, label=f'Mean sample temperature evolution', color='r')
axes[j].set_xlabel('Time, ps')
axes[j].set_ylabel('Temperature, K')
axes[j].legend(loc='upper left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')
j += ncolumns


axes[j].plot(timethermo, energy1d, label=f'Total sample energy evolution', color='b')
axes[j].set_xlabel('Time, ps')
axes[j].set_ylabel('Energy, eV')
axes[j].legend(loc='upper left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')
j += ncolumns



axes[j].plot(timethermo, ke1d, label=f'Total sample KE evolution', color='b')
axes[j].set_xlabel('Time, ps')
axes[j].set_ylabel('Energy, eV')
axes[j].legend(loc='upper left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')
j += ncolumns



axes[j].plot(timethermo, pe1d, label=f'Total sample PE evolution', color='b')
axes[j].set_xlabel('Time, ps')
axes[j].set_ylabel('Energy, eV')
axes[j].legend(loc='upper left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')
j += ncolumns



iTS=int(timestart/(ts*dsteps)-1)
iTSmax=int(iTS+2*avpoints)
iavmax=int((iTSmax+1)/avpoints-1)
iav=int(iavmax-2)
j=1
print(f'iTSmax={iTSmax},iTS={iTS}, iavmax={iavmax}, iav={iav}')

while iav <= iavmax:
#	print(j)
#	print(i)
	axes[j].errorbar(avcoord2d[iav,:], avtemp2d[iav, :], yerr=stdtemp2d[iav, :], label=f'{time[iTS]} ps')
#	axes[j].plot(avcoord2d[iav,:], avtemp2d[iav, :], label=f'{time[iTS]} ps')
	axes[j].set_xlabel('Position, A')
	axes[j].set_ylabel('Temperature, K')
	axes[j].legend(loc='lower right')
	axes[j].adjustable='datalim'
	axes[j].set_aspect('auto')
	
	iav += 1
	iTS += avpoints
	
	axes[j+ncolumns].errorbar(avcoord2d[iav,:], avtemp2d[iav, :], yerr=stdtemp2d[iav, :], label=f'{time[iTS]} ps')
	axes[j+ncolumns].set_xlabel('Position, A')
	axes[j+ncolumns].set_ylabel('Temperature, K')
	axes[j+ncolumns].legend(loc='lower right')
	axes[j+ncolumns].adjustable='datalim'
	axes[j+ncolumns].set_aspect('auto')
	
	iav += 1
	iTS += avpoints

	axes[j+ncolumns*2].errorbar(avcoord2d[iav,:], avtemp2d[iav, :], yerr=stdtemp2d[iav, :], label=f'{time[iTS]} ps')
	axes[j+ncolumns*2].set_xlabel('Position, A')
	axes[j+ncolumns*2].set_ylabel('Temperature, K')
	axes[j+ncolumns*2].legend(loc='lower right')
	axes[j+ncolumns*2].adjustable='datalim'
	axes[j+ncolumns*2].set_aspect('auto')

	iav += 1
	iTS += avpoints
	j += ncolumns*3
	

#meanStdNetTempProfile=np.average(stdNetTempProfile)

	
axes[j].plot(timeconvergeplot, avTdiff, label=f'{avpoints*ts*dsteps} ps averaged $\Delta$ T',color='r')
axes[j].axhline(y=netTdiff, label=f'{timestart} to {timeend} ps averaged $\Delta$ T',color='orange')
axes[j].set_ylim(bottom=0,)
axes[j].set_xlabel('Time, ps')
axes[j].set_ylabel('Temperature, K')
axes[j].legend(loc='lower left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')

j=2

iTSmax=int(timeend/(ts*dsteps)-1)
iTS=int(iTSmax-2*avpoints)
iavmax=int((iTSmax+1)/avpoints-1)
iav=int(iavmax-2)

print(f'iTSmax={iTSmax},iTS={iTS}, iavmax={iavmax}, iav={iav}')

while iav <= iavmax:
#	print(j)
#	print(i)
	axes[j].errorbar(avcoord2d[iav,:], avtemp2d[iav, :], yerr=stdtemp2d[iav, :], label=f'{time[iTS]} ps')
#	axes[j].plot(avcoord2d[iav,:], avtemp2d[iav, :], label=f'{time[iTS]} ps')
	axes[j].set_xlabel('Position, A')
	axes[j].set_ylabel('Temperature, K')
	axes[j].legend(loc='lower right')
	axes[j].adjustable='datalim'
	axes[j].set_aspect('auto')
	
	iav += 1
	iTS += avpoints
	
	axes[j+ncolumns].errorbar(avcoord2d[iav,:], avtemp2d[iav, :], yerr=stdtemp2d[iav, :], label=f'{time[iTS]} ps')
	axes[j+ncolumns].set_xlabel('Position, A')
	axes[j+ncolumns].set_ylabel('Temperature, K')
	axes[j+ncolumns].legend(loc='lower right')
	axes[j+ncolumns].adjustable='datalim'
	axes[j+ncolumns].set_aspect('auto')
	
	iav += 1
	iTS += avpoints

	axes[j+ncolumns*2].errorbar(avcoord2d[iav,:], avtemp2d[iav, :], yerr=stdtemp2d[iav, :], label=f'{time[iTS]} ps')
	axes[j+ncolumns*2].set_xlabel('Position, A')
	axes[j+ncolumns*2].set_ylabel('Temperature, K')
	axes[j+ncolumns*2].legend(loc='lower right')
	axes[j+ncolumns*2].adjustable='datalim'
	axes[j+ncolumns*2].set_aspect('auto')

	iav += 1
	iTS += avpoints
	j += ncolumns*3
	

	
axes[j].plot(timeconvergeplot, spatialMeanStdAvT, label=f'Mean Std. of slab T\'s\nacross {avpoints*ts*dsteps} ps',color='r')
axes[j].axhline(y=meanStdNetTempProfile, label=f'Mean Std. of {timestart} to {timeend} ps\naveraged slab T',color='orange')
#axes[j].set_ylim(bottom=0,)
axes[j].set_xlabel('Time, ps')
axes[j].set_ylabel('Temperature, K')
axes[j].legend(loc='lower left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')

j=3

axes[j].plot(timethermo100, energy100, label=f'{avpointsthermo1*ts*dstepsthermo} ps averaged E', color='b')
#axes[j].set_ylim(bottom=0,)
axes[j].set_xlabel('Time, ps')
axes[j].set_ylabel('Energy, eV')
axes[j].legend(loc='upper center')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')

j += ncolumns
axes[j].plot(timethermo100, energy100std, label=f'Std. in {avpointsthermo1*ts*dstepsthermo} ps averaged E', color='b')
#axes[j].set_ylim(bottom=0,)
axes[j].set_xlabel('Time, ps')
axes[j].set_ylabel('Energy, eV')
axes[j].legend(loc='upper center')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')

j += ncolumns
axes[j].plot(timeconvergeplot, avslopeforward, label=f'{avpoints*ts*dsteps} ps averaged forward slope',color='r')
axes[j].plot(timeconvergeplot, -avslopereverse, label=f'{avpoints*ts*dsteps} ps averaged reverse slope',color='orange')
axes[j].axhline(y=slopenetforward[0], label=f'{timestart} to {timeend} ps averaged forward slope',color='r', linestyle=':')
axes[j].axhline(y=-slopenetreverse[0], label=f'{timestart} to {timeend} ps averaged reverse slope',color='orange', linestyle=':')
axes[j].set_ylim(bottom=0,)
axes[j].set_xlabel('Time, ps')
axes[j].set_ylabel('dT/dx, K/A')
axes[j].legend(loc='lower left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')

j += ncolumns
axes[j].plot(timeconvergeplot, stdslopeforward, label=f'Std. of forward\ndT/dx across {avpoints*ts*dsteps} ps',color='r')
axes[j].axhline(y=stdnetslopeforward, label=f'Std. of forward\ndT/dx across {timestart} to {timeend} ps',color='r',linestyle=':')
axes[j].plot(timeconvergeplot, stdslopereverse, label=f'Std. of reverse\ndT/dx across {avpoints*ts*dsteps} ps',color='orange')
axes[j].axhline(y=stdnetslopereverse, label=f'Std. of reverse\ndT/dx across {timestart} to {timeend} ps',color='orange',linestyle=':')
#axes[j].set_ylim(bottom=0,)
axes[j].set_xlabel('Time, ps')
axes[j].set_ylabel('Std in dT/dx, K/A')
axes[j].legend(loc='lower left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')


plt.suptitle(f"Convergence of thermodynamic parameters to measure conductivity")

plt.show()
fig.savefig('Conductivity-temperature energy evolution16.jpg',dpi=300, bbox_inches='tight')
plt.close()	

####################################### Temperature profiles throughout averaging duration

nrows=2
ncolumns=4

fig, axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(16,25))
axes=axes.flatten()
#print(np.shape(axes))

dtGap=(timeend-timestart-av2points*ts*dsteps)/5

iTSmax=int(timeend/(ts*dsteps)-1)
iTS=int((timestart+av2points*ts*dsteps)/(ts*dsteps)-1)
diTS=int(dtGap/(ts*dsteps))

iavmax=int((iTSmax+1)/av2points-1)
iav=int((iTS+1)/av2points-1)
diav=int((dtGap/(ts*dsteps))/av2points)

j=0
print(f'\niTSmax2={iTSmax},iTSmin2={iTS}, dtGap={dtGap}, iavmax2={iavmax}, iavmin2={iav}, diav2={diav}')

while iav <= iavmax:
#	print(j)
#	print(i)
#	axes[j].errorbar(avcoord2d[iav,:], avtemp2d[iav, :], yerr=stdtemp2d[iav, :], label=f'{time[iTS]} ps')
	print(f'j={j},iav={iav}')
	axes[j].plot(av2coord2d[iav,:], av2temp2d[iav, :], label=f'{time[iTS]} ps')
	axes[j].set_xlabel('Position, A')
	axes[j].set_ylabel('Temperature, K')
	axes[j].legend(loc='lower right')
	axes[j].adjustable='datalim'
	axes[j].set_aspect('auto')
	
	j += 1
	iav += diav
	iTS += diTS
	
#netTempProfile=np.average(netTemp2d,axis=0)
#netCoordProfile=np.average(netCoord2d,axis=0)
#stdNetTempProfile=np.std(netTemp2d,axis=0)
#stdNetCoordProfile=np.std(netCoord2d,axis=0)


axes[j].plot(netCoordProfile, netTempProfile, label=f'{timestart} to {timeend} ps averaged')
axes[j].set_xlabel('Position, A')
axes[j].set_ylabel('Temperature, K')
axes[j].legend(loc='lower right')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')
		
plt.suptitle(f"Temperature profiles (100 ps averaged) through the k calculation window")

plt.show()
fig.savefig('Temp profile through averaging16.jpg',dpi=300, bbox_inches='tight')
plt.close()	

######################################### With Error Bars

nrows=2
ncolumns=4

fig, axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(18,25))
axes=axes.flatten()
#print(np.shape(axes))

dtGap=(timeend-timestart-av2points*ts*dsteps)/5

iTSmax=int(timeend/(ts*dsteps)-1)
iTS=int((timestart+av2points*ts*dsteps)/(ts*dsteps)-1)
diTS=int(dtGap/(ts*dsteps))

iavmax=int((iTSmax+1)/av2points-1)
iav=int((iTS+1)/av2points-1)
diav=int((dtGap/(ts*dsteps))/av2points)

j=0
print(f'\niTSmax2={iTSmax},iTSmin2={iTS}, dtGap={dtGap}, iavmax2={iavmax}, iavmin2={iav}, diav2={diav}')


while iav <= iavmax:
#	print(j)
#	print(i)
	axes[j].errorbar(av2coord2d[iav,:], av2temp2d[iav, :], yerr=std2temp2d[iav, :], label=f'{time[iTS]} ps')
#	axes[j].plot(av2coord2d[iav,:], av2temp2d[iav, :], label=f'{time[iTS]} ps')
	axes[j].set_xlabel('Position, A')
	axes[j].set_ylabel('Temperature, K')
	axes[j].legend(loc='lower right')
	axes[j].adjustable='datalim'
	axes[j].set_aspect('auto')
	
	j += 1
	iav += diav
	iTS += diTS
	
#stdNetTempProfile=np.std(netTemp2d,axis=0)
#stdNetCoordProfile=np.std(netCoord2d,axis=0)

axes[j].errorbar(netCoordProfile, netTempProfile, yerr=stdNetTempProfile, label=f'{timestart} to {timeend} ps averaged')
axes[j].set_xlabel('Position, A')
axes[j].set_ylabel('Temperature, K')
axes[j].legend(loc='lower right')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')

		
plt.suptitle(f"Temperature profiles (100 ps averaged) through the k calculation window")

plt.show()
fig.savefig('Temp profile with error through averaging16.jpg',dpi=300, bbox_inches='tight')
plt.close()	


################################################### Average conductivity calculation

heatstr=input("\nPlease input the Heating rate in eV/ps: ")
heat=float(heatstr)
power=heat*1.602176634		#converting to e-19 J/ps
heatflux=power/(ylength*zlength)

#timestartstr=input("Please give the start of the time range for averaging (in ps): ")
#timestart=float(timestartstr)

#timeendstr=input("Please give the end of the time range for averaging (in ps): ")
#timeend=float(timeendstr)

#TSindexstart=timestart/(ts*dsteps)				#this index is easier to compare with the timestep
#avindexstart=int(TSindexstart/avpoints-1)		#corresponding index of averaged quantities
#TSindexend=timeend/(ts*dsteps)				#this index is easier to compare with the timestep
#avindexend=int(TSindexend/avpoints-1)

T1=np.max(netTempProfile)
T2=np.min(netTempProfile)
delT=abs(T1-T2)

print(f'\naverage Tmax between {timestart} ps to {timeend} ps = {np.max(netTempProfile): 0.2f} K')# +/- {np.average(stdTmiddle1d[avindexstart:avindexend])} which is at slab 11')

print(f'average Tmin between {timestart} ps to {timeend} ps = {np.min(netTempProfile): 0.2f} K')# +/- {stdTfirst1d[avindex]} which is at slab 1')

print(f'\nTemp difference = {netTdiff: 0.2f} K over {(xlength/2)/10: 0.2f} nm\nHeat flux = {heatflux: 0.6f} e+13 W/m2')# +/- {stdTinstantdiff1d[avindex]} K')	
	
k = heatflux*0.5*xlength*1000/netTdiff

kslopeforward = heatflux*1000/slopenetforward[0]

kslopereverse = -heatflux*1000/slopenetreverse[0]

kslopeaverage = (kslopeforward + kslopereverse)/2

G = k/(0.5*xlength*0.1)

print(f'k = {k: 0.4f} W/m.K \n')

de, deVar = np.polyfit(timethermo500, energy500, 1, cov=True)

print(f'Rate of change of system total energy is {de[0]: 0.8f} +/- {np.sqrt(deVar[0][0]): 0.8f} eV/ps') 

############################################################ WRITE RESULTS

fa=open('k_log.txt','w')
fa.write(f'Average heat exchange rate in one direction of symmetric graded sample is {heat} eV/ps.\nx length is {xlength: 0.2f} A, y length (width) is {ylength: 0.2f} A, and thickness taken as {zlength: 0.2f} A\nAverage Tmin between {timestart} ps to {timeend} ps = {T2: 0.2f} K\nTemp difference = {delT: 0.2f} K over {(xlength/2)/10: 0.2f} nm\nMean Std. of slab T\'s from {timestart} to {timeend} ps={meanStdNetTempProfile} K\nwhich is {100*meanStdNetTempProfile/delT: 0.2f} % of the temperature difference\nHeat flux = {heatflux: 0.6f} e+13 W/m2\nk = {k: 0.4f} W/m.K		G = {G: 0.4f} e9 W/m2.K\nconsidering {timestart} ps to {timeend} ps averaged temperature profile atleast {distanceFromSource} nm away from source and sink,\nk from forward region slope = {kslopeforward: 0.4f} W/m.K\nk from reverse region slope = {kslopereverse: 0.4f} W/m.K\naverage k from slope =  {kslopeaverage: 0.4f} W/m.K\n\nRate of change of system total energy is {de[0]: 0.8f} +/- {np.sqrt(deVar[0][0]): 0.8f} eV/ps')
fa.close()
print("\n\n#############################################################\n\n")
print(f'Average heat exchange rate in one direction of symmetric graded sample is {heat} eV/ps.\nx length is {xlength: 0.2f} A, y length (width) is {ylength: 0.2f} A, and thickness taken as {zlength: 0.2f} A\nAverage Tmin between {timestart} ps to {timeend} ps = {T2: 0.2f} K\nTemp difference = {delT: 0.2f} K over {(xlength/2)/10: 0.2f} nm\nMean Std. of slab T\'s from {timestart} to {timeend} ps={meanStdNetTempProfile} K\nwhich is {100*meanStdNetTempProfile/delT: 0.2f} % of the temperature difference\nHeat flux = {heatflux: 0.6f} e+13 W/m2\nk = {k: 0.4f} W/m.K		G = {G: 0.4f} e9 W/m2.K\nconsidering {timestart} ps to {timeend} ps averaged temperature profile atleast {distanceFromSource} nm away from source and sink,\nk from forward region slope = {kslopeforward: 0.4f} W/m.K\nk from reverse region slope = {kslopereverse: 0.4f} W/m.K\naverage k from slope =  {kslopeaverage: 0.4f} W/m.K\n\nRate of change of system total energy is {de[0]: 0.8f} +/- {np.sqrt(deVar[0][0]): 0.8f} eV/ps')

print("\n\n#############################################################\n\n")



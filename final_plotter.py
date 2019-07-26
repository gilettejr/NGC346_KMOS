import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from astropy import units as u
from astropy import constants as const
from astropy.visualization import quantity_support
quantity_support()


def flux_to_ap_mag(flux):
	
	zpoint = 2
	
	m=zpoint-2.5*np.log10(flux)
	
	return m
	
#function takes fluxes in janskys and distance in lightyears

def flux_to_lum(flux,distance):
	
	#unit conversions
	
	bgammawave = 2166 * u.um
	
	bgammafreq = (const.c)/bgammawave
	
	bgammafreq=bgammafreq.to(u.Hz)
	
	flux = flux * u.Jy
	
	distance = distance * u.lyr
	
	#basic formula for converting flux density into luminosity
	
	lum=flux*bgammafreq*4*np.pi*distance**2
	
	lum = lum.to(u.L_sun)
	
	
	
	return lum
	
	
#function for converting fluxes in different bands 
	
def make_bolometric(flux_list,band_nus):
	
	#bandnames=[U,B,V,I,J,H,K,3.6,4.5,5.8,8.20,24.0]
	#Y553,Y535,Y545
	#band wavelengths in microns
	fbolist=[]
	for i in range(len(flux_list)-1):
		
		fbolist.append((band_nus[i]-band_nus[i+1])*(flux_list[i]+(flux_list[i+1])))
	
	fbo=0.5*sum(fbolist)
	lbol=4*np.pi*(210000*u.lyr)**2*fbo
	lbol=lbol.to(u.L_sun)
	#bands taken from https://www.eso.org/public/images/eso9934b/, should be at least roughly correct
	#integration technique from https://arxiv.org/pdf/1608.08631.pdf
	
	return lbol


bandfluxes=np.array([[0,0,0,0.018,0.345,0.353,0.432,4.669,1.32,5.4,42.886,409.8],
	[4.597,3.85,0.325,1.982,2.556,4.301,9.342,25.6,30.6,37.6,81.9,2048],
	[1.697,1.409,0.897,0.662,1.077,0.785,0.82,1.882,0.643,1.87,13.791,110.4]])
bandlambdas=np.array([0.38,0.44,0.55,0.90,1.25,1.65,2.16,3.6,4.5,5.8,8,24])
bandlambdas=bandlambdas*u.um
bandnus=const.c/bandlambdas
bandnus=bandnus.to(u.Hz)
bandfluxes=bandfluxes*u.mJy

sns.set()


#Y553 = make_bolometric(bandfluxes[0],bandnus)
#Y535 = make_bolometric(bandfluxes[1],bandnus)
#Y545 = make_bolometric(bandfluxes[2],bandnus)

#print(Y553)
#print(Y535)
#print(Y545)

line_lists=[0.0007365052610700439,0.0019782488323760506,0.0005496097015496698,0.0053499198102429715,0.00010559136590432348,0.011703934704719818,0.000684790482743213]
x_data=np.array([5910,8160,3160,3240,6000,40280,5180])



y_data=np.array([0.13208818,0.35478809,0.09856949,0.95947881,0.01893723,2.09903657,0.12281342])


#for i in line_lists:
	
	#y_data.append(flux_to_lum(i,210000))


	
print(y_data)
print(x_data)

locus1= x_data * np.exp(-6.6775)
locus01= x_data*0.1 * np.exp(-6.6775)
locus001 = x_data*0.01* np.exp(-6.6775)

plt.scatter(x_data,y_data,label='KMOS data')
plt.plot(x_data,locus1,label='Lacc=Lbol')
plt.plot(x_data,locus01,label='Lacc=0.1Lbol')
plt.plot(x_data,locus001,label='Lacc=0.01Lbol')
plt.legend()

plt.xlabel('Lbol/L_sun')
plt.ylabel('Lbgamma/L_sun')
plt.yscale('log')
plt.xscale('log')
plt.show()







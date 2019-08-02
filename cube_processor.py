import numpy as np
import pandas as pd
import glob, os                          # Operating system module and pathnames
import statistics
import matplotlib
#matplotlib.use('nbagg')
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.io import ascii, fits
from astropy.table import Table, Column
from astropy.coordinates import Angle, Latitude, Longitude, SkyCoord
from astropy.convolution import convolve, Box1DKernel
from astropy.stats import sigma_clipped_stats, mad_std
from astropy.modeling import models
from photutils import DAOStarFinder
from astropy.wcs import WCS
from specutils import Spectrum1D,SpectralRegion
from specutils.manipulation import noise_region_uncertainty,extract_region
from specutils.fitting import find_lines_threshold,fit_lines,estimate_line_parameters,fit_generic_continuum
from specutils.analysis import centroid,fwhm



from cube_class import cube_proc



# To deal with spectralcubes
from spectral_cube import SpectralCube, BooleanArrayMask

# Set colour palette & plot layouts
import seaborn as sns     #; sns.set()
sns.set_context("paper")

#class for named object


	
	#function to show collapsed cubes of object at given waveband
	
	
		
#function for checking uncombined observing block images
	
def check_single_image():
	
	#while loop so that user can continuously carry out the process on different stars without having to rerun the script
	
	cont = 'yes'
	
	#input at end of script allows user to terminate program
	
	while cont != 'no':
	
		#initialise class
	
		star=cube_proc()
		
		#assign variables from class initialisaton
		
		star_name=star.star_name
		band=star.band
		
		#function for finding path to each combined cube is called
		
		object_path_list=star.single_file_list(star_name,band)
		
		
		
		#loop for filling list with single spectra
		
		image_list=[]
		
		for i in object_path_list:
			
			#function for collapsing each cube into image called
			
			[image,cube] = star.collapse_into_image(i)
			
			#image list created
			
			image_list.append(image)
		
		#input to choose which image to view
		
		print('Which out of the '+str(len(image_list)) + ' images would you like to view?')
		
		x=input()
		
		#chosen figure plotted
		
		fig = plt.figure()
		plt.imshow((image_list[int(x)-1]).value)
		plt.tight_layout()


		plt.show()
		
#function for checking images of combined data
			
def check_combined_image():
			
	#while loop so that user can continuously carry out the process on different stars without having to rerun the script
	
	cont = 'yes'
	
	#input at end of script allows user to terminate program
	
	while cont != 'no':
		
		#initialise class
		
		star=cube_proc()
		
		#variables assigned from class initialisation
		
		star_name=star.star_name
		band=star.band
		
		#path to file set
		
		path=star.path_to_combined_file(star_name,band)
		
		
		#function for creating image called
		
		[image,cube] = star.collapse_into_image(path)
		
		#image plotted
		
		fig=plt.figure()
		plt.imshow(image.value)
		plt.tight_layout()
		plt.show()
	

def combine_and_plot():
	
	
	
	#while loop so that user can continuously carry out the process on different stars without having to rerun the script
	
	cont = 'yes'
	
	#input at end of script allows user to terminate program
	
	while cont != 'no':
		
		star=cube_proc()
		
		star_name=star.star_name
		band=star.band
		
		#function for finding path to each combined cube is called
		
		object_path_list=star.single_file_list(star_name,band)
		
		#list for single spectra is initialised
		
		list_of_spectra=[]
		
		#loop for filling list with single spectra
		
		for i in object_path_list:
			
			#function for collapsing each cube into image called
			
			[image,cube] = star.collapse_into_image(i)
			
			#function for creating spectrum from selected stars called
			
			plotting_data=star.construct_spectrum(image,cube)
			
			#construct_spectrum returns masked cube object and background corrected spectrum object
		
			maskedcube=plotting_data[0]
			correctedSp_Jy=plotting_data[1]
			uncorrectedSp_Jy=plotting_data[2]
			
			#single spectrum added to list
			list_of_spectra.append((plotting_data[1]).value)
			list_of_unspectra.append((plotting_data[2]).value)
		
			#list for medians of each point in spectrum initialised
							

			
		median_spectra=[]
		
		#reshaper used to trim spectra to prevent NAN blocks at either end, which return an error. Each star has a unique number of NAN blocks, so reshaper has to be dependent on star_name variable
		
		reshaper = 0
		
		if star_name == 'Y535':
			
			reshaper = 1
			
		elif star_name=='Y533':
			
			reshaper = 18
			
		elif star_name=='Y532':
			
			reshaper = 3
		uncmed=[]
		
		#reshaper variable applied to trim NAN blocks in loop
		
		for j in range(len((list_of_spectra)[0])-reshaper):
			
			#median of each point from four single spectra found
			combined = [list_of_spectra[0][j],list_of_spectra[1][j],list_of_spectra[2][j],list_of_spectra[3][j]]
			uncmed.append((np.std(combined)/2) * 1.2533)
			median = statistics.median(combined)
			#median of each point added to list
			median_spectra.append(median)
		
		#median spectrum figure initialised
			
		fig = plt.figure(figsize=(14,6))
		
		#further reshaping required, again due to NAN blocks in positions which result in errors
		
		if star_name=='Y551' or star_name=='Y538' or star_name=='Y545':
			
			spectral_axis=(maskedcube.spectral_axis.value)[1:]
		
		#only applies to three objects, otherwise no change is required
		
		else:
			
			spectral_axis=(maskedcube.spectral_axis.value)
			
			
			median_unspectra=[]
		
		#reshaper used to trim spectra to prevent NAN blocks at either end, which return an error. Each star has a unique number of NAN blocks, so reshaper has to be dependent on star_name variable
		
		reshaper = 0
		
		if star_name == 'Y535':
			
			reshaper = 1
			
		elif star_name=='Y533':
			
			reshaper = 18
			
		elif star_name=='Y532':
			
			reshaper = 3
		
		
		#reshaper variable applied to trim NAN blocks in loop
		
		for j in range(len((list_of_unspectra)[0])-reshaper):
			
			#median of each point from four single spectra found
			combined = [list_of_unspectra[0][j],list_of_unspectra[1][j],list_of_unspectra[2][j],list_of_unspectra[3][j]]
			
			median = statistics.median(combined)
			#median of each point added to list
			median_unspectra.append(median)
		
		#median spectrum figure initialised
			
		fig = plt.figure(figsize=(14,6))
		
		#further reshaping required, again due to NAN blocks in positions which result in errors
		
		if star_name=='Y551' or star_name=='Y538' or star_name=='Y545':
			
			spectral_axis=(maskedcube.spectral_axis.value)[1:]
		
		#only applies to three objects, otherwise no change is required
		
		else:
			
			spectral_axis=(maskedcube.spectral_axis.value)
			
		#median spectrum plotted and labelled
		
		plt.plot(spectral_axis,median_spectra, label='Bkg corr Source') # background corrected source spectrum 

		plt.xlabel('Wavelength (microns)')
		plt.ylabel(correctedSp_Jy.unit)
			
			#xlim for h band
		
		#x limits set from spectral axis
		
		plt.xlim(min(maskedcube.spectral_axis).value, (max(maskedcube.spectral_axis).value))
		
		#y limits chosen that best frame the spectrum 
		
		plt.ylim(-0.1e-3, 0.15e-2)
			
			

		plt.legend(frameon=False, fontsize='medium')

		plt.tight_layout()
		plt.show()
		
		#figure saved if required
		
		#plt.savefig('Spectra/'+star_name + '_bkg corrected spectrum_' + band +'_band.png')
		
		#spectrum object created for line fitting
		
		spectrum =Spectrum1D(spectral_axis=spectral_axis*u.micron,flux=median_spectra*u.Jansky)
		unspectrum =Spectrum1D(spectral_axis=spectral_axis*u.micron,flux=median_unspectra*u.Jansky)
		
		
		
	
		
		noise_region=SpectralRegion(2.26*u.um,2.3*u.um)
		spectrum=noise_region_uncertainty(spectrum,noise_region)
		uncmed=np.array(uncmed)
		uncmed=uncmed *u.Jy
		#unc = spectrum.uncertainty[0]
		f=open(star_name+band+ '_bkg_corrected_spectrum.csv','a')
		f.write('Wavelength(um),Flux(Jy),Fluxunc(Jy)\n')
		for i in range(len(uncmed)):
			f.write(str(spectrum.spectral_axis[i])+','+str(spectrum.flux[i]) + ',' + str(uncmed[i])+'\n')
		#errors = np.sqrt((unc)**2 + uncmed**2) 
		#function to apply single gaussian fit applied. Limits for subregion containing BGamma line used.
		
		#star.gaussian_fit(spectrum,2.162,2.172,star_name)
		f.close()
		
		noise_region=SpectralRegion(2.26*u.um,2.3*u.um)
		unspectrum=noise_region_uncertainty(unspectrum,noise_region)
	
		#unc = spectrum.uncertainty[0]
		f=open(star_name+band+ '_uncorrected_spectrum.csv','a')
		f.write('Wavelength(um),Flux(Jy),Fluxunc(Jy)\n')
		for i in range(len(uncmed)):
			f.write(str(unspectrum.spectral_axis[i])+','+str(unspectrum.flux[i]) + ',' + str(uncmed[i])+'\n')
		#errors = np.sqrt((unc)**2 + uncmed**2) 
		#function to apply single gaussian fit applied. Limits for subregion containing BGamma line used.
		
		#star.gaussian_fit(spectrum,2.162,2.172,star_name)
		f.close()
		
		
def show_spectra():
	

	
	#while loop so that user can continuously carry out the process on different stars without having to rerun the script
	
	cont = 'yes'
	
	#input at end of script allows user to terminate program
	
	while cont != 'no':
		
		#initialise class
		
		star=cube_proc()
		
		#variables set from class initialisation
		
		star_name=star.star_name
		band=star.band
		
		#not necessary to work with K data, spectra created using combine_and_plot function
		
		if band!='DUH':
			
			#path to file found
			
			path=star.path_to_combined_file(star_name,band)
			
			#image created
		
			[image,cube] = star.collapse_into_image(path)
			
			#function for creating spectrum from selected stars called
			
			plotting_data=star.construct_spectrum(image,cube)
			
			#construct_spectrum returns masked cube object and background corrected spectrum object
		
			maskedcube=plotting_data[0]
			correctedSp_Jy=plotting_data[1]
			uncorrectedSp_Jy = plotting_data[2]
			
			#datasets assigned
			
			spectral_axis=maskedcube.spectral_axis
			flux = correctedSp_Jy.value
			
			spectrum =Spectrum1D(spectral_axis=maskedcube.spectral_axis.value*u.micron,flux=correctedSp_Jy.value*u.Jansky)
			
			if band=='H':
				lower=1.515
				upper=1.535
				
			elif band=='YJ':
				lower = 1.172
				upper = 1.182
			
			noise_region=SpectralRegion(lower*u.um,upper*u.um)
			spectrum=noise_region_uncertainty(spectrum,noise_region)
			
			f=open(star_name+band+ '_bkg_corrected_spectrum.csv','a')
			f.write('Wavelength(um),Flux(Jy),unc=' + str(spectrum.uncertainty[0]) +'\n')
			for i in range(len(spectrum.spectral_axis)):
				f.write(str(spectrum.spectral_axis[i])+','+str(spectrum.flux[i]) + '\n')
				
			
			spectral_axis=maskedcube.spectral_axis
			flux = correctedSp_Jy.value
			
			unspectrum =Spectrum1D(spectral_axis=maskedcube.spectral_axis.value*u.micron,flux=uncorrectedSp_Jy.value*u.Jansky)
			
			if band=='H':
				lower=1.114
				upper=1.124
				
			elif band=='YJ':
				lower = 1.172
				upper = 1.182
			
			noise_region=SpectralRegion(lower*u.um,upper*u.um)
			unspectrum=noise_region_uncertainty(unspectrum,noise_region)
			
			f=open(star_name+band+ '_uncorrected_spectrum.csv','a')
			f.write('Wavelength(um),Flux(Jy),unc=' + str(unspectrum.uncertainty[0]) +'\n')
			for i in range(len(unspectrum.spectral_axis)):
				f.write(str(unspectrum.spectral_axis[i])+','+str(unspectrum.flux[i]) + '\n')
			
			#spectrum plotted
			
			#H1.114-1.124
			#YJ1.172-1.182
			
			
			
			fig = plt.figure(figsize=(14,6))
			plt.plot(spectral_axis,flux)
			plt.title((star_name)+' spectrum, ' + band + ' Band')
			plt.show()
		
		else:
			
			print('Nah, this does not work for K. Use 4 at the start')
			
			
			
#main function. Interactive input for choice of which function to use			
		
def main():
	
	print('\nWhat would you like to do with your datacubes?\n')
	
	print('1: Check Images from Single Observing Blocks\n')
	print('2: Check Images from Combined Observing Blocks\n')
	print('3: Check Combined Spectra from H and YJ bands\n')
	print('4: Combine Single OB data into spectra and fit Gaussian to Bracket Gamma\n')
	
	
	x=input()
	
	if x == '1':
		
		check_single_image()
		
	elif x=='2':
		
		check_combined_image()
		
	elif x=='3':
		
		show_spectra()
		

		
	elif x=='4':
		
		combine_and_plot()
		
main()


		
		
		
		
		

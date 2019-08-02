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
from spectral_cube import SpectralCube, BooleanArrayMask

#function to create image object from reduced datacube

class cube_proc:
	
	def __init__(self):
		
		#catalogue of objects opened and printed to select from
		
		f = open('starcat.txt','r')
		
		lines = f.readlines()
				
		f.close()
		
		print('Enter the object name from the following list:')
		
		for i in lines:
			print(i)

		#object name from input assigned to variable

		x = input()
		
		print('Enter the waveband you would like to view. Type all to view all 3')
		
		#band from input assigned to variable
		
		y = input()
		
		self.star_name=x
		self.band=y
		
	
	def path_to_combined_file(self,star_name,band):
			
	#path to directory of KMOS combined output data
	
		root_path = 'reflex_data/reflex_end_products/2019-06-20T16:49:58/'
		
		#directories for each individual band
		
		k_path = 'KMOS.2018-09-29T02:15:14.819_combine_OBs/'
		h_path = 'KMOS.2018-09-30T02:08:21.672_combine_OBs/'
		yj_path = 'KMOS.2018-09-30T05:42:33.077_combine_OBs/'
		
		#lists to be used for looping to get correct files
		
		bands = ['K','H','YJ']
		
		bandfiles=[k_path, h_path, yj_path]
		
		#looped if statements to retrieve the correct path depending on band
		
		if band=='H':
			j=1
			
		elif band=='YJ':
			j=2
			
		elif band == 'K':
			j=0
		
			
			#function output constructed from paths and filename from the loop
		
		filename = root_path + bandfiles[j]+'NGC346_'+bands[j]+'_2c_COMBINED_CUBE_'+star_name+'.fits'
		
		return filename
				
#function returns list of paths to each combined datacube for each observing block

	def single_file_list(self,star_name,band):
				
		#path to directory containing individual observing block data

		root_path = 'singles/'
		
		#paths to subdirectories containing data from each band
		
		k_path = 'K_data_'
		h_path = 'H_data_'
		yj_path = 'YJ_data_'
		
		
		
		
		#list of paths to files initialised

		filepaths=[]
		
		#list of bands created for looping and constructing filenames
			
		bands = ['K','H','YJ']
		
		#list of filepaths to find appropriate path for given band
		
		bandfiles=[k_path, h_path, yj_path]
		
		#loop to find correct path to each cube
		
		for j in range(len(bands)):
			
			#if statements to determine if chosen band is K
			
			if j==0:
			
			#iteration skipped if not K
			
				if band == 'H' or band=='YJ':
					
					continue
				
				#loop to fill filepaths list with 4 entries for K band
				
				for i in range(4):
					
					#appropriate path to file constructed, and filename found depending on iteration of loop
					
					filename = root_path + bandfiles[j]+str(i+1)+'/NGC346_'+bands[j]+'_2c_COMBINED_CUBE_'+star_name+'.fits'
					
					#list appended with path to file
					
					filepaths.append(filename)
					
			#if statements to determine if chosen band is H
				
			elif j==1:
			
			#iteration skipped if not H
			
				if band=='K' or band == 'YJ':
					
					continue
				
				#loop to fill filepaths list with 2 entries for H band
				
				for i in range(2):
					
					if i == 0:
					
					#appropriate path to file constructed, and filename found depending on iteration of loop
					
						filename=root_path + bandfiles[j] + str(i+1)+'/NGC346_'+bands[j]+'_2c_COMBINED_CUBE_'+star_name+'.fits'
					#list appended with path to file	
						filepaths.append(filename)
					
					else:
						
					#same as before, taking the 1C files into account
						filename=root_path + bandfiles[j] + str(i+1)+'/NGC346_'+bands[j]+'_1c_COMBINED_CUBE_'+star_name+'.fits'
						filepaths.append(filename)
						
			#if this point is reached, band has to be YJ
					
			else:
			
			#iteration skippe if not YJ
			
				if band == 'K' or band == 'H':
					
					continue
			#loop to fill filepaths list with 3 entries for YJ band
			
				for i in range(3):
					
					if i==0:
						
						#appropriate path to file constructed, and filename found depending on iteration of loop
						
						filename=root_path + bandfiles[j] + str(i+1)+'/NGC346_'+bands[j]+'_2c_COMBINED_CUBE_'+star_name+'.fits'
						
						#list appended with path to file	
						
						filepaths.append(filename)
						
					else:
						
						#same as before, taking the 1C files into account
						
						filename=root_path + bandfiles[j] + str(i+1)+'/NGC346_'+bands[j]+'_1c_COMBINED_CUBE_'+star_name+'.fits'
						
						filepaths.append(filename)
						
		#list of paths to single files returned		
					
		return filepaths

	#function to return user input for chosen object and band			

		
		
		
		#starname and band returned
		

	def collapse_into_image(self,file_name):
		
		#cube read by SpectralCube package
		
		cube = SpectralCube.read(file_name, hdu=1)

		#dimensions of cube set

		wave, dec, ra, = cube.world[:]

		# Data order in cube is (n_spectral, n_y, n_x)
		image = cube[0, :, :]

		# To get the units of the spectral axis
		cube.spectral_axis

		# Make a continuum image
		cont_img = cube.median(axis=0)

		# Find its shape
		cont_img.shape
		
		# Chop out the NaN borders
		cmin = cube.minimal_subcube()
		cube = cmin

		# Make a continuum image
		collapsed_image = cmin.median(axis=0)

		# Find its shape

		# Plot the cube with the NaN borders removed 
		#fig = plt.figure()
		#plt.imshow(cont_img.value)
		#plt.tight_layout()


		#plt.show()

		#image and datacube returned for plotting or further processing if required
		
		return [collapsed_image,cube]

	#function to create a spectrum from image and datacube. User input required in case of multiple stars to choose which to use

	def construct_spectrum(self,cont_img,cube):

		#image plotted to view object(s) for selection
		
		fig = plt.figure()

		plt.imshow(cont_img.value)

		plt.tight_layout()
		
		plt.show()
		
		#statistics defined
		
		mean, median, std = sigma_clipped_stats(cont_img.value, sigma=3.0)

		# Get a list of sources. Threshold can be very low if required since objects are manually chosen from image
		daofind = DAOStarFinder(fwhm=3.0, threshold=1*std)  
		sources = daofind(cont_img.value - median) 
		print("Number of sources in field: ", len(sources))
		print()
		print(sources)
		
		
		#star positions mapped
		
		positions = Table([sources['xcentroid'], sources['ycentroid']])                # Table of positions in pixels
		w = WCS(cont_img.header)                                                       # Instantiate WCS object
		radec_lst = w.pixel_to_world(positions['xcentroid'], positions['ycentroid'])   # Convert to RA and Dec in ICRS

		#print(radec_lst[0].to_string('hmsdms'))


		#extracting spectra from the source, circular aperture extraction



		# Size of frame 
		ysize_pix = cube.shape[1]
		xsize_pix = cube.shape[2]

		#sources and x positions printed. User chooses which object to select for spectroscopy

		if len(sources['xcentroid'])>1:
			print('Which source do you want the spectrum of?')
			
			x = input()
			choice=int(x)
			index = choice-1
			ycent_pix=sources['ycentroid'][choice-1]
			xcent_pix=sources['xcentroid'][choice-1]
			
			
		#if only one object is present, no input is needed	
			
		else:
			ycent_pix = sources['ycentroid'][0]
			xcent_pix = sources['xcentroid'][0]
			
		#optional input to save co-ordinates to a file	

		#print('Write co-ordinates to file? y/n')

		#x = input()

		#if x == 'y':

			#f=open('sub__star_pos.txt','a')
			#f.write(star_name+band+'	x = ' + str(xcent_pix) + '	' + 'y = ' + str(ycent_pix) +'\n')
			#f.close()

		#continue

		# Make an aperture radius for source
		apertureRad_pix = 2

		# Make a masked array for the apeture
		yy, xx = np.indices([ysize_pix,xsize_pix], dtype='float')    # Check ycentpix, xcentpix are in correct order 
		radius = ((yy-ycent_pix)**2 + (xx-xcent_pix)**2)**0.5        # Make a circle in the frame

		mask = radius <= apertureRad_pix                # Select pixels within the aperture radius


		maskedcube = cube.with_mask(mask)               # Make a masked cube

		spectrum = maskedcube.sum(axis=(1,2))           # Extract the spectrum from only the anulus - use sum
		spectrum_max = maskedcube.max(axis=(1,2))           # Extract the spectrum from only the anulus - use max
		noisespectrum = maskedcube.std(axis=(1, 2))     # Extract the noise spectrum for the source 

		#fig = plt.figure(figsize=(14,6))

		#plt.plot(maskedcube.spectral_axis.value,spectrum.value, label='Source sum')                 # source spectrum 
		#plt.plot(maskedcube.spectral_axis.value,spectrum_max.value, label='Source max')             # source spectrum 


		#plt.xlabel('Wavelength (microns)')
		#plt.ylabel('erg / (Angstrom cm2 s)') #plt.ylabel(spectrum.unit)

		#plt.xlim(min(cube.spectral_axis).value, 1.83)
		#plt.ylim(-0.1e-15, 0.2e-15)

		#plt.legend(frameon=False, fontsize='medium')

		#plt.tight_layout()
		#plt.show()

		#background subtraction

		# Find out how many pixels are in the aperture - short version

		pixInAp = np.count_nonzero(mask == 1) # pix in apeture
		print(pixInAp)

		pixNotAp = np.count_nonzero(mask == 0) # pix not in apeture
		print(pixNotAp)

		# Try measuring a spectrum from the background -> Use an anulus arround the source

		an_mask = (radius > apertureRad_pix + 1) & (radius <= apertureRad_pix + 2)   # Select pixels within an anulus
		an_maskedcube = cube.with_mask(an_mask)                                      # Make a masked cube

		bkg_spectrum = an_maskedcube.median(axis=(1,2))

		# Try measuring a spectrum from the background -> Use a box awway from source 

		#bkgcube = cube[:, 1:4, 10:13]  
		#bkgbox_spectrum = bkgcube.median(axis=(1,2))  

		# Background corrected spectrum - box
		#correctedSp_box = spectrum - (bkgbox_spectrum * pixInAp)
		
		#background corrected spectrum-anulus
		
		correctedSp = spectrum - (bkg_spectrum * pixInAp)

		# Plot the spectrum extracted from cirular aperture 

		#fig = plt.figure(figsize=(14,6))

		#plt.plot(maskedcube.spectral_axis.value,spectrum.value, label='Source')             # source spectrum 
		#plt.plot(maskedcube.spectral_axis.value,correctedSp.value, label='Bkg corr Source') # background corrected source spectrum 
		#plt.plot(maskedcube.spectral_axis.value,noisespectrum.value, label='Noise')         # noise spectrum 
		#plt.plot(an_maskedcube.spectral_axis.value,bkg_spectrum.value, label='Background')  # background spectrum 

		#plt.xlabel('Wavelength (microns)')
		#plt.ylabel('erg / (Angstrom cm2 s)') #plt.ylabel(spectrum.unit)

		#plt.xlim(min(cube.spectral_axis).value, 1.83)
		#plt.ylim(-0.1e-15, 0.2e-15)

		#plt.legend(frameon=False, fontsize='medium')

		#plt.tight_layout()
		#plt.show()

		#converting to correct units

		# Convert flux from erg / (Angstrom cm2 s) to Jy 
		correctedSp_Jy = correctedSp.to(u.Jansky, equivalencies=u.spectral_density(maskedcube.spectral_axis))
		spectrum_Jy = spectrum.to(u.Jansky, equivalencies=u.spectral_density(maskedcube.spectral_axis))
		noiseSp_Jy   = noisespectrum.to(u.Jansky, equivalencies=u.spectral_density(maskedcube.spectral_axis))

		# Plot the spectrum extracted from cirular aperture in Jy

		#fig = plt.figure(figsize=(14,6))
		#plt.plot(maskedcube.spectral_axis.value,correctedSp_Jy.value, label='Bkg corr Source') # background corrected source spectrum 

		#plt.xlabel('Wavelength (microns)')
		#plt.ylabel(correctedSp_Jy.unit)

		#xlim for h band

		##plt.xlim(min(cube.spectral_axis).value, (max(cube.spectral_axis).value))
		#plt.ylim(-0.1e-3, 0.15e-2)

		#xlim for yj band

		#plt.xlim(min(cube.spectral_axis).value, max(cube.spectral_axis).value)

		#hlim for y532

		#plt.ylim(-0.1e-3,0.5e-1)

		#plt.legend(frameon=False, fontsize='medium')

		#plt.tight_layout()

		#plt.savefig('Spectra/'+star_name + identifier + '_bkg corrected spectrum_' + band +'_band.png') 
		#plt.show()
		
		#masked cube and background corrected spectrum in Janskys returned
		
		return[maskedcube,correctedSp_Jy,spectrum_Jy]
		
	def gaussian_fit(self,spectrum,submin,submax,star_name):	
		

		noise_region=SpectralRegion(2.26*u.um,2.3*u.um)
		spectrum=noise_region_uncertainty(spectrum,noise_region)
		#lines = find_lines_threshold(spectrum, noise_factor=4)
		
		#print(lines[lines['line_type'] == 'emission'])
		#print(lines[lines['line_type'] == 'absorption'])
		
		#amp=0.00075
		#centre=2.1675
		
		#sub region defined from function parameters
		
		sub_region=SpectralRegion(submin*u.um,submax*u.um)
		
		#spectrum extracted from subregion
		
		sub_spectrum=extract_region(spectrum,sub_region)
		
		#continuum fitted
		
		g1_fit=fit_generic_continuum(spectrum)
		
		y_continuum_fitted=g1_fit(sub_spectrum.spectral_axis)
		
		#continuum fit plotted with subregion spectrum
		
		plt.plot(sub_spectrum.spectral_axis,sub_spectrum.flux)
		plt.plot(sub_spectrum.spectral_axis,y_continuum_fitted)
		plt.title('Continuum Fitting')
		plt.grid(True)
		plt.show()
		
		#continuum substracted to show intensity
		
		sub_spectrum=sub_spectrum-y_continuum_fitted
		
		#initial gaussian fitted to estimate parameters for more accurate fitting
		
		estimates=estimate_line_parameters(sub_spectrum,models.Gaussian1D())
		
		#new gauss_axis variable created for a smooth fit
		
		gauss_axis=np.linspace(min(sub_spectrum.spectral_axis),max(sub_spectrum.spectral_axis),2000)
		
		#gaussian fit from previous estimates produced
		
		g_init=estimates
		
		g_fit=fit_lines(sub_spectrum,g_init)
		
		y_fit=g_fit(gauss_axis)
		
		#fit plotted
		
		plt.plot(sub_spectrum.spectral_axis,sub_spectrum.flux)
		plt.plot(gauss_axis,y_fit)
		
		#linestrength found
		
		strength=(max(y_fit))
		
		#f=open('brackett_strengths.txt','a')
		#f.write(star_name + ' - '+ str(strength) + 'Jy') 
		#f.close()
		
		plt.title('Single fit peak')
		plt.grid(True)
		plt.legend('Original Spectrum','Specutils Fit Result')
		
		plt.show()

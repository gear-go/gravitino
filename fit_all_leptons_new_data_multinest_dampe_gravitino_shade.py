#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Monday July 10 11:45 2018

@author: germangomez
"""



###################################################################################################
"""
Fit to AMS02 data and DAMPE electron + positron data using gravitino dark matter
"""



####################################################################################################

import time
tic = time.clock()



import json
#import sys
#from numpy import log, exp, pi
#import scipy.stats, scipy
import pymultinest
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d


import scipy.stats

class Shade(object):
	"""
	The Shade class allows plotting of model predictions as calculated from a chain.
	
	:param: x the independent variable
	
	call add to add predictions from the chain.
	
	Example:
	
	x = numpy.linspace(0, 1, 100)
	shade = Shade(x)
	for c in chain:
		shade.add(c[0] * x + c[1])
	shade.line()
	shade.shade()
	plt.show()
	
	"""
	def __init__(self, x, shadeargs={}, lineargs={}):
		self.x = x
		self.ys = []
		self.shadeargs = shadeargs
		self.lineargs = lineargs
	def add(self, y):
		""" add a possible prediction """
		self.ys.append(y)
	def set_shadeargs(self, **kwargs):
		self.shadeargs = kwargs
	def set_lineargs(self, **kwargs):
		self.lineargs = kwargs
	def get_line(self, q=0.5):
		assert len(self.ys) > 0, self.ys
		return scipy.stats.mstats.mquantiles(self.ys, q, axis=0)[0]
	def shade(self, q=0.341, **kwargs):
		""" Use the stored predictions to plot a shaded region
		using 0.5-q and 0.5+q as limits. """
		shadeargs = dict(self.shadeargs)
		shadeargs.update(kwargs)
		lo = self.get_line(0.5 - q)
		hi = self.get_line(0.5 + q)
		return plt.fill_between(self.x, lo, hi, **shadeargs)
	def line(self, **kwargs):
		""" Use the stored predictions to plot the median """
		lineargs = dict(self.lineargs)
		lineargs.update(kwargs)
		mid = self.get_line(0.5)
		return plt.plot(self.x, mid, **lineargs)

__doc__ = Shade.__doc__
__all__ = ['Shade']

#############################################3
folder_flux_propt = '/home/ctachile/Documents/Gravitino_PUC_2018/AMS/PropagatedFlux/'
folder_flux_input = '/home/ctachile/Documents/Gravitino_PUC_2018/AMS/'
folder_output_ptn = '/home/ctachile/Documents/Gravitino_PUC_2018/AMS/'

mass = np.linspace(1000,4000, 7)

spec01 = {}
spec02 = {}
spec03 = {}
spec04 = {}
spec05 = {}
spec06 = {}
spec07 = {}
spec08 = {}
spec09 = {}

spec10 = {}
spec11 = {}
spec12 = {}
spec13 = {}
spec14 = {}
spec15 = {}
spec16 = {}
spec17 = {}
spec18 = {}

spec19 = {}
spec20 = {}
spec21 = {}
spec22 = {}
spec23 = {}
spec24 = {}
spec25 = {}
spec26 = {}
spec27 = {}

#read data from the latest files produced by mathematica program in ElePosGamma-M0 folder
def readdata(file_name, mass_index):
    file_cols =  np.genfromtxt(file_name,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
                                                  19,20,21,22,23,24,25,26,27),
                               dtype=[('Ener', float), 
                                      ('Fc01', float), ('Fc02', float), ('Fc03', float), 
                                      ('Fc04', float), ('Fc05', float), ('Fc06', float), 
                                      ('Fc07', float), ('Fc08', float), ('Fc09', float),
                                      ('Fc10', float), ('Fc11', float), ('Fc12', float),
                                      ('Fc13', float), ('Fc14', float), ('Fc15', float),
                                      ('Fc16', float), ('Fc17', float), ('Fc18', float),
                                      ('Fc19', float), ('Fc20', float), ('Fc21', float),
                                      ('Fc22', float), ('Fc23', float), ('Fc24', float),
                                      ('Fc25', float), ('Fc26', float), ('Fc27', float),
                                      ])
    #plt.plot(name['Ener'], name['Flux'],alpha=0.7, label=labelt)
   
    #9-final state channels, wich are phenomenologically distinguishable
    #Channel 1: AEMuNue    : antielectron-muon-neutrino
    #Channel 2: AETauNue   : antielectron-tau-neutrino
    #Channel 3: EAMuNue    : electron-antimuon-neutrino
    #Channel 4: MuAMuNue   : muon-antimuon-neutrino
    #Channel 5: TauATauNue : tau-antitau-neutrino
    #Channel 6: TauAMuNue  : tau-antimuon-neutrino
    #Channel 7: ATauMuNue  : antitau-muon-neutrino
    #Channel 8: EAENue     : electron-antielectron-neutrino
    #Channel 9: EATauNue   : electron-antitau-neutrino
    
    #three final state prompt particles (no inverse compton effects considered)
    #01-09 prompt electrons
    #10-18 prompt positrons
    #19-27 prompt photons
    
    #9 final channels, electron final state (see trilinear_channels.txt)
    spec01[mass_index] = ([ file_cols['Ener'], file_cols['Fc01'] ])
    spec02[mass_index] = ([ file_cols['Ener'], file_cols['Fc02'] ])
    spec03[mass_index] = ([ file_cols['Ener'], file_cols['Fc03'] ])
    spec04[mass_index] = ([ file_cols['Ener'], file_cols['Fc04'] ])
    spec05[mass_index] = ([ file_cols['Ener'], file_cols['Fc05'] ])
    spec06[mass_index] = ([ file_cols['Ener'], file_cols['Fc06'] ])
    spec07[mass_index] = ([ file_cols['Ener'], file_cols['Fc07'] ])
    spec08[mass_index] = ([ file_cols['Ener'], file_cols['Fc08'] ])
    spec09[mass_index] = ([ file_cols['Ener'], file_cols['Fc09'] ])
    
    #9 decay channels, positron final state
    spec10[mass_index] = ([ file_cols['Ener'], file_cols['Fc10'] ])
    spec11[mass_index] = ([ file_cols['Ener'], file_cols['Fc11'] ])
    spec12[mass_index] = ([ file_cols['Ener'], file_cols['Fc12'] ])
    spec13[mass_index] = ([ file_cols['Ener'], file_cols['Fc13'] ])
    spec14[mass_index] = ([ file_cols['Ener'], file_cols['Fc14'] ])
    spec15[mass_index] = ([ file_cols['Ener'], file_cols['Fc15'] ])
    spec16[mass_index] = ([ file_cols['Ener'], file_cols['Fc16'] ])
    spec17[mass_index] = ([ file_cols['Ener'], file_cols['Fc17'] ])
    spec18[mass_index] = ([ file_cols['Ener'], file_cols['Fc18'] ])
    
    #9 decay channels, gamma final state    
    spec19[mass_index] = ([ file_cols['Ener'], file_cols['Fc19'] ])
    spec20[mass_index] = ([ file_cols['Ener'], file_cols['Fc20'] ])
    spec21[mass_index] = ([ file_cols['Ener'], file_cols['Fc21'] ])
    spec22[mass_index] = ([ file_cols['Ener'], file_cols['Fc22'] ])
    spec23[mass_index] = ([ file_cols['Ener'], file_cols['Fc23'] ])
    spec24[mass_index] = ([ file_cols['Ener'], file_cols['Fc24'] ])
    spec25[mass_index] = ([ file_cols['Ener'], file_cols['Fc25'] ])
    spec26[mass_index] = ([ file_cols['Ener'], file_cols['Fc26'] ])
    spec27[mass_index] = ([ file_cols['Ener'], file_cols['Fc27'] ])
    
    return file_cols

def interpolateSpec(Mdm, spec, final_state,x):
    
    #looking for the mass indices around Mdm for mass interpolation
    for i in range(len(spec)-1):
        m1=mass[i]  
        m2=mass[i+1] 
        
        if m1 <= Mdm < m2:
            break
        
    #exponential factor for logarithmic interpolation?    
    alpha = np.log(float(Mdm)/float(m2)) / np.log(float(m1)/float(m2))

    #reading the flux and energy corresponding to each mass index around Mdm value
    dE1 = np.array(spec[m1][0])#[:,0]
    dN1 = np.array(spec[m1][1])#[:,1]
    
    dE2 = np.array(spec[m2][0])#[:,0]
    dN2 = np.array(spec[m2][1])#[:,1]
    
    #manual energy-flux mass interpolattion array 
    (xp,yp)=(dE1**alpha*dE2**(1.-alpha),dN1**alpha*dN2**(1.-alpha))
   
    #automatic energy-flux interpolation array
    dNdE = interp1d(xp,yp,bounds_error=False,fill_value=0.)
      
    #electron final state, (hopefully) bins from electron AMS-02
    return dNdE(x) 
    
   
    
#This function that returns an array with the flux of different final state particles
#Here we consider each branching fraction as independent, but in fact this can be reduced
#by considering some coincidences and symmetries between the spectra (trilinear_channels.txt)        

#Thus, in practice we should use EffTotGravFluxLep for the fit of AMS and
#EffTotGravFluxPho for the computation of the photon spectrum
       

#This function returns the effective lepton flux [1/GeV m2 s] considering that 
#there are channels with equal spectrum and
#there are BRs which are dependent because of CC (for instance AEMuNu = EMuNu)
#see trilinear_channels.txt for more information about the derivation of the effective flux    
    
def EffTotGravFluxLep(Mdm, fluxNorm, a1, a2, x):
    
    # interpolateSpec(Mdm,spec01) returns an array binned as given by some input data
    # currently we are using the binning of flux_positrons_AMS02.dat but maybe we could
    # systematize this option in interpolateSpec(Mdm,spec01) to consider ele, pos or gamma data
    
    if ( a1 + a2 <= 1.0):
        a3 = 1.0 - (a1 + a2)    
    else:
        a1 = a1/(a1 + a2)
        a2 = a2/(a1 + a2)
        a3 = 0
    
    #since spec01-ele = spec04-ele = spec07-ele we can use only spec01-ele and
    #a1 = Br1 + Br4 + Br7
    
    #and similarly for 2-5-6 and 3-8-9
    #a2 = Br2 + Br5 + Br6
    #a3 = Br3 + Br8 + Br9
    
    #Because of CC symmetry applied to 
    #channel topology - final state lepton (ele-spec01 == pos-spec03, etc) and 
    #Brs (Br1 = Br3, Br2 = Br9 and Br6 = Br7) we get that
    #the positron flux is equal to the electron one, neat.
    
    s1 = interpolateSpec(Mdm,spec01,0,x)
    s2 = interpolateSpec(Mdm,spec02,0,x)
    s3 = interpolateSpec(Mdm,spec03,0,x)
       
    #total gravitino flux considering the lifetime factor and pheno branching ratios
    #the units depend on the final state particle (ep 1/m2, gamma 1/cm2), see above
    #pheno branching ratios should be linked or derived from values of lambdas (to do task)
    dNTdE = (fluxNorm)*(a1*s1 + a2*s2 + a3*s3)#*x*x
    return dNTdE



#This function returns the effective photon spectrum [1/GeV cm2] derived 
#from the parameters a1 and a2 of the fit to AMS
    
def EffTotGravFluxPho(Mdm, fluxNorm, a1, a2,x):
    
    # interpolateSpec(Mdm,spec01) returns an array binned as given by some input data
    # currently we are using the binning of flux_positrons_AMS02.dat but maybe we could
    # systematize this option in interpolateSpec(Mdm,spec01) to consider ele, pos or gamma data
    
    if ( a1 + a2 <= 1.0):
        a3 = 1.0 - (a1 + a2)
    else:
        a1 = a1/(a1 + a2)
        a2 = a2/(a1 + a2)
        a3 = 0
    
    s1 = interpolateSpec(Mdm,spec19,2,x)
    s2 = interpolateSpec(Mdm,spec20,2,x)
    s3 = interpolateSpec(Mdm,spec21,2,x)
    s4 = interpolateSpec(Mdm,spec22,2,x)
    s5 = interpolateSpec(Mdm,spec23,2,x)
    s6 = interpolateSpec(Mdm,spec24,2,x)
    s7 = interpolateSpec(Mdm,spec25,2,x)
    s8 = interpolateSpec(Mdm,spec26,2,x)
    s9 = interpolateSpec(Mdm,spec27,2,x)
    
    #From the fit to AMS (DAMPE) we just can fix a1, a2 and a3
    #then we have some freedom to choose the photon branching fractions
    #Thus, we randomly generate three sets of three numbers each based on a1, a2 and a3    
    #Such that we respect the constrains
    
    #a1 = Br1 + Br4 + Br7
    #a2 = Br2 + Br5 + Br6
    #a3 = Br3 + Br8 + Br9
    
    #c1 = np.random.uniform(0,1,3)
    #c1 /= c1.sum()
    #c1 = a1*c1
    
    #c2 = np.random.uniform(0,1,3)
    #c2 /= c2.sum()
    #c2 = a2*c2
    
    #c3 = np.random.uniform(0,1,3)
    #c3 /= c3.sum()
    #c3 = a3*c3

    c1 = [a1,0,0]
    c2 = [a2,0,0]
    c3 = [0,a3,0]    

    #total gravitino flux considering the lifetime factor and pheno branching ratios
    #the units depend on the final state particle (ep 1/m2, gamma 1/cm2), see above
    #pheno branching ratios should be linked or derived from values of lambdas (to do task)
    
    dNTdE = (fluxNorm)*(c1[0]*s1 + c2[0]*s2 + c3[0]*s3 + 
                           c1[1]*s4 + c2[1]*s5 + c2[2]*s6 + 
                           c1[2]*s7 + c3[1]*s8 + c3[2]*s9)*x*x
    return dNTdE


#############################################

nthreads = 2

ep_data = np.genfromtxt('electron_positron_fraction_AMS02_usine.dat',usecols=(3,4,5,6,7,8,9,10),
                     dtype=[('Ener', float),
                            ('EBIN_LOW', float),
                            ('EBIN_HIGH', float),
                            ('VALUE', float),
                            ('ERR_STATme', float),
                            ('ERR_STATma', float),
                            ('ERR_SYSTme', float),
                            ('ERR_SYSTma', float),]) # e+/(e-+e+) from Accardo et al., PRL 113, 121101 (2014)

e_data = np.genfromtxt('flux_electron_AMS02.dat',usecols=(3,4,5,6,7,8,9,10),
                     dtype=[('Ener', float),
                            ('EBIN_LOW', float),
                            ('EBIN_HIGH', float),
                            ('VALUE', float),
                            ('ERR_STATme', float),
                            ('ERR_STATma', float),
                            ('ERR_SYSTme', float),
                            ('ERR_SYSTma', float),]) # e- from Accardo et al., P2014PhRvL.113l1102A 

p_data = np.genfromtxt('flux_positrons_AMS02.dat',usecols=(3,4,5,6,7,8,9,10),
                     dtype=[('Ener', float),
                            ('EBIN_LOW', float),
                            ('EBIN_HIGH', float),
                            ('VALUE', float),
                            ('ERR_STATme', float),
                            ('ERR_STATma', float),
                            ('ERR_SYSTme', float),
                            ('ERR_SYSTma', float),]) # e- from Accardo et al., P2014PhRvL.113l1102A 


epp_dampe_data = np.genfromtxt('electron_positron_sum_DAMPE_China.dat',usecols=(2,0,1,9,10),
                     dtype=[('Ener', float),
                            ('EBIN_LOW', float),
                            ('EBIN_HIGH', float),
                            ('VALUE', float),
                            ('ERR_STATma', float),]) # (e-+e+) from DAMPE

epp_AMS02_data = np.genfromtxt('flux_electronANDpositrons_AMS02.dat',usecols=(3,4,5,6,7,8,9,10),
                     dtype=[('Ener', float),
                            ('EBIN_LOW', float),
                            ('EBIN_HIGH', float),
                            ('VALUE', float),
                            ('ERR_STATme', float),
                            ('ERR_STATma', float),
                            ('ERR_SYSTme', float),
                            ('ERR_SYSTma', float),]) # e- from Accardo et al., P2014PhRvL.113l1102A 


epp_CALET_data = np.genfromtxt('electron_positron_sum_CALET.dat',usecols=(2,0,1,3,4),
                     dtype=[('Ener', float),
                            ('EBIN_LOW', float),
                            ('EBIN_HIGH', float),
                            ('VALUE', float),
                            ('ERR_STATma', float),]) # e- from Accardo et al., P2014PhRvL.113l1102A 


readdata(folder_flux_propt + 'PropFlux_MED_dnde_9channels_epg_FinalStates_mass_1000_GeV.dat', 1000.0)
readdata(folder_flux_propt + 'PropFlux_MED_dnde_9channels_epg_FinalStates_mass_1500_GeV.dat', 1500.0)
readdata(folder_flux_propt + 'PropFlux_MED_dnde_9channels_epg_FinalStates_mass_2000_GeV.dat', 2000.0)
readdata(folder_flux_propt + 'PropFlux_MED_dnde_9channels_epg_FinalStates_mass_2500_GeV.dat', 2500.0)
readdata(folder_flux_propt + 'PropFlux_MED_dnde_9channels_epg_FinalStates_mass_3000_GeV.dat', 3000.0)
readdata(folder_flux_propt + 'PropFlux_MED_dnde_9channels_epg_FinalStates_mass_3500_GeV.dat', 3500.0)
readdata(folder_flux_propt + 'PropFlux_MED_dnde_9channels_epg_FinalStates_mass_4000_GeV.dat', 4000.0)

#def ep(CA,CC,gammaA,gammaC,Esi,x):
#    model = CA*(x**gammaA) * ((1+ (CC/CA)*(x**(gammaC-gammaA)) * np.exp(-x*Esi)) / (1 + CA*(x**gammaA) + (2*CC*(x**gammaC)*np.exp(-x*Esi)) ))
#    return model

#def p(Ce, gammae, Cs, gammas, Esi, x):
#    model = Ce*(x**(-gammae)) + Cs*(x**(-gammas))*np.exp(-x*Esi)
#    return model

x = p_data['Ener']
xpos = p_data['Ener']
xele = e_data['Ener']
xpho = xpos


# a more elaborate prior
# parameters are pos1, width, height1, [height2]
def prior(cube, ndim, nparams):
        cube[0] = cube[0]*40            # Cp uniform prior between 0:40
        cube[1] = cube[1]*5            # gammap uniform prior between 0:5
        cube[2] = cube[2]*300+300            # Ce uniform prior between 0:400
        cube[3] = cube[3]*5            # gammae uniform prior between 0:5
        cube[4] = (cube[4]*3000 + 1000 )         # Mdm uniform prior between 1000:4000
        cube[5] = (cube[5]*(3)-3)           # fluxNorm log-uniform prior between 0.001:1
        cube[6] = cube[6]*(1)  # a1 log-uniform prior between 10^-4 and 1
        cube[7] = cube[7]*(1)  # a2 log-uniform prior between 10^-4 and 1
        cube[8] = (cube[8]*(4) - 8) # g log-uniform prior between 10^-4 and 1
        cube[9] = (cube[9]*(4) - 8) # h log-uniform prior between 10^-4 and 1
        cube[10] = (cube[10]*(4) - 8) # i log-uniform prior between 10^-4 and 1
        cube[11] = (cube[11]*(4) - 4) # j log-uniform prior between 10^-4 and 1


def fraction(Ce,gammae,Cp,gammap,Mdm, fluxNorm, a1, a2):
    def spectrum(x):
        return (Cp*(x**(-gammap)) + EffTotGravFluxLep(Mdm, fluxNorm, a1, a2, x) )/ ( Cp*(x**(-gammap)) + 2*EffTotGravFluxLep(Mdm, fluxNorm, a1, a2, x) + Ce*(x**(-gammae))  )
        #return Cp*(x**(-gammap)) + Cs*(x**(-gammas))*np.exp(-x*Esi) / ( Cp*(x**(-gammap)) + Cs*(x**(-gammas))*np.exp(-x*Esi) + Ce*(x**(-gammae)) + Cs*(x**(-gammas))*np.exp(-x*Esi) )
#        return positrons(Ce,gammae,Cp, gammap, Cs, gammas, Esi) / ( positrons(Ce,gammae,Cp, gammap, Cs, gammas, Esi) + electrons(Ce, gammae,Cp, gammap, Cs, gammas, Esi) )
    return spectrum(ep_data['Ener'])

def suma(Ce,gammae,Cp,gammap,Mdm, fluxNorm, a1, a2):
    def spectrum(x):
        return  ( Cp*(x**(-gammap)) + Ce*(x**(-gammae)) + 2*EffTotGravFluxLep(Mdm, fluxNorm, a1, a2, x) )
#        return positrons(Ce,gammae,Cp, gammap, Cs, gammas, Esi) / ( positrons(Ce,gammae,Cp, gammap, Cs, gammas, Esi) + electrons(Ce, gammae,Cp, gammap, Cs, gammas, Esi) )
    return spectrum(epp_dampe_data['Ener'])

def positrons(Cp,gammap,Mdm, fluxNorm, a1, a2):
    def spectrum(x):
        return  Cp*(x**(-gammap)) + EffTotGravFluxLep(Mdm, fluxNorm, a1, a2, x)#Cs*(x**(-gammas))*np.exp(-x*Esi)
    return spectrum(p_data['Ener'])


def electrons(Ce,gammae,Mdm, fluxNorm, a1, a2):
    def spectrum(x):
        return  Ce*(x**(-gammae)) + EffTotGravFluxLep(Mdm, fluxNorm, a1, a2, x)#
    return spectrum(e_data['Ener'])



def sour(Cs, gammas, Esi, x):# Luego es remplazado por el gravitino
    model2 =  Cs*(x**(-gammas))*np.exp(-x*Esi)
    return model2

def back(Cp, gammap, x):
    model = Cp * (x**(-gammap)) 
    return model

#def loglikee(cube, ndim, nparams):
#        Ce, gammae, Cs, gammas, Esi, h  = cube[2], cube[3], cube[4],cube[5], cube[6], cube[8]
#        ymodel = electrons(Ce,gammae,Cs,gammas, Esi)
#        yerr = e_data['ERR_STATma']
#        y = e_data['VALUE']
#        inv_sigma2 = 1.0/(yerr**2 + ymodel**2*np.exp(2*np.log(h)))
#        return -0.5*(np.sum((y-ymodel)**2*inv_sigma2 - np.log(inv_sigma2)))
#
#
#def loglikep(cube, ndim, nparams):
#        Cp, gammap, Cs, gammas, Esi, g  = cube[0], cube[1], cube[4],cube[5], cube[6], cube[7]
#        ymodel = positrons(Cp,gammap,Cs,gammas, Esi)
#        yerr = p_data['ERR_STATma']
#        y = p_data['VALUE']
#        inv_sigma2 = 1.0/(yerr**2 + ymodel**2*np.exp(2*np.log(g)))
#        return -0.5*(np.sum((y-ymodel)**2*inv_sigma2 - np.log(inv_sigma2)))

def loglike(cube,ndim,nparams):
        Cp, gammap,Ce, gammae, Mdm, logfluxNorm, a1, a2, logg, logh, logi, logj =cube[0], cube[1], cube[2], cube[3], cube[4],cube[5], cube[6], cube[7], cube[8], cube[9], cube[10], cube[11]
        fluxNorm = 10**logfluxNorm
        #print(Mdm,fluxNorm, logfluxNorm, Ce)
        g = 10**logg
        h = 10**logh
        i = 10**logi
        j = 10**logj
        ymodelp = positrons(Cp,gammap,Mdm, fluxNorm, a1, a2)
        yerrp = p_data['ERR_STATma']
        yp = p_data['VALUE']
        inv_sigma2p = 1.0/(yerrp**2 + ymodelp**2*np.exp(2)*g)
        loglikep = -0.5*(np.sum((yp-ymodelp)**2*inv_sigma2p - np.log(inv_sigma2p)))
        ymodele = electrons(Ce,gammae,Mdm, fluxNorm, a1, a2)
        yerre = e_data['ERR_STATma']
        ye = e_data['VALUE']
        inv_sigma2e = 1.0/(yerre**2 + ymodele**2*np.exp(2)*h)
        loglikee = -0.5*(np.sum((ye-ymodele)**2*inv_sigma2e - np.log(inv_sigma2e)))
        ymodelep = fraction(Ce,gammae,Cp,gammap,Mdm, fluxNorm, a1, a2)
        yerrep = ep_data['ERR_STATma']
        yep = ep_data['VALUE']
        inv_sigma2ep = 1.0/(yerrep**2 + ymodelep**2*np.exp(2)*i)
        loglikeep = -0.5*(np.sum((yep-ymodelep)**2*inv_sigma2ep - np.log(inv_sigma2ep)))
        ymodeleppc = suma(Cp,gammap,Ce,gammae,Mdm, fluxNorm, a1, a2)
        yerreppc = epp_dampe_data['ERR_STATma']
        yeppc = epp_dampe_data['VALUE']
        inv_sigma2eppc = 1.0/(yerreppc**2 + ymodeleppc**2*np.exp(2)*j)
        loglikeeppc = -0.5*(np.sum((yeppc-ymodeleppc)**2*inv_sigma2eppc - np.log(inv_sigma2eppc)))
        return loglikep + loglikee + loglikeep + loglikeeppc

ypdata = p_data['VALUE']
yedata = e_data['VALUE']
datafile = "fit_multines_C_Gravitino_DAMPE_shade"
parameters = ["Cp", "gammap", "Ce", "gammae", "Mdm", "logfluxNorm", "a1", "a2", "g", "h", "i","j"]
n_params = len(parameters)

# run MultiNest
pymultinest.run(loglike, prior, n_params, outputfiles_basename=datafile + '_4_', resume = False, verbose = True)
json.dump(parameters, open(datafile + '_4_params.json', 'w')) # save parameter names
print('it is safe here')
# plot the distribution of a posteriori possible models


Cp_true =  28.1 #18.0#26.9
gammap_true = 3.666#3.5#3.79 

Mdm_true = 1200#3.0# 4.21
fluxNorm_true = 1.7#2.647#2.707
a1_true = 0.2#2.4e-4#1.5e-4
a2_true = 0.3
#f_true = 0.001
g_true = 2.15e-2#1e-5#1.85e-2

#plt.plot(x, ydata, '+ ', color='red', label='data')

             
a = pymultinest.Analyzer(outputfiles_basename=datafile + '_4_', n_params = n_params)

###############PLOT Positron flux FROM FIT TO all#############





fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

plt.errorbar(p_data['Ener'], (p_data['Ener']**3)*p_data['VALUE'], yerr=(p_data['Ener']**3)*p_data['ERR_STATma'],color='#3288bd',label='positrons', fmt='o',alpha=0.7)

x = p_data['Ener']
shade = Shade(x)
for (Cp, gammap, Ce, gammae, Mdm, logfluxNorm, a1, a2, g, h, i, j) in a.get_equal_weighted_posterior()[::10,:-1]:
	shade.add((p_data['Ener']**3)*positrons(Cp, gammap,  Mdm, 10**logfluxNorm, a1, a2))
shade.line(color='k')
shade.shade(q=0.25, color='#008837', alpha=0.3)
shade.shade(q=0.5, color='#7b3294', alpha=0.3)


#for (Cp, gammap, Ce, gammae, Mdm, logfluxNorm, a1, a2, logg, h, i, j) in a.get_equal_weighted_posterior()[::100,:-1]:
#    plt.plot(p_data['Ener'], (p_data['Ener']**3)*positrons(Cp, gammap,  Mdm, 10**logfluxNorm, a1, a2),color='red',alpha=0.3)
#        plt.plot(p_data['Ener'], (p_data['Ener']**3)*sour(Cs, gammas, Esi,p_data['Ener']),color='blue',ls='--', alpha=0.2)
#        plt.plot(x, model(pos1, width, height1, 0), '-', color='blue', alpha=0.3, label='data')
#plt.plot(p_data['Ener'], (p_data['Ener']**3)*positrons(Cp_true, gammap_true, Cs_true, gammas_true, Esi_true),color='black',alpha=0.7, label='marginal')
#plt.plot(p_data['Ener'], (p_data['Ener']**3)*sour(Cs_true, gammas_true, Esi_true,p_data['Ener']),color='green',ls='--', alpha=0.7)

ax.legend()
ax.legend(fontsize=5)
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylim((5, 30))
ax.set_xlim(5, 1e3)
ax.minorticks_on()
ax.legend(loc='best')
ax.set_xlabel('Energy [GeV]',fontsize=15)
ax.set_ylabel('E^3 Positron flux',fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
ax.legend(loc='best')

plt.savefig(datafile + '_4_positron_flux.pdf')
###############PLOT Electron flux FROM FIT TO ELECTRONS AND POSITRONS#############c1 = np.array(a1,0,0)


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)


x = e_data['Ener']
shade = Shade(x)
for (Cp, gammap, Ce, gammae, Mdm, logfluxNorm, a1, a2, g, h, i, j) in a.get_equal_weighted_posterior()[::10,:-1]:
	shade.add((e_data['Ener']**3)*electrons(Ce, gammae, Mdm, 10**logfluxNorm, a1, a2))
shade.line(color='k')
shade.shade(q=0.25, color='#008837', alpha=0.3)
shade.shade(q=0.5, color='#7b3294', alpha=0.3)


#for (Cp, gammap, Ce, gammae, Mdm, logfluxNorm, a1, a2, g, h, i, j) in a.get_equal_weighted_posterior()[::100,:-1]:
#        plt.plot(e_data['Ener'], (e_data['Ener']**3)*electrons(Ce, gammae, Mdm, 10**logfluxNorm, a1, a2),color='red',alpha=0.3)

plt.errorbar(e_data['Ener'], (e_data['Ener']**3)*e_data['VALUE'], yerr=(e_data['Ener']**3)*e_data['ERR_STATma'],color='#3288bd',label='Data', fmt='o',alpha=0.5)


ax.legend()
ax.legend(loc='best')
ax.legend(fontsize=9)
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylim((50, 250))
ax.set_xlim(5, 1e3)
ax.minorticks_on()
ax.legend(loc='best')
ax.set_xlabel('Energy [GeV]',fontsize=15)
ax.set_ylabel('Electron flux',fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig(datafile + '_4_electron_flux.pdf')


###############PLOT POSITRON FRACTION FROM FIT TO ELECTRONS AND POSITRONS#############
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
#plt.plot(ep_data['Ener'], fraction(Ce_true,gammae_true,Cp_true,gammap_true,Cs_true,gammas_true,Esi_true),color='purple',label='model', alpha=0.7)

x = ep_data['Ener']
shade = Shade(x)
for (Cp, gammap, Ce, gammae, Mdm, logfluxNorm, a1, a2, g, h, i, j) in a.get_equal_weighted_posterior()[::10,:-1]:
	shade.add(fraction(Ce,gammae,Cp,gammap, Mdm, 10**logfluxNorm, a1, a2))
shade.line(color='k')
shade.shade(q=0.25, color='#008837', alpha=0.3)
shade.shade(q=0.5, color='#7b3294', alpha=0.3)



#for (Cp, gammap, Ce, gammae, Mdm, logfluxNorm, a1, a2, g, h, i, j) in a.get_equal_weighted_posterior()[::100,:-1]:
#           plt.plot(ep_data['Ener'], fraction(Ce,gammae,Cp,gammap, Mdm, 10**logfluxNorm, a1, a2),color='red', alpha=0.3)

plt.errorbar(ep_data['Ener'], ep_data['VALUE'], yerr=ep_data['ERR_STATma'],color='#3288bd',label='Data', fmt='o',alpha=0.5)


ax.legend()
ax.legend(loc='best')
ax.legend(fontsize=9)
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylim((0, 2e-1))
ax.set_xlim(5, 1e3)
ax.minorticks_on()
ax.legend(loc='best')
ax.set_xlabel('Energy [GeV]',fontsize=15)
ax.set_ylabel('Positron fraction',fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig(datafile + '_4_ep_fraction.pdf')


###############PLOT POSITRON+ELECTRON FROM FIT TO EVERYTHING#############

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

x = epp_dampe_data['Ener']
shade = Shade(x)
for (Cp, gammap, Ce, gammae, Mdm, logfluxNorm, a1, a2, g, h, i, j) in a.get_equal_weighted_posterior()[::10,:-1]:
	shade.add((epp_dampe_data['Ener']**3)*suma(Cp,gammap,Ce,gammae, Mdm, 10**logfluxNorm, a1, a2))
shade.line(color='k')
shade.shade(q=0.25, color='#008837', alpha=0.3)
shade.shade(q=0.5, color='#7b3294', alpha=0.3)


#for (Cp, gammap, Ce, gammae, Mdm, logfluxNorm, a1, a2, g, h, i, j) in a.get_equal_weighted_posterior()[::100,:-1]:
#           plt.plot(epp_CALET_data['Ener'], (epp_CALET_data['Ener']**3)*suma(Cp,gammap,Ce,gammae, Mdm, 10**logfluxNorm, a1, a2),color='red', alpha=0.3)
#           plt.plot(epp_CALET_data['Ener'], (epp_CALET_data['Ener']**3)*EffTotGravFluxLep(Mdm, 10**logfluxNorm, a1, a2,epp_CALET_data['Ener']),color='green',ls='--', alpha=0.3)


#plt.plot(e_data['Ener'], (e_data['Ener']**3)*electrons(ml_Ce,ml_gammae,ml_Cp,ml_gammap,ml_Cs,ml_gammas,ml_Esi),color='black',label='model fitted')
#plt.plot(e_data['Ener'], (e_data['Ener']**3)*sour(ml_Cs, ml_gammas, ml_Esi,e_data['Ener']),color='red',label='source fitted')
plt.errorbar(epp_AMS02_data['Ener'], (epp_AMS02_data['Ener']**3)*epp_AMS02_data['VALUE'], yerr=(epp_AMS02_data['Ener']**3)*epp_AMS02_data['ERR_STATma'],color='#3288bd',label='AMS02', fmt='o',alpha=0.5)
plt.errorbar(epp_dampe_data['Ener'], (epp_dampe_data['Ener']**3)*epp_dampe_data['VALUE'], yerr=(epp_dampe_data['Ener']**3)*epp_dampe_data['ERR_STATma'],color='blue',label='DAMPE', fmt='o',alpha=0.5)
plt.errorbar(epp_CALET_data['Ener'], (epp_CALET_data['Ener']**3)*epp_CALET_data['VALUE'], yerr=(epp_CALET_data['Ener']**3)*epp_CALET_data['ERR_STATma'],color='green',label='CALET', fmt='o',alpha=0.5)
#plt.errorbar(p_data['Ener'], (p_data['Ener']**3)*(p_data['VALUE']), yerr=(p_data['Ener']**3)*p_data['ERR_STATma'],color='red',label='Positron AMS02', fmt='s',alpha=0.5)
#plt.errorbar(e_data['Ener'], (e_data['Ener']**3)*e_data['VALUE'], yerr=(e_data['Ener']**3)*e_data['ERR_STATma'],color='yellow',label='Electron AMS02', fmt='^',alpha=0.5)
print(Mdm, 10**logfluxNorm, a1, a2)

ax.legend()
ax.legend(loc='best')
ax.legend(fontsize=9)
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylim((0, 250))
ax.set_xlim(5, 1e4)
ax.minorticks_on()
ax.legend(loc='best')
ax.set_xlabel('Energy [GeV]',fontsize=15)
ax.set_ylabel('Electron + Positron ',fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig(datafile + '_4_electron_plus_positron.pdf')



###############PLOT GAMMA RAY FROM FIT TO EVERYTHING#############

EGB_upper_lower_limits_file = np.loadtxt('EGB_summed1_Ajello.dat')# This file contains all the background components added considering the uncertainty in each of them
EGB_data_file = np.loadtxt('EGB_data_comparison.dat')

EGB_upper_lower_limits_energy=[]
EGB_upper_lower_limits_upper=[]
EGB_upper_lower_limits_lower=[]
EGB_upper_lower_limits_mid=[]

for i in range(EGB_upper_lower_limits_file.shape[0]):
    EGB_upper_lower_limits_energy.append(EGB_upper_lower_limits_file[i][0])
    EGB_upper_lower_limits_upper.append(EGB_upper_lower_limits_file[i][2])
    EGB_upper_lower_limits_lower.append(EGB_upper_lower_limits_file[i][1])
    EGB_upper_lower_limits_mid.append(np.power((EGB_upper_lower_limits_file[i][1]+EGB_upper_lower_limits_file[i][2]),0.5))


EGB_data_energy_mid=[]
EGB_data_energy_lower=[]
EGB_data_energy_upper=[]
EGB_data_flux=[]
EGB_data_yerror=[]
EGB_data_UL_DM=[]
EGB_data_UL_DM_L=[]



for i in range(EGB_data_file.shape[0]):
    EGB_data_energy_mid.append(np.power((EGB_data_file[i][0]*EGB_data_file[i][1]),0.5) )
    EGB_data_UL_DM.append( (EGB_data_file[i][0]*EGB_data_file[i][1]) * EGB_data_file[i][8]/(EGB_data_file[i][1] - EGB_data_file[i][0])  )
    EGB_data_UL_DM_L.append( (EGB_data_file[i][0]*EGB_data_file[i][1]) * 0.2* EGB_data_file[i][8]/(EGB_data_file[i][1] - EGB_data_file[i][0])  )
    EGB_data_flux.append( (EGB_data_file[i][0]*EGB_data_file[i][1]) * EGB_data_file[i][2]/(EGB_data_file[i][1] - EGB_data_file[i][0])  )
    EGB_data_yerror.append( (EGB_data_file[i][0]*EGB_data_file[i][1]) * (2*EGB_data_file[i][3])/(EGB_data_file[i][1] - EGB_data_file[i][0])  )




##############################################################################


#for (Cp, gammap, Ce, gammae, Mdm, logfluxNorm, a1, a2, g, h, i, j) in a.get_equal_weighted_posterior()[::100,:-1]:
#           plt.plot(EGB_data_energy_mid, EffTotGravFluxPho(Mdm, 10**logfluxNorm, a1, a2,EGB_data_energy_mid),color='green', alpha=0.3)


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

x = epp_dampe_data['Ener']
shade = Shade(x)
for (Cp, gammap, Ce, gammae, Mdm, logfluxNorm, a1, a2, g, h, i, j) in a.get_equal_weighted_posterior()[::10,:-1]:
	shade.add(EffTotGravFluxPho(Mdm, 10**logfluxNorm, a1, a2,epp_dampe_data['Ener']))
shade.line(color='k')
shade.shade(q=0.25, color='#008837', alpha=0.3)
shade.shade(q=0.5, color='#7b3294', alpha=0.3)
plt.errorbar(EGB_data_energy_mid, EGB_data_flux, yerr=EGB_data_yerror,color='#3288bd',label=' EGB data', fmt='o',alpha=.7)
plt.errorbar(EGB_data_energy_mid,EGB_data_UL_DM, yerr=EGB_data_UL_DM_L, uplims=True,color='red', label='DM upper limit', linestyle='solid', alpha=0.7)
#plt.fill_between(excess['Energy'], maxfluxgFB['GCexcess_Sample'], minfluxgFB['GCexcess_Sample'], color='#e34a33',alpha=0.7,label='Syst. modeling bck')
plt.fill_between(EGB_upper_lower_limits_energy, EGB_upper_lower_limits_upper, EGB_upper_lower_limits_lower, color='#fdbb84',label='Total bckg models',alpha=0.5)
ax.legend()
ax.legend(loc='best')
ax.legend(fontsize=9)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim((3e-9, 12e-6))
ax.set_xlim((1e1, 2e3))
ax.set_xlabel('Energy [GeV]',fontsize=15)
ax.set_ylabel('$E^2 dN/dE$ [GeV/cm$^2$s]',fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig(datafile + '_4_gamma_rays_EGB.pdf')


#############################################################################





fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
values = a.get_equal_weighted_posterior()
bestfit=a.get_best_fit()
cm = plt.cm.get_cmap('RdYlBu')
sc = plt.scatter(values[:,4],1e26/(10**values[:,5]),c=values[:,12], vmin=min(values[:,12]), vmax=max(values[:,12]), marker="o",alpha=0.5, label="all points")
#plt.scatter(values[:,4],10**values[:,5],alpha=0.9, label="all points")
plt.colorbar(sc,label='LogLikelihood')
plt.scatter(bestfit['parameters'][4],1e26/(10**bestfit['parameters'][5]), alpha=1, edgecolors="black",marker=(5, 1), label="best fit")
ax.legend()
ax.legend(loc='best')
ax.legend(fontsize=9)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='c', scilimits=(0,0))
#plt.ticklabel_format(axis='y',style='sci',scilimits=(26,27),useOffset=True)
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim((1e26, 1e27))
ax.set_xlim((1e3, 4.5e3))
ax.set_xlabel('gravitino mass [GeV]',fontsize=13)
ax.set_ylabel('Gravitino lifetime [s]',fontsize=13)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig(datafile + '_4_mass_tau_LogL.pdf')


a_lnZ = a.get_stats()['global evidence']
print ()
print( '************************')
print( 'MAIN RESULT: Evidence Z ')
print( '************************')
print( '  log Z for model = %.1f' % (a_lnZ / np.log(10)))
print()
toc = time.clock()
elapse=toc - tic
print("elapse time in min=",elapse/60)
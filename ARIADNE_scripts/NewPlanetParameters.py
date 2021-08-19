import numpy as np
import pickle
import pandas as pd
from astropy import constants as c
from astropy import units as u
import matplotlib.pyplot as plt

#####getting TOI name from TOI list#####
toi_list = pd.read_csv('ARIADNE_TOI_Results.csv', dtype=dict(GAIA=str))

#############create empty dataframe to store results###############
df = pd.DataFrame(columns=['TOI Number','TIC_Rstar', 
                'TIC_Mstar', 'Rplanet_old', 'Rplanet_old_e', 'Mplanet_old', 
                'semi_amp_old', 'ARIADNE_Rstar', 'ARIADNE_Mstar', 
                'ARIADNE_lum', 'Rplanet_new','Rplanet_new_e',
                'Mplanet_new', 'semi_amp_new'])

###################################################
########Find planet parameters for each TOI########
###################################################
teff = []
teff_e_old = [] #tic
teff_e_new = [] #ariadne
rad_diff = []
rstar_e_old = [] #tic
rstar_e_new = [] #ariadne

for m in range(len(toi_list)):
    word = toi_list['TOI_Number'][m]
    starname = 'TOI-'+ word.split('.')[0]
    toi_number = word.split('.')[0]

    if np.isnan(toi_list['ARIADNE_Radius'][m]) == False:  #if ariadne output exists for star, calculate new planet parameters
        
        #pickle file for corresponding TOI
        ariadne_results = pickle.load(open('../ARIADNE_FitResults/'+ starname +'/BMA.pkl', 'rb'))

        ########getting star parameters from pickle file#########
        rstar_new = ariadne_results['best_fit_averaged']['rad']
        teff_new = ariadne_results['best_fit_averaged']['teff']
        logg_new = ariadne_results['best_fit_averaged']['logg']
        lum_new = ariadne_results['best_fit_averaged']['lum']
        iso_mass_new = ariadne_results['best_fit_averaged']['iso_mass']


        ########original values from “TOIList_SG4.csv” file########
        dtype_dict={'TIC':str,'TOI':str,'SG1Disposition':str,'Rplanet':np.float64,
                    'Period':np.float64,'CK2016_Mass':np.float64,'PredictedK':np.float64,
                    'Mstar_TIC8':np.float64,'Rstar_TIC8':np.float64,'Teff_TIC8':np.float64}

        old_values=pd.read_csv('TOI_SG4List.csv',header=0,dtype=dtype_dict)

        #new subset dataframe with info for specific TOI & its planets
        TOI_sg4_values = old_values[old_values['TOI'].str.startswith(toi_number+'.')]

        #import csv of full toi list to get old radius errors
        full_toi_list = pd.read_csv('FullTOI_List.csv', dtype={'TOI':str, 'Planet Radius (R_Earth) err':np.float64,'Depth (ppm)': np.float64, 'Depth (ppm) err': np.float64})

        ###############################################################
        ## iterate over all rows of returned data frame to calculate ## 
        ###    parameters for all planets orbiting a single TOI     ###
        ###############################################################
        for i in TOI_sg4_values.index:
            TOI_name = TOI_sg4_values['TOI'][i]
            rstar_old = TOI_sg4_values['Rstar_TIC8'][i]
            Rp_old = TOI_sg4_values['Rplanet'][i]
            Rp_old_e = full_toi_list[full_toi_list['TOI']==TOI_name]['Planet Radius (R_Earth) err']
            depth = full_toi_list[full_toi_list['TOI']==TOI_name]['Depth (ppm)'].item()
            depth_e = full_toi_list[full_toi_list['TOI']==TOI_name]['Depth (ppm) err'].item()
            mstar_old = TOI_sg4_values['Mstar_TIC8'][i]
            Mplanet_old= TOI_sg4_values['CK2016_Mass'][i]
            rv_old = TOI_sg4_values['PredictedK'][i]
            Period = TOI_sg4_values['Period'][i]

            teff.append(toi_list['ARIADNE_Teff'][m])
            teff_e_old.append(toi_list['TICv8_Teff_err'][m])
            teff_e_new.append((toi_list['ARIADNE_teff_upper'][m]+toi_list['ARIADNE_teff_lower'][m])/2)
            rstar_e_old.append(toi_list['TICv8_Radius_err'][m])

            e_new_upper = float(toi_list['ARIADNE_rad_upper'][m])
            e_new_lower = float(toi_list['ARIADNE_rad_lower'][m])
            e_new_avg = (e_new_upper + e_new_lower)/2
            rstar_e_new.append(e_new_avg)

            ########find updated planet radius########
            Rp_new = (rstar_new/rstar_old) * Rp_old

            k = depth / 1e6 #fractional depth (original in ppm)
            k_e = depth_e / 1e6
            Rp_new_e = np.sqrt( (k * (e_new_avg**2)) + (.25 * (k_e**2) * (rstar_new**2)/k))

            ##use updated radius to find estimated planet mass using CK2016 relations##
            if Rp_new < 1.23:
                C = .00346
                S = .279

            if 1.23 < Rp_new < 11.1:
                C = -.0925
                S = .589

            if Rp_new >= 14.3:
                C = -2.85
                S = .881

            log_m = (np.log10(Rp_new) - C ) / S
            Mplanet_new = 10**(log_m)

            ########rv semi-amplitude value########

            # Normalization.
            # RV m/s of a 1.0 Jupiter mass planet tugging on a 1.0
            # solar mass star on a 1.0 year orbital period
            K_0 = 28.4329

            """Compute Doppler semi-amplitude

            Args:
                Msini (float): mass of planet [Mearth]
                P (float): Orbital period [days]
                Msun (float): Mass of star [Msun]
                e (float): eccentricity

            Returns:
                Doppler semi-amplitude [m/s]

            """

            # convert inputs to array so they work with units
            P = np.array(Period)  #keep array for TOIs with multiple planets
            Msini = np.array(Mplanet_new) #new planet mass (earth masses)
            Msun = np.array(iso_mass_new) #new stellar mass (solar masses)
            e = np.array(0) #assume 0

            P = (P * u.d).to(u.year).value
            Mtotal= (Msun * u.M_sun).value + (Msini * u.M_earth).to(u.M_sun).value
            Msini = (Msini * u.M_earth).to(u.M_jup).value

            new_semi_amplitude = K_0*(1 - e**2)**-0.5*Msini*P**(-1.0/3.0)*Mtotal**(-2.0 / 3.0)

            #if the period is unknown, assign rv semi-amplitude of -1 so we know to ignore it
            if Period == 9999:
                new_semi_amplitude = -1

            if 11.1 <= Rp_new <= 14.3:
                Mplanet_new = -1
                new_semi_amplitude = -1

            #append rows to empty dataframe after gathering all relevant information
            df = df.append({'TOI Number': TOI_name,'TIC_Rstar': rstar_old, 'TIC_Mstar': mstar_old, 
                            'Rplanet_old': Rp_old, 'Rplanet_old_e' : Rp_old_e,
                            'Mplanet_old': Mplanet_old, 'semi_amp_old': rv_old, 'ARIADNE_Rstar': rstar_new, 
                            'ARIADNE_Mstar': iso_mass_new, 'ARIADNE_lum': lum_new, 
                            'Rplanet_new': Rp_new, 'Rplanet_new_e': Rp_new_e, 'Mplanet_new': Mplanet_new, 
                            'semi_amp_new': new_semi_amplitude}, 
                            ignore_index = True)
            rad_diff.append(np.abs(Rp_old - Rp_new))

print(df) #final dataframe for planet/star information for all TOIs
df['RadDiff'] = rad_diff
df['teff_e_old'] = teff_e_old
df['teff_e_new'] = teff_e_new
df['rstar_e_old'] = rstar_e_old
df['rstar_e_new'] = rstar_e_new

###########################################################
#######plotting histograms for each planet parameter#######
###########################################################

#radius
fig = plt.figure()
ax = fig.add_subplot(111)
plt.hist(df['Rplanet_new'], bins=np.logspace(np.log10(0.01),np.log10(15), 50), label='New', alpha=.5)
plt.hist(df['Rplanet_old'], bins=np.logspace(np.log10(0.01),np.log10(15), 50), label='Old', alpha=.5)
plt.xlabel('Radius (R$_{\oplus}$)')
plt.ylabel('Count')
plt.title('Histogram of Radii for TOI Planets')
plt.legend()
plt.xscale('log')
fig.savefig('../TargetListOverview/PlanetParameterResults/radius_hist.png', dpi=300)

#Mass
fig = plt.figure()
ax = fig.add_subplot(111)
plt.hist(df['Mplanet_new'][df['Mplanet_new']>0], bins=np.logspace(np.log10(0.01),np.log10(200), 50), label='New', alpha=.5)
plt.hist(df['Mplanet_old'], bins=np.logspace(np.log10(0.01),np.log10(200), 50), label='Old', alpha=.5)
plt.xlabel('Mass (M$_{\oplus}$)')
plt.ylabel('Count')
plt.title('Histogram of Masses for TOI Planets')
plt.legend()
plt.xscale('log')
fig.savefig('../TargetListOverview/PlanetParameterResults/mass_hist.png', dpi=300)


#Semi-amplitude
fig = plt.figure()
ax = fig.add_subplot(111)
plt.hist(df['semi_amp_new'][df['semi_amp_new']>0], bins=np.logspace(np.log10(0.1),np.log10(200), 50), label='New', alpha=.5)
plt.hist(df['semi_amp_old'], bins=np.logspace(np.log10(0.1),np.log10(200), 50), label='Old', alpha=.5)
plt.xlabel('RV Semi-Amplitude (m/s)')
plt.ylabel('Count')
plt.title('Histogram of RV Semi-Amplitudes for TOI Planets')
plt.legend()
plt.xscale('log')
fig.savefig('../TargetListOverview/PlanetParameterResults/semi_amp_hist.png', dpi=300)


####################################################################
####################### Comparison Plots ###########################
####################################################################

#comparing new to old mass
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot([0,400],[0,400], '--')
plt.scatter(df['Mplanet_old'][df['Mplanet_new']>0], df['Mplanet_new'][df['Mplanet_new']>0],alpha=.8)
plt.xlabel('Old Values (M$_{\oplus}$)')
plt.ylabel('New Values (M$_{\oplus}$)')
plt.title('New vs Old Planet Mass values')
plt.xscale('log')
plt.yscale('log')
fig.savefig('../TargetListOverview/PlanetParameterResults/mass_comparison.png', dpi=300)

#comparing new to old radius
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot([0,400],[0,400], '--')
plt.scatter(df['Rplanet_old'], df['Rplanet_new'],alpha=.8)
plt.xlabel('Old Values (R$_{\oplus}$)')
plt.ylabel('New Values (R$_{\oplus}$)')
plt.title('New vs Old Planet Radius values')
plt.xscale('log')
plt.yscale('log')
plt.xlim([0,15])
plt.ylim([0,15])
fig.savefig('../TargetListOverview/PlanetParameterResults/radius_comparison.png', dpi=300)


####################################################################
#################### New - Old Value Plots #########################
####################################################################

#new - old radius
fig = plt.figure()
ax = fig.add_subplot(111)
plt.scatter(teff,df['Rplanet_new']-df['Rplanet_old'])
plt.xlabel('Teff (K)')
plt.ylabel('New-Old Radius Values (R$_{\oplus}$)')
plt.title('New - Old Planet Radius values')
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim([0,15])
plt.ylim([-10,10])
fig.savefig('../TargetListOverview/PlanetParameterResults/radius_diff_comparison.png', dpi=300)

#new-old mass
fig = plt.figure()
ax = fig.add_subplot(111)
plt.scatter(teff,(df['Mplanet_new']-df['Mplanet_old']))
plt.xlabel('Teff (K)')
plt.ylabel('New-Old Mass Values (M$_{\oplus}$)')
plt.title('New - Old Planet Mass values')
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim([0,15])
plt.ylim([-10,10])
fig.savefig('../TargetListOverview/PlanetParameterResults/mass_diff_comparison.png', dpi=300)

#new-old semi-amp
fig = plt.figure()
ax = fig.add_subplot(111)
plt.scatter(teff,(df['semi_amp_new']-df['semi_amp_old']))
plt.xlabel('Teff (K)')
plt.ylabel('New-Old Semi-Amplitude Values (m/s)')
plt.title('New - Old Planet Semi-Amplitude values')
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim([0,15])
plt.ylim([-10,10])
fig.savefig('../TargetListOverview/PlanetParameterResults/semiamp_diff_comparison.png', dpi=300)

########################################################
######################find outliers#####################
########################################################
pd.set_option('max_columns', None)
#print(df.loc[df['TOI Number'] == '924.01'])
#print(df.loc[df['TOI Number'] == '2666.01'])


#plot to see impact on radii
fig = plt.figure()
ax = fig.add_subplot(111)
plt.scatter(teff, ((df['Rplanet_new'] - df['Rplanet_old'])/df['Rplanet_old_e']))
plt.ylabel('(New Radius - Old Radius) / Old Radius Error')
plt.xlabel('Teff (K)')
plt.title('Difference in Radius / Original Radius Error')
#plt.ylim([-1,5])
fig.savefig('../TargetListOverview/PlanetParameterResults/radius_impact_plot', dpi=300)

#uncertainties of stellar radii plot as function of teff
fig, (ax1,ax2) = plt.subplots(1,2, figsize=(14,10))
ax1.plot(teff, [i/j for i, j in zip(teff_e_old,teff_e_new)],'o',c='orange')
ax1.set_title('Effective Temperature',size=18)
ax1.set_xlabel('Teff (K)',size=16)
ax1.set_ylabel('TIC / ARIADNE Teff Uncertainty',size=16)
ax1.set_ylim([0,8])
ax1.tick_params(axis='both', which='major', labelsize=16)
#ax1.set_xlabel('Planet radius (R$_{\oplus}$)',fontsize=30,fontweight='bold')
#ax1.set_ylabel('Planet mass (M$_{\oplus}$)',fontsize=30,fontweight='bold')
plt.title('Ratio of Uncertainties as a Function of Spectral Type')

ax2.plot(teff, [i / j for i, j in zip(rstar_e_old,rstar_e_new)], 'o')
ax2.set_title('Radius',size=18)
ax2.set_xlabel('Teff (K)',size=16)
ax2.set_ylabel('TIC / ARIADNE Radius Uncertainty',size=16)
ax2.set_ylim([0,8])
ax2.tick_params(axis='both', which='major', labelsize=16)
fig.savefig('../TargetListOverview/PlanetParameterResults/uncertainty_ratio_plot', dpi=300)

#uncertainties of planet radii (old & new) as function of new planet radius
fig, (ax1,ax2) = plt.subplots(1,2, figsize=(14,10))
ax1.plot(df['Rplanet_new'], df['Rplanet_new_e']/df['Rplanet_new'],'o',c='orange')
ax1.set_title('New Planet Radius Values',size=18)
ax1.set_xlabel('New Planet Radius (R$_{\oplus}$)',size=16)
ax1.set_ylabel('New Radius Error / New Radius',size=16)
#ax1.set_ylim([0,8])
ax1.tick_params(axis='both', which='major', labelsize=16)

ax2.plot(df['Rplanet_new'], df['Rplanet_old_e']/df['Rplanet_old'], 'o')
ax2.set_title('Old Planet Radius Values',size=18)
ax2.set_xlabel('New Planet Radius (R$_{\oplus}$)',size=16)
ax2.set_ylabel('Old Radius Error / Old Radius',size=16)
#ax2.set_ylim([0,8])
ax2.tick_params(axis='both', which='major', labelsize=16)
fig.savefig('../TargetListOverview/PlanetParameterResults/uncertainty_planet_ratio_plot', dpi=300)

#ratio of uncertainties of planet radii (old / new) & diff in radii as function of new planet radius
fig, (ax1,ax2) = plt.subplots(1,2, figsize=(14,10))
ax1.plot(df['Rplanet_new'], df['Rplanet_old_e']/df['Rplanet_new_e'],'o',c='orange')
ax1.set_title('Ratio of Radius Error Values',size=18)
ax1.set_xlabel('New Planet Radius (R$_{\oplus}$)',size=16)
ax1.set_ylabel('Old Radius Error / New Radius Error',size=16)
#ax1.set_ylim([0,.4])
ax1.tick_params(axis='both', which='major', labelsize=16)

ax2.plot(df['Rplanet_new'], df['Rplanet_new'] - df['Rplanet_old'], 'o')
ax2.set_title('Difference between Planet Radius Values',size=18)
ax2.set_xlabel('New Planet Radius (R$_{\oplus}$)',size=16)
ax2.set_ylabel('New Radius - Old Radius',size=16)
#ax2.set_ylim([0,1])
ax2.tick_params(axis='both', which='major', labelsize=16)
fig.savefig('../TargetListOverview/PlanetParameterResults/planet_radius_diff_ratio_plot', dpi=300)

print('new errors=',df['Rplanet_new_e'])
print('old errors=',df['Rplanet_old_e'])

#print('median teff error (old)=', np.median(teff_e_old),df.loc[df['teff_e_old']==np.median(teff_e_old),'TOI Number'])
#print('median teff error (new)=', np.median(teff_e_new),df.loc[df['teff_e_new']==np.median(teff_e_new),'TOI Number'])
#print('median radius error (old)=', np.median(rstar_e_old),df.loc[df['rstar_e_old']==np.median(rstar_e_old),'TOI Number'])
#print('median radius error (new)=', np.median(rstar_e_new),df.loc[df['rstar_e_new']==np.median(rstar_e_new),'TOI Number'])
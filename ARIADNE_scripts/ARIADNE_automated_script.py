from astroARIADNE.star import Star
from astroARIADNE.fitter import Fitter
from astroARIADNE.plotter import SEDPlotter
import pandas as pd
import numpy as np


toi_list = pd.read_csv('Cycle3_TOI_List.csv')


#Definition for string syntax

tick = "'"



###################################################################
## The lines below will need to be changed on a per-target basis ##
###################################################################

for i in range(len(toi_list)):
    word = toi_list['TOI_Number'][i]
    starname = 'TOI-'+ word.split('.')[0]
    ra = toi_list['RA'][i]
    dec = toi_list['Dec'][i]
    gaia_id = toi_list['GAIA'][i]


    #checking for priors

    #teff
    if np.isnan(toi_list['Teff'][i]) == True:  #if star has no priors, set to default
        teff_default_str = tick + 'teff' + tick + ': (' + tick + 'default' + tick + ')'
    else:
        teff= toi_list['Teff'][i]
        teff_unc= toi_list['e_Teff'][i]

        teff_default_str = tick + 'teff' + tick + ': (' + tick + 'normal' + tick + ', ' + 'teff' + ', ' + 'teff_unc' + ')'
    

    #logg
    if np.isnan(toi_list['logg'][i]) == True:
        logg_default_str = tick + 'logg' + tick + ': (' + tick + 'default' + tick + ')'
    else:
        logg = toi_list['logg'][i]
        logg_unc = toi_list['e_logg'][i]

        logg_default_str = tick + 'logg' + tick + ': (' + tick + 'normal' + tick + ', ' + 'logg' + ', ' + 'logg_unc' + ')' 


    #Fe/H
    if np.isnan(toi_list['metallicity'][i]) == True:
        feh_default_str = tick + 'z' + tick + ': (' + tick + 'default' + tick + ')'
    else:
        feh= toi_list['metallicity'][i]
        feh_unc= toi_list['e_metallicity'][i]

        feh_default_str = tick + 'z' + tick + ': (' + tick + 'normal' + tick + ', ' + 'feh' + ', ' + 'feh_unc' + ')'


    #Av
    if np.isnan(toi_list['av'][i]) == True:
        av_default_str = tick + 'Av' + tick + ': (' + tick + 'default' + tick + ')'
    else:
        av = toi_list['av'][i]
        av_unc = toi_list['e_av'][i]

        av_default_str = tick + 'Av' + tick + ': (' + tick + 'normal' + tick + ', ' + 'av' + ', ' + 'av_unc' + ')'
    

    ###############################################################
    ##         From here on, everything should be static         ##
    ###############################################################

    out_folder='../ARIADNE_FitResults/'+starname
    in_file = out_folder + 'BMA.pkl'
    plots_out_folder=out_folder+'plots/'

    s = Star(starname, ra, dec, g_id=gaia_id)

    #s.estimate_logg()

    engine = 'dynesty'
    nlive = 500
    dlogz = 0.5
    bound = 'multi'
    sample = 'rwalk'
    threads = 4
    dynamic = False

    setup = [engine, nlive, dlogz, bound, sample, threads, dynamic]

    # Feel free to uncomment any unneeded/unwanted models
    models = ['phoenix','btsettl','btnextgen','btcond','kurucz','ck04']

    f = Fitter()
    f.star = s
    f.setup = setup
    f.av_law = 'fitzpatrick'
    f.out_folder = out_folder
    f.bma = True
    f.models = models
    f.n_samples = 100000

    f.prior_setup = {teff_default_str,logg_default_str,feh_default_str,'dist': ('default'),'rad': ('default'),av_default_str}

    f.initialize()
    f.fit_bma()

    artist = SEDPlotter(in_file, plots_out_folder)

    artist.plot_SED_no_model()
    artist.plot_bma_hist()
    artist.plot_bma_HR(10)
    artist.plot_corner()
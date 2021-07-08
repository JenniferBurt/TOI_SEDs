from astroARIADNE.star import Star
from astroARIADNE.fitter import Fitter, multinest_log_like
from astroARIADNE.plotter import SEDPlotter
import pandas as pd
import numpy as np


toi_list = pd.read_csv('../TargetListOverview/Cycle3_TOI_List.csv')


#Definition for string syntax

tick = "'"



###################################################################
## The lines below will need to be changed on a per-target basis ##
###################################################################

for i in range(0,len(toi_list)):
    word = toi_list['TOI_Number'][i]
    starname = 'TOI-'+ word.split('.')[0]
    ra = toi_list['RA'][i]
    dec = toi_list['Dec'][i]
    gaia_id = toi_list['GAIA'][i]


    #checking for priors by building dictionary one parameter at a time
    my_dict = {}

    #teff
    if np.isnan(toi_list['Teff'][i]) == True:  #if star has no priors, set to default
        teff_dict = {'teff': 'default'}
        my_dict.update(teff_dict)  #update dictionary with teff info
    else:
        teff= toi_list['Teff'][i]      #if info is listed on csv, use those values as priors
        teff_unc= toi_list['e_Teff'][i]

        teff_dict = {'teff': ('normal', teff, teff_unc)}
        my_dict.update(teff_dict)
    

    #logg
    if np.isnan(toi_list['logg'][i]) == True:
        logg_dict = {'logg': 'default'}
        my_dict.update(logg_dict)

    else:
        logg = toi_list['logg'][i]
        logg_unc = toi_list['e_logg'][i]

        logg_dict = {'logg': ('normal', logg, logg_unc)}
        my_dict.update(logg_dict)


    #Fe/H
    if np.isnan(toi_list['metallicity'][i]) == True:
        feh_dict = {'z': 'default'}
        my_dict.update(feh_dict)

    else:
        feh= toi_list['metallicity'][i]
        feh_unc= toi_list['e_metallicity'][i]

        feh_dict = {'z': ('normal', feh, feh_unc)}
        my_dict.update(feh_dict)


    #distance - no priors
    dist_dict = {'dist': 'default'}
    my_dict.update(dist_dict)

    #radius - no priors
    rad_dict = {'rad': 'default'}
    my_dict.update(rad_dict)


    #Av - priors from Stilism
    if np.isnan(toi_list['av'][i]) == True:
        av_dict = {'Av': 'default'}
        my_dict.update(av_dict)

    else:
        av = toi_list['av'][i]
        av_unc = toi_list['e_av'][i]

        av_dict = {'Av': ('normal', av, av_unc)}
        my_dict.update(av_dict)
    

    ###############################################################
    ##         From here on, everything should be static         ##
    ###############################################################

    out_folder='../ARIADNE_FitResults/'+starname
    in_file = out_folder + '/BMA.pkl'
    plots_out_folder=out_folder+'/plots/'

    s = Star(starname, ra, dec, g_id=gaia_id)

    s.estimate_logg()

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

    f.prior_setup = my_dict

    f.initialize()
    f.fit_bma()

    artist = SEDPlotter(in_file, plots_out_folder)

    artist.plot_SED_no_model()
    artist.plot_bma_hist()
    artist.plot_bma_HR(10)
    artist.plot_corner()
from astroARIADNE.star import Star
from astroARIADNE.fitter import Fitter
from astroARIADNE.plotter import SEDPlotter
import pandas as pd
import numpy as np


star_list = pd.read_csv('../TargetListOverview/Boyajian_Star_List.csv')


#Definition for string syntax

tick = "'"



###################################################################
## The lines below will need to be changed on a per-target basis ##
###################################################################

for i in range(0,len(star_list)):
    starname = star_list['ID'][i]
    ra = star_list['RA'][i]
    dec = star_list['Dec'][i]
    gaia_id = star_list['GAIA'][i]


    #checking for priors
    my_dict = {"teff": [],"logg": [],"z": [], "dist": 'default', "rad": 'default', 'Av': []}

    #teff
    if np.isnan(star_list['Teff'][i]) == True:  #if star has no priors, set to default
        my_dict["teff"].append("default")
    else:
        teff= star_list['Teff'][i]      #if info is listed on csv, use those values as priors
        teff_unc= star_list['e_Teff'][i]

        my_dict["teff"].append('normal')
        my_dict["teff"].append(teff)
        my_dict["teff"].append(teff_unc)
    

    #logg
    if np.isnan(star_list['logg'][i]) == True:
        my_dict["logg"].append('default')
    else:
        logg = star_list['logg'][i]
        logg_unc = star_list['e_logg'][i]

        my_dict["logg"].append('normal')
        my_dict["logg"].append(logg)
        my_dict["logg"].append(logg_unc)


    #Fe/H
    if np.isnan(star_list['metallicity'][i]) == True:
        my_dict["z"].append('default')
    else:
        feh= star_list['metallicity'][i]
        feh_unc= star_list['e_metallicity'][i]

        my_dict["z"].append('normal')
        my_dict["z"].append(feh)
        my_dict["z"].append(feh_unc)


    #Av
    if np.isnan(star_list['av'][i]) == True:
        my_dict['Av'].append('default')
    else:
        av = star_list['av'][i]
        av_unc = star_list['e_av'][i]

        my_dict["Av"].append('normal')
        my_dict["Av"].append(av)
        my_dict["Av"].append(av_unc)
    

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
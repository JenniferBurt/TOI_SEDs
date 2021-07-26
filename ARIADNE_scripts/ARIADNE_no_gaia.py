from astroARIADNE.star import Star
from astroARIADNE.fitter import Fitter
from astroARIADNE.plotter import SEDPlotter
import pandas as pd
import numpy as np


star_list = pd.read_csv('BoyajianStars.csv', dtype=dict(GAIA=str))


###################################################################
## The lines below will need to be changed on a per-target basis ##
###################################################################

for i in range(4,5):
    starname = star_list['ID'][i]
    ra = star_list['RA'][i]
    dec = star_list['Dec'][i]

    gaia_str = star_list['GAIA'][i]  #read in string value that preserves all digits of ID
    gaia_id = int(gaia_str)          #convert string to integer for ARIADNE to use


    #checking for priors, only want to use stilism Av values in this case
    my_dict = {}

    #teff
    teff_dict = {'teff': 'default'}
    my_dict.update(teff_dict)

    #logg
    logg_dict = {'logg': 'default'}
    my_dict.update(logg_dict)

    #feh
    feh_dict = {'z': 'default'}
    my_dict.update(feh_dict)

    #dist
    dist_dict = {'dist': 'default'}
    my_dict.update(dist_dict)

    #rad
    rad_dict = {'rad': 'default'}
    my_dict.update(rad_dict)

    #Av
    #if np.isnan(star_list['av'][i]) == True:
    #    av_dict = {'Av': 'default'}
    #    my_dict.update(av_dict)
    #else:
    #    av = star_list['av'][i]
    #    av_unc = star_list['e_av'][i]

    #   av_dict = {'Av': ('normal', av, av_unc)}
    #   my_dict.update(av_dict)

    av_dict = {'Av': ('fixed', 0)}
    my_dict.update(av_dict)

        
    

    ###############################################################
    ##         From here on, everything should be static         ##
    ###############################################################

    out_folder='../ARIADNE_FitResults/'+starname+'_no_gaia'
    in_file = out_folder + '/BMA.pkl'
    plots_out_folder=out_folder+'/plots/'

    s = Star(starname, ra, dec, g_id=gaia_id)

    #remove TESS magnitude if it exists
    s.remove_mag('TESS')
    #s.remove_mag('GALEX_NUV')
    s.remove_mag('GALEX_FUV')
    
    #remove gaia magnitudes
    s.remove_mag('GaiaDR2v2_G')
    s.remove_mag('GaiaDR2v2_BP')
    s.remove_mag('GaiaDR2v2_RP')

    #establish minimum value for all photometric errors of .03
    for i in range(len(s.filter_names)): #go through each filter ARIADNE could possibly use

        if s.used_filters[i] == 1.0:   #if s.used_filters = 1, that means ARIADNE is using that particular filter

            if s.mag_errs[i] < .03:  #if the error for that particular filter is less than .03...

                s.mag_errs[i] = .03  #set it to .03
            else:
                pass
        else:   #otherwise, leave the error as it is
            pass

    s.print_mags() #prints updated photometry just to check

    s.estimate_logg()

    engine = 'dynesty'
    nlive = 150
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
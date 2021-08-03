from astroARIADNE.star import Star
from astroARIADNE.fitter import Fitter, multinest_log_like
from astroARIADNE.plotter import SEDPlotter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


#read in cycle 3 TOI list containing star info
toi_list = pd.read_csv('Cycle3_TOI_List.csv', dtype=dict(GAIA=str))


###################################################################
## The lines below will need to be changed on a per-target basis ##
###################################################################

for i in range(110,115):   #iterates through csv to run fit for each star
    word = toi_list['TOI_Number'][i]
    starname = 'TOI-'+ word.split('.')[0]    #omits everything after the period (ex. TOI-2447.01 becomes TOI-2447)
    ra = toi_list['RA'][i]
    dec = toi_list['Dec'][i]

    gaia_str = toi_list['GAIA'][i]  #gaia number must be read in as string to preserve all digits
    gaia_id = int(gaia_str)    #convert string to integer for ARIADNE to read


    #checking for priors by building dictionary one parameter at a time
    my_dict = {}

    #teff
    teff_dict = {'teff': ('uniform', 3000, 7000)} #uniform priors from to-do list
    my_dict.update(teff_dict)
    

    #logg
    logg_dict = {'logg': ('uniform', 3.5, 5.5)}
    my_dict.update(logg_dict)


    #Fe/H
    feh_dict = {'z': ('uniform', -1, .6)}
    my_dict.update(feh_dict)


    #distance - no priors
    dist_dict = {'dist': ('default')}
    my_dict.update(dist_dict)

    #radius - no priors
    rad_dict = {'rad': ('default')}
    my_dict.update(rad_dict)


    #Av
    if toi_list['Distance'][i] < 50: #if distance is < 50pc, av = 0
        av_dict = {'Av': ('fixed', 0)}
        my_dict.update(av_dict)

    else: #use priors from stilism
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

    #remove GALEX & TESS mags
    s.remove_mag('TESS')
    s.remove_mag('GALEX_NUV')
    s.remove_mag('GALEX_FUV')

    g_corrected = None
    bp_corrected = None
    rp_corrected = None

    # if gaia g magnitude is within certain range, we must make a correction (Evans 2018)
    #gaia G band
    if 2 < s.mags[s.filter_names=='GaiaDR2v2_G'][0] < 6.5:

        #save gaia dr2 magnitudes + uncertainties for G, BP, and RP bandpasses
        gaia_g = s.mags[s.filter_names=='GaiaDR2v2_G'][0]
        g_uncertainty = s.mag_errs[s.filter_names=='GaiaDR2v2_G'][0]

        #perform correction from Evans (2018)
        g_corrected = -.047344 + 1.16405*gaia_g - 0.046799 * (gaia_g**2) +0.0035015 * (gaia_g**3)
    

    #range for gaia BP band correction
    if 2 < s.mags[s.filter_names=='GaiaDR2v2_G'][0] < 4:

        #get gaia magnitude + uncertainty
        gaia_bp = s.mags[s.filter_names=='GaiaDR2v2_BP'][0]
        bp_uncertainty = s.mag_errs[s.filter_names=='GaiaDR2v2_BP'][0]

        #correction
        bp_corrected = 1.95282 * gaia_bp - 0.11018* (gaia_bp**2) - 2.0384


    #gaia RP band correction
    if 2 < s.mags[s.filter_names=='GaiaDR2v2_G'][0] < 3.5:

        #gaia red-band magnitude + uncertainty
        gaia_rp = s.mags[s.filter_names=='GaiaDR2v2_RP'][0]
        rp_uncertainty = s.mag_errs[s.filter_names=='GaiaDR2v2_RP'][0]
    
        #correction
        rp_corrected = -13.946 + 14.239*gaia_rp - 4.23*(gaia_rp**2) + 0.4532*(gaia_rp**3)


    #if gaia mags were corrected, replace old photometry with corrected ones (use original uncertainties)
    if g_corrected is not None:
        s.remove_mag('GaiaDR2v2_G')
        s.add_mag(g_corrected, g_uncertainty, 'GaiaDR2v2_G')

    if bp_corrected is not None:
        s.remove_mag('GaiaDR2v2_BP')
        s.add_mag(bp_corrected, bp_uncertainty, 'GaiaDR2v2_BP')

    if rp_corrected is not None:
        s.remove_mag('GaiaDR2v2_RP')
        s.add_mag(rp_corrected, rp_uncertainty, 'GaiaDR2v2_RP')

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

################## Manual SED Plotting ####################
    """Plot raw photometry."""

    self = s
    
    # Get plot ylims
    yy=self.flux * self.wave
    ymax = np.max(yy[np.nonzero(yy)])
    ymin = np.min(yy[np.nonzero(yy)])

    f, ax = plt.subplots(figsize=(12,8))

    # Model plot
    used_f = self.filter_names[self.filter_mask]
    colors = np.array([
    'tomato', 'indianred', 'tab:red',
    'salmon', 'coral',
    'mediumorchid', 'mediumslateblue', 'tab:blue',
    'darkslateblue', 'darkblue',
    'olivedrab', 'yellowgreen', 'greenyellow', 'yellow',
    'orangered', 'chocolate', 'khaki',
    'limegreen', 'darkgreen', 'lime', 'seagreen', 'lawngreen', 'green',
    'aquamarine', 'turquoise', 'lightseagreen', 'teal', 'cadetblue',
    'firebrick', 'darkred',
    'blueviolet', 'darkviolet',
    'midnightblue', 'blue',
    'deeppink', 'fuchsia', 'mediumslateblue'])

    zipped_data=zip(colors[self.filter_mask],self.wave[self.filter_mask], self.flux[self.filter_mask], self.flux_er[self.filter_mask],self.bandpass[self.filter_mask], used_f)

    for c, w, fl, fe, bp, fi in zipped_data:    

        ax.errorbar(w, fl * w,
            xerr=bp, yerr=fe,
            fmt='',
            ecolor=c,
            marker=None)

        ax.scatter(w, fl * w,
            edgecolors='black',
            marker='o',
            c=c,
            s=60,
            alpha=0.85, label=fi)
    
    #    print(w,fl*w)

    ax.set_ylim([ymin * .8, ymax * 1.25])
    ax.set_xscale('log', nonposx='clip')
    ax.set_yscale('log', nonposy='clip')
    ax.set_ylabel(r'$\lambda$F$_\lambda$ (erg cm$^{-2}$s$^{-1}$)',
        fontsize=26,
        fontname='serif')
    ax.set_xlabel(r'Wavelength ($\mu$m)',
        fontsize=26,
        fontname='serif') 
    ax.legend(loc=0)

    ax.tick_params(axis='both', which='major',labelsize=22)
    ax.set_xticks(np.linspace(1, 10, 10))
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.set_xlim([0.1, 6])

    for tick in ax.get_yticklabels():
        tick.set_fontname('serif')

    plt.savefig(out_folder +'/'+starname+'_SED.png',bbox_inches='tight')
    plt.savefig('../ARIADNE_FitResults/ProblemStarSEDs/'+starname+'_SED.png',bbox_inches='tight')

###################### End Manual SED Plotting ########################

    #s.estimate_logg()

    #engine = 'dynesty'
    #nlive = 150
    #dlogz = 0.5
    #bound = 'multi'
    #sample = 'rwalk'
    #threads = 4
    #dynamic = False

    #setup = [engine, nlive, dlogz, bound, sample, threads, dynamic]

    # Feel free to uncomment any unneeded/unwanted models
    #models = ['phoenix','btsettl','btnextgen','btcond','kurucz','ck04']

    #f = Fitter()
    #f.star = s
    #f.setup = setup
    #f.av_law = 'fitzpatrick'
    #f.out_folder = out_folder
    #f.bma = True
    #f.models = models
    #f.n_samples = 100000

    #f.prior_setup = my_dict

    #f.initialize()
    #f.fit_bma()

    #artist = SEDPlotter(in_file, plots_out_folder)

    #artist.plot_SED_no_model()
    #artist.plot_bma_hist()
    #artist.plot_bma_HR(10)
    #artist.plot_corner()
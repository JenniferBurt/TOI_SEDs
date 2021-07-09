from astroARIADNE.star import Star
from astroARIADNE.fitter import Fitter
from astroARIADNE.plotter import SEDPlotter

###################################################################
## The lines below will need to be changed on a per-target basis ##
###################################################################
starname = 'HD 146233'
ra = 243.9063361
dec = -8.371641162
gaia_id = 4345775217221821312

#leave all priors as default except for av
av=0
av_unc=0.04965

###############################################################
##         From here on, everything should be static         ##
###############################################################

out_folder='../ARIADNE_FitResults/'+starname+'_default/'
in_file = out_folder + 'BMA.pkl'
plots_out_folder=out_folder+'plots/'

s = Star(starname, ra, dec, g_id=gaia_id)

#add 2MASS JK bands and WISE photometry
#s.add_mag(4.667, .26, '2MASS_J') #mag, uncertainty, filter name
#s.add_mag(4.186, 0.292, '2MASS_K')
#s.add_mag(3.988, 0.234, 'WISE_RSR_W1') #1st and 2nd WISE bands
#s.add_mag(3.642, 0.173, 'WISE_RSR_W2')

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

f.prior_setup = {'teff': ('default'),'logg': ('default'),'z': ('default'),'dist': ('default'),'rad': ('default'),'Av': ('normal', av, av_unc)}

f.initialize()
f.fit_bma()

artist = SEDPlotter(in_file, plots_out_folder)

artist.plot_SED_no_model()
artist.plot_bma_hist()
artist.plot_bma_HR(10)
artist.plot_corner()
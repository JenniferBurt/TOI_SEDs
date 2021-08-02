from astroARIADNE.star import Star
from astroARIADNE.fitter import Fitter
from astroARIADNE.plotter import SEDPlotter

###################################################################
## The lines below will need to be changed on a per-target basis ##
###################################################################
starname = 'NGTS-6'
ra = 75.795
dec = -30.399
gaia_id = 4875693023844840448

teff=4700
teff_unc=100

feh=0.11
feh_unc=0.09

av=0.02
av_unc=0.01
###############################################################
##         From here on, everything should be static         ##
###############################################################

out_folder='../ARIADNE_FitResults/'+starname+'/'
in_file = out_folder + 'BMA.pkl'
plots_out_folder=out_folder+'plots/'

s = Star(starname, ra, dec, g_id=gaia_id)

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

f.prior_setup = {'teff': ('normal', teff, teff_unc),'logg': ('default'),'z': ('normal', feh, feh_unc),'dist': ('default'),'rad': ('default'),'Av': ('normal', av, av_unc)}

f.initialize()
f.fit_bma()

artist = SEDPlotter(in_file, plots_out_folder)

artist.plot_SED_no_model()
artist.plot_bma_hist()
artist.plot_bma_HR(10)
artist.plot_corner()
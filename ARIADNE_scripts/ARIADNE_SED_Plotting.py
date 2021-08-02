from astroARIADNE.star import Star
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

ra = 75.795
dec = -30.399
starname = 'NGTS-6'
gaia_id = 4875693023844840448
out_folder='/Users/jenburt/Dropbox/Research/PythonPackages/astroARIADNE/SED_Tests/'

s = Star(starname, ra, dec, g_id=gaia_id)

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

plt.savefig(out_folder + starname+'_SED.png',bbox_inches='tight')
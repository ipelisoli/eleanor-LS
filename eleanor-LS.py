import eleanor
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import optimize
from astropy.io import ascii
from astropy.convolution import convolve, Box1DKernel
from astropy.stats import LombScargle,sigma_clip
import sys

## Function definitions ##

def periodogram(bjd, flux, flux_err):
    t = (bjd - bjd[0])*24.0
    mean = np.mean(flux)
    flux = flux - mean
    dt = [ t[i+1] - t[i-1] for i in range(1,len(t)-1)]
    fmax = 1.0/np.median(dt)
    fmin = 2.0/(max(t))
    ls = LombScargle(t, flux, flux_err)
    #Oversampling a factor of 10 to achieve frequency resolution
    freq, power = ls.autopower(minimum_frequency=fmin,
                               maximum_frequency=fmax,
                               samples_per_peak=10)
    best_f = freq[np.argmax(power)]
    period = 1.0/best_f #period from the LS periodogram
    fap_p = ls.false_alarm_probability(power.max())

    amp = np.sqrt(power)/mean
    return np.array(freq), np.array(amp), period, fap_p

def sine_func(x, a, b):
    return 1.0+a*np.sin(2.*np.pi*x + b)

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def avg_array(ar, n):
    ar2 = np.nanmedian(np.pad(ar, (0, n - ar.size%n), mode='constant',
                              constant_values=np.NaN).reshape(-1, n), axis=1)
    return ar2

def phase_data(bjd, flux, flux_err, period):
    phase = ((bjd - bjd[0])*24.0 / period) % 1.0
    flux_phased = [flux for phase, flux in sorted(zip(phase, flux))]
    flux_err_phased = [flux_err for phase, flux_err in sorted(zip(phase, flux_err))]
    phase = np.sort(phase)
    # Fit the data
    params, params_covariance = optimize.curve_fit(sine_func, phase,
                                                   flux_phased,
                                                   p0=[np.mean(flux), 0.0])
    flux_fit = 1.0 + params[0] * np.sin(2.*np.pi*phase + params[1])
    return np.array(phase), np.array(flux_phased), np.array(flux_err_phased), np.array(flux_fit), params[0]

#######

## User inputs ##
gaia_id = np.int(sys.argv[1])
aperture = input("Which aperture type would you like?\n(Check elanor documentation to see what this implies)\nOptions are normal, small, large: ")
clip = int(input("How many standard deviations should be clipped?\n"))
#######

# Download the data for this star
star = eleanor.multi_sectors(gaia=gaia_id, sectors='all')

# Do the photometry for all sectors using the selected aperture size

data = []
for s in star:
    datum = eleanor.TargetData(s, height=15, width=15, bkg_size=31,
                               aperture_mode=aperture, regressors='corner',
                               do_psf=False, do_pca=False)
    data.append(datum)

# Create flux vector

time = []
flux = []
err = []

for sector, datum in enumerate(data):
    # Only data not flagged for bad quality
    q = datum.quality == 0
    # Time
    time_i = datum.time[q].flatten()
    time = np.hstack([time, time_i])
    # Flux
    flux_i = (datum.corr_flux[q]/np.median(datum.corr_flux[q])).flatten()
    flux = np.hstack([flux, flux_i])
    # Uncertainty
    err_i = (datum.flux_err[q]/np.median(datum.corr_flux[q])).flatten()
    err = np.hstack([err, err_i])

# sigma-clipping
if clip > 0:
    filtered_data = sigma_clip(flux, sigma=clip, maxiters=None)
    clip = ~(filtered_data.mask)
    time = time[clip]
    flux = flux[clip]
    err = err[clip]

# calculate FT
freq, amp, period, fap = periodogram(time, flux, err)
# phase data
phase, flux_phased, flux_err_phased, flux_fit, p_amp = phase_data(time, flux, err, period)

fig = plt.figure(figsize=(12,10))
plt.rcParams.update({'font.size': 15})

gridspec.GridSpec(2,2)

plt.subplot2grid((2,2),(0,0),colspan=2,rowspan=1)
plt.xlim(np.min(time), np.max(time))
plt.xlabel("BJD - 2457000")
plt.ylabel("Relative flux")
plt.scatter(time, flux, marker='.', color='k', zorder = 1)

plt.subplot2grid((2,2),(1,0),colspan=1,rowspan=1)
plt.xscale("log")
plt.xlim(np.min(1/freq), np.max(1/freq))
plt.xlabel("Period (hours)")
plt.ylabel("Amplitude")
plt.plot(1/freq, amp, c = 'k')

plt.subplot2grid((2,2),(1,1),colspan=1,rowspan=1)
plt.ylabel("Relative flux")
plt.xlabel("Phase")
plt.xlim(0,2)
plt.scatter(phase, flux_phased, marker='.', color='0.75', zorder = 4)
plt.scatter(running_mean(phase, 50), running_mean(flux_phased, 50), marker='.',
            color='k', zorder = 5)
plt.plot(phase, flux_fit, 'r--', zorder = 6)
plt.scatter(1+phase, flux_phased, marker='.', color='0.75', zorder = 4)
plt.scatter(1+running_mean(phase, 50), running_mean(flux_phased, 50), marker='.',
            color='k', zorder = 5)
plt.plot(1+phase, flux_fit, 'r--', zorder = 6)

plt.tight_layout()

fig.savefig("GaiaDR2_%s_eleanor.png" % gaia_id)

print("I have a found a period of %7.3f hours with a FAP of %7.5g and amplitude of %5.2f per cent." % (period, fap, 100*abs(p_amp)))
print("Output has been saved to GaiaDR2_%s_eleanor.png." % gaia_id)

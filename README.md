# eleanor-LS

Given an object's Gaia DR2 source_id, this code uses eleanor (http://adina.feinste.in/eleanor/) to retrieve TESS FFI data for the object and perform its photometry.
It will then calculate a Lomb-Scargle periodogram and output the best period, with its false alarm probability (using the Baluev 2008, MNRAS 385, 1279 approximation) and amplitude. A plot showing the light curve, periodogram, and phased data will be created.

To run it, simply do 'python eleanor-LS.py [source_id]' and follow prompt instructions.

There are three types of possible apertures: normal, default, and large. In short, 'small' will only consider apertures 8 pixels in size or smaller, 'large' only considers apertures larger than 8 pixels. The former can be useful for faint stars, whereas the latter is recommended for bright stars. For more info, check the eleanor documentation.

IMPORTANT CAVEATS:
1. Take the reported amplitudes with a huge grain of salt. Due to the large TESS pixel size, contamination from other stars is almost inevitable, so the estimated amplitude largely depends on the chosen aperture.
2. Keep in mind the long exposure time of sectors 1-26 (30-minutes). If the period of your star is also of this order, you're likely to see two peaks in the light curve: the real period, and a combination between the period and the exposure time.

Requires numpy, matplotlib, eleanor, and astropy.

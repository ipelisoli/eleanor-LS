# eleanor-LS

Given an object's Gaia DR2 source_id, this code uses eleanor (http://adina.feinste.in/eleanor/) to retrieve TESS FFI data for the object and perform its photometry.
It will then calculate a Lomb-Scargle periodogram and output the best period, with its false alarm probability (using the Baluev 2008, MNRAS 385, 1279 approximation) and amplitude. A plot showing the light curve, periodogram, and phased data will be created.

To run it, simply do 'python eleanor-LS.py [source_id]'

Requires numpy, matplotlib, eleanor, and astropy.

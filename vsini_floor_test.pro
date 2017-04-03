; Writing this up in order to test the suspicious vsini_0 vs. teff
; relationship we found before.  When I plotted up vsini vs. Teff
; (just straight plotted it), it looked like the minimum measurement
; was a function of Teff (ie there were no measurements below some
; minimum vsini at a given Teff, and the minimum had a slope to it as
; a function of Teff)

; Here, my plan is to select a few spectra from each bin and see if I
; can either duplicate the measurement, or simulate what the detection
; floor should be.  Why does it appear to be Teff-dependent?

pro vsini_floor_test

; Read in allparam_dr13
; Choose some spectra
; Read in the chosen spectra
; CONSTRUCT mini grid of PHOENIX spectra centered on reported
; parameters
; Run my vsini rotation code

end

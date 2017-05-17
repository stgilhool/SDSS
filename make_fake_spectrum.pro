; Make a fake SDSS spectrum from a PHOENIX template

function make_fake_spectrum, TEFF=teff, $
                             FEH=feh, $
                             LOGG=logg, $
                             VSINI=vsini, $
                             EPSILON=epsilon, $
                             SNR=snr

; Readin the PHOENIX template corresponding to TEFF, FEH and LOGG

; Interpolate onto SDSS wl grid
                        
; Broaden with VSINI and EPSILON

; Add Poisson noise according to SNR

; Convolve with APOGEE LSF (need to make a fiducial LSF)

; Output the spectrum, along with relevant header info
; input variables
; bmin, bmax etc?



end

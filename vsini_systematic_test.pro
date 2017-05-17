; Procedure to test the PHOENIX template spectra for systematic vsini
; errors at cool temperatures/check that the 5km/s limit is reasonable
; at all temperatures

pro vsini_systematic_test

; Read in PHOENIX spectra at T=2500-4000

; Read in SDSS parameters, etc.
; LSF
; WL_SOLN
; DELTA_WL



; LOOP
;;; Make fake observation

; Rebin PHOENIX spectrum at the APOGEE resolution
; Add 1% noise
; Convolve with APOGEE LSF

;;; Measure Vsini
;; Model Observation
; Same process, but also convolve with broadening kernel


end

; Testing the DR13 ASPCAP parameter accuracy 
pro dr13_test

data_path = '/home/stgilhool/APOGEE/APOGEE_data/'

allstar_file = '/home/stgilhool/APOGEE/APOGEE_data/allStar-l30e.2.fits'

mstar_file = data_path + 'allparams_dr13.fits'


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read in the all Stars file and keep only our 1350

;star_info = mrdfits(allstar_file,1)

;mdwarf = where((star_info.apogee_target1 and 2L^19) ne 0L, nmdwarf)

;m_info = star_info[mdwarf]

; Write it out

;mwrfits, m_info, mstar_file, /create
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Read in the fits file we made
m_info = mrdfits(mstar_file,1)

; Read in my table
master_info = mrdfits('/home/stgilhool/APOGEE/master_info5.fits',1)

; Get the indices for IRTF entries
;irtf_idx = where(finite(master_info.teff_irtf) and m_info.teff gt 0,nirtf)
irtf_idx = where(finite(master_info.teff_irtf) and m_info.teff gt 0 and m_info.teff lt 4200,nirtf)
hot_idx = where(finite(master_info.teff_irtf) and m_info.teff gt 4000,nhot)
;irtf_idx = where(finite(master_info.teff_irtf), nirtf)

print, nirtf

teff_irtf = master_info[irtf_idx].teff_irtf
teff_sdss = m_info[irtf_idx].teff
e_teff_sdss = m_info[irtf_idx].teff_err


err_irtf = replicate(60, nirtf)
; Plot
line_x = dindgen(100)*500
line_y = line_x

res = teff_irtf-teff_sdss
rmse = stddev(res)

set_plot, 'ps'
device, filename = 'TeffDR13_test.eps'
device, /color, bits=8
device, xs = 13, ys= 8, /inches
loadct, 13


plot, teff_irtf, teff_sdss, ps = 4, yr=[2500,4500], /ys, xtitle = 'Teff_IRTF', ytit = 'Teff_ASPCAP', title = "Test of ASPCAP parameter reliability. RMSE = " + strtrim(rmse,2) + " K"
oploterror, teff_irtf, teff_sdss, err_irtf, e_teff_sdss, ps=4
 
oplot, line_x, line_y, linest=2

stop




end

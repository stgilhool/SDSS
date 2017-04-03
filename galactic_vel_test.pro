; See if I can get space velocities for some of my stars
pro galactic_vel_test

; Globals
skip_load=1
napo = 1350

; Master Info Table
mi_file = '/home/stgilhool/APOGEE/master_info6.fits'
mi_str = mrdfits(mi_file, 1)

; pmra, pmdec
ra_all = mi_str.ra
dec_all = mi_str.dec
pmra_all = mi_str.pmra
pmde_all = mi_str.pmde
pm = mi_str.pm
rpmra = mi_str.rpmra
rpmde = mi_str.rpmde
rpm = mi_str.rpm

; plx
plx = mi_str.plx_dtm

;vsini
vsini_vfit = mi_str.vsini_vfit

; Teff from aspcap
aspcap_str = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/allparams_dr13.fits',1)
teff_asp = aspcap_str.teff
teff_apg = mi_str.teff_apg
teff_vfit = mi_str.teff_vfit


if skip_load then goto, skip_load_mark
; Read in the apg headers to get galactic latitude
apg_files = readin_apo(nfiles = napo)
;apg_files = readin_apo(nfiles = 100)

save, /variables, filename = 'galactic_vel_test.sav'

skip_load_mark:
restore, 'galactic_vel_test.sav'

head_ptr = apg_files.header

; spec 17 has plx
test = *head_ptr[17]
; galactic lat, vrad
glat = fxpar(test, 'GLAT')
glon = fxpar(test, 'GLON')
vrad = fxpar(test, 'VRAD')

height = []
glat = []
dist = []
glon = []
parallax = plx/1d3
rv = []

for snum = 0, napo-1 do begin

    plxi = plx[snum]/1d3
    disti = 1d0/plxi
    dist = [dist, disti]

    ; get head
    head = *head_ptr[snum]
    
    glati = fxpar(head, 'GLAT')
    glati = glati * !const.DtoR
    glat = [glat, glati]

    gloni = fxpar(head, 'GLON')
    glati = glati * !const.DtoR
    glat = [glat, glati]

    rvi = fxpar(head,'VRAD')
    rv = [rv, rvi]

    if finite(plxi) and finite(glati) then begin
        heighti = disti*sin(glati)
    endif else heighti = !values.d_nan

    height = [height, heighti]

endfor

fin_idx = where(finite(height), nfin)


; all data needed for transformation to galactic velocities
distance1 = dist[fin_idx]
parallax1 = parallax[fin_idx]
ra1 = ra_all[fin_idx]
dec1 = dec_all[fin_idx]
pmra1 = pmra_all[fin_idx]
pmde1 = pmde_all[fin_idx]
vrad1 = rv[fin_idx]
; At this point, pm and plx are in arcsec/year and arcsec
; (Dittman gives plx in mas, and pm in arcsec)

; Data on the stars themselves
teff = teff_vfit[fin_idx]
vsini = vsini_vfit[fin_idx]

; Calculate UVW values
; if using DISTANCE instead of PLX, pm's must be in mas/year
; if using PLX instead of DISTANCE, pm's and plx must be in same arc units
gal_uvw, uarr, varr, warr, ra=ra1, dec=dec1, pmra=pmra1, pmdec=pmde1, vrad=vrad1, plx=parallax1, /lsr

wabs = abs(warr)

!p.multi=[0,1,2]

plot, uarr, warr, ps=6, yr=[-100,100]
plot, uarr, varr, ps=6, yr=[-100,100]

vtot = sqrt(uarr^2 + varr^2 + warr^2)
dxstop
!p.multi=0

plot, wabs, vsini, ps=6


outstr = create_struct('U',uarr,'V',varr,'W',warr,'VSINI',vsini,'TEFF',teff)

help, outstr
dxstop
mwrfits, outstr, 'spacevelocity.fits', /create
print, "File written"
dxstop

; Bin by vtot
vtvec = [0,20,40,60,1000]
bin_idx = value_locate(vtvec, vtot)

b0idx = where(bin_idx eq 0)
b1idx = where(bin_idx eq 1)
b2idx = where(bin_idx eq 2)
b3idx = where(bin_idx eq 3)

; binned vtot
vt0 = vtot[b0idx]
vt1 = vtot[b1idx]
vt2 = vtot[b2idx]
vt3 = vtot[b3idx]
; binned vsini
vs0 = vsini[b0idx]
nvs0 = n_elements(vs0) 
vs1 = vsini[b1idx]
nvs1 = n_elements(vs1) 
vs2 = vsini[b2idx]
nvs2 = n_elements(vs2) 
vs3 = vsini[b3idx]
nvs3 = n_elements(vs3) 


rotidx0 = where(vs0 gt 5, nr0)
rotidx1 = where(vs1 gt 5, nr1)
rotidx2 = where(vs2 gt 5, nr2)
rotidx3 = where(vs3 gt 5, nr3)

rfrac0 = double(nr0)/double(nvs0)
rfrac1 = double(nr1)/double(nvs1)
rfrac2 = double(nr2)/double(nvs2)
rfrac3 = double(nr3)/double(nvs3)

lim0 = binomial_stat_errors(nvs0, nr0, 68d0)
lim1 = binomial_stat_errors(nvs1, nr1, 68d0)
lim2 = binomial_stat_errors(nvs2, nr2, 68d0)
lim3 = binomial_stat_errors(nvs3, nr3, 68d0)

 set_plot, 'ps'
 device, filename = 'galactic_vel_test.eps'
 device, /color, bits=8
 device, xs = 13, ys= 8, /inches
 loadct, 13


plot, minmax(vtot), [0,1], /nodata, xtit = "Total Velocity w.r.t. LSR (km/s)", $
  ytit = "Fraction of rapid rotators in bin", charsize = 1.5
oplot, [mean(vt0)], [rfrac0], ps=6
oploterror, [mean(vt0)], [rfrac0], [mean(vt0)-min(vt0)], [rfrac0-lim0[0]], /lobar, ps=6
oploterror, [mean(vt0)], [rfrac0], [max(vt0)-mean(vt0)], [lim0[1]-rfrac0], /hibar, ps=6
oplot, [mean(vt1)], [rfrac1], ps=6
oploterror, [mean(vt1)], [rfrac1], [mean(vt1)-min(vt1)], [rfrac1-lim1[0]], /lobar, ps=6
oploterror, [mean(vt1)], [rfrac1], [max(vt1)-mean(vt1)], [lim1[1]-rfrac1], /hibar, ps=6
oplot, [mean(vt2)], [rfrac2], ps=6
oploterror, [mean(vt2)], [rfrac2], [mean(vt2)-min(vt2)], [rfrac2-lim2[0]], /lobar, ps=6
oploterror, [mean(vt2)], [rfrac2], [max(vt2)-mean(vt2)], [lim2[1]-rfrac2], /hibar, ps=6
oplot, [mean(vt3)], [rfrac3], ps=6
oploterror, [mean(vt3)], [rfrac3], [mean(vt3)-min(vt3)], [rfrac3-lim3[0]], /lobar, ps=6
oploterror, [mean(vt3)], [rfrac3], [max(vt3)-mean(vt3)], [lim3[1]-rfrac3], /hibar, ps=6

xyouts, 1.1*[mean(vt0)], 1.1*[rfrac0], strtrim(nvs0,2)
xyouts, 1.1*[mean(vt1)], 1.1*[rfrac1], strtrim(nvs1,2)
xyouts, 1.1*[mean(vt2)], 1.1*[rfrac2], strtrim(nvs2,2)
xyouts, 1.1*[mean(vt3)], 1.1*[rfrac3], strtrim(nvs3,2)
device, /close_file

dxstop




; Bin by wabs
wvec = [0,10,20,30,1000]
bin_idx = value_locate(wvec, wabs)

b0idx = where(bin_idx eq 0)
b1idx = where(bin_idx eq 1)
b2idx = where(bin_idx eq 2)
b3idx = where(bin_idx eq 3)

; binned vtot
vt0 = wabs[b0idx]
vt1 = wabs[b1idx]
vt2 = wabs[b2idx]
vt3 = wabs[b3idx]
; binned vsini
vs0 = vsini[b0idx]
nvs0 = n_elements(vs0) 
vs1 = vsini[b1idx]
nvs1 = n_elements(vs1) 
vs2 = vsini[b2idx]
nvs2 = n_elements(vs2) 
vs3 = vsini[b3idx]
nvs3 = n_elements(vs3) 


rotidx0 = where(vs0 gt 5, nr0)
rotidx1 = where(vs1 gt 5, nr1)
rotidx2 = where(vs2 gt 5, nr2)
rotidx3 = where(vs3 gt 5, nr3)

rfrac0 = double(nr0)/double(nvs0)
rfrac1 = double(nr1)/double(nvs1)
rfrac2 = double(nr2)/double(nvs2)
rfrac3 = double(nr3)/double(nvs3)

lim0 = binomial_stat_errors(nvs0, nr0, 68d0)
lim1 = binomial_stat_errors(nvs1, nr1, 68d0)
lim2 = binomial_stat_errors(nvs2, nr2, 68d0)
lim3 = binomial_stat_errors(nvs3, nr3, 68d0)

 set_plot, 'ps'
 device, filename = 'galactic_wvel_test.eps'
 device, /color, bits=8
 device, xs = 13, ys= 8, /inches
 loadct, 13


plot, minmax(wabs), [0,1], /nodata, xtit = "Magnitude of W w.r.t. LSR (km/s)", $
  ytit = "Fraction of rapid rotators in bin", charsize = 1.5
oplot, [mean(vt0)], [rfrac0], ps=6
oploterror, [mean(vt0)], [rfrac0], [mean(vt0)-min(vt0)], [rfrac0-lim0[0]], /lobar, ps=6
oploterror, [mean(vt0)], [rfrac0], [max(vt0)-mean(vt0)], [lim0[1]-rfrac0], /hibar, ps=6
oplot, [mean(vt1)], [rfrac1], ps=6
oploterror, [mean(vt1)], [rfrac1], [mean(vt1)-min(vt1)], [rfrac1-lim1[0]], /lobar, ps=6
oploterror, [mean(vt1)], [rfrac1], [max(vt1)-mean(vt1)], [lim1[1]-rfrac1], /hibar, ps=6
oplot, [mean(vt2)], [rfrac2], ps=6
oploterror, [mean(vt2)], [rfrac2], [mean(vt2)-min(vt2)], [rfrac2-lim2[0]], /lobar, ps=6
oploterror, [mean(vt2)], [rfrac2], [max(vt2)-mean(vt2)], [lim2[1]-rfrac2], /hibar, ps=6
oplot, [mean(vt3)], [rfrac3], ps=6
oploterror, [mean(vt3)], [rfrac3], [mean(vt3)-min(vt3)], [rfrac3-lim3[0]], /lobar, ps=6
oploterror, [mean(vt3)], [rfrac3], [max(vt3)-mean(vt3)], [lim3[1]-rfrac3], /hibar, ps=6

xyouts, 1.1*[mean(vt0)], 1.1*[rfrac0], strtrim(nvs0,2)
xyouts, 1.1*[mean(vt1)], 1.1*[rfrac1], strtrim(nvs1,2)
xyouts, 1.1*[mean(vt2)], 1.1*[rfrac2], strtrim(nvs2,2)
xyouts, 1.1*[mean(vt3)], 1.1*[rfrac3], strtrim(nvs3,2)
device, /close_file


dxstop

;;; Plot V vs. U
 set_plot, 'ps'
 device, filename = 'galactic_vu_test.eps'
 device, /color, bits=8
 device, xs = 13, ys= 8, /inches
 loadct, 13

plot, uarr, varr, ps=6, xtit="U (km/s)", ytit="V (km/s)"

device, /close_file

;;; Plot Vtot vs Vsini
 set_plot, 'ps'
 device, filename = 'galactic_vtot_vsini_test.eps'
 device, /color, bits=8
 device, xs = 13, ys= 8, /inches
 loadct, 13

plot, vsini, vtot, ps=6, xtit='Vsini (km/s)',ytit='Vtot (km/s)'

device, /close_file

;;; Plot each component vs Vsini
 set_plot, 'ps'
 device, filename = 'galactic_vels_vsini_test.eps'
 device, /color, bits=8
 device, xs = 13, ys= 8, /inches
 loadct, 13

!p.multi = [0,1,3]
plot, vsini, uarr, ps=6, xtit='Vsini (km/s)', ytit='U (km/s)'

plot, vsini, varr, ps=6, xtit='Vsini (km/s)', ytit='V (km/s)'

plot, vsini, warr, ps=6, xtit='Vsini (km/s)', ytit='W (km/s)'

device, /close_file

print, "Plots written"
dxstop

; Characteristic velocity dispersions
;Bensby
; x_d = 0.94d0
; x_td = 0.06d0
; x_h = 0.0015d0

;Juric 2008 (used by Newton 2016)
x_d = 0.89d0
x_td = 0.106d0
x_h = 0.004d0

udisp_d = 35.
udisp_td = 67.
udisp_h = 160.

vdisp_d = 20.
vdisp_td = 38.
vdisp_h = 90.

wdisp_d = 16.
wdisp_td = 35.
wdisp_h = 90.

vasym_d = -15.
vasym_td = -46.
vasym_h = -220.

k_d = 1d0/(((2d0*!pi)^1.5d0) * udisp_d * vdisp_d * wdisp_d)
k_td = 1d0/(((2d0*!pi)^1.5d0) * udisp_td * vdisp_td * wdisp_td)
k_h = 1d0/(((2d0*!pi)^1.5d0) * udisp_h * vdisp_h * wdisp_h)

f_d = k_d * exp(-1d0*((uarr^2)/(2d0*udisp_d^2)) - $
                (((varr-vasym_d)^2)/(2d0*udisp_d^2)) - $
                ((warr^2)/(2d0*wdisp_d^2)))

f_td = k_td * exp(-1d0*((uarr^2)/(2d0*udisp_td^2)) - $
                (((varr-vasym_td)^2)/(2d0*udisp_td^2)) - $
                ((warr^2)/(2d0*wdisp_td^2)))

f_h = k_h * exp(-1d0*((uarr^2)/(2d0*udisp_h^2)) - $
                (((varr-vasym_h)^2)/(2d0*udisp_h^2)) - $
                ((warr^2)/(2d0*wdisp_h^2)))

td_d = (x_td/x_d)*(f_td/f_d)


 set_plot, 'ps'
 device, filename = 'galactic_thin_thick.eps'
 device, /color, bits=8
; device, xs = 13, ys= 8, /inches
 loadct, 13

!p.multi=0

plot, varr, td_d, /ylog, ps=6, xtit = 'V_LSR (km/s)', ytit = 'P(thick disk)/P(thin disk)'
oplot, [-200,200], [10,10], linest=2
xyouts, [-70], [11], ['Likely Thick Disk Members']
oplot, [-200,200], [0.1,0.1], linest=2
xyouts, [-70], [0.07], ['Likely Thin Disk Members']

device, /close_file

dxstop



 set_plot, 'ps'
 device, filename = 'galactic_thin_thick_vsini.eps'
 device, /color, bits=8
; device, xs = 13, ys= 8, /inches
 loadct, 13

!p.multi=0

plot, vsini, td_d, /ylog, ps=6, xtit = 'Vsini (km/s)', ytit = 'P(thick disk)/P(thin disk)'
oplot, [-200,200], [10,10], linest=2
xyouts, [30], [11], ['Likely Thick Disk Members']
oplot, [-200,200], [0.1,0.1], linest=2
xyouts, [30], [0.07], ['Likely Thin Disk Members']

device, /close_file

dxstop

  set_plot, 'ps'
  device, filename = 'galactic_toomre.eps'
  device, /color, bits=8
  device, xs = 5, ys= 5, /inches
  loadct, 13

theta = dindgen(200)*!pi/200d0

r25 = replicate(25d0, 200)
r50 = replicate(50d0, 200)
r75 = replicate(75d0, 200)
r100 = replicate(100d0, 200)

; bin by prob
d_idx = where(td_d le 0.1, nd)
td_idx = where(td_d ge 10, ntd)
int_idx = where(td_d gt 0.1 and td_d lt 10, nint)



plot, varr, sqrt(uarr^2 + warr^2), ps=6, xtit='V_LSR (km/s)', ytit='(U_LSR^2 + W_LSR^2)^1/2', title='Toomre diagram'
oplot, varr[td_idx], sqrt(uarr[td_idx]^2 + warr[td_idx]^2), ps=6, color=250, thick=4
oplot, varr[int_idx], sqrt(uarr[int_idx]^2 + warr[int_idx]^2), ps=6, color=125, thick=4


;oplot, r25, theta, /polar
oplot, r50, theta, /polar, linest=2
;oplot, r75, theta, /polar
oplot, r100, theta, /polar, linest=2

device, /close_file


  set_plot, 'ps'
  device, filename = 'galactic_colortest.eps'
  device, /color, bits=8
  device, xs = 5, ys= 5, /inches
  loadct, 13

plot, dindgen(1000), dindgen(1000), /nodata
for i = 0,999 do begin
oplot, [i], [i], color = i, ps=6
endfor

device, /close_file



dxstop

; Bin by teff
vsiniy = vsini[d_idx]
teffy = teff[d_idx]

cool_idx = where(teffy lt 3250, ncool)
hot_idx = where(teffy ge 3250, nhot)

vcool = vsiniy[cool_idx]
vhot = vsiniy[hot_idx]

tcool = teffy[cool_idx]
thot = teffy[hot_idx]

fcool_idx = where(vcool ge 5, ncoolrot)
fhot_idx = where(vhot ge 5, nhotrot)

frac_cool = double(ncoolrot)/double(ncool)
frac_hot = double(nhotrot)/double(nhot)

;set_plot, 'x'

;plot, teffy, vsiniy, ps=6, /nodata, yr=[0,1]

;oplot, [mean(tcool)], [frac_cool], ps=6
;oplot, [mean(thot)], [frac_hot], ps=6




lim0 = binomial_stat_errors(ncool, ncoolrot, 68d0)
lim1 = binomial_stat_errors(nhot, nhotrot, 68d0)

 set_plot, 'ps'
 device, filename = 'galactic_young_rotfrac.eps'
 device, /color, bits=8
 device, xs = 13, ys= 8, /inches
 loadct, 13


plot, minmax(teffy), [0,1], /nodata, xtit = "Teff", $
  ytit = "Fraction of rapid rotators in bin", charsize = 1.5, title="Fractions of rotators for likely young stars"
oplot, [mean(tcool)], [frac_cool], ps=6
oploterror, [mean(tcool)], [frac_cool], [mean(tcool)-min(tcool)], [frac_cool-lim0[0]], /lobar, ps=6
oploterror, [mean(tcool)], [frac_cool], [max(tcool)-mean(tcool)], [lim0[1]-frac_cool], /hibar, ps=6
oplot, [mean(thot)], [frac_hot], ps=6
oploterror, [mean(thot)], [frac_hot], [mean(thot)-min(thot)], [frac_hot-lim1[0]], /lobar, ps=6
oploterror, [mean(thot)], [frac_hot], [max(thot)-mean(thot)], [lim1[1]-frac_hot], /hibar, ps=6
; oplot, [mean(vt2)], [rfrac2], ps=6
; oploterror, [mean(vt2)], [rfrac2], [mean(vt2)-min(vt2)], [rfrac2-lim2[0]], /lobar, ps=6
; oploterror, [mean(vt2)], [rfrac2], [max(vt2)-mean(vt2)], [lim2[1]-rfrac2], /hibar, ps=6
; oplot, [mean(vt3)], [rfrac3], ps=6
; oploterror, [mean(vt3)], [rfrac3], [mean(vt3)-min(vt3)], [rfrac3-lim3[0]], /lobar, ps=6
; oploterror, [mean(vt3)], [rfrac3], [max(vt3)-mean(vt3)], [lim3[1]-rfrac3], /hibar, ps=6

xyouts, [mean(tcool)]+100, [frac_cool]+0.05, strtrim(ncool,2)
xyouts, [mean(thot)]+100, [frac_hot]+0.05, strtrim(nhot,2)
;xyouts, 1.1*[mean(vt2)], 1.1*[rfrac2], strtrim(nvs2,2)
;xyouts, 1.1*[mean(vt3)], 1.1*[rfrac3], strtrim(nvs3,2)
device, /close_file










stop

wvec = [0,10,20,30,1000]
bin_idx = value_locate(wvec, wabs)

b0idx = where(bin_idx eq 0)
b1idx = where(bin_idx eq 1)
b2idx = where(bin_idx eq 2)
b3idx = where(bin_idx eq 3)

; binned vtot
vt0 = wabs[b0idx]
vt1 = wabs[b1idx]
vt2 = wabs[b2idx]
vt3 = wabs[b3idx]
; binned vsini
vs0 = vsini[b0idx]
nvs0 = n_elements(vs0) 
vs1 = vsini[b1idx]
nvs1 = n_elements(vs1) 
vs2 = vsini[b2idx]
nvs2 = n_elements(vs2) 
vs3 = vsini[b3idx]
nvs3 = n_elements(vs3) 


rotidx0 = where(vs0 gt 5, nr0)
rotidx1 = where(vs1 gt 5, nr1)
rotidx2 = where(vs2 gt 5, nr2)
rotidx3 = where(vs3 gt 5, nr3)

rfrac0 = double(nr0)/double(nvs0)
rfrac1 = double(nr1)/double(nvs1)
rfrac2 = double(nr2)/double(nvs2)
rfrac3 = double(nr3)/double(nvs3)

lim0 = binomial_stat_errors(nvs0, nr0, 68d0)
lim1 = binomial_stat_errors(nvs1, nr1, 68d0)
lim2 = binomial_stat_errors(nvs2, nr2, 68d0)
lim3 = binomial_stat_errors(nvs3, nr3, 68d0)

 set_plot, 'ps'
 device, filename = 'galactic_wvel_test.eps'
 device, /color, bits=8
 device, xs = 13, ys= 8, /inches
 loadct, 13


plot, minmax(wabs), [0,1], /nodata, xtit = "Magnitude of W w.r.t. LSR (km/s)", $
  ytit = "Fraction of rapid rotators in bin", charsize = 1.5
oplot, [mean(vt0)], [rfrac0], ps=6
oploterror, [mean(vt0)], [rfrac0], [mean(vt0)-min(vt0)], [rfrac0-lim0[0]], /lobar, ps=6
oploterror, [mean(vt0)], [rfrac0], [max(vt0)-mean(vt0)], [lim0[1]-rfrac0], /hibar, ps=6
oplot, [mean(vt1)], [rfrac1], ps=6
oploterror, [mean(vt1)], [rfrac1], [mean(vt1)-min(vt1)], [rfrac1-lim1[0]], /lobar, ps=6
oploterror, [mean(vt1)], [rfrac1], [max(vt1)-mean(vt1)], [lim1[1]-rfrac1], /hibar, ps=6
oplot, [mean(vt2)], [rfrac2], ps=6
oploterror, [mean(vt2)], [rfrac2], [mean(vt2)-min(vt2)], [rfrac2-lim2[0]], /lobar, ps=6
oploterror, [mean(vt2)], [rfrac2], [max(vt2)-mean(vt2)], [lim2[1]-rfrac2], /hibar, ps=6
oplot, [mean(vt3)], [rfrac3], ps=6
oploterror, [mean(vt3)], [rfrac3], [mean(vt3)-min(vt3)], [rfrac3-lim3[0]], /lobar, ps=6
oploterror, [mean(vt3)], [rfrac3], [max(vt3)-mean(vt3)], [lim3[1]-rfrac3], /hibar, ps=6

xyouts, 1.1*[mean(vt0)], 1.1*[rfrac0], strtrim(nvs0,2)
xyouts, 1.1*[mean(vt1)], 1.1*[rfrac1], strtrim(nvs1,2)
xyouts, 1.1*[mean(vt2)], 1.1*[rfrac2], strtrim(nvs2,2)
xyouts, 1.1*[mean(vt3)], 1.1*[rfrac3], strtrim(nvs3,2)
device, /close_file


stop



end

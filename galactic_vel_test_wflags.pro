; See if I can get space velocities for ALL of my stars
; UPDATE to remove stars flagged as bad
pro galactic_vel_test_wflags

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



; Calculate distances from parallaxes and add those estimated from absk_model
plx_idx = where(finite(plx),nfin, complement=nan_idx, ncomplement=nnan)

dist = 1d0/parallax

; Calculate a distance for stars with just absolute magnitude
if nnan ne 0 then begin
    absk_file = 'absk_fitall_model.fits'
    absk_info = mrdfits(absk_file,1)
    absk_measured = absk_info.absk
    absk_model = absk_info.absk_model
    kmag = aspcap_str.k

    dist_modulus_model = 10d0^((kmag-absk_model+5d0)/5d0)
    dist_modulus_measured = 10d0^((kmag-absk_measured+5d0)/5d0)
endif

dist[nan_idx] = dist_modulus_model[nan_idx]
;;;

; Load ASPCAP PM's for stars without Dittman PM's
nopm_idx = where(not finite(pmra_all), nnopm)
pmra = pmra_all
pmde = pmde_all
if nnopm gt 0 then begin
    pmra[nopm_idx] = aspcap_str[nopm_idx].pmra/1d3
    pmde[nopm_idx] = aspcap_str[nopm_idx].pmdec/1d3
endif


;;; Loop through all spectra and get galactic coords, rv
for snum = 0, napo-1 do begin

    ; get head
    head = *head_ptr[snum]
    
    glati = fxpar(head, 'GLAT')
    glati = glati * !const.DtoR
    glat = [glat, glati]

    gloni = fxpar(head, 'GLON')
    gloni = gloni * !const.DtoR
    glon = [glon, gloni]

    rvi = fxpar(head,'VRAD')
    rv = [rv, rvi]

endfor

; all data needed for transformation to galactic velocities
height = dist*sin(glat)
distance1 = dist
parallax1 = parallax
ra1 = ra_all
dec1 = dec_all
pmra_mas = pmra*1d3
pmde_mas = pmde*1d3
vrad1 = rv
; distance1 = dist[fin_idx]
; parallax1 = parallax[fin_idx]
; ra1 = ra_all[fin_idx]
; dec1 = dec_all[fin_idx]
; pmra1 = pmra_all[fin_idx]
; pmde1 = pmde_all[fin_idx]
; vrad1 = rv[fin_idx]
; At this point, pm and plx are in arcsec/year and arcsec
; (Dittman gives plx in mas, and pm in arcsec)

; Data on the stars themselves
teff = teff_vfit
vsini = vsini_vfit

; Calculate UVW values
; if using DISTANCE instead of PLX, pm's must be in mas/year
; if using PLX instead of DISTANCE, pm's and plx must be in same arc units

; This one should have everything in arcsec
;gal_uvw, uarr, varr, warr, ra=ra1, dec=dec1, pmra=pmra, pmdec=pmde,
;vrad=vrad1, plx=parallax1, /lsr

; Here, we have PM in mas/year, and distance in parsecs
gal_uvw, uarr, varr, warr, ra=ra1, dec=dec1, pmra=pmra_mas, pmdec=pmde_mas, vrad=vrad1, distance=distance1, /lsr

; Remove appropriate entries
uarr_all = uarr
varr_all = varr
warr_all = warr


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


;flags =
;mrdfits('/home/stgilhool/APOGEE/APOGEE_data/spectral_flags_sg.fits',0)
flags = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/gold_sample.fits',0)
flag_idx = where(flags eq 1, complement=okay_idx)

;Output DATA
 outstr = create_struct('U',uarr,'V',varr,'W',warr,'Z',height,'PLX',parallax,'DIST',distance1,'GLAT',glat,'GLON',glon,'VSINI',vsini,'TEFF',teff,'THICKTHIN_PROB',td_d, 'GS_FLAGS', flags, absk_info)

 help, outstr
 dxstop
 mwrfits, outstr, 'spacevelocity_wflags.fits', /create
 print, "File written"
 dxstop


; Keep only gold sample entries
uarr = uarr[okay_idx]
varr = varr[okay_idx]
warr = warr[okay_idx]
height = height[okay_idx]
parallax = parallax[okay_idx]
distance1 = distance1[okay_idx]
glat = glat[okay_idx]
glon = glon[okay_idx]
vsini = vsini[okay_idx] 
teff = teff[okay_idx]

wabs = abs(warr)

!p.multi=[0,1,2]

plot, uarr, warr, ps=6, yr=[-100,100]
plot, uarr, varr, ps=6, yr=[-100,100]

vtot = sqrt(uarr^2 + varr^2 + warr^2)
dxstop
!p.multi=0

plot, vtot, vsini, ps=6

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
 device, filename = 'galactic_vel_test_wflag.eps'
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


; Plot abs(height) vs abs(W), color-coded by vsini
COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr
set_plot, 'ps'
;  device, filename = 'galactic_height_zvel_vsini_wflags.eps'

device, /color, decomposed = 0
loadct, 3
loadct, 3, rgb_table=rgb_tab
; Mess with colors for color-coding
gamma_ct, 0.3, /current

rgb_tab = [[R_curr], [G_curr], [B_curr]]


minvsini = 0
maxvsini = max(vsini)
binsize = 5
vsini_hist = histogram(vsini, binsize=binsize, min=minvsini, max=maxvsini,$
                       reverse_indices=ri)
nbins = n_elements(vsini_hist) 
null = where(vsini_hist gt 0, nbins_filled)
vsini_node_vec = lindgen(nbins+1)*binsize + minvsini
vsini_ctr_vec = vsini_node_vec[0:-2]+(binsize/2.)

color_vec_orig = lindgen(nbins_filled+2)*(256/(nbins_filled+2))

color_vec = byte(color_vec_orig[1:-2])
;color_vec = byte(color_vec_orig)

vdist = abs(height)
vspeed = abs(warr)

p1 = plot(vspeed, vdist, /ylog,/xlog, xtit = 'Vertical Speed', $
          ytit = 'Distance from Galactic midplane',$
          title = "Vertical Position vs. Vertical Speed, color-coded by Vsini", $
          /nodata, /buffer)

color_idx = 0
for bin=0, nbins-1 do begin
    ; Recover the data from the bin
    if ri[bin] ne ri[bin+1] then begin
        bin_idx = ri[ri[bin]:ri[bin+1]-1]
        vdist_bin = vdist[bin_idx]
        vspeed_bin = vspeed[bin_idx]
        color_i = color_vec[color_idx]
        color_i_vec = replicate(color_i, n_elements(vspeed_bin)) 
        
        if bin eq 0 then begin
            p1 = plot(vspeed_bin, vdist_bin, rgb_table = rgb_tab, symbol='plus', $
                      sym_size = 0.5, linestyle=6, $
                      vert_colors=color_i_vec,/overplot, $
                      position=[0.2, 0.16, 0.89, 0.91], /nodata)
        endif else begin
            p1 = plot(vspeed_bin, vdist_bin, rgb_table = rgb_tab, symbol='square', $
                      sym_size = 0.5, linestyle=6, $
                      vert_colors=color_i_vec,/overplot, $
                      position=[0.2, 0.2, 0.89, 0.91])
        endelse
        
        color_idx++
    endif
endfor
    

cb = colorbar(target=p1, orientation=0, position=[0.25, 0.08, 0.76, 0.1], $
                            text_orientation=90, $
              tickname = ['ND','5-10','10-15','15-20','20-25','25-30','35-40','40-45',$
                         '45-50','50-55','55-60'])

p1.save, "galactic_height_zvel_vsini_nonds_wflags.eps"

device, /close_file

print, "Check out the file"
stop
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
 device, filename = 'galactic_wvel_test_wflags.eps'
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
 device, filename = 'galactic_vu_test_wflags.eps'
 device, /color, bits=8
 device, xs = 13, ys= 8, /inches
 loadct, 13

plot, uarr, varr, ps=6, xtit="U (km/s)", ytit="V (km/s)"

device, /close_file

;;; Plot Vtot vs Vsini
 set_plot, 'ps'
 device, filename = 'galactic_vtot_vsini_test_wflags.eps'
 device, /color, bits=8
 device, xs = 13, ys= 8, /inches
 loadct, 13

plot, vsini, vtot, ps=6, xtit='Vsini (km/s)',ytit='Vtot (km/s)'

device, /close_file


;;; Plot each component vs Vsini
 set_plot, 'ps'
 device, filename = 'galactic_vels_vsini_test_wflags.eps'
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
 device, filename = 'galactic_thin_thick_wflags.eps'
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

; OUTPUT DATA

; outstr = create_struct('U',uarr,'V',varr,'W',warr,'Z',height,'PLX',parallax,'DIST',distance1,'GLAT',glat,'GLON',glon,'VSINI',vsini,'TEFF',teff,'THICKTHIN_PROB',td_d)

; help, outstr
; dxstop
; mwrfits, outstr, 'spacevelocity_wflags.fits', /create
; print, "File written"
; dxstop




 set_plot, 'ps'
 device, filename = 'galactic_thin_thick_vsini_wflags.eps'
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
  device, filename = 'galactic_toomre_wflags.eps'
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

;;; Plot V_LSR vs Vsini, color coded by membership
 set_plot, 'ps'
 device, filename = 'galactic_vlsr_vsini_test_wflags.eps'
 device, /color, bits=8
 device, xs = 13, ys= 8, /inches
 loadct, 13

plot, vsini, varr, ps=6, xtit='Vsini (km/s)',ytit='V_LSR (km/s)'
oplot, vsini[int_idx], varr[int_idx], ps=6, color=125, thick=4
oplot, vsini[td_idx], varr[td_idx], ps=6, color=250, thick=4

device, /close_file

;;; Plot Wabs_LSR vs Vsini, color-coded by membership
 set_plot, 'ps'
 device, filename = 'galactic_wvel_vsini_test_wflags.eps'
 device, /color, bits=8
 device, xs = 13, ys= 8, /inches
 loadct, 13

plot, vsini, wabs, ps=6, xtit='Vsini (km/s)',ytit='Magnitude of W_LSR (km/s)'
oplot, vsini[int_idx], wabs[int_idx], ps=6, color=125, thick=4
oplot, vsini[td_idx], wabs[td_idx], ps=6, color=250, thick=4

device, /close_file


;;; Plot W_LSR vs Vsini, color-coded by membership
 set_plot, 'ps'
 device, filename = 'galactic_wlsr_vsini_test_wflags.eps'
 device, /color, bits=8
 device, xs = 13, ys= 8, /inches
 loadct, 13

plot, vsini, warr, ps=6, xtit='Vsini (km/s)',ytit='W_LSR (km/s)'
oplot, vsini[int_idx], warr[int_idx], ps=6, color=125, thick=4
oplot, vsini[td_idx], warr[td_idx], ps=6, color=250, thick=4

device, /close_file

;;; Plot RMS(W) in 10km/s bins of Vsini
set_plot, 'ps'
 device, filename = 'galactic_rmsw_vsini_test_wflags.eps'
 device, /color, bits=8
 device, xs = 13, ys= 8, /inches
 loadct, 13

binsize = 10
binmin = 0
w_hist = histogram(vsini, binsize=binsize, min=binmin, reverse_indices=ri)
nbins = n_elements(w_hist) 
bin0 = lindgen(nbins)*binsize+binmin
binmid = lindgen(nbins)*binsize+binmin+(binsize/2)
binfull = lindgen(nbins+1)*binsize+binmin

rms_vec = []
std_vec = []

for bin = 0, nbins-1 do begin

    ; If the bin has data
    if ri[bin] ne ri[bin+1] then begin
        bin_idx = ri[ri[bin]:ri[bin+1]-1]
        rms_bin = sqrt(mean(warr[bin_idx]^2, /double))
        ;std_bin = stddev(warr[bin_idx])
        std_bin = rms_error(rms_bin, n_elements(bin_idx))
        
        if std_bin eq -1 then message, "rms error error"

        rms_vec = [rms_vec, rms_bin]
        std_vec = [std_vec, std_bin]
    endif else begin
        rms_vec = [rms_vec, 0]
        std_vec = [std_vec, 0]
    endelse

endfor

plot, binmid, rms_vec, ps=6, xr = minmax(binfull), xtitle = 'Vsini (km/s)', $
  ytitle = 'RMS of W_LSR'
oploterror, binmid, rms_vec, std_vec, ps=6
xyouts, binmid+2, rms_vec, strtrim(w_hist,2)

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




lim0 = binomial_stat_errors(ncool, ncoolrot, 90d0)
;lim1 = binomial_stat_errors(nhot, nhotrot, 90d0)
lim1 = [double(nhotrot)/nhot*0.8, double(nhotrot)/nhot*1.2]

 set_plot, 'ps'
 device, filename = 'galactic_young_rotfrac_wflags.eps'
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



end

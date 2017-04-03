; See if I can get galactic height for some of my stars
pro galactic_height_test

; Globals
skip_load=1
napo = 1350

; Master Info Table
mi_file = '/home/stgilhool/APOGEE/master_info6.fits'
mi_str = mrdfits(mi_file, 1)

; pmra, pmdec
pmra = mi_str.pmra
pmde = mi_str.pmde
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

save, /variables, filename = 'galactic_height_test.sav'

skip_load_mark:
restore, 'galactic_height_test.sav'

head_ptr = apg_files.header

; spec 17 has plx
test = *head_ptr[17]
; galactic lat, vrad
glat = fxpar(test, 'GLAT')
vrad = fxpar(test, 'VRAD')

height = []
glat = []
dist = []

for snum = 0, napo-1 do begin

    plxi = plx[snum]/1d3
    disti = 1d0/plxi
    dist = [dist, disti]

    ; get head
    head = *head_ptr[snum]
    
    glati = fxpar(head, 'GLAT')
    glati = glati * !const.DtoR
    glat = [glat, glati]
    
    if finite(plxi) and finite(glati) then begin
        heighti = disti*sin(glati)
    endif else heighti = !values.d_nan

    height = [height, heighti]

endfor

plot, abs(height), xtit = 'Spec Num', ytit = 'Gal Plane Height (pc)', ps=6
    
height = abs(height)

fin_idx = where(finite(height))

h = height[fin_idx]
teff = teff_vfit[fin_idx]
vsini = vsini_vfit[fin_idx]


; Bin by height
hvec = [0,5]

ri_vec = value_locate(hvec, h)

loidx = where(ri_vec eq 0, nlow)
hiidx = where(ri_vec eq 1, nhigh)

help, nlow
help, nhigh

hlo = h[loidx]
hhi = h[hiidx]

tlo = teff[loidx]
thi = teff[hiidx]

vlo = vsini[loidx]
vhi = vsini[hiidx]

loloidx = where(tlo lt 3250, nlolo)
lohiidx = where(tlo ge 3250, nlohi)

hiloidx = where(thi lt 3250, nhilo)
hihiidx = where(thi ge 3250, nhihi)

; teff and vsini for each part of each plot
tlolo = tlo[loloidx]
vlolo = vlo[loloidx]

tlohi = tlo[lohiidx]
vlohi = vlo[lohiidx]

thilo = thi[hiloidx]
vhilo = vhi[hiloidx]

thihi = thi[hihiidx]
vhihi = vhi[hihiidx]

; compute fractions of detections
llvidx = where(vlolo gt 5, nll)
lhvidx = where(vlohi gt 5, nlh)
hlvidx = where(vhilo gt 5, nhl)
hhvidx = where(vhihi gt 5, nhh)

llfrac = double(nll)/double(nlolo)
lllim = binomial_stat_errors(nlolo, nll, 68d0)

lhfrac = double(nlh)/double(nlohi)
lhlim = binomial_stat_errors(nlohi, nlh, 68d0)

hlfrac = double(nhl)/double(nhilo)
hllim = binomial_stat_errors(nhilo, nhl, 68d0)

hhfrac = double(nhh)/double(nhihi)
hhlim = binomial_stat_errors(nhihi, nhh, 68d0)


set_plot, 'ps'
device, filename = 'galactic_height_test.eps'
device, /color, bits=8
device, xs = 13, ys= 8, /inches
loadct, 13

!p.multi = [0,2,1]

tvec = lindgen(15)*100 + 2500

plot, tvec, [0,1], /nodata, xtit='Teff', ytit='Fraction of fast rotators', charsize=1.5, title = 'Galactic Height < 5 pc (Young?)'
oplot, [mean(tlolo)], [llfrac], ps=6
oploterror, [mean(tlolo)], [llfrac], [mean(tlolo)-min(tlolo)], [llfrac]-lllim[0], /lobar, ps=6
oploterror, [mean(tlolo)], [llfrac], [max(tlolo)-mean(tlolo)], lllim[1]-[llfrac], /hibar, ps=6
xyouts, [2800], [0.70], strtrim(nlolo,2)
xyouts, [3500], [0.1], strtrim(nlohi,2)


oplot, [mean(tlohi)], [lhfrac], ps=6
oploterror, [mean(tlohi)], [lhfrac], [mean(tlohi)-min(tlohi)], [lhfrac]-lhlim[0], /lobar, ps=6
oploterror, [mean(tlohi)], [lhfrac], [max(tlohi)-mean(tlohi)], lhlim[1]-[lhfrac], /hibar, ps=6

plot, tvec, [0,1], /nodata, xtit='Teff', ytit='Fraction of fast rotators', charsize=1.5, title = 'Galactic Height > 5 pc (Old?)'
oplot, [mean(thilo)], [hlfrac], ps=6
oploterror, [mean(thilo)], [hlfrac], [mean(thilo)-min(thilo)], [hlfrac]-hllim[0], /lobar, ps=6
oploterror, [mean(thilo)], [hlfrac], [max(thilo)-mean(thilo)], hllim[1]-[hlfrac], /hibar, ps=6

oplot, [mean(thihi)], [hhfrac], ps=6
oploterror, [mean(thihi)], [hhfrac], [mean(thihi)-min(thihi)], [hhfrac]-hhlim[0], /lobar, ps=6
oploterror, [mean(thihi)], [hhfrac], [max(thihi)-mean(thihi)], hhlim[1]-[hhfrac], /hibar, ps=6
xyouts, [2800], [0.8], strtrim(nhilo,2)
xyouts, [3600], [0.25], strtrim(nhihi,2)

device, /close_file


stop

outstr = {HEIGHT:h, TEFF:teff, VSINI:vsini}

mwrfits, outstr, 'gal_height.fits',/create


stop


end

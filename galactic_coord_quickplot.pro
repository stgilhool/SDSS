; Quickly plot the galactic coords of all stars
pro galactic_coord_quickplot

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

if not skip_load then message, "must edit save/restore code before proceeding"
;save, /variables, filename = 'galactic_vel_test.sav'

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
    ;glati = glati * !const.DtoR
    glat = [glat, glati]

    gloni = fxpar(head, 'GLON')
    ;gloni = gloni * !const.DtoR
    glon = [glon, gloni]

    rvi = fxpar(head,'VRAD')
    rv = [rv, rvi]

    if finite(plxi) and finite(glati) then begin
        heighti = disti*sin(glati)
    endif else heighti = !values.d_nan

    height = [height, heighti]

endfor


plot, glon, glat, ps=6, xr = [0,360], yr=[-180,180], /xs,/ys, $
  xtit='Galactic Longitude (deg)', ytit='Galactic Latitude (deg)'

plane_idx = where(abs(glat) le 10, nplane)
print, nplane
print, n_elements(glat) 

stop


fin_idx = where(finite(height), nfin)

glat_fin = glat[fin_idx]
glon_fin = glon[fin_idx]

plot, glon_fin, glat_fin, ps=6, xr = [0,360], yr=[-180,180], /xs,/ys, $
  xtit='Galactic Longitude (deg)', ytit='Galactic Latitude (deg)'

plane_idx = where(abs(glat_fin) le 10, nplane)
print, nplane
print, n_elements(glat_fin) 


stop


end

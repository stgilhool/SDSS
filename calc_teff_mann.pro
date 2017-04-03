; Calculate Teff based on Mann 2015 empirical relations
pro calc_teff_mann

; Load spectral info
infofile = '/home/stgilhool/APOGEE/master_info5.fits'
info = mrdfits(infofile,1)

; These are the coeffs of the polynomial fits
; 1 refers to 
; Teff/3500 = a + bx + cx^2 + dx^3 + ex^4 + f[Fe/H]
; 2 refers to 
; Teff/3500 = a + bx + cx^2 + dx^3 + ex^4 + fz + gz^2
; where x = (r-J), and z = (J-H)
a1 = 2.532d0
b1 = -1.319d0
c1 = 0.4449d0
d1 = -0.07151d0
e1 = 0.004333d0
f1 = 0.05629d0

a2 = 2.151d0
b2 = -1.092d0
c2 = 0.3767d0
d2 = -0.06292d0
e2 = 0.003950d0
f2 = 0.1697d0
g2 = 0.03106d0

rj_feh_coeff = [a1,b1,c1,d1,e1]
feh_coeff    = [f1]
rj_jh_coeff  = [a2,b2,c2,d2,e2]
jh_coeff     = [f2,g2]

; Get full r-J, J-H and [Fe/H], and M_K
r_j = info.r_mag - info.jmag
j_h = info.jmag - info.hmag
feh = info.feh_irtf

kmag = info.kmag
plx  = (info.plx_dtm)/1d3
eplx = (info.eplx_dtm)/1d3
absk = kmag + 5d0 + 5d0*alog10(plx)
eabsk = 5d0*0.434*(eplx/plx)

; Relation 1
; only using irtf feh's
rel1_idx = where(finite(feh) and finite(r_j), rel1cnt)
if rel1cnt eq 0 then message, "No spectra suitable for Relation 1"

r_j1 = r_j[rel1_idx]
feh1 = feh[rel1_idx]

r_jpoly = poly(r_j1, rj_feh_coeff)
fehpoly = poly(feh1, feh_coeff)

teff1_3500 = r_jpoly + fehpoly

teff1 = 3500d0 * teff1_3500
absk1 = absk[rel1_idx]

; compare with other Teffs
teff1_irtf = info[rel1_idx].teff_irtf
teff1_apg  = info[rel1_idx].teff_apg
teff1_vfit = info[rel1_idx].teff_vfit

; Relation 2

rel2_idx = where(finite(r_j), rel2cnt)
if rel2cnt eq 0 then message, "No spectra suitable for Relation 2"

r_j2 = r_j[rel2_idx]
j_h2 = j_h[rel2_idx]

r_jpoly = poly(r_j2, rj_jh_coeff)
j_hpoly = poly(j_h2, jh_coeff)

teff2_3500 = r_jpoly + j_hpoly

teff2 = 3500d0 * teff2_3500
absk2 = absk[rel2_idx]

; compare with other Teffs
teff2_irtf = info[rel2_idx].teff_irtf
teff2_apg  = info[rel2_idx].teff_apg
teff2_vfit = info[rel2_idx].teff_vfit


; Plot
window, 0, xs = 1500, ys = 800
plot, teff1, absk1, ps = 6, xtit = 'Teff', ytit = 'M_K', yr=[12,3],xr =[4500,2500], $
  title = "M_K vs. Teff for r-J,[Fe/H] relation", charsize = 1.5
dxstop
window, 1, xs = 1500, ys = 800
plot, teff2, absk2, ps = 6, xtit = 'Teff', ytit = 'M_K', yr=[12,3],xr =[4500,2500], $
  title = "M_K vs. Teff for r-J,J-H relation", charsize = 1.5
dxstop 
window, 2, xs=1500, ys=800
t1sort = sort(teff1)
t2sort = sort(teff2)

teff1 = teff1[t1sort]
teff1_irtf = teff1_irtf[t1sort]
teff1_vfit = teff1_vfit[t1sort]
teff1_apg = teff1_apg[t1sort]

teff2 = teff2[t2sort]
teff2_irtf = teff2_irtf[t2sort]
teff2_vfit = teff2_vfit[t2sort]
teff2_apg = teff2_apg[t2sort]

plot, teff1, ps=6, yr=[2000,5000]
oplot, teff1_irtf, ps=1
dxstop
plot, teff1, ps=6, yr=[2000,5000]
oplot, teff1_apg, ps=2
dxstop
plot, teff1, ps=6, yr=[2000,5000]
oplot, teff1_vfit, ps=4
dxstop
plot, teff2, ps=6, yr=[2000,5000]
oplot, teff2_irtf, ps=1
dxstop
plot, teff2, ps=6, yr=[2000,5000]
oplot, teff2_apg, ps=2
dxstop
plot, teff2, ps=6, yr=[2000,5000]
oplot, teff2_vfit, ps=4

stop

end

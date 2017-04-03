; Readin vsini from literature
; and compare with DR13 aspcap values
pro vsini_lit_all_aspcap_src

skip = 1
if skip then goto, reiners
;;; Paths
apg_path = '/home/stgilhool/APOGEE/vsini_results/'
apg_vsini_file = apg_path+'initial_run/rfile_master.fits'

aspcap_file = '/home/stgilhool/APOGEE/APOGEE_data/allparams_dr13.fits'

;;; APG info
napgspec = 1350L

; Readin APG headers
apg_str = readin_apo(nfiles = napgspec)
; Initialize array to store APG 2mass id's
apg_2mid = strarr(napgspec)
apg_2mid2 = strarr(napgspec)
; Readin my results
apg_vsini_str = mrdfits(apg_vsini_file,1)

; Readin aspcap results
aspcap_str = mrdfits(aspcap_file,1)



;;; Find matches between published tables and my targets

; Info for Newton's vsini's
nwt_id = ['J02170993+3526330', $
          'J02204625+0258375', $
          'J03205965+1854233', $
          'J03425325+2326495', $
          'J06000351+0242236', $
          'J07444018+0333089', $
          'J08294949+2646348', $
          'J08505062+5253462', $
          'J09301445+2630250', $
          'J09535523+2056460', $
          'J10163470+2751497', $
          'J10521423+0555098', $
          'J11474074+0015201', $
          'J12185939+1107338', $
          'J13003350+0541081', $
          'J13564148+4342587', $
          'J17195422+2630030', $
          'J18024624+3731048', $
          'J18452147+0711584', $
          'J19173151+2833147', $
          'J22245593+5200190', $
          'J03132299+0446293', $
          'J05011802+2237015', $
          'J06022918+4951561', $
          'J09002359+2150054', $
          'J11005043+1204108', $
          'J12265737+2700536', $
          'J16370146+3535456', $
          'J22081254+1036420', $
          'J23134727+2117294', $
          'J23354132+0611205' $
          ]

nwt_vsini = [ 28.2, $
              23.3, $
              8.0, $
              12.7, $
              5.8, $
              4.5, $
              8.1, $
              13.1, $
              6.7, $
              16.5, $
              4.0, $
              19.1, $
              5.6, $
              9.2, $
              16.8, $
              14.0, $
              1.9, $
              4.5, $
              16.1, $
              13.2, $
              4.5, $
              2.2, $
              8.8, $
              4.3, $
              20.0, $
              26.5, $
              13.5, $
              7.0, $
              18.6, $
              16.0, $
              9.8 $
              ]

nwt_upperlimit = [0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  1, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  1, $
                  0, $
                  0, $
                  1, $
                  1, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0 $
                  ]

nwt_vsini_err = [0.7, $
                 0.7, $
                 !values.f_nan, $
                 0.5, $
                 0.3, $
                 !values.f_nan, $
                 1.1, $
                 0.7, $
                 1.5, $
                 0.4, $
                 !values.f_nan, $
                 0.2, $
                 1.4, $
                 1.9, $
                 2.1, $
                 2.0, $
                 0.3, $
                 !values.f_nan, $
                 0.1, $
                 0.5, $
                 !values.f_nan, $
                 !values.f_nan, $
                 0.3, $
                 1.2, $
                 0.6, $
                 0.8, $
                 0.6, $
                 1.8, $
                 2.0, $
                 4.0, $
                 1.1 $
                 ]
                 
nnwt = n_elements(nwt_id) 


nwtapg_str = {id:nwt_id[0], $
              nwt_idx:0L, $
              nwt_vsini:0d0, $
              nwt_vsini_ul:0d0, $
              nwt_vsini_ll:0d0, $
              apg_idx:0L, $
              apg_vsini:0d0, $
              aspcap_vsini:0d0 $
             }
nwtapg = replicate(nwtapg_str, nnwt)
nmatch = 0

; initialize
nwt_aspcap_idx = []

for i = 0L, napgspec - 1 do begin
    ; Read in APG header for a given file
    head = *(apg_str[i].header)
    apgid = fxpar(head, 'OBJID')

    ; Get aspcap id
    aspcap_id = aspcap_str[i].apogee_id
    ; make sure they match
    if apgid ne aspcap_id then message, "ID's don't match"

    ; Replace the '2M' with 'J' in the OBJID field
    id = mg_streplace(apgid, '2M', 'J')
    id = strtrim(id,2)
    ; Store the APG id
    apg_2mid[i] = id
    apg_2mid2[i] = strtrim(apgid,2)
    ; Check if there's a match in the Newton data
    nwt_idx = where(nwt_id eq id, match)
    if match eq 1 then begin
        
        ; Save aspcap index
        nwt_aspcap_idx = [nwt_aspcap_idx, long(i)]

        ; Retrieve my result
        apg_vsini = apg_vsini_str[i].vsini_best

        ; Retrieve aspcap result
        aspcap_vsini = aspcap_str[i].vsini

        ; Store nwt_idx, apg_idx, vsini's and id
        nwtapg[nmatch].id = id
        nwtapg[nmatch].nwt_idx = nwt_idx[0]
        nwtapg[nmatch].apg_idx = long(i)
        nwtapg[nmatch].nwt_vsini = nwt_vsini[nwt_idx[0]]
        nwtapg[nmatch].apg_vsini = apg_vsini
        nwtapg[nmatch].aspcap_vsini = aspcap_vsini
        
        ; Deal with error bars of Newton data
        if nwt_upperlimit[nwt_idx[0]] eq 1 then begin
            nwtapg[nmatch].nwt_vsini_ul = nwt_vsini[nwt_idx[0]]
            nwtapg[nmatch].nwt_vsini_ll = 0d0
        endif else if nwt_vsini_err[nwt_idx[0]] eq !values.f_nan then begin
            nwtapg[nmatch].nwt_vsini_ul = nwt_vsini[nwt_idx[0]]
            nwtapg[nmatch].nwt_vsini_ll = nwt_vsini[nwt_idx[0]]
        endif else begin
            nwtapg[nmatch].nwt_vsini_ul = nwt_vsini[nwt_idx[0]] + $
              nwt_vsini_err[nwt_idx[0]]

            nwtapg[nmatch].nwt_vsini_ll = nwt_vsini[nwt_idx[0]] - $
                          nwt_vsini_err[nwt_idx[0]]
        endelse
        

        nmatch++
    endif else if match gt 1 then message, "Too many matches" 
        
endfor

nwtapg = nwtapg[0:nmatch-1]

fx = lindgen(nmatch)
window, 0, xs = 1300, ys = 800
plot, fx, nwtapg.aspcap_vsini, ps = 6, $
  title = "ASPCAP vsini's (boxes) vs. Newton's vsini's (diamonds w/ error bars)", $
  charsize=2, xtitle = "File", ytitle = "vsini (km/s)", xs = 2
oplot, fx, nwtapg.nwt_vsini, ps = 4
errplot, fx, nwtapg.nwt_vsini_ll, nwtapg.nwt_vsini_ul

;save, /variables, filename= 'vsini_lit.sav'
save, /variables, filename= 'vsini_lit_aspcap.sav'
stop

reiners:
restore, 'vsini_lit_aspcap.sav'

;;; Match with Reiner's data

;;; Reiner's info
rnr_path = '/home/stgilhool/APOGEE/vsini_literature/reiners/'
rnr_file = rnr_path + 'vsini_reiners.fits'

rnr_str = mrdfits(rnr_file,1)
rnr_ra = rnr_str.raj2000
rnr_de = rnr_str.dej2000
rnr_ll = rnr_str.l_vsini
rnr_vsini = rnr_str.vsini


apg_ra_arr = dblarr(napgspec)
apg_de_arr = dblarr(napgspec)

for i = 0, napgspec-1 do begin
    ;Get APG ra and dec
    head = *(apg_str[i].header)
    apg_ra = fxpar(head, 'RA')
    apg_de = fxpar(head, 'DEC')
    apg_ra_arr[i] = apg_ra
    apg_de_arr[i] = apg_de
endfor
    
;See if there's a match in Reiner's data
srcor, rnr_ra, rnr_de, apg_ra_arr, apg_de_arr, 4000d0, out1, out2, /option, /spherical
    
; Get overlapping Reiner's info
r_ra = rnr_str[out1].raj2000
r_de = rnr_str[out1].dej2000


r_vsini = rnr_str[out1].vsini
r_ul = r_vsini + 3.0
r_ll = r_vsini - 3.0
r_ll_idx = where(strtrim(rnr_str[out1].l_vsini,2) eq '<', rllcnt, complement=r_det_idx)

if rllcnt gt 0 then begin
    r_ll[r_ll_idx] = 0.0
    r_ul[r_ll_idx] = r_vsini[r_ll_idx]
endif

; Get overlapping APG info
a_ra = apg_ra_arr[out2]
a_de = apg_de_arr[out2]
;a_vsini = apg_vsini_str[out2].vsini_best
a_vsini = aspcap_str[out2].vsini

; Check if coord match is reasonable
;delta_ra = fix(abs(r_ra-a_ra)*3600d0*1000)/1000.
;delta_dec = fix(abs(r_de-a_de)*3600d0*1000)/1000.
;delta_ra = abs(r_ra-a_ra)*3600d0
;delta_dec = abs(r_de-a_de)*3600d0

;;; Match with Davison
dav_id = ['LHS 1785', 'G108-021A', 'GJ 1105', 'GJ1182', 'LHS 3075', $
          'GJ 333.2B', 'GJ 363', 'GJ403']
dav_2M = ['2M05470907-0512106', $
          '2M06421118+0334527', $
          '2M07581269+4118134', $
          '2M14153253+0439312', $
          '2M15294392+4252498', $
          '2M09005033+0514293', $
          '2M09422327+5559015', $
          '2M10520440+1359509' $
          ]

dav_vsini = [4.5, $
             0.9, $
             2.0, $
             6.8, $
             2.5, $
             3.0, $
             3.0, $
             3.0 $
             ]

dav_vsini_ul = dav_vsini + [0.6, $
                            0.16, $
                            0.0, $
                            0.5, $
                            0.0, $
                            0.0, $
                            0.0, $
                            0.0 $
                            ]

dav_vsini_ll = dav_vsini - [0.6, $
                            0.22, $
                            dav_vsini[2], $
                            0.5, $
                            dav_vsini[4], $
                            dav_vsini[5], $
                            dav_vsini[6], $
                            dav_vsini[7] $
                            ]

apg_dav_vsini = fltarr(n_elements(dav_vsini)) 
; Loop through the Davison matches and retrieve the APG info
for i = 0, n_elements(dav_2M)-1 do begin
    dav_match_idx = where(apg_2mid2 eq dav_2M[i], davcnt)
    if davcnt ne 1 then message, "Perhaps there's a typo in your dav_2M"
    ; Save apg_vsini
    ;apg_dav_vsini[i] = apg_vsini_str[dav_match_idx].vsini_best
    apg_dav_vsini[i] = aspcap_str[dav_match_idx].vsini
endfor



;;;;;;;;;;;;;;;;;;;;;:HOUDEBINE:;;;;;;;;;;;;;;;;;;;;;
houd_id = ['HIP 53175', $
           'GJ 3643', $
           'LHS 2795', $
           'LHS 2942', $
           'GL 570.2', $
           'GL 629.1' ]
houd_2M = ['2M10523960+0029019', $
          '2M11054316+1014093', $
          '2M13455527+2723131', $
          '2M14333985+0920094', $
          '2M14573227+3123446', $
          '2M16325791-1235296']

houd_vsini = [3.7, $
             1.0, $
             6.4, $
             5.3, $
             2.63, $
             2.3]
; houd_vsini[0] is alternately 2.48 +/- 0.45
; houd_vsini[5] is alternately 2.01 + 0.27 - 0.37
; houd_vsini[5] is alternately 2.7

houd_vsini_ul = houd_vsini + [0.0, $
                              0.0, $
                              0.0, $
                              0.0, $
                              0.0, $
                              0.0]
    
houd_vsini_ll = houd_vsini - [0.0, $
                              0.0, $
                              0.0, $
                              0.0, $
                              0.0, $
                              0.0]
    
nhoud = n_elements(houd_vsini) 


apg_houd_vsini = fltarr(nhoud) 
; Loop through the Houdebine matches and retrieve the APG info
for i = 0, nhoud-1 do begin
    houd_match_idx = where(apg_2mid2 eq houd_2M[i], houdcnt)
    if houdcnt ne 1 then message, "Perhaps there's a typo in your houd_2M"
    ; Save apg_vsini
    
    apg_houd_vsini[i] = aspcap_str[houd_match_idx].vsini
endfor










fx = lindgen(n_elements(out1)+nmatch+n_elements(dav_vsini)+nhoud)
nwtx = lindgen(nmatch)
rnrx = lindgen(n_elements(out1)) + nmatch
davx = lindgen(n_elements(dav_vsini)) + n_elements(rnrx) + nmatch  
houdx = lindgen(nhoud) + n_elements(rnrx) + n_elements(davx) + nmatch  

window, 0, xs = 1850, ys = 1050
plot, nwtx, nwtapg.aspcap_vsini, ps = 6, $
  title = "My vsini's (boxes) vs. Literature (diamonds)", $
  charsize=2, xtitle = "File", ytitle = "vsini (km/s)", xs = 2, yr = [0,35], xr = [0, n_elements(fx)], /ys 
oplot, nwtx, nwtapg.nwt_vsini, ps = 4
errplot, nwtx, nwtapg.nwt_vsini_ll, nwtapg.nwt_vsini_ul

oplot, rnrx, a_vsini, ps = 6
oplot, rnrx, r_vsini, ps = 4, color=200
errplot, rnrx, r_ll, r_ul, color=200

oplot, davx, apg_dav_vsini, ps = 6
oplot, davx, dav_vsini, ps = 4, color=99999
errplot, davx, dav_vsini_ll, dav_vsini_ul, color=99999

oplot, houdx, apg_houd_vsini, ps = 6
oplot, houdx, houd_vsini, ps = 4, color=9999999
errplot, houdx, houd_vsini_ll, houd_vsini_ul, color=9999999

xyouts, davx[3]+0.1, dav_vsini[3], "Highly uncertain according to Jenkins et al 2009"
xyouts, mean(nwtx), 30, "Newton data", charsize=1.5
xyouts, mean(rnrx), 30, "Reiner data", charsize=1.5
xyouts, mean(davx), 30, "Davison and others data", charsize=1.5
xyouts, mean(houdx), 30, "Houdebine and others data", charsize=1.5

;stop

window, 0, xs = 1850, ys = 1000
vsini_lit = [nwtapg.nwt_vsini, r_vsini, dav_vsini, houd_vsini]
vsini_lit_ll = [nwtapg.nwt_vsini_ll, r_ll, dav_vsini_ll, houd_vsini_ll]
vsini_lit_ul = [nwtapg.nwt_vsini_ul, r_ul, dav_vsini_ul, houd_vsini_ul]
vsini_apg = [nwtapg.aspcap_vsini, a_vsini, apg_dav_vsini, apg_houd_vsini]

vsortidx = sort(vsini_lit)
vsini_lit_sort = vsini_lit[vsortidx]
vsini_apg_sort = vsini_apg[vsortidx]
vsini_lit_ll_sort = vsini_lit_ll[vsortidx]
vsini_lit_ul_sort = vsini_lit_ul[vsortidx]
;rot_idx = where(vsini_apg_sort gt 4, rotcnt)

; rot_idx = where(vsini_lit_ll_sort gt 0 and vsini_apg_sort gt 4, rotcnt)
rotcnt = n_elements(vsini_lit_sort) 
rot_idx = lindgen(rotcnt)
litx = vsini_lit_sort[rot_idx]
apgy = vsini_apg_sort[rot_idx]

litxul = vsini_lit_ul_sort[rot_idx]
litxll = vsini_lit_ll_sort[rot_idx]
sigxul = litxul-litx
sigxll = litx-litxll
sigyul = replicate(2,rotcnt)
sigyll = replicate(2,rotcnt)

fitco = poly_fit(litx, apgy, 1, yfit=fit)
;fitx = dindgen(fix(max(vsini_lit)*1.1-4))+4
fitx = dindgen(40)
fity = poly(fitx, fitco)
apgfitx = dindgen(fix(max(vsini_apg)*1.1-4))+4
perfecty = fitx

fitco2 = poly_fit(apgy, litx, 1, yfit=fit2, measure_errors=sigxll)
;plot, vsini_lit, vsini_apg, ps = 6, xtit="vsini from literature", ytit = "vsini fit", title = "APOGEE vsini vs literature vsini", charsize = 2.0
;plot, litx, apgy, ps = 6, xtit="vsini from literature", ytit = "vsini
;fit", title = "APOGEE vsini vs literature vsini", charsize = 2.0,
;xr=[0,30], yr=[0,35], /xs, /ys

entry_device = !d.name
; set_plot, 'ps'
; device, filename = 'vsini_lit_all_aspcap.ps'
; device, /color, bits=8
; device, xs = 13, ys= 8, /inches
; loadct, 13



plot, litx, apgy, ps = 6, xtit="vsini from literature", ytit = "vsini ASPCAP", title = "ASPCAP vsini vs literature vsini", charsize = 2.0, xr=[0,30], yr=[0,35], /xs, /ys
oploterror, litx, apgy, sigxul, sigyul, /hibar, ps = 3
oploterror, litx, apgy, sigxll, sigyll, /lobar, ps = 3
;oplot, fitx, fity
;oplot, poly(apgfitx, fitco2), apgfitx

;oplot, fitx, perfecty, linestyle=2, color=200
oplot, fitx, perfecty, linestyle=2

device, /close_file
set_plot, entry_device

;xyouts, mean(fitx), mean(fity)*1.5, "Best-fit slope = "+strtrim(1./fitco2[1],2), charsize = 2.0
    

stop

end

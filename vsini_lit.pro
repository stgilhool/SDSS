; Readin vsini from literature
pro vsini_lit

skip = 1
if skip then goto, reiners
;;; Paths
apg_path = '/home/stgilhool/APOGEE/vsini_results/'
apg_vsini_file = apg_path+'rfile_master.fits'

;;; APG info
napgspec = 1350L

; Readin APG headers
apg_str = readin_apo(nfiles = napgspec)
; Initialize array to store APG 2mass id's
apg_2mid = strarr(napgspec)
; Readin my results
apg_vsini_str = mrdfits(apg_vsini_file,1)


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
              apg_vsini:0d0 $
             }
nwtapg = replicate(nwtapg_str, nnwt)
nmatch = 0

for i = 0L, napgspec - 1 do begin
    ; Read in APG header for a given file
    head = *(apg_str[i].header)
    apgid = fxpar(head, 'OBJID')
    ; Replace the '2M' with 'J' in the OBJID field
    id = mg_streplace(apgid, '2M', 'J')
    id = strtrim(id,2)
    ; Store the APG id
    apg_2mid[i] = id
    ; Check if there's a match in the Newton data
    nwt_idx = where(nwt_id eq id, match)
    if match eq 1 then begin
        ; Retrieve my result
        apg_vsini = apg_vsini_str[i].vsini_best
        ; Store nwt_idx, apg_idx, vsini's and id
        nwtapg[nmatch].id = id
        nwtapg[nmatch].nwt_idx = nwt_idx[0]
        nwtapg[nmatch].apg_idx = long(i)
        nwtapg[nmatch].nwt_vsini = nwt_vsini[nwt_idx[0]]
        nwtapg[nmatch].apg_vsini = apg_vsini
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
plot, fx, nwtapg.apg_vsini, ps = 6, $
  title = "My vsini's (boxes) vs. Newton's vsini's (diamonds w/ error bars)", $
  charsize=2, xtitle = "File", ytitle = "vsini (km/s)", xs = 2
oplot, fx, nwtapg.nwt_vsini, ps = 4
errplot, fx, nwtapg.nwt_vsini_ll, nwtapg.nwt_vsini_ul

save, /variables, filename= 'vsini_lit.sav'
stop

reiners:
restore, 'vsini_lit.sav'

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
r_ll = rnr_str[out1].l_vsini
r_vsini = rnr_str[out1].vsini
; Get overlapping APG info
a_ra = apg_ra_arr[out2]
a_de = apg_de_arr[out2]
a_vsini = apg_vsini_str[out2].vsini_best

; Check if coord match is reasonable
;delta_ra = fix(abs(r_ra-a_ra)*3600d0*1000)/1000.
;delta_dec = fix(abs(r_de-a_de)*3600d0*1000)/1000.
delta_ra = abs(r_ra-a_ra)*3600d0
delta_dec = abs(r_de-a_de)*3600d0


fx = lindgen(n_elements(out1)) 

window, 0, xs = 1700, ys=900

plot, fx, a_vsini, ps = 6, title = "My vsini (boxes) vs. Reiner's vsini (diamonds)", $
  xtitle = "File", ytitle = "vsini (km/s)", charsize = 2, xs = 2, yr = [0, max(a_vsini)>max(r_vsini)+5]
oplot, fx, r_vsini, ps = 4

;xyouts, fx+0.1, a_vsini, string(delta_ra,delta_dec, format='("(",F0.1,",",F0.1,")")')
xyouts, fx+0.1, a_vsini, string(delta_ra, format='("(",F0.1,",")')
xyouts, fx+0.1, a_vsini-0.5, string(delta_dec, format='(F0.1,")")')
xyouts, fx, a_vsini+2, rnr_str[out1].name
    

stop

end

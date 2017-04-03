; Read in VALD line list
function readin_lines

file_path = '/home/stgilhool/APOGEE/absk_corr/'
file_name = 'line_list.txt'

file = file_path + file_name

; Read it in
readcol, file, element_ion, wl_ctr, excite, vmic, log_gf, rad, $
  stark, waals, lande_factor, depth, reference, skipline=3, $
  delimiter=',', format="A,D,D,F,D,D,D,D,D,D,A"

; Initialize
nlines = n_elements(element_ion) 

fiducial_structure = {SPECIES:'', $
                      WL_CTR:0d0, $
                      DEPTH:0d0 $
                      }

line_list = replicate(fiducial_structure, nlines)

; Populate the line list
element = strarr(nlines)
for line = 0, nlines-1 do begin
    ; Take just the element name
    el = stregex(element_ion[line], '([A-Za-z][A-Za-z]*)', $
                  /extract)
    element[line] = el
    line_list[line].species = el

    ; Add wl's
    line_list[line].wl_ctr = wl_ctr[line]
    
    ; Add depths
    line_list[line].depth = depth[line]

endfor

return, line_list

end

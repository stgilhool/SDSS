; Get APG header info and write to a structure
pro apg_head

; Read in header file
hfile = '/home/stgilhool/APOGEE/apg_headers.fits'
head_byte = mrdfits(hfile,0)
head = string(head_byte)

;;; Make structure
nspec = n_elements(head[0,*]) 
; Get field names for each line
; And create the current structure
final_struct = []
for spec = 0, nspec-1 do begin

    current_header = head[*,spec]
    nlines = n_elements(current_header) 

    tags = []
    current_struct = []

    for line = 0, nlines-1 do begin
        header_line = current_header[line]
        tag = stregex(header_line, '^([^ =][^ ]*).*=', /extract, /subexpr)
        tag = tag[1]
        if tag ne '' then begin
            tags = [tags, tag]
            val = fxpar(current_header,tag)
            
            if current_struct eq !null then current_struct = create_struct(tag, val) $
            else current_struct = create_struct(current_struct, tag, val)
        endif
        
    endfor
    ; Add the full header as a field
    current_struct = create_struct(current_struct, 'HEADER', current_header)

    ; Save completed structure
    final_struct = [final_struct, current_struct]
endfor

; Write the final_structure
mwrfits, final_struct, '/home/stgilhool/APOGEE/apg_head_info.fits', /create

end

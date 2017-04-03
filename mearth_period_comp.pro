pro mearth_period_comp

; Read in Newton's MEarth names and periods
;template = mrdfits('/home/stgilhool/APOGEE/vsini_literature/newton_mearth_ascii_template.fits',1)

mearth_file = '/home/stgilhool/APOGEE/vsini_literature/NewtonMEarth.txt'

;mearth_data = read_ascii(mearth_file, template=template)

openr, lun, mearth_file, /get_lun

fmt_str = '(A1,X,A17,X,A11,X,F10.6,X,F11.8,X,F6.4,X,F6.4,X,A19,X,F7.4,X,F7.4,X,F6.1,X,F4.1,X,A19,X,F6.1,X,F6.1,X,F6.1,X,F8.3,X,F6.4,X,F6.4,X,F5.3,X,F5.3,X,F5.1,X,I1,X,I4,X,F6.4,X,I6)'

lnum=0
type = []
id2mass=[]
idlspm =[]
ra = []
de = []
plx = []
e_plx = []
r_plx = []
pmra = []
pmde = []
rvel = []
e_rvel = []
r_rvel = []
U = []
V = []
W = []
P = []
a = []
e_a = []
mass = []
radius = []
vrot = []
contflag = []
num = []
error = []
fstat = []

; define the strings and ints
line=''
type_i=''
id2mass_i=''
idlspm_i = ''
r_plx_i = ''
r_rvel_i=''
contflag_i = 0s
num_i = 0s
fstat_i = 0s

while not eof(lun) do begin
    
    if lnum lt 49 then begin
        readf, lun, line
        print, lnum
        print, line
        lnum++
        continue
    endif else begin    
        
        readf, lun, type_i, $
          id2mass_i, $
          idlspm_i, $
          ra_i, $
          de_i, $
          plx_i, $
          e_plx_i, $
          r_plx_i, $
          pmra_i, $
          pmde_i, $
          rvel_i, $
          e_rvel_i, $
          r_rvel_i, $
          U_i, $
          V_i, $
          W_i, $
          P_i, $
          a_i, $
          e_a_i, $
          mass_i, $
          radius_i, $
          vrot_i, $
          contflag_i, $
          num_i, $
          error_i, $
          fstat_i, $
          format=fmt_str
        
        ;Save
        type = [type, type_i]
        id2mass = [id2mass, id2mass_i]
        idlspm = [idlspm, idlspm_i]
        ra = [ra, ra_i]
        de = [de, de_i]
        plx = [plx, plx_i]
        e_plx = [e_plx, e_plx_i]
        r_plx = [r_plx, r_plx_i]
        pmra = [pmra, pmra_i]
        pmde = [pmde, pmde_i]
        rvel = [rvel, rvel_i]
        e_rvel = [e_rvel, e_rvel_i]
        r_rvel = [r_rvel, r_rvel_i]
        U = [U, U_i]
        V = [V, V_i]
        W = [W, W_i]
        P = [P, P_i]
        a = [a, a_i]
        e_a = [e_a, e_a_i]
        mass = [mass, mass_i]
        radius = [radius, radius_i]
        vrot = [vrot, vrot_i]
        contflag = [contflag, contflag_i]
        num = [num, num_i]
        error = [error, error_i]
        fstat = [fstat, fstat_i]



        lnum++
    endelse
endwhile

free_lun, lun
    

; Read in my table, just in case

; Read in ASPCAP table
aspcap_file = '/home/stgilhool/APOGEE/APOGEE_data/allparams_dr13.fits'
a_str = mrdfits(aspcap_file,1)
a_objid = a_str.apogee_id

; Find matches
n_nwt = n_elements(id2mass) 

match_idx_apg = []
match_idx_nwt = []
match_id = []
match_per = []
match_vrot = []
match_vsini = []
match_teff = []
match_type = []

for nwt_i = 0, n_nwt-1 do begin
    
    nwt_id = id2mass[nwt_i]
    nwt_id = mg_streplace(nwt_id, 'J','2M')
        
    match_idx_i = where(a_objid eq nwt_id, match_check)

    if match_check then begin
        match_idx_apg = [match_idx_apg, match_idx_i]
        match_idx_nwt = [match_idx_nwt, nwt_i]
        match_id = [match_id, nwt_id]
        match_per = [match_per, P[nwt_i]]
        match_vrot = [match_vrot, vrot[nwt_i]]
        match_vsini = [match_vsini, a_str[match_idx_i].vsini]
        match_teff = [match_teff, a_str[match_idx_i].teff]
        match_type = [match_type, type[nwt_i]]
    endif
endfor

outstr = {apg_idx:match_idx_apg, $
          nwt_idx:match_idx_nwt, $
          objid:match_id, $
          type:match_type, $
          period:match_per, $
          vrot:match_vrot, $
          vsini:match_vsini, $
          teff:match_teff}

mwrfits, outstr, '/home/stgilhool/APOGEE/vsini_literature/mearth_period_comp.fits', /create

print, "N matches: " + strtrim(n_elements(match_idx_apg),2)

stop

end

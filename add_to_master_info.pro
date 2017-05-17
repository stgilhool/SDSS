; Add fields to master_info.fits

pro add_to_master_info, keys, values, master_info_file=master_info_file

if n_elements(master_info_file) eq 0 then $
  master_info_file='/home/stgilhool/APOGEE/master_info6.fits'

; Read in master_info
info = mrdfits(master_info_file,1)
n_elem = n_elements(info)

; Loop through and append data

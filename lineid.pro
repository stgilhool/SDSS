pro init_line_id
;----------------
; label a plot with line indentifcations
; name: string: linelist to read
; thres: threshold strength for reading
; a line into the arrays and, here, for
; plotting it in the plot
;
common linelist,nlines,maxpoint,wl,fl,bb,iext,icnt,$
wll,idl,iml,wlm,idm,imm,wlj,idj,imj
nlines=0L
maxpoint=200000L
wl=fltarr(maxpoint)     ; wavelength of spectrum
fl=fltarr(maxpoint)     ; emitted flux
bb=fltarr(maxpoint)     ; BB of Teff
iext=intarr(maxpoint)   ; index where tau(total)=1
icnt=intarr(maxpoint)   ; index where tau(cont)=1
wll=fltarr(maxpoint)    ; lambda of strongest atomic line
idl=intarr(maxpoint)    ; it's ID (Z*100+ion)
iml=intarr(maxpoint)   ; index where it's tau=1
wlm=fltarr(maxpoint)    ; lambda of strongest molecular line
idm=intarr(maxpoint)    ; it's ID (PHOENIX internal code)
imm=intarr(maxpoint)   ; index where it's tau=1
wlj=fltarr(maxpoint)    ; lambda of strongest JOLA band
idj=intarr(maxpoint)    ; it's ID (France's internal code)
imj=intarr(maxpoint)   ; index where it's tau=1

;
; give the element name and the ioization stages:
common species,elnames,ionstages
elnames=['??',' H','He','Li','Be',' B',' C',' N',' O',$
' F','Ne','Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca',$
'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga',$
'Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc',$
'Ru','Rh','??']
ionstages=['I','II','III','IV','V','VI','VII','VIII','IX','X',$
'XI','XII','XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX',$
'XXI','XXII','XXIII','XXIV','XXV','XXVI','XXVII','XXVIII','XXIX','XXX']
;
;
; give the element name and the ioization stages:
common mspecies,molnames,maxmol,ifacod
maxmol = 151
molnames = strarr(maxmol+1)
ifacod = intarr(maxmol+1)
   molnames(01) = ' ?????????'
   molnames(01) = ' H2O      '
   molnames(02) = ' 46TiO    '
   molnames(03) = ' 47TiO    '
   molnames(04) = ' 48TiO    '
   molnames(05) = ' 49TiO    '
   molnames(06) = ' 50TiO    '
   molnames(07) = ' H2       '
   molnames(08) = ' 12CH     '
   molnames(09) = ' 13CH     '
   molnames(10) = ' 14NH     '
   molnames(11) = ' 15NH     '
   molnames(12) = ' 16OH     '
   molnames(13) = ' 17OH     '
   molnames(14) = ' 18OH     '
   molnames(15) = ' 24MgH    '
   molnames(16) = ' 25MgH    '
   molnames(17) = ' 26MgH    '
   molnames(18) = ' 28SiH    '
   molnames(19) = ' 29SiH    '
   molnames(20) = ' 30SiH    '
   molnames(21) = ' 12C12C   '
   molnames(22) = ' 12C13C   '
   molnames(23) = ' 13C13C   '
   molnames(24) = ' 12C14N   '
   molnames(25) = ' 13C14N   '
   molnames(26) = ' 12C15N   '
   molnames(27) = ' 13C15N   '
   molnames(28) = ' 12C16O   '
   molnames(29) = ' 13C16O   '
   molnames(30) = ' 12C17O   '
   molnames(31) = ' 13C17O   '
   molnames(32) = ' 12C18O   '
   molnames(33) = ' 13C18O   '
   molnames(34) = ' 28Si16O  '
   molnames(35) = ' 29Si16O  '
   molnames(36) = ' 30Si16O  '
   molnames(37) = ' 28Si18O  '
   molnames(38) = ' 46TiO (J)'
   molnames(39) = ' 47TiO (J)'
   molnames(40) = ' 48TiO (J)'
   molnames(41) = ' 49TiO (J)'
   molnames(42) = ' 50TiO (J)'
   molnames(43) = ' H2-16O(H)'
   molnames(44) = ' H2-18O(H)'
   molnames(45) = ' H2-17O(H)'
   molnames(46) = ' H16OD (H)'
   molnames(47) = ' 12C16O2  '
   molnames(48) = ' 13C16O2  '
   molnames(49) = ' 12C16O18O'
   molnames(50) = ' 12C16O17O'
   molnames(51) = ' 13C16O18O'
   molnames(52) = ' 13C16O17O'
   molnames(53) = ' 12C18O2  '
   molnames(54) = ' 12C18O17O'
   molnames(55) = ' 16O3     '
   molnames(56) = ' 16O16O18O'
   molnames(57) = ' 16O18O16O'
   molnames(58) = ' 14N2-16O '
   molnames(59) = ' 14N15N16O'
   molnames(60) = ' 15N14N16O'
   molnames(61) = ' 14N2-18O '
   molnames(62) = ' 14N2-17O '
   molnames(63) = ' 12C16O(H)'
   molnames(64) = ' 13C16O(H)'
   molnames(65) = ' 12C18O(H)'
   molnames(66) = ' 12C17O(H)'
   molnames(67) = ' 13C18O(H)'
   molnames(68) = ' 12CH4    '
   molnames(69) = ' 13CH4    '
   molnames(70) = ' 12CH3D   '
   molnames(71) = ' 16O2     '
   molnames(72) = ' 16O18O   '
   molnames(73) = ' 16O17O   '
   molnames(74) = ' 14N16O   '
   molnames(75) = ' 15N16O   '
   molnames(76) = ' 14N18O   '
   molnames(77) = ' 32S16O2  '
   molnames(78) = ' 34S16O2  '
   molnames(79) = ' 14N16O2  '
   molnames(80) = ' 14NH3    '
   molnames(81) = ' 15NH3    '
   molnames(82) = ' H14N16O3 '
   molnames(83) = ' 16OH (H) '
   molnames(84) = ' 18OH (H) '
   molnames(85) = ' 16OD (H) '
   molnames(86) = ' H19F     '
   molnames(87) = ' H35Cl    '
   molnames(88) = ' H37Cl    '
   molnames(89) = ' H79Br    '
   molnames(90) = ' H81Br    '
   molnames(91) = ' H127I    '
   molnames(92) = ' 35Cl16O  '
   molnames(93) = ' 37Cl16O  '
   molnames(94) = ' 16O12C32S'
   molnames(95) = ' 16O12C34S'
   molnames(96) = ' 16O13C32S'
   molnames(97) = ' 18O12C32S'
   molnames(98) = ' H2-12C16O'
   molnames(99) = ' H2-13C16O'
  molnames(100) = ' H2-12C18O'
  molnames(101) = ' H16O35Cl '
  molnames(102) = ' H16O37Cl '
  molnames(103) = ' 14N2     '
  molnames(104) = ' H12C14N  '
  molnames(105) = ' H13C14N  '
  molnames(106) = ' H12C15N  '
  molnames(107) = '12CH3-35Cl'
  molnames(108) = '12CH3-37Cl'
  molnames(109) = ' H2-16O2  '
  molnames(110) = ' 12C2H2   '
  molnames(111) = ' 12C13CH2 '
  molnames(112) = ' 12C2H6   '
  molnames(113) = ' 31PH3    '
  molnames(114) = '12C16O19F2'
  molnames(115) = ' 32S19F6  '
  molnames(116) = ' H2-32S   '
  molnames(117) = ' H2-34S   '
  molnames(118) = ' H2-33S   '
  molnames(119) = ' 12C16O(G)'
  molnames(120) = ' 13C16O(G)'
  molnames(121) = ' 12C18O(G)'
  molnames(122) = ' 12C17O(G)'
  molnames(123) = ' 13C18O(G)'
  molnames(124) = ' 13C17O(G)'
  molnames(125) = ' 14C16O(G)'
  molnames(126) = ' 12C14N(J)'
  molnames(127) = ' 13C14N(J)'
  molnames(128) = ' 28SiO (N)'
  molnames(129) = ' 29SiO (N)'
  molnames(130) = ' 30SiO (N)'
  molnames(131) = ' H3+      '
  molnames(132) = ' YO       '
  molnames(133) = ' 90ZrO    '
  molnames(134) = ' 91ZrO    '
  molnames(135) = ' 92ZrO    '
  molnames(136) = ' 93ZrO    '
  molnames(137) = ' 94ZrO    '
  molnames(138) = ' 95ZrO    '
  molnames(139) = ' 96ZrO    '
  molnames(140) = ' H2O(AMES)'
  molnames(141) = '12CH4(GH) '
  molnames(142) = '13CH4(GH) '
  molnames(143) = '12CH3D(GH)'
  molnames(144) = '12CH4(GE) '
  molnames(145) = '13CH4(GE) '
  molnames(146) = '12CH3D(GE)'
  molnames(147) = ' 46TiO (S)'
  molnames(148) = ' 47TiO (S)'
  molnames(149) = ' 48TiO (S)'
  molnames(150) = ' 49TiO (S)'
  molnames(151) = ' 50TiO (S)'

;
; translation table to France's codes:

;-- H2O:
      ifacod(1) = 112 
;
;-- TiO:
      ifacod(2) = 111 
      ifacod(3) = 111 
      ifacod(4) = 111 
      ifacod(5) = 111 
      ifacod(6) = 111 
;
;-- TiO (Jorgensen):
      ifacod(38) = 111 
      ifacod(39) = 111 
      ifacod(40) = 111 
      ifacod(41) = 111 
      ifacod(42) = 111 
;
;-- H2:
      ifacod(7) = 113 
;
;-- CH:
      ifacod(8) = 103 
      ifacod(9) = 103 
;
;-- NH:
      ifacod(10) = 104 
      ifacod(11) = 104 
;
;-- OH:
      ifacod(12) = 102 
      ifacod(13) = 102 
      ifacod(14) = 102 
;
;-- MgH:
      ifacod(15) = 108 
      ifacod(16) = 108 
      ifacod(17) = 108 
;
;-- SiH:
      ifacod(18) = 110 
      ifacod(19) = 110 
      ifacod(20) = 110 
;
;-- C2:
      ifacod(21) = 105 
      ifacod(22) = 105 
      ifacod(23) = 105 
;
;-- CN:
      ifacod(24) = 106 
      ifacod(25) = 106 
      ifacod(26) = 106 
      ifacod(27) = 106 
;
;-- CO:
      ifacod(28) = 107 
      ifacod(29) = 107 
      ifacod(30) = 107 
      ifacod(31) = 107 
      ifacod(32) = 107 
      ifacod(33) = 107 
;
;-- SiO:
      ifacod(34) = 121 
      ifacod(35) = 121 
      ifacod(36) = 121 
      ifacod(37) = 121 
;
;-- H2O HITRAN :
      ifacod(43) = 112 
      ifacod(44) = 112 
      ifacod(45) = 112 
      ifacod(46) = 112 
;
;-- CO2 HITRAN :
      ifacod(47) = 116 
      ifacod(48) = 116 
      ifacod(49) = 116 
      ifacod(50) = 116 
      ifacod(51) = 116 
      ifacod(52) = 116 
      ifacod(53) = 116 
      ifacod(54) = 116 
;
;-- O3 HITRAN :
      ifacod(55) = 20000 
      ifacod(56) = 20000 
      ifacod(57) = 20000 
;
;-- N2O HITRAN :
      ifacod(58) = 20001 
      ifacod(59) = 20001 
      ifacod(60) = 20001 
      ifacod(61) = 20001 
      ifacod(62) = 20001 
;
;-- CO HITRAN :
      ifacod(63) = 107 
      ifacod(64) = 107 
      ifacod(65) = 107 
      ifacod(66) = 107 
      ifacod(67) = 107 
;
;-- CH4 HITRAN :
      ifacod(68) = 139 
      ifacod(69) = 139 
      ifacod(70) = 139 
;
;-- O2 HITRAN :
      ifacod(71) = 117 
      ifacod(72) = 117 
      ifacod(73) = 117 
;-- NO HITRAN :
      ifacod(74) = 115 
      ifacod(75) = 115 
      ifacod(76) = 115 
;
;-- SO2 HITRAN :
      ifacod(77) = 20002 
      ifacod(78) = 20002 
;
;-- NO2 Hitran :
      ifacod(79) = 20003 
;
;-- NH3 HITRAN :
      ifacod(80) = 170 
      ifacod(81) = 170 
;
;-- HNO3 HITRAN :
      ifacod(82) = 20004 
;
;-- OH HITRAN :
      ifacod(83) = 102 
      ifacod(84) = 102 
      ifacod(85) = 102 
;
;-- HF HITRAN :
      ifacod(86) = 124 
;
;-- HCl HITRAN :
      ifacod(87) = 123 
      ifacod(88) = 123 
;
;-- HBr HITRAN :
      ifacod(89) = 20005 
      ifacod(90) = 20005 
;
;-- HI HITRAN :
      ifacod(91) = 20006 
;
;-- ClO HITRAN :
      ifacod(92) = 20007 
      ifacod(93) = 20007 
;
;-- OCS HITRAN :
      ifacod(94) = 155 
      ifacod(95) = 155 
      ifacod(96) = 155 
      ifacod(97) = 155 
;
;-- H2CO HITRAN :
      ifacod(98) = 20008 
      ifacod(99) = 20008 
      ifacod(100) = 20008 
;
;-- HOCl HITRAN :
      ifacod(101) = 20009 
      ifacod(102) = 20009 
;-- N2 HITRAN :
      ifacod(103) = 114 
;
;-- HCN HITRAN :
      ifacod(104) = 137 
      ifacod(105) = 137 
      ifacod(106) = 137 
;
;-- CH3Cl HITRAN :
      ifacod(107) = 20010 
      ifacod(108) = 20010 
;
;-- H2O2 HITRAN :
      ifacod(109) = 20011 
;
;-- C2H2 HITRAN :
      ifacod(110) = 138 
      ifacod(111) = 138 
;
;-- C2H6 HITRAN :
      ifacod(112) = 20012 
;
;-- PH3 HITRAN :
      ifacod(113) = 20013 
;
;-- COF2 HITRAN :
      ifacod(114) = 20014 
;
;-- SF6 HITRAN :
      ifacod(115) = 20015 
;
;-- H2S HITRAN :
      ifacod(116) = 154 
      ifacod(117) = 154 
      ifacod(118) = 154 
;
;-- CO NASA-Goorvitch :
      ifacod(119) = 107 
      ifacod(120) = 107 
      ifacod(121) = 107 
      ifacod(122) = 107 
      ifacod(123) = 107 
      ifacod(124) = 107 
      ifacod(125) = 107 
;
;-- CN Joergensen
      ifacod(126) = 106 
      ifacod(127) = 106 
;
; SiO
      ifacod(128) = 121
      ifacod(129) = 121    
      ifacod(130) = 121    
;
; H3+
      ifacod(131) = 21301  
;
; YO
      ifacod(132) = 133    
;
; ZrO
      ifacod(133) = 118    
      ifacod(134) = 118    
      ifacod(135) = 118    
      ifacod(136) = 118    
      ifacod(137) = 118    
      ifacod(138) = 118    
      ifacod(139) = 118    
;
; H2O (AMES)
      ifacod(140) = 112    
;
; CH4 (GEISA)
      ifacod(141) = 139    
      ifacod(142) = 139    
      ifacod(143) = 139    
      ifacod(144) = 139    
      ifacod(145) = 139    
      ifacod(146) = 139    
;
;-- TiO (Schwenke):
      ifacod(147) = 111 
      ifacod(148) = 111 
      ifacod(149) = 111 
      ifacod(150) = 111 
      ifacod(151) = 111 

;
;
; translation: France's code to string ID
common jspecies,coden,fname,nmols
nmols = 99
coden = intarr(nmols+1)
fname = strarr(nmols+1)

coden=[$
102, 103, 104, 105, 106, 107, 108, 109, 110, 111,$
112, 113, 114, 115, 116, 117, 118, 119, 120, 121,$
122, 123, 124, 125, 126, 127, 128, 129, 130, 131,$
132, 133, 134, 135, 136, 137, 138, 139, 140, 141,$
142, 143, 144, 145, 146, 147, 148, 149, 150, 151,$
152, 153, 154, 155, 156, 157, 158, 159, 160, 161,$
162, 163, 164, 165, 166, 167, 168, 169, 170, 171,$
172, 173, 174, 175, 176, 177, 178, 179, 180, 181,$
182, 183, 184, 185, 186, 187, 188, 189, 190, 191]
coden=[coden,$
192, 193, 194, 195, 196, 197, 198, 199]

fname=[$
'  OH ','  CH ','  NH ','  C2 ','  CN ','  CO ','  MgH',$
'  CaH','  SiH','  TiO','  H2O','  H2 ','  N2 ','  NO ',$
'  CO2','  O2 ','  ZrO','  VO ','  MgS','  SiO','  AlH',$
'  HCl','  HF ','  HS ','  TiH','  AlO','  BO ','  CrO',$
'  LaO','  MgO','  ScO','  YO ','  SiF',' NaCl',' CaOH',$
'  HCN',' C2H2','  CH4','  CH2','  C2H','  HCO','  NH2',$
' LiOH','  C2O',' AlOF',' NaOH',' MgOH',' AlO2',' Al2O',$
' AlOH',' SiH2',' SiO2','  H2S','  OCS','  KOH',' TiO2',$
'TiOCl','  VO2',' FeF2','  YO2',' ZrO2',' BaOH',' LaO2',$
' C2H4','  C3 ',' SiC2','  CH3','  C3H','  NH3',' C2N2',$
'  C2N',' CaF2','AlOCl',' Si2C','  CS2','CaCl2','  AlF',$
'  CaF','  Si2','  SiS','  CS ',' AlCl','  KCl',' CaCl']
fname=[fname,$
'  TiS',' TiCl','  SiN','  AlS','  AL2','  FeO','  SiC',$
' TiF2','  FeH',' LiCl','   NS','  NaH','   SO','   S2']


end

;
pro rdid,name,w,f,b
;------------------
; label a plot with line indentifcations
; name: string: linelist to read
; thres: threshold strength for reading
; a line into the arrays and, here, for
; plotting it in the plot
;
common linelist,nlines,maxpoint,wl,fl,bb,iext,icnt,$
wll,idl,iml,wlm,idm,imm,wlj,idj,imj
print,'reading spectrum+ID file:',name
openr,/get_lun,lun,name
;
i = -1L
while not eof(lun) do begin
 i = i+1L
 readf,lun,wldum,fldum,bbdum,iextdum,icntdum,$
  wlldum,idldum,imldum,wlmdum,idmdum,immdum,$
  wljdum,idjdum,imjdum,format='(1x,f11.3,2d12.4,2i3,3(f11.4,i5,i3))'
 wl(i) = wldum
 wll(i) = wlldum
 wlm(i) = wlmdum
 wlj(i) = wljdum
 fl(i) = fldum
 bb(i) = bbdum
 iext(i) = iextdum
 icnt(i) = icntdum
 idl(i) = idldum
 iml(i) = imldum
 idm(i) = idmdum
 idj(i) = idjdum
 imm(i) = immdum
 imj(i) = imjdum
endwhile
close,lun
nlines = i
wl = wl(0:nlines)
wll = wll(0:nlines)
wlm = wlm(0:nlines)
fl = fl(0:nlines)
bb = bb(0:nlines)
iext = iext(0:nlines)
icnt = icnt(0:nlines)
idl = idl(0:nlines)
iml = iml(0:nlines)
idm = idm(0:nlines)
imm = imm(0:nlines)
fl=10.^(fl-40.)
bb=10.^(bb-40.)
w = wl
f = fl
b = bb
print,'number of lines stored:',nlines
end
;
;
;
pro rdidg,name,w,f,b
;-------------------
; label a plot with line indentifcations
; name: string: linelist to read
; thres: threshold strength for reading
; a line into the arrays and, here, for
; plotting it in the plot
;
common linelist,nlines,maxpoint,wl,fl,bb,iext,icnt,$
wll,idl,iml,wlm,idm,imm,wlj,idj,imj
nlines=0L
;
print,'reading spectrum+ID file:',name
spawn,string("gunzip -c "+name+' >rdyng_tmp'),ierr
openr,/get_lun,lun,"rdyng_tmp"
;
i = -1L
while not eof(lun) do begin
 i = i+1L
 readf,lun,wldum,fldum,bbdum,iextdum,icntdum,$
  wlldum,idldum,imldum,wlmdum,idmdum,immdum,$
  wljdum,idjdum,imjdum
 wl(i) = wldum
 wll(i) = wlldum
 wlm(i) = wlmdum
 wlj(i) = wljdum
 fl(i) = fldum
 bb(i) = bbdum
 iext(i) = iextdum
 icnt(i) = icntdum
 idl(i) = idldum
 iml(i) = imldum
 idm(i) = idmdum
 idj(i) = idjdum
 imm(i) = immdum
 imj(i) = imjdum
endwhile
close,lun
spawn,string("rm rdyng_tmp"),ierr
nlines = i
wl = wl(0:nlines)
wll = wll(0:nlines)
wlm = wlm(0:nlines)
fl = fl(0:nlines)
bb = bb(0:nlines)
iext = iext(0:nlines)
icnt = icnt(0:nlines)
idl = idl(0:nlines)
iml = iml(0:nlines)
idm = idm(0:nlines)
imm = imm(0:nlines)
fl=10.^(fl-40.)
bb=10.^(bb-40.)
w = wl
f = fl
b = bb
print,'number of lines stored:',nlines
end
;
;
;
pro rdidgf,name,w,f,b
;--------------------
; read a NEW format spectrum:
; label a plot with line indentifcations
; name: string: linelist to read
; thres: threshold strength for reading
; a line into the arrays and, here, for
; plotting it in the plot
;
common linelist,nlines,maxpoint,wl,fl,bb,iext,icnt,$
wll,idl,iml,wlm,idm,imm,wlj,idj,imj
nlines=0L
;
print,'reading spectrum+ID file:',name
spawn,string("gunzip -c "+name+' >rdyng_tmp'),ierr
openr,/get_lun,lun,"rdyng_tmp"
;
i = -1L
while not eof(lun) do begin
 i = i+1L
 readf,lun,wldum,fldum,bbdum,iextdum,icntdum,$
  wlldum,idldum,imldum,wlmdum,idmdum,immdum,$
  wljdum,idjdum,imjdum,format='(f12.3,2e12.4,2i3,3(f11.4,i5,i3))'
 wl(i) = wldum
 wll(i) = wlldum
 wlm(i) = wlmdum
 wlj(i) = wljdum
 fl(i) = fldum
 bb(i) = bbdum
 iext(i) = iextdum
 icnt(i) = icntdum
 idl(i) = idldum
 iml(i) = imldum
 idm(i) = idmdum
 idj(i) = idjdum
 imm(i) = immdum
 imj(i) = imjdum
endwhile
close,lun
free_lun,lun
spawn,string("rm rdyng_tmp"),ierr
nlines = i
wl = wl(0:nlines)
wll = wll(0:nlines)
wlm = wlm(0:nlines)
fl = fl(0:nlines)
bb = bb(0:nlines)
iext = iext(0:nlines)
icnt = icnt(0:nlines)
idl = idl(0:nlines)
iml = iml(0:nlines)
idm = idm(0:nlines)
imm = imm(0:nlines)
fl=10.^fl
bb=10.^bb
w = wl
f = fl
b = bb
print,'number of lines stored:',nlines
end
;
;
;
pro plotid,labst,txtst,flag,molflg,delatom,delmol
;------------------------------------------------
; lable a plot with line identifications
; take from common block linelist.
;-- input:
; thres: threshold for line strength for a line
;        to be labeled.
; labst: relative y-coordiate to start label (line), e.g. 0.75
; txtst: relative y-coordiate to start text , e.g. 0.80
; flag: 0: print numeric identification
;       1: print text identification
; molflg: 0: do not molecular lines, 1: plot 'em
; delatom:
; delmol: parameter for molecular print threshold. gives
;         difference between imm and iext for lines to plot
common linelist,nlines,maxpoint,wl,fl,bb,iext,icnt,$
wll,idl,iml,wlm,idm,imm
common species,elnames,ionstages
common mspecies,molnames,maxmol,ifacod
;
; get the y axis minimum and maximum:
ymin=!y.crange(0)
ymax=!y.crange(1)
;
; same for x axis
xmin=!x.crange(0)
xmax=!x.crange(1)
;
; set the position for the label: in the upper 75% 
; of the y-axis range
ystart = labst*ymax
ytext =  txtst*ymax
;
; go through the list and label the plot, 
; if a line above the threshold is found
wllo = 0.d0
idlo = 0
wlmo = 0.d0
idmo = 0
i=0L
for i=0L,nlines do begin
 if iml(i) gt icnt(i)+delatom or iml(i) gt imm(i) $
   or wll(i) lt xmin  or wll(i) gt xmax then goto,l666
;
   if wllo eq wll(i) and idlo eq idl(i) then goto,l666
;
; atomic lines
 if flag eq 1 then begin
  el = fix(idl(i)/100)
  ion = fix(idl(i)-100*el)
  if el le 46 then begin
   sid = elnames(min([el,46]))
   sion = ionstages(ion)
   swl = strtrim(string(wll(i),format='(f7.1)'),2)
   labstr = '!B '+sid+' '+sion+' '+swl
  endif else  begin
   sid = string(fix(idl(i)))
   swl = strtrim(string(wll(i),format='(f7.1)'),2)
   labstr = '!B '+sid+' '+swl
  endelse
 endif
 if flag eq 0 then begin
  sid = string(fix(idl(i)))
  swl = strtrim(string(wll(i),format='(f7.1)'),2)
  labstr = '!B '+sid+' '+swl
 endif
 plots, [wll(i),wll(i)],[ystart,ytext]
 xyouts, wll(i),ytext,labstr,orientation=90,/normal
 wllo = wll(i)
 idlo = idl(i)
;
; molecular lines
l666:
 if molflg eq 0 then goto,l667
 if imm(i) gt iext(i)+delmol or imm(i) gt iml(i) $
   or wlm(i) lt xmin  or wlm(i) gt xmax then goto,l667
 if wlmo eq wlm(i) and idmo eq idm(i) then goto,l667
  if flag eq 1 then begin
   el = fix(idm(i))
   if el le maxmol then begin
    sid = molnames(min([el,maxmol,ifacod+1]))
    swl = strtrim(string(wlm(i),format='(f7.1)'),2)
    labstr = '!B '+sid+' '+swl
   endif else  begin
    sid = string(fix(idm(i)))
    swl = strtrim(string(wlm(i),format='(f7.1)'),2)
    labstr = '!B '+sid+' '+swl
   endelse
  endif
  if flag eq 0 then begin
   sid = string(fix(ident(i)))
   swl = strtrim(string(wavel(i),format='(f7.1)'),2)
   labstr = '!B '+sid+' '+swl
  endif
 plots, [wlm(i),wlm(i)],[ystart,ytext]
 xyouts, wlm(i),ytext,labstr,orientation=90
 wlmo = wlm(i)
 idmo = idm(i)
l667:
endfor
end
;
;
;
pro molid,labst,txtst,flag,delmol
;--------------------------------
; lable a plot with line identifications
; take from common block linelist.
;-- input:
; thres: threshold for line strength for a line
;        to be labeled.
; labst: relative y-coordiate to start label (line), e.g. 0.75
; txtst: relative y-coordiate to start text , e.g. 0.80
; flag: 0: print numeric identification
;       1: print text identification
; delmol: parameter for molecular print threshold. gives
;         difference between imm and iext for lines to plot
common linelist,nlines,maxpoint,wl,fl,bb,iext,icnt,$
wll,idl,iml,wlm,idm,imm,wlj,idj,imj
common species,elnames,ionstages
common mspecies,molnames,maxmol,ifacod
;
; get the y axis minimum and maximum:
ymin=!y.crange(0)
ymax=!y.crange(1)
;
; same for x axis
xmin=!x.crange(0)
xmax=!x.crange(1)
;
; set the position for the label: in the upper 75% 
; of the y-axis range
ystart = labst*ymax
ytext =  txtst*ymax
;
; go through the list and label the plot, 
; if a line above the threshold is found
wllo = 0.d0
idlo = 0
wlmo = 0.d0
idmo = 0
idmcurr = 0
i=0L
endflg = 0
for i=0L,nlines do begin
 if endflg eq 1 then return
 if((i eq nlines) or (xmax le wlm(i))) then begin
  endflg = 1
  if(wlmo eq 0) then begin
   wlmo = wlm(i)
   wlmmin = wlmo
  endif
  goto,l666
 endif
 if imm(i) gt iext(i)+delmol $;or imm(i) gt iml(i) $
   or wlm(i) lt xmin  or wlm(i) gt xmax then goto,l667
 if wlmo eq wlm(i) and idmo eq idm(i) then goto,l667
 plots, [wlm(i),wlm(i)],[ystart,ytext]
l666:
 if(endflg eq 1 or ((idm(i) ne idmcurr) and (idmcurr ne 0))) then begin
; we need to plot the horizontal bar and the label
; for the last series of molecular species to plot 
  el = fix(idmcurr)
;  if el le maxmol then begin
;   sid = molnames(min([el,maxmol,ifacod+1]))
;   labstr = '!B '+sid
;  endif else  begin
;   sid = string(fix(idmcurr))
;   labstr = '!B '+sid
;  endelse
 if(flag ne 0) then begin
;  id = min(where(idmcurr eq coden))
  id = idmcurr
  if(id ge 0) then labstr = '!B'+molnames(id)
endif
  if flag eq 0 then begin
   sid = string(fix(idmcurr))
   labstr = '!B '+sid
  endif
  plots, [wlmmin,wlmo],[ytext,ytext]
  xyouts, 0.5*(wlmo+wlmmin),ytext,labstr,orientation=90
  wlmmin = wlm(i)
  idmcurr = idm(i)
 endif
 wlmo = wlm(i)
 idmo = idm(i)
 if(idmcurr eq 0) then begin
  wlmmin = wlmo
  idmcurr = idmo
 endif
l667:
endfor
end
;
;
;
pro atomid,labst,txtst,flag,delatom
;----------------------------------
; lable a plot with line identifications
; take from common block linelist.
;-- input:
; thres: threshold for line strength for a line
;        to be labeled.
; labst: relative y-coordiate to start label (line), e.g. 0.75
; txtst: relative y-coordiate to start text , e.g. 0.80
; flag: 0: print numeric identification
;       1: print text identification
; molflg: 0: do not molecular lines, 1: plot 'em
; delatom:
;  parameter for molecular print threshold. gives
;         difference between imm and iext for lines to plot
common linelist,nlines,maxpoint,wl,fl,bb,iext,icnt,$
wll,idl,iml,wlm,idm,imm
common species,elnames,ionstages
common mspecies,molnames,maxmol,ifacod
;
; get the y axis minimum and maximum:
ymin=!y.crange(0)
ymax=!y.crange(1)
;
; same for x axis
xmin=!x.crange(0)
xmax=!x.crange(1)
;
; set the position for the label: in the upper 75% 
; of the y-axis range
ystart = labst*ymax
ytext =  txtst*ymax
;
; go through the list and label the plot, 
; if a line above the threshold is found
wllo = 0.d0
idlo = 0
wlmo = 0.d0
idmo = 0
i=0L
for i=0L,nlines do begin
 if iml(i) gt iext(i)+delatom $
   or wll(i) lt xmin  or wll(i) gt xmax then goto,l666
;
   if wllo eq wll(i) and idlo eq idl(i) then goto,l666
;
; atomic lines
 if flag eq 1 then begin
  el = fix(idl(i)/100)
  ion = fix(idl(i)-100*el)
  if el le 46 then begin
   sid = elnames(min([el,46]))
   sion = ionstages(ion)
   swl = strtrim(string(wll(i),format='(f7.1)'),2)
   labstr = '!B '+sid+' '+sion+' '+swl
  endif else  begin
   sid = string(fix(idl(i)))
   swl = strtrim(string(wll(i),format='(f7.1)'),2)
   labstr = '!B '+sid+' '+swl
  endelse
 endif
 if flag eq 0 then begin
  sid = string(fix(idl(i)))
  swl = strtrim(string(wll(i),format='(f7.1)'),2)
  labstr = '!B '+sid+' '+swl
 endif
 plots, [wll(i),wll(i)],[ystart,ytext]
 xyouts, wll(i),ytext,labstr,orientation=90
 wllo = wll(i)
 idlo = idl(i)
;
l666:
l667:
endfor
end
;
;
;
pro atomidw,labst,txtst,flag,delatom,xmin,xmax
;----------------------------------
; lable a plot with line identifications
; take from common block linelist.
;-- input:
; thres: threshold for line strength for a line
;        to be labeled.
; labst: relative y-coordiate to start label (line), e.g. 0.75
; txtst: relative y-coordiate to start text , e.g. 0.80
; flag: 0: print numeric identification
;       1: print text identification
; molflg: 0: do not molecular lines, 1: plot 'em
; delatom:
;  parameter for molecular print threshold. gives
;         difference between imm and iext for lines to plot
common linelist,nlines,maxpoint,wl,fl,bb,iext,icnt,$
wll,idl,iml,wlm,idm,imm
common species,elnames,ionstages
common mspecies,molnames,maxmol,ifacod
;
; get the y axis minimum and maximum:
ymin=!y.crange(0)
ymax=!y.crange(1)
;
; same for x axis
;xmin=!x.crange(0)
;xmax=!x.crange(1)
;
; set the position for the label: in the upper 75% 
; of the y-axis range
ystart = labst*ymax
ytext =  txtst*ymax
;
; go through the list and label the plot, 
; if a line above the threshold is found
wllo = 0.d0
idlo = 0
wlmo = 0.d0
idmo = 0
i=0L
for i=0L,nlines do begin
 if iml(i) gt iext(i)+delatom $
   or wll(i) lt xmin  or wll(i) gt xmax then goto,l666
;
   if wllo eq wll(i) and idlo eq idl(i) then goto,l666
;
; atomic lines
 if flag eq 1 then begin
  el = fix(idl(i)/100)
  ion = fix(idl(i)-100*el)
  if el le 46 then begin
   sid = elnames(min([el,46]))
   sion = ionstages(ion)
   swl = strtrim(string(wll(i),format='(f7.1)'),2)
   labstr = '!B '+sid+' '+sion+' '+swl
  endif else  begin
   sid = string(fix(idl(i)))
   swl = strtrim(string(wll(i),format='(f7.1)'),2)
   labstr = '!B '+sid+' '+swl
  endelse
 endif
 if flag eq 0 then begin
  sid = string(fix(idl(i)))
  swl = strtrim(string(wll(i),format='(f7.1)'),2)
  labstr = '!B '+sid+' '+swl
 endif
 plots, [wll(i),wll(i)],[ystart,ytext]
 xyouts, wll(i),ytext,labstr,orientation=90
 wllo = wll(i)
 idlo = idl(i)
;
l666:
l667:
endfor
end
;
;
;
pro jolaid,labst,txtst,flag,delmol
;---------------------------------
; lable a plot with JOLA identifications
; take from common block linelist.
;-- input:
; thres: threshold for line strength for a line
;        to be labeled.
; labst: relative y-coordiate to start label (line), e.g. 0.75
; txtst: relative y-coordiate to start text , e.g. 0.80
; flag: 0: print numeric identification
;       1: print text identification
; delmol: parameter for molecular print threshold. gives
;         difference between imm and iext for lines to plot
common linelist,nlines,maxpoint,wl,fl,bb,iext,icnt,$
wll,idl,iml,wlm,idm,imm,wlj,idj,imj
common jspecies,coden,fname,nmols
;
; get the y axis minimum and maximum:
ymin=!y.crange(0)
ymax=!y.crange(1)
;
; same for x axis
xmin=!x.crange(0)
xmax=!x.crange(1)
;
; set the position for the label: in the upper 75% 
; of the y-axis range
ystart = labst*ymax
ytext =  txtst*ymax
;
; go through the list and label the plot, 
; if a line above the threshold is found
wllo = 0.d0
idlo = 0
wlmo = 0.d0
idmo = 0
idmcurr = 0
i=0L
endflg = 0
for i=0L,nlines do begin
 if endflg eq 1 then return
 if((i eq nlines) or (xmax le wlj(i))) then begin
  endflg = 1
  if(wlmo eq 0) then begin
   wlmo = wlj(i)
   wlmmin = wlmo
  endif
  goto,l666
 endif
 if imj(i) gt iext(i)+delmol $; or imj(i) gt iml(i) $
  or wlj(i) lt xmin  or wlj(i) gt xmax then goto,l667
 if wlmo eq wlj(i) and idmo eq idj(i) then goto,l667
 plots, [wlj(i),wlj(i)],[ystart,ytext]
l666:
 if((endflg eq 1 and idmcurr ne 0) or ((idj(i) ne idmcurr) and (idmcurr ne 0))) then begin
; we need to plot the horizontal bar and the label
; for the last series of molecular species to plot 
  if(flag ne 0) then begin
   for id=0,nmols-1 do begin
    if(idmcurr eq coden(id)) then begin
     sid = fname(id)
     labstr = '!B'+sid+' (JOLA)'
     goto,l888
    endif
   endfor
  endif
l888:
  if flag eq 0 then begin
   sid = string(fix(idmcurr))
   labstr = '!B'+sid+' (JOLA)'
  endif
  plots, [wlmmin,wlmo],[ytext,ytext]
  xyouts, 0.5*(wlmo+wlmmin),ytext,labstr,orientation=90
  wlmmin = wlj(i)
  idmcurr = idj(i)
 endif
 wlmo = wlj(i)
 idmo = idj(i)
 if(idmcurr eq 0) then begin
  wlmmin = wlmo
  idmcurr = idmo
 endif
l667:
endfor
end
;
;
;
pro molidc,labst,txtst,flag,delmol,docodes
;---------------------------------
; lable a plot with line identifications
; take from common block linelist.
;-- input:
; thres: threshold for line strength for a line
;        to be labeled.
; labst: relative y-coordiate to start label (line), e.g. 0.75
; txtst: relative y-coordiate to start text , e.g. 0.80
; flag: 0: print numeric identification
;       1: print text identification
; delmol: parameter for molecular print threshold. gives
;         difference between imm and iext for lines to plot
common linelist,nlines,maxpoint,wl,fl,bb,iext,icnt,$
wll,idl,iml,wlm,idm,imm
common species,elnames,ionstages
common mspecies,molnames,maxmol,ifacod
common jspecies,coden,fname,nmols
;
; get the y axis minimum and maximum:
ymin=!y.crange(0)
ymax=!y.crange(1)
;
; same for x axis
xmin=!x.crange(0)
xmax=!x.crange(1)
;
; set the position for the label: in the upper 75% 
; of the y-axis range
ystart = labst*ymax
ytext =  txtst*ymax
;
; go through the list and label the plot, 
; if a line above the threshold is found
wllo = 0.d0
idlo = 0
wlmo = 0.d0
idmo = 0
idmcurr = 0
i=0L
icount = 0
endflg = 0
for i=0L,nlines do begin
 if endflg eq 1 then return
 if((i eq nlines) or (xmax le wlm(i))) then begin
  endflg = 1
  if(wlmo le 0.) then begin
   wlmo = wlm(i)
   wlmmin = wlmo
   idmcurr = ifacod(idm(i))
  endif
  goto,l666
 endif
 if imm(i) gt iext(i)+delmol $; or imm(i) gt iml(i) $
   or wlm(i) lt xmin  or wlm(i) gt xmax then goto,l667
 if wlmo eq wlm(i) and idmo eq ifacod(idm(i)) then goto,l667
;
; check if this code is in the list of 'doit' code, ignore it
; if not:
 if(min(where(ifacod(idm(i)) eq docodes)) eq -1) then begin
  ;print,'ignoring code: ',ifacod(idm(i))
  goto,l667
 endif
 plots, [wlm(i),wlm(i)],[ystart,ytext]
l666:
 if(endflg eq 1 or ((ifacod(idm(i)) ne idmcurr) and (idmcurr ne 0))) then begin
; we need to plot the horizontal bar and the label
; for the last series of molecular species to plot 
;
 if(flag ne 0) then begin
  id = min(where(idmcurr eq coden))
  if(id eq -1) then begin
   print,idmcurr,' is not known'
   goto,l889
  endif
  if(id ge 0) then labstr = '!B'+fname(id)
;  for id=0,nmols-1 do begin
;   if(idmcurr eq coden(id)) then begin
;    sid = fname(id)
;    labstr = '!B'+sid
;    ;if(2*(icount/2) eq icount) then labstr = '!B    '+sid
;    goto,l888
;   endif
;  endfor
 endif
  if flag eq 0 then begin
   sid = string(fix(idmcurr))
   labstr = '!B'+sid
   ;if(2*(icount/2) eq icount) then labstr = '!B    '+sid
  endif
l888:
  plots, [wlmmin,wlmo],[ytext,ytext]
  xyouts, 0.5*(wlmo+wlmmin),ytext,labstr,orientation=90
  icount = icount+1
l889:
  wlmmin = wlm(i)
  idmcurr = ifacod(idm(i))
 endif
 wlmo = wlm(i)
 idmo = ifacod(idm(i))
 if(idmcurr eq 0) then begin
  wlmmin = wlmo
  idmcurr = idmo
 endif
l667:
endfor
end
;
;
;
pro molidi,labst,txtst,flag,delmol,docodes
;---------------------------------
; lable a plot with line identifications
; take from common block linelist.
;-- input:
; thres: threshold for line strength for a line
;        to be labeled.
; labst: relative y-coordiate to start label (line), e.g. 0.75
; txtst: relative y-coordiate to start text , e.g. 0.80
; flag: 0: print numeric identification
;       1: print text identification
; delmol: parameter for molecular print threshold. gives
;         difference between imm and iext for lines to plot
common linelist,nlines,maxpoint,wl,fl,bb,iext,icnt,$
wll,idl,iml,wlm,idm,imm
common species,elnames,ionstages
common mspecies,molnames,maxmol,ifacod
common jspecies,coden,fname,nmols
;
; get the y axis minimum and maximum:
ymin=!y.crange(0)
ymax=!y.crange(1)
;
; same for x axis
xmin=!x.crange(0)
xmax=!x.crange(1)
;
; set the position for the label: in the upper 75% 
; of the y-axis range
ystart = labst*ymax
ytext =  txtst*ymax
;
; go through the list and label the plot, 
; if a line above the threshold is found
wllo = 0.d0
idlo = 0
wlmo = 0.d0
idmo = 0
idmcurr = 0
i=0L
icount = 0
endflg = 0
for i=0L,nlines do begin
 if endflg eq 1 then return
 if((i eq nlines) or (xmax le wlm(i))) then begin
  endflg = 1
  if(wlmo le 0.) then begin
   wlmo = wlm(i)
   wlmmin = wlmo
   idmcurr = (idm(i))
  endif
  goto,l666
 endif
 if imm(i) gt iext(i)+delmol $; or imm(i) gt iml(i) $
   or wlm(i) lt xmin  or wlm(i) gt xmax then goto,l667
 if wlmo eq wlm(i) and idmo eq (idm(i)) then goto,l667
;
; check if this code is in the list of 'doit' code, ignore it
; if not:
 if(min(where((idm(i)) eq docodes)) eq -1) then begin
  ;print,'ignoring code: ',(idm(i))
  goto,l667
 endif
 plots, [wlm(i),wlm(i)],[ystart,ytext]
l666:
 if(endflg eq 1 or (((idm(i)) ne idmcurr) and (idmcurr ne 0))) then begin
; we need to plot the horizontal bar and the label
; for the last series of molecular species to plot 
;
 if(flag ne 0) then begin
  ;id = min(where(idmcurr eq coden))
  ;if(id eq -1) then begin
  ; print,idmcurr,' is not known'
  ; goto,l889
  ;endif
  ;if(id ge 0) then labstr = '!B'+fname(id)
  labstr = '!B'+molnames(idmcurr)
;  for id=0,nmols-1 do begin
;   if(idmcurr eq coden(id)) then begin
;    sid = fname(id)
;    labstr = '!B'+sid
;    ;if(2*(icount/2) eq icount) then labstr = '!B    '+sid
;    goto,l888
;   endif
;  endfor
 endif
  if flag eq 0 then begin
   sid = string(fix(idmcurr))
   labstr = '!B'+sid
   ;if(2*(icount/2) eq icount) then labstr = '!B    '+sid
  endif
l888:
  plots, [wlmmin,wlmo],[ytext,ytext]
  xyouts, 0.5*(wlmo+wlmmin),ytext,labstr,orientation=90
  icount = icount+1
l889:
  wlmmin = wlm(i)
  idmcurr = (idm(i))
 endif
 wlmo = wlm(i)
 idmo = (idm(i))
 if(idmcurr eq 0) then begin
  wlmmin = wlmo
  idmcurr = idmo
 endif
l667:
endfor
end

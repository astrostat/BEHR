;+
;behr_test.pro
;	an example program that reads in data from an ascii file containing
;	source and background counts and backscale factors, calls the IDL
;	wrapper routine to BEHR to generate hardness ratios and plots them.
;
;usage
;	infile='behr_test.inp'
;	BEHRdir='./'	;REQUIRED full path to location of BEHR executable
;	psfile='behr_test.ps'	;if you want the output plots in postscript
;	outfile='behr_test.asc'	;if you want the output in an ascii table format
;	algo='quad'		;algo='gibbs' for using the Gibbs sampler,
;				;algo='quad' for Quadrature
;				;algo='auto=20' to choose automatically
;	.run behr_test
;
;	NOTE: additional named input keyword variables for BEHR_HUG()
;	may also be defined prior to running this script.  Type
;		help,behr_hug()
;	to obtain a list of accepted keywords.  See the documentation
;	for BEHR_HUG or the BEHR Tutorial Guide for a description of
;	the variables.
;
;restrictions
;	it is assumed that the input file contains data as an ASCII table,
;	with columns
;	-- soft-band counts in source region
;	-- hard-band counts in source region
;	-- soft-band counts in background region
;	-- hard-band counts in background region
;	-- ratio of background-to-source area for soft band
;	-- ratio of background-to-source area for hard band
;
;requires
;	PINTofALE function BEHR_HUG
;	IDLAstro routine READCOL (and associated)
;
;history
;	vinay kashyap (2005-dec-17)
;	added option to print output to ascii file (VK; 2006-may-04)
;-

;	check if input file is defined
if n_elements(infile) eq 0 then begin
  message,'INFILE is not defined; assuming default of "behr_test.inp"',/informational
  infile='behr_test.inp'
endif

;	check where the BEHR executable is
if n_elements(BEHRdir) eq 0 then begin
  message,'BEHRdir is not defined; please specify directory where executable is'
endif else begin
  fil=findfile(BEHRdir[0]+'/BEHR',count=nfil)
  if nfil eq 0 then begin
    message,BEHRdir[0]+': does not contain BEHR executable'
  endif
endelse

;	check if input file exists
filename=findfile(infile[0],count=nfil)
if nfil eq 0 then message,infile[0]+': file not found'

;	some initializations
if n_elements(verbose) eq 0 then verbose=1

;	read from input file
readcol,filename[0],softsrc,hardsrc,softbkg,hardbkg,softarea,hardarea

;	call the wrapper
ok='ok'
if n_tags(behr) eq 0 then ok='must run behr_hug' else $
 if (tag_names(behr))[0] ne 'SOFTSRC' then ok='structure not in right format' else $
  if n_elements(behr.SOFTSRC) ne n_elements(softsrc) then ok='new file read in'

if ok ne 'ok' then behr=behr_hug(softsrc,hardsrc,softbkg,hardbkg,softarea,hardarea,$
	softeff=softeff,hardeff=hardeff,softidx=softidx,hardidx=hardidx,$
	softscl=softscl,hardscl=hardscl,post=post,$
	level=level,algo=algo,details=details,nsim=nsim,nburnin=nburnin,$
	nbin=nbins,outputf=outputf,outputR=outputR,outputHR=outputHR,$
	outputC=outputC,outputMC=outputMC,outputPr=outputPr,BEHRdir=BEHRdir[0],$
	verbose=verbose)

;	make plots

;	set up display
pmulti=!p.multi & if !d.NAME eq 'X' then !p.multi=[0,2,2]
if keyword_set(psfile) then begin
  !p.multi=pmulti
  setplot=!d.NAME & set_plot,'ps' & device,file=psfile[0],/landscape,/color
endif else begin
  message,'output will be directed to screen, not to PSFILE',/informational
  if !d.N_COLORS gt 256 then device,decomposed=0
endelse
loadct,3

;	logS v/s logH
plot,behr.S.mode,behr.H.mode,psym=4,symsize=2,charsize=1.4,thick=2,$
	xtitle='S',ytitle='H',/xlog,/ylog,xr=[1,3e3],yr=[1,1e3]
for i=0,n_elements(softsrc)-1L do $
	oplot,[behr.S.lowerbound[i],behr.S.upperbound[i]]>1e-10,behr.H.mode[i]*[1,1],color=175
for i=0,n_elements(softsrc)-1L do $
	oplot,behr.S.mode[i]*[1,1],[behr.H.lowerbound[i],behr.H.upperbound[i]]>1e-10,color=175

;	R v/s H+S
plot,behr.R.mode,behr.HpS.mode,psym=4,symsize=2,charsize=1.4,thick=2,$
	xtitle='R=S/H',ytitle='H+S',/ylog,xr=[0,10],yr=[1,1e3]
for i=0,n_elements(softsrc)-1L do $
	oplot,[behr.R.lowerbound[i],behr.R.upperbound[i]]>1e-10,behr.HpS.mode[i]*[1,1],color=175
for i=0,n_elements(softsrc)-1L do $
	oplot,behr.R.mode[i]*[1,1],[behr.HpS.lowerbound[i],behr.HpS.upperbound[i]]>1e-10,color=175

;	C v/s H
plot,behr.C.mode,behr.H.mode,psym=4,symsize=2,charsize=1.4,thick=2,$
	xtitle='C=log!D10!N(S/H)',ytitle='H',/ylog,xr=[-1.5,1.5],yr=[1,1e3]
for i=0,n_elements(softsrc)-1L do $
	oplot,[behr.C.lowerbound[i],behr.C.upperbound[i]],behr.H.mode[i]*[1,1],color=175
for i=0,n_elements(softsrc)-1L do $
	oplot,behr.C.mode[i]*[1,1],[behr.H.lowerbound[i],behr.H.upperbound[i]]>1e-10,color=175

;	HR v/s H+S
plot,behr.HR.mean,behr.HpS.mode,psym=4,symsize=2,charsize=1.4,thick=2,$
	xtitle='HR=(H-S)/(H+S)',ytitle='H+S',/ylog,xr=[-1,1],yr=[1,1e3]
for i=0,n_elements(softsrc)-1L do $
	oplot,[behr.HR.lowerbound[i],behr.HR.upperbound[i]],behr.HpS.mode[i]*[1,1],color=175
for i=0,n_elements(softsrc)-1L do $
	oplot,behr.HR.mean[i]*[1,1],[behr.HpS.lowerbound[i],behr.HpS.upperbound[i]]>1e-10,color=175

!p.multi=pmulti
if keyword_set(psfile) then begin
  device,/close & set_plot,setplot
endif

;	print output to ascii table
if keyword_set(outfile) then begin
  print,'writing S,H,R,C,HR,H+S to file '+outfile[0]+' ...'
  openw,uout,outfile[0],/get_lun

  printf,uout,'S	S.lowerbound	S.upperbound'
  for i=0,n_elements(softsrc)-1L do printf,uout,(behr.S.mode)[i],$
	'	',(behr.S.lowerbound)[i],'	',(behr.S.upperbound)[i]

  printf,uout,''
  printf,uout,'H	H.lowerbound	H.upperbound'
  for i=0,n_elements(softsrc)-1L do printf,uout,(behr.H.mode)[i],$
	'	',(behr.H.lowerbound)[i],'	',(behr.H.upperbound)[i]

  printf,uout,''
  printf,uout,'R	R.lowerbound	R.upperbound'
  for i=0,n_elements(softsrc)-1L do printf,uout,(behr.R.mode)[i],$
	'	',(behr.R.lowerbound)[i],'	',(behr.R.upperbound)[i]

  printf,uout,''
  printf,uout,'C	C.lowerbound	C.upperbound'
  for i=0,n_elements(softsrc)-1L do printf,uout,(behr.C.mode)[i],$
	'	',(behr.C.lowerbound)[i],'	',(behr.C.upperbound)[i]

  printf,uout,''
  printf,uout,'HR	HR.lowerbound	HR.upperbound'
  for i=0,n_elements(softsrc)-1L do printf,uout,(behr.HR.mean)[i],$
	'	',(behr.HR.lowerbound)[i],'	',(behr.HR.upperbound)[i]

  printf,uout,''
  printf,uout,'H+S	H+S.lowerbound	H+S.upperbound'
  for i=0,n_elements(softsrc)-1L do printf,uout,(behr.HpS.mean)[i],$
	'	',(behr.HpS.lowerbound)[i],'	',(behr.HpS.upperbound)[i]

  close,uout & free_lun,uout
  print,'... done.'
endif else begin
  message,'Set OUTFILE to name of file to contain ascii output',/informational
endelse

save,file='behr_test.save'

end


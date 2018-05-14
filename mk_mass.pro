;+
; NAME: mk_mass
;
; PURPOSE: Calculate a mass posterior given an input K magnitude and distance
;          Based on The M_{K_S}-Mass relation from Mann et al. (2018b).
;          
;          http://adsabs.harvard.edu/abs/PAPER
;
; CALLING SEQUENCE:  mass = mk_mass(K,dist,eK,edist,feh=feh,efeh=efeh,
;                                   /oned,/silent)
;
; INPUTS:
;          All input must be scalars
;
;          K:     2MASS Ks-band (magnitudes)
;          dist:  distance (parsecs)
;          ek:    error on K (magnitudes). If not provided assumed to be 0
;          edist: error on distance (parcsecs). If not provided assumed to be 0
;
; OPTIONAL INPUTS:  
;          feh:  Iron abundance. If provided the code will use the
;                Mk-Mass-[Fe/H] relation instead
;          efeh: error on iron abundance (dex). If not provided, but feh
;                is provied, assumed to be 0
;          post: MCMC posterior. Will read it in if not provided. This is
;                useful if you need to run a lot of values (so it does not
;                need to be read in many times). 
;
; KEYWORD PARAMETERS:
;          oned:  Returns a simple 1D error instead of a posterior
;
;          silent: surpress the many (annoying) warnings
;
; OUTPUTS:
;          mass: if /oned is set this is a 2-element array containing the
;                median and standard deviation of the mass posterior. Else
;                this contains the posterior on mass of 1600000
;                elements (the length of the trimmed posterior). 
;
; EXAMPLE (see tester below for more examples):
;    Typical case
;          IDL> k = 8.782
;          IDL> ek = 0.02
;          IDL> dist = 14.55
;          IDL> edist = 0.13
;          IDL> mass = mk_mass(k,dist,ek,edist)
;          IDL> print,median(mass),stdev(mass)
;    If you have a posterior on K and distance instead of 1D errors
;          IDL> post = mrdfits('',/silent)
;          IDL> mass = []
;          IDL> for i = 0,n_elements(k)-1 do mass = [mass,mk_mass(k[i],dist[i],0d0,0d0,post=post)]
;          IDL> cghistoplot,mass,/outline
;
; MODIFICATION HISTORY:
;          May 10 2018: Ported from scratch code by A. Mann
;          May 14 2018: Added testing modules. A. Mann
;         
;
;If you use this code, please cite the relevant paper:
;
;-

function mk_mass,k,dist,ek,edist,feh=feh,efeh=efeh,post=post,silent=silent,oned=oned

  if n_elements(oned) eq 0 then oned = 0
  if n_elements(silent) eq 0 then silent = 0
  if n_elements(ek) eq 0 then begin
     if silent eq 0 then print,'Warning, assuming no error on K, mass errors underestimated'
     ek = 0d0
  endif
  if n_elements(edist) eq 0 then begin
     if silent eq 0 then print,'Warning, assuming no error on distance, mass errors underestimated'
     edist = 0d0
  endif
  if n_elements(efeh) eq 0 and n_elements(feh) gt 0 then begin
     if silent eq 0 then print,'Warning, assuming no error on [Fe/H], mass errors underestimated'
     efeh = 0d0
  endif
  if n_elements(k) ne 1 or n_elements(ek) ne 1 or n_elements(dist) ne 1 or n_Elements(edist) ne 1 then begin
     if silent eq 0 then print,'mk_mass,kmag,ekmag,distance,dist_err; all input must be the same size'
     return,-1
  endif
  if n_elements(feh) gt 0 and (n_elements(feh) ne 1 or n_elements(efeh) ne 1) then begin
     if silent eq 0 then print,'If [Fe/H] is provided as an array, it must have the same size as other parameters'
     return,-1
  endif

  if n_elements(post) eq 0 then begin
     if n_elements(feh) ne 0 then post = mrdfits('resources/Mk-M_5_feh_trim.fits',/silent) else $
        post = mrdfits('resources/Mk-M_5_trim.fits',/silent)
  endif
  ntot = n_elements(post[0,*])
  a = (post[0,*])[*]
  b = (post[1,*])[*]
  c = (post[2,*])[*]
  d = (post[3,*])[*]
  if n_elements(feh) eq 0 then begin
     e = (post[4,*])[*]
     f = 0d0*e
     feh = 0d0
  endif else begin
     e = 0d0*d
     f = (post[4,*])[*]
  endelse
  m = dblarr(n_elements(k))
  kmag = k+ek*randomn(seed,n_elements(a))
  distance = dist + edist*randomn(seed,n_elements(a))
  mk = kmag-5d0*(alog10(distance)-1d0)
  zp = 7.5d0
  mass = (10d0^(a+b*(mk-zp)+c*(mk-zp)^2d0+d*(mk-zp)^3d0+e*(mk-zp)^4d0))*(1d0+feh*f)
  if oned eq 1 then m = [mean(mass),stdev(mass)] else $
     m = mass
  l = where(finite(m) eq 0)
  if l[0] ne -1 then print,'warning, some output is NaN'
  return,m[*]
  
end


;; This is a small code to help test the output of the mk_mass code
;; It just runs a few well-known systems (mostly planet hosts) and
;; checkes how the masses compare to semi-independent determinations
;; (e.g., using models, transit-fit density, etc.)
PRO tester

;; trappist-1
  k = 10.296
  ek = 0.023
  ;;dist = 12.136
  ;;edist = 0.118
  dist = 12.42989539 ;; Gaia DR2 parallax
  edist = 0.018710230
  mass = mk_mass(k,dist,ek,edist)
  print,'Trappist-1:'
  print,'Our mass:'+String(median(mass),format="(D6.3)")+'+/-'+string(stdev(mass),format="(D6.3)")
  print,'Van Grootel et al: 0.089+/-0.006'
  mVG = 0.089+0.006*randomn(seed,n_elements(mass))
  cghistoplot,mass,/outline,thick=3,datacolorname='red',binsize=0.0001,xtitle='Mass (Solar)'
  cghistoplot,mVG,/outline,thick=3,datacolorname='blue',/oplot,binsize=0.0001
  legend,['Our Mass','Van Grootel et al.'],color=cgcolor(['red','blue']),linestyle=0,thick=3,/top,/left,textcolor=cgcolor(['red','blue'])
  print,'Diff = '+string(median(mass-mVG),format="(D6.3)")
  print,'Sig = '+string(mean(mass-mVG)/stdev(mass-mVG),format="(D6.3)")

  ;;stop
  
  ;;GJ1214
  k = 8.782
  ek = 0.02
  dist = 14.55
  edist = 0.13
  feh = 0.30
  efeh = 0.10
  mass = mk_mass(k,dist,ek,edist)
  mass_feh = mk_mass(k,dist,ek,edist,feh=feh,efeh=efeh)
  print,'GJ1214:'
  print,'Our mass:'+String(median(mass),format="(D6.3)")+'+/-'+string(stdev(mass),format="(D6.3)")
  print,'Our mass (with [Fe/H]):'+String(median(mass_feh),format="(D6.3)")+'+/-'+string(stdev(mass_feh),format="(D6.3)")
  print,'Anglada-Escudé: 0.176+/-0.009'
  mAE = 0.176+0.009*randomn(seed,n_elementS(mass))
  cghistoplot,mass,/outline,thick=3,datacolorname='red',binsize=0.0001,xtitle='Mass (Solar)'
  cghistoplot,mass_feh,/outline,thick=3,datacolorname='orange',binsize=0.0001,/oplot
  cghistoplot,mAE,/outline,thick=3,datacolorname='blue',/oplot,binsize=0.0001
  legend,['Our Mass','Our Mass with [Fe/H]','Anglada-Escude et al.'],color=cgcolor(['red','orange','blue']),linestyle=0,thick=3,/top,/right,textcolor=cgcolor(['red','orange','blue'])
  print,'Diff = '+string(median(mass-mAE),format="(D6.3)")
  print,'Sig = '+string(mean(mass-mAE)/stdev(mass-mAE),format="(D6.3)")
  print,'Diff(feh) = '+string(median(mass_feh-mAE),format="(D6.3)")
  print,'Sig(feh = '+string(mean(mass_feh-mAE)/stdev(mass_Feh-mAE),format="(D6.3)")

  ;;stop
  ;; GU Boo
  ;; this is an EB with near equal-mass components, so we can roughly assume DeltaK ~ 0.1
  ;; has a gaia plx
  dist = 1000d0/6.1468
  edist = (0.0159/6.1468)*dist
  ktot = 10.222
  delk = 0.1
  fluxratio = 10d0^(0.05/2.5)
  del_eps = 2.5*alog10(1d0+1d0/fluxratio)
  ka = del_eps+ktot
  kb = ka + delk
  ek = 0.021
  mass_a = mk_mass(ka,dist,ek,edist)
  mass_b = mk_mass(kb,dist,ek,edist)
  print,'GU Boo:'
  print,'Our mass (A):'+String(median(mass_a),format="(D6.3)")+'+/-'+string(stdev(mass_a),format="(D6.3)")
  print,'Our mass (B):'+String(median(mass_b),format="(D6.3)")+'+/-'+string(stdev(mass_b),format="(D6.3)")
  print,'M_A = 0.616+/-0.006 (Lopez-Morales & Ribas (2005))'
  print,'M_B = 0.600+/-0.006 (Lopez-Morales & Ribas (2005))'
  Ma = 0.616+0.006*randomn(seed,n_elementS(mass))
  Mb = 0.600+0.006*randomn(seed,n_elementS(mass))
  cghistoplot,mass_a,/outline,thick=3,datacolorname='red',binsize=0.0001,xtitle='Mass (Solar)'
  cghistoplot,mass_b,/outline,thick=3,datacolorname='orange',binsize=0.0001,/oplot
  cghistoplot,Ma,/outline,thick=3,datacolorname='blue',/oplot,binsize=0.0001
  cghistoplot,Mb,/outline,thick=3,datacolorname='violet',/oplot,binsize=0.0001
  legend,['Our A','Our B','LR05 A','LR05 B'],color=cgcolor(['red','orange','blue','violet']),linestyle=0,thick=3,/top,/right,textcolor=cgcolor(['red','orange','blue','violet'])
  print,'Diff = '+string(median(mass_a-ma),format="(D6.3)")+', '+string(median(mass_b-mb),format="(D6.3)")
  print,'Sig = '+string(mean(mass_a-ma)/stdev(mass_a-ma),format="(D6.3)")+', '+string(mean(mass_b-mb)/stdev(mass_b-mb),format="(D6.3)")


  
  ;; Let's say we have some assymetric posterior on distance, and want a posterior on mass
  ;; note this runs really slowly. Maybe not so practical!
  post = mrdfits('resources/Mk-M_5_trim.fits',/silent)
  tmp = findgen(n_elements(post[0,*]))
  l = wherE(tmp mod 100 eq 1)
  post = post[*,l] ;; runs a bit faster if you trim this down.
  num = 1d3
  logdist = generatearray(1.15,1.3,num)
  dist = 10.0^logdist
  k = 8.0+0.01*randomn(seed,num)
  mass = []
  start = systime(/seconds)
  for i = 0,num-1 do begin
     m = mk_mass(k[i],dist[i],0d0,0d0,post=post)
     ;; select a random 1 out of 100
     m = shuffle(m)
     tmp = findgen(n_elements(m))
     l = where(tmp mod 100 eq 5) ;; just trimming for speed. 
     mass = [mass,m[l]]
     ;;mass = [mass,m]
     ;;if i mod 1000 eq 500 then begin
     ;;   print,(systime(/seconds)-start)/i,n_elements(mass)
     ;;endif
  endfor
  cghistoplot,mass,/outline,thick=3,xtitle='Mass (Solar)' ;; assymetric distribution in mass, as expected
  
END


Function generatearray,min,max,elements
  min *= 1d0
  max *= 1d0
  elements *= 1d0
  array = findgen(elements)
  array = (array*((max-min)/(elements-1d0))) + min
  return,array
END

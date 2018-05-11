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
; EXAMPLE:
;    Create a color table of only greens and plot to the screen
;          IDL>
;    If you have a posterior on K and distance instead of 1D errors
;          IDL> post = mrdfits('',/silent)
;          IDL> mass = []
;          IDL> for i = 0,n_elements(k)-1 do mass = [mass,mk_mass(k[i],dist[i],0d0,0d0,post=post)]
;          IDL> cghistoplot,mass,/outline
;
; MODIFICATION HISTORY:
;          May 10 2018: Ported from scratch code by A. Mann
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
     if n_elements(feh) ne 0 then post = mrdfits('Mk-M_5_feh_trim.fits',/silent) else $
        post = mrdfits('Mk-M_5_trim.fits',/silent)
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
  dist = 12.136
  edist = 0.118
  mass = mk_mass(k,dist,ek,edist)
  cghistoplot,mass,/outline,thick=3
  print,'Trappist-1:'
  print,'Our mass:'+String(median(mass),format="(D6.3)")+'+/-'+string(stdev(mass),format="(D6.3)")
  print,'Van Grootel et al: 0.089+/-0.006'
  mVG = 0.089+0.006*randomn(seed,n_elements(mass))
  cghistoplot,mass-mVG,/outline,thick=3
  print,median(mass-mVG),stdev(mass-mVG),median(mass-mVG)/stdev(mass-mVG)

  ;;GJ1214
  k = 8.782
  ek = 0.02
  dist = 14.55
  edist = 0.13
  mass = mk_mass(k,dist,ek,edist)
  print,median(mass),stdev(mass)
  mAE = 0.176+0.009*randomn(seed,n_elementS(mass))
  cghistoplot,mass-mAE,/outline,thick=3
  print,median(mass-mAE),stdev(mass-mAE),median(mass-mAE)/stdev(mass-mAE)

  ;; Let's say we have some assymetric posterior on distance, and want a posterior on mass
  post = mrdfits('Mk-M_5_trim.fits',/silent)
  num = 1d5
  logdist = generatearray(1.15,1.3,num)
  dist = 10.0^logdist
  k = 8.0+0.01*randomn(seed,num)
  mass = []
  for i = 0,num-1 do mass = [mass,mk_mass(k[i],dist[i],0d0,0d0,post=post)]
  cghistoplot,mass,/outline,thick=3
  
  stop

END


Function generatearray,min,max,elements
  min *= 1d0
  max *= 1d0
  elements *= 1d0
  array = findgen(elements)
  array = (array*((max-min)/(elements-1d0))) + min
  return,array
END

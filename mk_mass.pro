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
;          sige_u: This is for adjusting the value of sigma_e. Probably
;                don't adjust this unless you know what you are
;                doing.
;          fast: Doesn't read in posteriors, just uses 1D errors. Note
;                that this method produces slightly less accurate uncertainties
;
; KEYWORD PARAMETERS:
;          oned:  Returns a simple 1D error instead of a posterior
;
;          silent: surpress the many (annoying) warnings
;
; OUTPUTS:
;          mass: if /oned is set this is a 2-element array containing the
;                median and standard deviation of the mass posterior. Else
;                this contains the posterior on mass of 800000
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
;          Sep 06 2018: Fast mode by A. Mann
;          Aug 20 2018: Allow for user defined sige by A. Mann
;          Jun 21 2018: Added sigma_e, edge issues by A. Mann
;          May 10 2018: Ported from scratch code by A. Mann
;          May 14 2018: Added testing modules. A. Mann
;         
;
;If you use this code, please cite the relevant paper:
; 
;-

function mk_mass,k,dist,ek,edist,feh=feh,efeh=efeh,post=post,silent=silent,oned=oned,sige_u=sige_u,fast=fast

  ;; you might need to adjust this depending on where you install
  ;; this, and if you want to call the code from another location.
  ;;path_to_posteriors = './resources/'
  path_to_posteriors = '~/Dropbox/MMK/resources/'

  if n_elements(oned) eq 0 then oned = 0
  if n_elements(silent) eq 0 then silent = 0
  if n_elements(fast) eq 0 then fast = 0
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
  mk_tmp = k-5d0*(alog10(dist)-1d0)
  if mk_tmp gt 11.6 or mk_tmp lt 3.75 then begin
     if silent eq 0 then print,'Beyond edge of relation'
     return,!Values.F_NAN
  endif

  ;; fast mode!
  if fast eq 1 then begin
     ;; read in the table;; should this just be hard-coded for speed?
     ;; reading files over and over if someone runs this on 100,000
     ;; stars... hrmmm
     readcol,path_to_posteriors+'table.txt',mkinterpol,errinterpol,format='d,d',/silent
     a0 = -0.63770
     a1 = -0.21125
     a2 = -5.0137d-3
     a3 = 8.6686d-3
     a4 = 5.7932d-4
     a5 = -2.7138d-4
     f = 0.0
     feh = 0d0
     tmp = k-5d0*(alog10(dist)-1d0)
     err2 = interpol(errinterpol,mkinterpol,tmp)
     nmonte = 1d5 ;; we should swap this to derivatives from MC errors. It's needlessly slow to generate 1d5 instances!
  endif else begin
     
     if n_elements(post) eq 0 then begin
        if n_elements(feh) ne 0 then post = mrdfits(path_to_posteriors+'Mk-M_8_feh_trim.fits',/silent) else $
           post = mrdfits(path_to_posteriors+'Mk-M_7_trim.fits',/silent)
     endif
     ntot = n_elements(post[0,*])
     a0 = (post[0,*])[*]
     a1 = (post[1,*])[*]
     a2 = (post[2,*])[*]
     a3 = (post[3,*])[*]
     a4 = (post[4,*])[*]
     a5 = (post[5,*])[*]
     nmonte = ntot
     if n_elements(feh) eq 0 then begin
        sige = exp((post[6,*])[*])
        f = 0d0*a0
        feh = 0d0
     endif else begin
        f = (post[6,*])[*]
        sige = exp((post[7,*])[*]) ;; 
     endelse
     if n_elements(sige_u) eq 1 then sige = sige_u ;; user defined sig_e value. This is for testing what happens when we fiddle with this (e.g., how does chi^2 change). Maybe if you are skeptical of our ~2% errors or you are working on young stars, etc.
  endelse

  m = dblarr(n_elements(k))
  kmag = k+ek*randomn(seed,nmonte)
  distance = dist + edist*randomn(seed,nmonte)
  mk = kmag-5d0*(alog10(distance)-1d0)

  ;; output a warning if posterior over the edges of the relation
  ll = where(mk lt 4.0 or mk gt 11)
  if n_elements(ll) gt 1 then print,'Warning, '+strtrim(string(100*(n_elements(ll)*1d0)/(n_elements(mk)*1d0),format="(D6.2)"),2)+'% of posterior beyond relation edge.'
  ;;if (mean(mk) gt 10.5 or mean(mk) lt 4.5) and silent eq 0 then print,'Warning, near the edges of the relation'
  
  zp = 7.5d0
  mass = (10d0^(a0+a1*(mk-zp)+a2*(mk-zp)^2d0+a3*(mk-zp)^3d0+a4*(mk-zp)^4d0+a5*(mk-zp)^5d0))*(1d0+feh*f)
  if fast eq 1 then begin
     mass += mass*(err2/100d0)*randomn(seed,nmonte)
  endif else begin
     mass += median(sige)*mass*randomn(seed,nmonte)
  endelse
  if oned eq 1 then m = [mean(mass),stdev(mass)] else $
     m = mass
  l = where(finite(m) eq 0)
  bad = where(mk gt 11.5 or mk lt 4.0)
  if bad[0] ne -1 then mass[bad] = !Values.F_NAN
  if l[0] ne -1 and silent eq 0 then print,'warning, some output is NaN (points beyond edge of relation?)'
  return,m[*]
  
end


;; This is a small code to help test the output of the mk_mass code
;; It just runs a few well-known systems and
;; checkes how the masses compare to semi-independent determinations
;; (e.g., using models, transit-fit density, etc.)
PRO tester

  ;; trappist-1
  k = 10.296
  ek = 0.023
  dist = 12.42989539 ;; Gaia DR2 parallax
  edist = 0.018710230
  mass = mk_mass(k,dist,ek,edist)
  print,'Trappist-1:'
  print,'Our mass:'+String(median(mass),format="(D6.3)")+'+/-'+string(stdev(mass),format="(D6.3)")+' ('+strtrim(string(100*stdev(mass)/median(mass),format="(D4.1)"),2)+'% error)'
  print,'Van Grootel et al: 0.089+/-0.006'
  mVG = 0.089+0.006*randomn(seed,n_elements(mass))
  cghistoplot,mass,/outline,thick=3,datacolorname='red',binsize=0.0001,xtitle='Mass (Solar)'
  cghistoplot,mVG,/outline,thick=3,datacolorname='blue',/oplot,binsize=0.0001
  legend,['Our Mass','Van Grootel et al.'],color=cgcolor(['red','blue']),linestyle=0,thick=3,/top,/left,textcolor=cgcolor(['red','blue'])
  print,'Diff = '+string(median(mass-mVG),format="(D6.3)")
  print,'Sig = '+string(mean(mass-mVG)/stdev(mass-mVG),format="(D6.3)")

  
  print,'--------'
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
  print,'Diff(with feh) = '+string(median(mass_feh-mAE),format="(D6.3)")
  print,'Sig(with feh) = '+string(mean(mass_feh-mAE)/stdev(mass_Feh-mAE),format="(D6.3)")

  print,'--------'
  ;;stop
  ;; GU Boo
  ;; this is an EB with near equal-mass components, so we can roughly
  ;; assume DeltaK = 0.13 (adopted from optical flux ratio)
  ;; has a gaia plx
  dist = 1000d0/6.1468
  edist = (0.0159/6.1468)*dist
  ktot = 10.222
  delk = 0.13
  fluxratio = 10d0^(delk/2.5)
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
  
  print,'--------'
  ;; Let's say we have some assymetric posterior on distance, and want a posterior on mass
  ;; note this runs really slowly. Maybe not so practical!
  post = mrdfits('resources/Mk-M_7_trim.fits',/silent)
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
  cghistoplot,mass,/outline,thick=3,xtitle='Mass (Solar)' ;; assymetric (non-Gaussian) distribution in mass, as expected
  ;; faster way is to do this with the /fast option, but I wrote the example before writing that fix...
  
END


Function generatearray,min,max,elements
  min *= 1d0
  max *= 1d0
  elements *= 1d0
  array = findgen(elements)
  array = (array*((max-min)/(elements-1d0))) + min
  return,array
END

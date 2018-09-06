PRO mmk_tester

  k = 9.0
  ek = 0.02
  dist = 25.
  edist = 1.0

  noerr = mk_mass(k,dist,/oned)

  mass = mk_mass(k,dist,ek,edist,/oned)
  print,mass

  arr = [-0.64661,-0.21246,-2.6534d-3,7.9485d-3,3.6899d-4,-1.9226d-4]
  mk = (k-5.0*(alog10(dist)-1.0))-7.5
  emk = sqrt(((2.17157/dist)*edist)^2. + ek^2.)

  ;tmp1 = k+ek*randomn(seed,1d5)
  ;tmp2 = dist+edist*randomn(seed,1d5)
  ;tmp3 = tmp1-5.0*(alog10(tmp2)-1.0)-7.5
  ;print,emk,stdev(tmp3)
  
  ms = 10.0^(arr[0]+arr[1]*mk+arr[2]*mk^2.+arr[3]*mk^3.+arr[4]*mk^4.+arr[5]*mk^5.)
  ems = sqrt((alog(10)*(arr[1]+arr[2]*2*mk+arr[3]*3*mk^2.+arr[4]*4*mk^3.+arr[5]*5*mk^4.)*ms*emk)^2.)

  mk = mk+emk*randomn(Seed,1d5)
  ms_arr = 10.0^(arr[0]+arr[1]*mk+arr[2]*mk^2.+arr[3]*mk^3.+arr[4]*mk^4.+arr[5]*mk^5.)
  print,stdev(ms_arr)

  ;; final error is the two errors in quad:
  finerr = sqrt(stdev(ms_arr)^2.+(0.022*ms)^2.)
  print,mass[1] ,finerr
  stop

END

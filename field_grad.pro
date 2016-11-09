

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; A routine to determine the gradient of a scalar.
; If grid sizes are supplied, they are used everywhere.
; If a coordinate grid is supplied, the grid size is determined
; locally.
; velocities are used to find the upwind direction.
; uxf is the xcomponent of velocity, divided by the local courant
; speed. 
; uwf is the upwind fraction, it defines how heavily to weight toward
; the upwind end. 

pro field_grad, source, dx_source, dy_source, dx=dx, dy=dy, xx=xx, yy=yy, uxb=uxb, uyb=uyb, uw_diff=uw_diff, rbc=rbc, lbc=lbc, tbc=tbc, bbc=bbc

; We'll default to smooth gradients.
if keyword_set(tbc) eq 0 then tbc = 'open'
if keyword_set(bbc) eq 0 then bbc = 'open'
if keyword_set(rbc) eq 0 then rbc = 'open'
if keyword_set(lbc) eq 0 then lbc = 'open'

; First we note that the source may be either 2D or 3D. To accomodate
; the shift action, we need two cases for the argument.

   sl = [-1,0]
   sr = [ 1,0]
   sd = [0,-1]
   su = [0, 1]

if n_elements(size(source, /dim)) ne 2 then begin
   print, 'input array must be of 2 dimensions'
   stop
endif

; First we'll do the simple center difference method

if keyword_set(uw_diff) eq 0 then begin

   ; First we'll set up the coordinate grid

   if keyword_set(dx) eq 0 then begin
      if keyword_set(xx) eq 1 then begin
         if keyword_set(dx) then print, 'dx cannot be defined simultaneously with the coordinate grid - coordinate definition used instead.'
         dx = (shift(xx, sl) - shift(xx, sr))/2.
      endif else dx = 1.
   endif

   if keyword_set(dy) eq 0 then begin
      if keyword_set(yy) eq 1 then begin
         if keyword_set(dy) then print, 'dy cannot be defined simultaneously with the coordinate grid - coordinate definition used instead.'
         dy = (shift(yy, sd) - shift(yy, su))/2.
      endif else dy = 1.
   endif

   ; And now we'll find the delta values

   dx_source = (shift(source, sl) - shift(source, sr)) / (2. * dx)
   dy_source = (shift(source, sd) - shift(source, su)) / (2. * dy)

endif else begin
; Now we perform our upwind differencing.
  
   ; First we check for the necessary velocity

   if n_elements(uxb) eq 0 or n_elements(uyb) eq 0 then begin
      print, 'cannot perform upwind difference without a supplied velocity direction'
      stop
   endif

   ; Then we'll ensure that the supplied velocity is unitary
   if min(abs(uxb)) ne 1 or min(abs(uyb)) ne 1 then begin
      print, 'velocity input boolean must be unitary'
      stop
   endif

   ; Now we determine the grid deltas on each side

   if keyword_set(xx) eq 1 then begin
      dxl = xx - shift(xx, sr)
      dxr = shift(xx, sl) - xx
   endif else begin
      if keyword_set(dx) eq 1 then begin
         dxl = dx
         dxr = dx
      endif else begin
         dxl = 1.
         dxr = 1.
      endelse
   endelse

   if keyword_set(yy) eq 1 then begin
      dyb = yy - shift(yy, su)
      dyt = shift(yy, sd) - yy
   endif else begin
      if keyword_set(dy) eq 1 then begin
         dyb = dy
         dyt = dy
      endif else begin
         dyb = 1.
         dyt = 1.
      endelse
   endelse

   ; Now we'll determine the four derivatives
   
   dxl_source = ( source - shift(source, sr) ) / dxl
   dxr_source = ( shift(source, sl) - source ) / dxr
   dyb_source = ( source - shift(source, su) ) / dyb
   dyt_source = ( shift(source, sd) - source ) / dyt

   ; Now we need our 'choice boolean'

                                ; If the wind blows from the left then
                                ; the x velocity is positive and we
                                ; choose the left derivative. 

 
                                ; If we're doing full upwind
                                ; differencing then we only
                                ; care about the direction of the
                                ; velocity


   dx_source = ( (1 + uxb) * dxl_source + (1 - uxb) * dxr_source ) / 2.
   dy_source = ( (1 + uyb) * dyb_source + (1 - uyb) * dyt_source ) / 2.

endelse

if n_elements(size(source, /dim)) eq 3 then begin

   if lbc eq 'open' then dx_source[ 0,*,*] = dx_source[ 1,*,*]
   if rbc eq 'open' then dx_source[-1,*,*] = dx_source[-2,*,*]
   if bbc eq 'open' then dy_source[*, 0,*] = dy_source[*, 1,*]
   if tbc eq 'open' then dy_source[*,-1,*] = dy_source[*,-2,*]

   if lbc eq 'zero' then dx_source[ 0,*,*] = 0
   if rbc eq 'zero' then dx_source[-1,*,*] = 0
   if bbc eq 'zero' then dy_source[*, 0,*] = 0
   if tbc eq 'zero' then dy_source[*,-1,*] = 0

endif else begin

   if lbc eq 'open' then dx_source[ 0,*] = dx_source[ 1,*]
   if rbc eq 'open' then dx_source[-1,*] = dx_source[-2,*]
   if bbc eq 'open' then dy_source[*, 0] = dy_source[*, 1]
   if tbc eq 'open' then dy_source[*,-1] = dy_source[*,-2]

   if lbc eq 'zero' then dx_source[ 0,*] = 0
   if rbc eq 'zero' then dx_source[-1,*] = 0
   if bbc eq 'zero' then dy_source[*, 0] = 0
   if tbc eq 'zero' then dy_source[*,-1] = 0

endelse

end

   

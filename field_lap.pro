

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; This routine finds the laplacian.
; It is inherently centered
; We'll naively assume uniform grid spacing. If this becomes a
; problem we can fix it later.

pro field_lap, source, lap_source, dx=dx, dy=dy, xx=xx, yy=yy, rbc=rbc, lbc=lbc, bbc=bbc, tbc=tbc, x_only=x_only, y_only=y_only

if keyword_set(tbc) eq 0 then tbc = 'open'
if keyword_set(bbc) eq 0 then bbc = 'open'
if keyword_set(lbc) eq 0 then lbc = 'open'
if keyword_set(rbc) eq 0 then rbc = 'open'

if size(source, /n_dimensions) ne 2 then begin
   print, 'input array must be of 2 dimensions'
   stop
endif

sl = intarr(size(source, /n_dimensions))
sr = sl
sd = sl
su = sl

sl[0] = -1
sr[0] =  1
su[1] =  1
sd[1] = -1

if keyword_set(xx) then begin
   if keyword_set(dx) then print, 'Cannot specify coordinate array and differential simultaneously. Using coordinate information.'
   dxl = xx - shift(xx, sr)
   dxr = shift(xx, sl) - xx
endif else begin
   if keyword_set(dx) eq 0 then dx = 1.
   dxr = dx
   dxl = dx
endelse

if keyword_set(yy) then begin
   if keyword_set(dy) then print, 'Cannot specify coordinate array and differential simultaneously. Using coordinate information.'
   dyt = shift(yy, sd) - yy
   dyb = yy - shift(yy, su)
endif else begin
   if keyword_set(dy) eq 0 then dy = 1.
   dyt = dy
   dyb = dy
endelse

if keyword_set(x_only) then x_only = 1. else x_only = 0.
if keyword_set(y_only) then y_only = 1. else y_only = 0.

xlap = ( ( ( shift(source, sl) - source ) / dxr ) - ( ( source - shift(source, sr) ) / dxl ) ) / ( 0.5 * ( dxl + dxr ) )
ylap = ( ( ( shift(source, sd) - source ) / dyt ) - ( ( source - shift(source, su) ) / dyb ) ) / ( 0.5 * ( dyt + dyb ) ) 

lap_source = (1. - y_only) * xlap + (1. - x_only) * ylap

if tbc eq 'open' then lap_source[*,-1,*] = lap_source[*,-2,*]
if bbc eq 'open' then lap_source[*, 0,*] = lap_source[*, 1,*]
if lbc eq 'open' then lap_source[ 0,*,*] = lap_source[ 1,*,*]
if rbc eq 'open' then lap_source[-1,*,*] = lap_source[-2,*,*]

if tbc eq 'zero' then lap_source[*,-1,*] = 0
if bbc eq 'zero' then lap_source[*, 0,*] = 0
if lbc eq 'zero' then lap_source[ 0,*,*] = 0
if rbc eq 'zero' then lap_source[-1,*,*] = 0

end





















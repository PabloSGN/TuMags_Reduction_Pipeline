FUNCTION Shc, d1, d2, FILTER=filt, INTERPOLATE=int_max, N2=n2, $
              POLY=pfit, RANGE=range
;+
; NAME:
;       SHC
; PURPOSE:
;       Find the linear shift between two 2-d images using the
;       forier crosscorrelation method and optional interpolation for
;       sub-pixel shifts.
; CATEGORY:
;
; CALLING SEQUENCE:
;       RESULT = SHC ( IMG1, IMG2 )
; INPUTS:
;       IMG1,IMG2 : 2-d arrays of same size
; KEYWORDS:
;       FILTER : (Flag) if set and not zero, an edge filter is used
;                for the fft.
;   INTERPOLATE: (Flag) If set and not zero, interpolate the position
;                of the maximum of the cross correlation function to
;                sub-pixel accuracy.
;            N2: (Flag) Only use a subarea of 2^n1+2^n2 for the FFT to
;                speed things up
;          POLY: (int) subtract a polynomial fit of degree POLY from
;                the data before doing the correlation.
;         RANGE: (int) constraing the possible solution to be lower
;                than the given value.
; OUTPUTS:
;       RESULT : 2-element vector with the shift of IMG2 relative to
;                IMG1. RESULT(0) is shift in X-direction, RESULT(1) is
;                shift in Y-direction.
; PROCEDURE:
;       Compute the crosscorrelation function using FFT and locate the
;       maximum.
; MODIFICATION HISTORY:
;       13-Aug-1992  P.Suetterlin, KIS
;       29-Aug-1995  PS Added Edge filter
;       30-Aug-1995  PS Add Subpix interpolation. Slight rewrite of
;                       normal maximum finding (Use shift)
;       11-Sep-1995  PS Forgot the 1-d case. Re-implemented.
;       10-Sep-1996  PS Auto-Crop datasets to same dimensions.
;                       Also trim uneven dimensions, it's ugly for FFT
;       08-Jan-2001  PS Keyword N2 to use good powers of 2
;       22-Nov-2001  Subtract of polynomial surface.  Some formatting.
;                    add on_error return.
;-

on_error, 2

p1 = reform(d1-avg(d1)) &  s1=size(p1)
p2 = reform(d2-avg(d2)) &  s2=size(p2)
  ;;; find common size
sx = s1(1) < s2(1) & sy=s1(2) < s2(2)
  ;;; make them even
sx = 2*(sx/2)
sy = 2*(sy/2)

IF (s1(0) GT 2) OR (s1(0) EQ 0) OR (s2(0) EQ 0) OR (s2(0) EQ 0) THEN $
  message, 'Only 1-d or 2-d Data!'

  ;;; use a fft-friendly 2^n size?
IF keyword_set(n2) THEN BEGIN
      ;;; use 2^n1+2^n2, but n2 only if it contributes a noticeable
      ;;; part to the size.  First for X-axis
    l = 2^(fix(alog10(sx)/0.30103)-1)
    IF (sx-2*l) NE 0 THEN BEGIN
        ll = 2^(fix(alog10(sx-2*l)/0.30103)-1)
        IF (ll NE 0) THEN IF l/ll LT 16 THEN $
          sx=2*(l+ll) $
        ELSE $
          sx=2*l
    ENDIF
      ;;; same for 2nd dimension
    l = 2^(fix(alog10(sy)/0.30103)-1)
    IF (sy-2*l) NE 0 THEN BEGIN
        ll = 2^(fix(alog10(sy-2*l)/0.30103)-1)
        IF (ll NE 0) THEN IF l/ll LT 16 THEN $
          sy=2*(l+ll) $
        ELSE $
          sy=l
    ENDIF
ENDIF
  ;;; clip to common size
p1 = p1(0:sx-1, *)
p2 = p2(0:sx-1, *)
IF s1(0) NE 1 THEN BEGIN
    p1 = p1(*, 0:sy-1)
    p2 = p2(*, 0:sy-1)
ENDIF
  ;;; is a removal of polynomial fits requested?
IF keyword_set(pfit) THEN BEGIN
    IF s1(0) EQ 1 THEN BEGIN
        x = findgen(sx)
        c = poly_fit(x, p1, pfit, fit)
        p1 = p1-fit
        c = poly_fit(x, p2, pfit, fit)
        p2 = p2-fit
    ENDIF ELSE BEGIN
        p1 = p1-sfit(p1, pfit)
        p2 = p2-sfit(p2, pfit)
    ENDELSE
ENDIF

IF s1(0) EQ 1 THEN GOTO, onedim

  ;;; Filtering the data?  Use exponential filter
IF keyword_set(filt) THEN BEGIN
    x = findgen(sx)
    x = exp(-(shift(x < (sx-x), sx/2)/(sx/2))^2)
    y = findgen(sy)
    y = exp(-(shift(y < (sy-y), sy/2)/(sy/2))^2)
    mm = x#y
    p1 = p1*mm
    p2 = p2*mm
ENDIF

  ;;; compute cros correlation using FFT, center it

cc = shift(abs(fft(fft(p1, -1)*conj(fft(p2, -1)), 1)), sx/2, sy/2)

  ;;; range clipping?  Only use central square

IF keyword_set(range) THEN cc = cc((sx/2-range) > 0:((sx/2+range) < sx)-1, $
                                   (sy/2-range) > 0:((sy/2+range) < sy)-1)

mx = max(cc, loc)
 ;;; Simple Maximum location
ccsz = size (cc)
xmax = loc mod ccsz(1)
ymax = loc/ccsz(1)

IF NOT keyword_set(int_max) THEN GOTO, Ende

  ;;; if requested: linear iterpolation of the exact position
IF (xmax*ymax GT 0) AND (xmax LT (ccsz(1)-1)) $
  AND (ymax LT (ccsz(2)-1)) THEN BEGIN
      ;;; Sep 91 phw try including more points in interpolations
    denom = mx*2 - cc(xmax-1, ymax) - cc(xmax+1, ymax)
    xfra = (xmax-.5) + (mx-cc(xmax-1, ymax))/denom
    denom = mx*2 - cc(xmax, ymax-1) - cc(xmax, ymax+1)
    yfra = (ymax-.5) + (mx-cc(xmax, ymax-1))/denom
    xmax = xfra
    ymax = yfra
ENDIF

Ende:

return, [xmax-ccsz(1)/2, ymax-ccsz(2)/2]

Onedim:
IF keyword_set(filt) THEN BEGIN
    x = findgen(sx)
    x = exp(-(shift(x < (sx-x), sx/2)/(sx/2))^2)
    p1 = p1*x
    p2 = p2*x
ENDIF

cc = shift(abs(fft(fft(p1, -1)*conj(fft(p2, -1)), 1)), sx/2)
mx = max(cc(1:sx-2), xmax)
xmax = xmax+1

IF keyword_set(int_max) THEN BEGIN
     ;;; Polyfit of degree 2 for three points, extremum
    c1 = (cc(xmax+1)-cc(xmax-1))/2.
    c2 = cc(xmax+1)-c1-cc(xmax)
    xmax = xmax-c1/c2/2.
ENDIF

return, xmax-sx/2

END

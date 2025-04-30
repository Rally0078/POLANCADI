function w_avg, inarr, yrange, thresh     ;weighted averaging

sumw = 0
sumwh = 0
sumi = 0
for i = 0,yrange - 1 do begin
  wt = abs(inarr(i))
  if wt gt thresh then begin       ;only accept data above threshold
    wtt = 10.0^(wt/70.0)         ;convert back from dB to power (db values were x7 for use by colorscale)
    wtt = wtt^3               ;weight peaks more heavily
    sumw = sumw + wtt
    sumwh = sumwh + i*wtt
    sumi = sumi + 1
  endif
endfor
;stop
if sumi ne 0 then return, sumwh/sumw else return, 0
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro realhtanalysis, dir, site, hhmm_str, hour, minute, farr, harr, freqs, heights
   dxarr = farr
   dyarr = harr
   close, /all       ;make sure file #2 is close, could still be open due to previous error

 ;write temporary file to polan directory for use by polrun
      openw,2,'c:\polan\in.dat'
      printf,2, 'OUTPUT MODE ==>          -9.00  0.0  0.0  0.0    0'
     ; printf,2, 'Date = '+ strmid(dir,3,8)+ '           1.60 85.0  0.0 0.00    0'
      printf,2, 'Date = '+ strmid(dir,3,8)+ '          '+$
             string(format = '(f5.2,f5.1)',site.gyrofreq,site.magdip) +'  0.0 0.00    0'
      line = hhmm_str + '                    0.0  '
      printf,2,line
      line = ''
        for n = 0, n_elements(dxarr)-1 do begin
          line = line + string(format='(f6.1,(", ",f8.2))',dxarr(n),dyarr(n))
          printf,2,line
          line = ''
        endfor
        if strlen(line) gt 0 then printf,2,line
      close,2      ;data in temporary file
	  cd, "c:\polan", CURRENT=owd
      spawn, 'c:\polan\polan.exe'         ;run polan
	  cd, owd
            ;next read back the real height calculated by polan
      textin = '  '                         ;establish textin as a text variable
     ; openr, 2, 'c:\polan\polout.t', ERROR = err     ;open file to which polan writes its output

      openr, 2, 'c:\polan\out.dat', ERROR = err     ;open file to which polan writes its output
       ;openr, 2, 'c:cdata\polout.t'
       if (err ne 0) then begin
          print, 'unable to open file c:\polan\out.dat'
       goto, nofile
       endif
           ;if file then read through the file until find line 'Real Heights'
      repeat begin
        readf,2, textin
        ;print, textin
      endrep until textin eq 'Real Heights' or eof(2)
        if eof(2) then begin
          print, 'no real heights'
          close, 2
          return
        endif
           ;ok, so now read real heights from the following lines
      n = 0                            ;counter for real height arrays
      freqs = fltarr(200)                ;array to hold real height frequency
      heights = fltarr(200)                  ;array to hold real height height
      repeat begin
        readf,2,textin                    ;read a line from the file
        ;stop
        lengthtextin = strlen(textin)
        noofvalues= floor(lengthtextin/14)       ;values in file are in format f8.2,f6.1??
        for i = 0,noofvalues-1 do begin
          freqt = strmid(textin,i*14+1,6)
          heightt = strmid(textin,i*14+7,6)
          freq = float(freqt)
          height = float(heightt)
          if (( n eq 0) or ((freq gt 0.8) and (height gt 80))) then begin     ;polan writes first freq as 0.5 MHz and puts low freq at end
            freqs(n) = freq
            heights(n) = height
            n = n+1
          endif ;else goto, endofrthts
         ;print, freq,height
        endfor
        ;stop
      endrep until eof(2) eq !true or textin eq ''  ;eof or blank line
      endofrthts:
      ;heights = heights(0:n-1)               ;trim real height array for plotting use
      ;freqs = freqs(0:n-1)                  ;trim freq array for plotting use
      close, 2                         ;all done reading from file so close
      ;so have suitable polout file -copy it to cdata dir with a suitable name
      months = ['a','b','c','d','e','f','g','h','i','j','k','l']
      if hour lt 10 then hrstr = '0' + string(format = '(i1)',hour) else hrstr = string(format = '(i2)',hour)
      if minute lt 10 then minstr = '0' + string(format = '(i1)',minute) else minstr = string(format = '(i2)',minute)
     ; ccommand = 'copy c:\rsi\idl52\polout.t c:\cdata\'+strmid(dir,4,1) + months(fix(strmid(dir,5,2))-1)+$
      ccommand = 'copy polout.t c:\cdata\'+strmid(dir,4,1) + months(fix(strmid(dir,5,2))-1)+$
                    strmid(dir,7,2) + hrstr + minstr + '.pol'
      spawn, ccommand
    nofile:                             ;jump here if unable to open file polout.t
end


;+
; NAME:
;       CADI_IGWSCALE
; PURPOSE:
;       Generate summary plots of ionograms from CADI data files.
;       Designed for MD and FFI data formats.
; CALLING SEQUENCE:
;       CADI_IGWSUMM, DIR=DIR, DATESTR=DATESTR, EXT=EXT
; INPUTS:
;       DIR = Scalar string containing directory in which data is located.
;          Syntax ready to prepend onto filename e.g. 'D:\CADIDATA\'.
;          Default is ''.
;       DATESTR = string holding date of data as 6 ASCII digits YYMMDD.
;          Default is to take date from header of first .ext file in directory.
;       EXT = scalar string containing extension of data files, without period.
;          Default is 'MD0'.
;     THOUR = gives starting hour for the ionograms (if cnatimein it is derivedfrom the cantime)
;     HEIGHTRANGE = 2 element array giving min height and max height
;     ELDENSITYRANGE = 2 element array giving range if x scale is electron density (not implemented in cewidgets)
;     NOOFPANELS, NOOFCOLUMNS< NOOFROWS: gives geometry of the plotting
;     CANTIMEIN = 2 element array of type Cantime giving beginning and end of plottimes
;     COLOUR = if set, plots in color
; KEYWORD PARAMETERS:
;       INFORM: If present and nonzero, progress messages are printed.
;       ERROR = integer scalar. Returned as true if an error occurred, else
;          false.
; SIDE EFFECTS:
;       Generates a plot file.
; MODIFICATION HISTORY:
;       JWM 4 April 2004 program is based on Ian Grant's cadi_igsumm.
;-

pro cadi_igwscale, dir=dir, datestr=datestr_in, ext=file_ext, $
   inform=inform, error=error, heightrange = heightrange,$
    cantimein = cantimein,  colour = colour, heightfrequencyscale = heightfrequencyscale

on_error, 0 ;2

psav = !p
xsav = !x
ysav = !y

;if keyword_set(heightfrequencyscale) then DEVICE, CURSOR_STANDARD = 32645 else DEVICE, CURSOR_STANDARD = 32515

autoscaledone = !false
if not keyword_set(file_ext) then print, 'no file extension spaecified'
if not keyword_set(dir) then dir = ''

; Get header information (CADI_READ_DIR get it from first file)

h = cadi_read_dir(dir=dir, ext=file_ext, /onlyheader, error=error )
if error then begin
   error_message = 'Error reading first file for header information'
   goto, error_return
   endif

; Input date from keyword or from first file

if not keyword_set(cantimein) then begin
  if not keyword_set(datestr_in) then $
   candate = canmidnight(h.time) $   ; Date from time in first file
  else begin
   datestr = datestr_in ; Unnecessary?
   reads,datestr,year,month,day,format='(3i2)'       ; Date from keyword
   if year lt 54 then year = year + 100        ;y2k adjustment
   ;print, 'doing timecalc'
   candate=cantime(1900+year,month,day,0,0,0)
   ;print, 'done'
   endelse
; Get DAYNUM for plot title

  civil_time, year, month, day, daynum, hour, minute, second, can_time = candate

endif else begin   ;yes there is cantimein
     cantimeset = cantimein
     civil_time, year, month, day, daynum, starthour, minute, second, cantimein(0)
     sthr = starthour
     noofhours = (cantimeset(1) - cantimeset(0))/3600L
;print, 'noofhours = ',noofhours
     ;endhr = sthr + noofhours - 1
endelse
;stop
rx = 0L         ; Index of receiver to be plotted
mhzrange = [1.0,20.0]
cursorfreq = 4.00
cursorfreq1 = 8.0  ; For fbEs
cursorheight = 250
cursorheight1 = 300; For range spread F
hpheight = 250
spreadf=0
;if keyword_set(spf3) then spreadf=1 print, 'spreadf value',spreadf
;if not keyword_set(spf3) then spreadf=0
tdiff = h.dtime/60   ;tdiff is  (for scaling) just the time increment (in minutes) between ionograms
;print, 'tdiff = ',tdiff
!p.noerase = 1
!p.ticklen = 0.04
charsize = 0.8
;stop
;reads,esf,format='(i2)'       ; spread F from keyboard
;print, 'captured spred F'

def_ticknames = replicate('',20)

   !p.background = 255b     ;white
   back_color = 255b
   !p.color = 0          ;black
   fore_color = 0b
erase

prev_timex = -1
prev_panelx = -1

htmin = 0
if keyword_set(heightrange) then htmin = heightrange(0)
if dir eq "md4" then htmax = 1020 else htmax = 510
if keyword_set(heightrange) then htmax = heightrange(1)
;stop

cadi = cadi_read_dir(dir=dir, timerange = cantimeset,$ ;timerange=timerange_hour, $
    ext=file_ext, tfh= 0, error=read_error, inf=keyword_set(inform), htlimitmax=htmax, htlimitmin = htmin )
if read_error then begin
    message, /inform, $
      'Error reading data in time range '+strcanrange(timerange_hour)
    goto, end_hourx_loop
endif

if not cadi.header_exists then goto, end_hourx_loop      ; No soundings for this hour

   h = cadi.header
   d = cadi.data

site = cadi_site(h.site, h.time)
interferometerconstEW = 3.0e8/(2.0*!pi*site.separation(0))
interferometerconstNS = 3.0e8/(2.0*!pi*site.separation(1))

nrcvrs = n_elements(d(0).IQ(0,*))
print, 'no of receivers = ',nrcvrs
  ; stop
   dh = median(difference( h.heights )) ; Typical height step; pixel height
dxarr = fltarr(1)     ;define arrays which will be used for scaled data
dyarr = fltarr(1)
;start plotting
selOmode = !true
selXmode = !true
enhanceoverhead = !false
datafileopen = !false
datascanned = !false
usersym, [-1.0,1.0,1.0,-1.0,-1.0],[-1.0,-1.0,1.0,1.0,-1.0],/fill
htimex = 0
boxcolor = 230

startplotpanels:
  erase
  civil_time, year, month, day, daynum, hour, minute, second, h.times(htimex)
  timex = round( (60*(hour-sthr)+minute) / tdiff )   ;get index of ionogram plot
  hrange = ends(h.heights)
  if keyword_set(heightrange) then hrange = heightrange
  if hrange(1) gt htmax then hrange(1) = htmax

;put on db scale
  yticklabels = ['0','6','12','18','24','30','36']
  plot, [0,1],[0,1],/nodata, yrange = [0,256], ystyle = 1, xrange = [0,1],$
     position = [0.97,0.18,0.99,0.85], xstyle = 4, yticks = 6, ytickname = yticklabels, title = 'dB'
  for c = 0,255 do oplot, [0,1],[c,c],color = 255 - c, thick = 3

  ;put on boxes for use by mouse
  plot, [0],[0],/nodata,/noerase, position = [0.05,0.05,0.99,0.15], xrange = [0,100], xstyle = 5,$
       yrange = [0,10], ystyle = 5
  oplot,[5],[2],psym = 8, symsize = 3, color = 240
  xyouts,  0,2, 'Done', color = 240
  oplot,[16, 25],[2, 2],psym = 8, symsize = 3, color = 50
  xyouts, 8,4, 'step:', color = 50
  xyouts, 8,2, 'forward 10', color = 50
  xyouts, 18,2, 'forward 1', color = 50
  oplot,[32.5, 40.5],[2, 2],psym = 8, symsize = 3, color = 130
  xyouts, 27,2, 'back 1', color = 130
  xyouts, 34,2, 'back 10', color = 130
  oplot,[32.5, 40.5],[2, 2],psym = 8, symsize = 3, color = 240
  ;if keyword_set(spf3) then
  ;xyouts, 43,2, 'Spread F (Yes)', color = 240
  ;endif
  if nrcvrs gt 1 then begin
    oplot,[80.5, 88],[2, 2],psym = 8, symsize = 3, color = 200
    xyouts, 74,2, 'O mode', color = 200
    xyouts, 82,2, 'X mode', color = 200
    if (selOmode eq !true) then xyouts, 821, 40, 'X', /device
    if (selXmode eq !true) then xyouts, 892, 40, 'X', /device
  endif
  if nrcvrs eq 4 then begin
    oplot, [97],[2], psym =8, symsize = 3,color = 180
    xyouts, 90,4, 'enhance' , color = 180
    xyouts, 90,2, 'overhead', color = 180
    if enhanceoverhead eq !true then xyouts, 978,40, 'X', /device
  endif
  oplot, [55],[2], psym = 8, symsize = 3, color = boxcolor
if not keyword_set(heightfrequencyscale) then xyouts, 55,2, 'clear data', color = 230 else $
      xyouts, 52.5,5, 'no ESF', color = 230

      oplot, [60],[2], psym = 8, symsize = 3, color = boxcolor
  if not keyword_set(heightfrequencyscale) then xyouts, 60,2, 'clear data', color = 230 else $
      xyouts, 57.5,5, 'ESF(R)', color = 230

      ;oplot, [45],[2],psym = 8, symsize = 5, color = 200
  oplot, [65],[2], psym = 8, symsize = 3, color = boxcolor
  if not keyword_set(heightfrequencyscale) then xyouts, 65,2, 'clear data', color = 230 else $
      xyouts, 62.5,5, 'ESF(F)', color = 230

 oplot, [70],[2], psym = 8, symsize = 3, color = boxcolor
  if not keyword_set(heightfrequencyscale) then xyouts, 70,2, 'clear data', color = 230 else $
      xyouts, 67.5,5, 'ESF(M)', color = 230


;set up 2d array to hold image of the ionogram
nfreqs = h.nfreqs
nhghts = h.nheights
igrmarr = fltarr(nfreqs, nhghts)
if not keyword_set(heightfrequencyscale) then title='CADI ' + h.site + " Ionograms    " + strcantime(/date,cantimeset(0))$
  else title='CADI ' + h.site + " Ionograms    " + strcantime(/date,cantimeset(0))+ "                      fof2 = "+$
        string(format = '(f4.1)',cursorfreq) + "  h'f1 = "+ string(format = '(f6.1)',cursorheight) + $
           " h'f2 = "+string(format = '(f6.1)',cursorheight1)+"  hpf2 = " + string(format = '(f6.1)',hpheight) + $
           "  fbEs = " + string(format = '(f4.1)',cursorfreq1)

  plot_oi, /nodata,/noerase, [0], [0], position = [0.05, 0.15, 0.95, 0.94], xtitle = 'MHz', ytitle = 'km',$
     xrange=mhzrange, xstyle=1, yrange=hrange, ystyle=1, xtickname = def_ticknames, ytickname = def_ticknames,$
        title = title, ymargin = [100,100],/noclip
if keyword_set(heightfrequencyscale) then begin
   oplot,[cursorfreq,cursorfreq], hrange,linestyle = 5
   oplot,[cursorfreq1,cursorfreq1], hrange,linestyle = 2, color = 240  ; for fbEs
   oplot, [cursorfreq+site.gyrofreq/2.0,cursorfreq+site.gyrofreq/2.0], hrange, linestyle = 2
   oplot, [0.834*cursorfreq,0.834*cursorfreq],hrange, linestyle = 1
   oplot, [1,20],[cursorheight,cursorheight], linestyle = 3
   oplot, [1,20],[cursorheight1,cursorheight1], linestyle = 3,color = 240  ; for range spread
   oplot, [1,20], [hpheight,hpheight], linestyle = 1
endif
  maxpowdb = 0
  if cadi.ndata gt 0 then begin          ; Any echoes in this hour?
    w = where(d.x.t eq htimex, wcount)  ; Any echoes for this time?
    if wcount gt 0 then begin
      d_ionogram = d(w)
      prevxpos = 0
      prevypos = 0
   harr = [0.0]
   warr = [0.0]
   t_harr = [0.0]
   t_farr = [0.0]
   elno = 1
      ;sumweightedheights= 0
      ;sumweights = 0
      for p = 0,n_elements(d_ionogram)-1 do begin
        rept1 = fix(d_ionogram(p).IQ(0,0))
        if (rept1 gt 127) then rept1 = rept1 - 256
        rept1 = float(site.polarity(0))*rept1
        impt1 = fix(d_ionogram(p).IQ(1,0))
        if (impt1 gt 127) then impt1 = impt1 - 256
        impt1 = float(site.polarity(0))*impt1
        rcvr1cmplx = complex(rept1,impt1)
        rcvr1phase = atan(impt1,rept1)
        if (nrcvrs gt 1) then begin
          rept2 = fix(d_ionogram(p).IQ(0,1))
          if (rept2 gt 127) then rept2 = rept2 - 256
          rept2 = float(site.polarity(1))*rept2
          impt2 = fix(d_ionogram(p).IQ(1,1))
          if (impt2 gt 127) then impt2 = impt2 - 256
          impt2 = float(site.polarity(1))*impt2
          rcvr2cmplx = complex(rept2,impt2)
          rcvr2phase = atan(impt2,rept2)
        endif
        if (nrcvrs gt 2) then begin
          rept3 = fix(d_ionogram(p).IQ(0,2))
          if (rept3 gt 127) then rept3 = rept3 - 256
          rept3 = float(site.polarity(2))*rept3
          impt3 = fix(d_ionogram(p).IQ(1,2))
          if (impt3 gt 127) then impt3 = impt3 - 256
          impt3 = float(site.polarity(2))*impt3
          rcvr3cmplx = complex(rept3,impt3)
          rcvr3phase = atan(impt3,rept3)
        endif
        if (nrcvrs gt 3) then begin
          rept4 = fix(d_ionogram(p).IQ(0,3))
          if (rept4 gt 127) then rept4 = rept4 - 256
          rept4 = float(site.polarity(3))*rept4
          impt4 = fix(d_ionogram(p).IQ(1,3))
          if (impt4 gt 127) then impt4 = impt4 - 256
          impt4 = float(site.polarity(3))*impt4
          rcvr4cmplx = complex(rept4,impt4)
          rcvr4phase = atan(impt4,rept4)
        endif
        if (nrcvrs eq 4) then begin
          rcvr12 = (rcvr1cmplx + rcvr2cmplx)/2
          rcvr34 = (rcvr3cmplx + rcvr4cmplx)/2
        endif else if nrcvrs eq 2 or nrcvrs eq 3 then begin
          rcvr12 = rcvr1cmplx
          rcvr34 = rcvr2cmplx
        endif ;nrcvrs 2 or 4
        sig = rcvr1cmplx
        if nrcvrs gt 1 then begin
        Omodeamp = abs((rcvr12 + rcvr34*complex(0,1))/2)
        Xmodeamp = abs((rcvr12 - rcvr34*complex(0,1))/2)
if (selOmode eq !true) and (selXmode eq !false) then if Omodeamp gt Xmodeamp then sig = 2.0*Omodeamp else sig = 1.5
if (selOmode eq !false) and (selXmode eq !true) then if Xmodeamp gt Omodeamp then sig = 2.0*Xmodeamp else sig = 1.5
        ;if (selOmode eq !true) and (selXmode eq !false) then sig = (rcvr12 + rcvr34*complex(0,1))/2
        ;if (selOmode eq !false) and (selXmode eq !true) then sig = (rcvr12 - rcvr34*complex(0,1))/2
        if (selOmode eq !false) and (selXmode eq !false) then sig = complex(0,0)
        ;if (p eq n_elements(d_ionogram)-1) then stop
        endif
        mag = abs(sig)
        if (mag lt 1) then mag = 1.5
        if nrcvrs eq 4 and enhanceoverhead eq !true then begin
          ;determine angular location of echo, eed phase difference between receiver pairs
          phase12 = rcvr1phase - rcvr2phase + site.phasecorrectionEW    ;phase difference rcvr1 and 2
          if phase12 gt !pi then phase12 = phase12 - 2.0*!pi       ;check if in range -pi to pi
          if phase12 lt -!pi then phase12 = phase12 + 2.0*!pi
          phase34 = rcvr3phase - rcvr4phase + site.phasecorrectionNS    ;phase difference rcvr 3 and 4
          if phase34 gt !pi then phase34 = phase34 - 2.0*!pi       ;check if in range -pi to pi
          if phase34 lt -!pi then phase34 = phase34 + 2.0*!pi
          ;ok have phase differences so relate these to distance from overhead using
                                                       ;sin(theta) = phase*wavelength/2*pi*d
          sinthetaEW = phase12*interferometerconstEW/h.freqs(d_ionogram(p).x.f)
          sinthetaNS = phase12*interferometerconstNS/h.freqs(d_ionogram(p).x.f)
          ;next checkif distance from overhead is small (sintheta = 0.25 is ~15 deg)
          if ((abs(sinthetaEW) gt 0.25) or (abs(sinthetaNS) gt 0.25)) then mag = 1.5 else mag = 2.0*mag ;smallest value for a color
        endif
        powdb = byte(7*20*alog10(mag))             ;convert power to db and then scale to 255 for color
        ;if (powdb gt 254) then powdb = 254
        xpos = h.freqs(d_ionogram(p).x.f)/1.0e6               ;get x
        ypos = 3*d_ionogram(p).x.h                      ;get y
;if ((xpos - 4) le 0.1) and (enhanceoverhead eq !true) then stop
        if (xpos eq prevxpos) and (ypos eq prevypos) then begin     ;check if this is the max value for this x,y
            if powdb gt maxpowdb then maxpowdb = powdb
        endif else begin  ;different height, freq
          oldmax = maxpowdb
          maxpowdb = powdb
        endelse

            ;next get weighted mean height,first make uparrays of heights and weights
        if xpos eq prevxpos then begin
          if ypos ne prevypos then begin
            if oldmax gt 100 then harr = [harr, float(prevypos)]
            if oldmax gt 100 then warr = [warr, 10.0^(float(oldmax)/140)]
          endif  ;ypos
        endif else begin    ;xpos ne prevxpos: gone to next frequency


   harr = [0.0]
   warr = [0.0]
   elno = 1

endelse ;next frequency, or height
        prevxpos = xpos
        prevypos = ypos

        oplot, [xpos], [ypos],psym = 8, $                       ;do the plotting
               color = 255 - maxpowdb, symsize = float(maxpowdb)/255 + 0.2

      ;put value into igrmarr
      igrmarr( d_ionogram(p).x.f,fix((float(ypos)-hrange(0))/3)) = maxpowdb
    endfor


    endif ; wcount gt 0
  endif ; cadi.ndata gt 0
;have finished plotting the ionogram

            ; Write time on ionogram
  hhmm_str = string(hour,':', minute, format='(i2.2,a,i2.2)')
  hhmm_str1 = string(hour,' ',minute, format='(i2.2,a,i2.2)')
  xyouts, 1.3, hrange(0) + 70, '!c'+hhmm_str, charsize = 2
;xyouts, 1.5, hrange(0) + 70, 'Spread F (yes(1)/no(0))', charsize = 2
if not keyword_set(heightfrequencyscale) then begin

if autoscaledone eq !false then begin      ;start the autoscaling
;next filter the ionogram to remove noise
transfigrm =fft(igrmarr,-1)               ;fft the ionogram
;filter by setting high freq elements to 0
for f = fix(0.2*float(h.nfreqs)),fix(0.8*float(h.nfreqs)) do transfigrm(f,*) = complex(0,0)
for hgt = fix(0.2*float(h.nheights)),fix(0.8*float(h.nheights)) do transfigrm(*,hgt) = complex(0,0)
igrmarr = fft(transfigrm,+1)              ;fft in reverse dirn


; for freq = 0,h.nfreqs -1 do begin         ;plot the filtered ionogram
;   for height = 0,h.nheights-1 do $
;     if abs(igrmarr(freq,height)) gt 30 then $
;     oplot,[h.freqs(freq)/1.0e6],[h.heights(height)], color = 255 - abs(igrmarr(freq,height)),$
;                                   psym = 8, symsize = abs(igrmarr(freq,height))/255 + 0.2
;
; endfor


 maxigrm = max(igrmarr,index)     ;find maximum of the ionogram
 IndexX = index MOD h.nfreqs
 IndexY = index/h.nfreqs

 thresh = 30
;PRINT, IndexX, IndexY,h.freqs(IndexX),h.heights(IndexY)
;oplot, [h.freqs(IndexX)/1.0e6],[h.heights(IndexY)],psym = 2,symsize = 2, color = 0 ;mark max

;start scaling at freq of the max as determined above, but use weighted means
ymindex = round(w_avg(igrmarr(IndexX,*),h.nheights, thresh))    ;find weighted mean for the max
if ymindex eq 0 then goto, noautoscaledata          ;no datapoints on ionogram above threshold
xwmeanarr = [IndexX]              ;start an array
ywmeanarr = [ymindex]             ;of the weighted means

;get datapoints at higher and lower points using weighted means
;need to decide whether to step in freq or in height (depends on slope)
;so do freqs above the max and then find the slope
xindex = IndexX
ylow = ymindex
deltay = 0

repeat begin
  if n_elements(igrmarr) lt 1 then goto, noautoscaledata
 ;go to next higher frequency and get weighted mean
  xindex = xindex + 1
  if xindex eq nfreqs then begin
    print, 'can''t do autoscale'
    goto, noautoscaledata
  endif
  ;if xindex gt n_elements(igrmarr) then goto, noautoscaledata
  ;print, 'index x = ',xindex
  ylow = ylow + deltay
  if ylow lt 10 then begin
    print, 'ylow too near border of array, ylow =',ylow
    goto, noautoscaledata
  endif else if n_elements(igrmarr(xindex,*)) lt 21 then begin
     print, 'not enough datapoints for autoscale'
    goto, noautoscaledata
  endif else if ylow ge nhghts-10 then begin
    print, 'ylow too near upper border of array, ylow =',ylow
    goto, noautoscaledata

   endif else $
    ymindex = round(w_avg(igrmarr(xindex,(ylow - 10):(ylow+10)),21,thresh))


    ;ymindex = round(w_avg(igrmarr(xindex,*),h.nheights, thresh))

  if ymindex gt 0 then begin
    deltay = ymindex - 10
    ylow = ylow - 10 + ymindex
    xwmeanarr = [xwmeanarr,xindex]
    ywmeanarr = [ywmeanarr,ylow]
  endif
endrep until deltay ge 2         ;steep slope so switch to incrementing height

;take height as the known (will increment height) and calculate weighted mean freq
;xmax = xindex
htptr = ylow
htptr_s = htptr
xptr = xindex+1
repeat begin
  xmindex = fix(w_avg(igrmarr((xptr-4):(xptr+4),htptr),9,30))   ;use 9 freq centred on last weighted freq mean
  if xmindex gt 0 then begin          ;a 0 value is returned if no suitable values for weighted mean calc
    xfr = xmindex+xptr-4              ;adjust freq value by lowest freq used for weighted mean calc
    if xfr lt xptr then begin
      ;calculate weighted mean using both heights and choose best value
      Carray = igrmarr((xptr-4):(xptr+4),htptr) + igrmarr((xptr-4):(xptr+4),htptr_s)       ;get combined values
      xmindex = fix(w_avg(Carray,9,30)) ;use 9 freq centred on last weighted freq mean   ;and recalculate weighted mean
      xfr = xmindex + xptr -4
      arrmax = max(xwmeanarr(0:n_elements(xwmeanarr)-2))
      if xfr lt arrmax then xfr = arrmax
        xwmeanarr(n_elements(xwmeanarr)-1) = xfr
     endif
    ;xfr = xptr          ;avoid decreses in frequency: causes error in polan
    xwmeanarr = [xwmeanarr,xfr]          ;add values to top of arrays
    ywmeanarr = [ywmeanarr,htptr]
    xptr = xfr
    htptr_s = htptr
  endif
  ;print, 'htptr, xptr, xmindex =',htptr,xptr,xmindex
  htptr = htptr + 1
endrep until htptr eq h.nheights or htptr gt htptr_s + 17 ;stop if at top or 51 km above last point saved

;ok have gone to the top of the ionogram, next start at max and go down
htptr = ylow
xptr = IndexX-1
xptr_s = xptr
repeat begin
  if htptr lt 15 then goto, cantgolower
  if htptr ge nhghts - 15 then goto, cantgohigher
  ymindex = fix(w_avg(igrmarr(xptr,(htptr-15):(htptr+15)),21,30))     ;use 31 heights centred on weighted mean
  if ymindex gt 0 then begin
    xwmeanarr = [xptr,xwmeanarr]
    ywmeanarr = [ymindex +htptr -15, ywmeanarr]
    htptr = ymindex+htptr-10
    xptr_s = xptr
  endif
  ;print, 'htptr, xptr, xmindex =',htptr,xptr,xmindex
  xptr = xptr - 1
endrep until xptr eq 0 or xptr lt xptr_s - 6    ;stop if at beginning or too large a gap from last saved value
cantgolower:
cantgohigher:

;smooth the x and y arrays if enough data points
if n_elements(xwmeanarr) gt 11 then begin
xsmoothed = fix(smooth(float(xwmeanarr),11))
ysmoothed = fix(smooth(float(ywmeanarr),11))
endif else begin
xsmoothed = xwmeanarr
ysmoothed = ywmeanarr
endelse

if n_elements(xwmeanarr) lt 5 then goto, noautoscaledata

;have autoscaled values, next interpolate to get values at 0.1 MHz or 0.2 MHz
xfreqsmooth = h.freqs(xsmoothed)/1.0e6
minf = min(xfreqsmooth)
maxf = max(xfreqsmooth)
freqspan = maxf - minf
if freqspan gt 8 then finc = 0.2 else finc = 0.1
fsarr10 = finc*findgen(fix(freqspan/finc)+1)+minf    ;make array with freq each 0.1 MHz

interpoly = interpol( ysmoothed,xfreqsmooth,fsarr10)       ;interpolate to get data at each 0.05 MHz
;oplot, fsarr10, h.heights(interpoly), linestyle = 3, thick =2, color = 240 ;plot the autoscaled values
hsa = h.heights(interpoly) + hrange(0)
;need to add end values for polan
fsa = [fsarr10,0.0]           ;polan end of array is 0.0
;hsa = [h.heights(interpoly),0.0]   ;polan end of array  is 0.0
hsa = [hsa,0.0]
hsa(n_elements(hsa)-2) = 0.0       ;polan critical frequency marker is to set height to 0.0
 realhtanalysis, dir, site, hhmm_str, hour, minute, fsa,hsa,autofreqs,autoheights
 ;oplot, autofreqs,autoheights, color = 240, thick = 3, linestyle = 1     ;plot the real ht analysis values
;if hhmm_str eq '00:25' then stop
autoscaledone = !true
endif ;autoscaledone !false

;plot/replot the autoscaled values
;oplot, xfreqsmooth,h.heights(ysmoothed), color = 130, thick = 4         ;show smoothed with heavy green line
;print,'fsarr10 = ',fsarr10,' h.heights(interpoly) = ', h.heights(interpoly)
if autoscaledone eq !true then begin
  if n_elements(fsarr10) gt 1 then oplot, fsarr10, h.heights(interpoly)+hrange(0), linestyle = 3, thick =2, color = 240 ;plot the autoscaled values
  if n_elements(autofreqs) gt 1 then oplot, autofreqs,autoheights, color = 240, thick = 3, linestyle = 1       ;plot the real ht analysis values
endif
noautoscaledata:
endif ;not keyword_set(heightfrequencyscale)
  ;if (htimex ge h.ntimes -1) then goto, endtimex   ;check if at end of data, if true then end plotting

if datascanned eq !true then begin
     oplot, pxarr, pyarr, linestyle = 2, thick = 2, color = 0                    ;superimpose the scanned values
  if n_elements(freqs) gt 0 then oplot, freqs, heights, linestyle = 1, color = 0, thick = 3 ;and the real heights
  datascanned = !false
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;MOUSE CLICKS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;next get mouse click for what to do

  CURSOR, x, y, /DEVICE, /wait
  print, 'cursor clicked at ',x,',',y

if keyword_set(heightfrequencyscale) and y gt 90 then begin
   cursor, xs,ys,/DATA, /wait
   if ((xs gt 2.5) and (xs lt 18) and (ys lt 500)) then cursorfreq = xs else begin
     if ((xs gt 2.5) and (xs lt 18) and (ys gt 500))then cursorfreq1 = xs
     if (xs lt 1.25) then cursorheight1 = ys
     if (xs gt 1.50) and (xs lt 2.5) then cursorheight = ys  else if (xs gt 19) then hpheight = ys
   endelse
   goto, startplotpanels
endif


  if (abs(x - 98) le 10) and (abs(y - 44) le 10) then begin     ;Done clicked
     close, /all
     datascanned = !false
     ;datafileopen = !false
     ;
     return
  endif

  if (abs(x - 204) le 10) and (abs(y - 44) le 10) then begin    ;step forward 10
    if datafileopen eq !true then close, 1
    if (htimex ge h.ntimes -1) then goto, endtimex  ;check if at end of data, if true then end plotting
    datascanned = !false
    autoscaledone = !false
    htimex = htimex + 10
    if (htimex gt h.ntimes-1) then htimex = h.ntimes -1
    boxcolor = 230
  endif

  if ((abs(x - 290) le 10) and (abs(y - 44) le 10)) then begin   ;step forward 1
    forward1:
    ;stop
     if datafileopen eq !true then close, 1
     if (htimex ge h.ntimes -1) then goto, endtimex     ;check if at end of data, if true then end plotting
     autoscaledone = !false
     datascanned = !false
     htimex = htimex + 1
     boxcolor = 230
     wait, 0.2
     goto, startplotpanels
  endif

  if ((abs(x - 362) le 10) and (abs(y - 44) le 10)) then begin      ;step back 1
     if datafileopen eq !true then close, 1
     ;datafileopen = !false
     autoscaledone = !false
     datascanned = !false
     htimex = htimex - 1
     boxcolor = 230
  endif

  if ((abs(x - 438) le 10) and (abs(y - 44) le 10)) then begin          ;step back 10
    if datafileopen eq !true then close, 1
    ;datafileopen = !false
    datascanned = !false
    autoscaledone = !false
  htimex = htimex -10
    if (htimex lt 0) then htimex = 0
    boxcolor = 230
  endif


  if nrcvrs gt 1 then begin
    if (abs(x - 823) le 10) and (abs(y - 44) le 10) then begin   ;user clicked on O mode box
     if (selOmode eq !true) then begin
       selOmode = !false
       xyouts, 821, 40, 'X', /device, color = !p.background
     endif else begin
       selOmode = !true
       xyouts, 821, 40, 'X', /device
     endelse
    autoscaledone = !false   ;need to redo autoscale
    endif
  endif

  if nrcvrs gt 1 then begin
    if (abs(x - 894) le 10) and (abs(y - 44) le 10) then begin   ;user clicked on X mode box
      if (selXmode eq !true) then begin
       selXmode = !false
       xyouts, 892, 40, 'X', /device, color = !p.background
     endif else begin
       selXmode = !true
       xyouts, 892, 40, 'X', /device
     endelse
     autoscaledone = !false     ;need to redo autoscale
   endif ;(x)
 endif ;nrcvrs gt 1

  if nrcvrs eq 4 then begin
    if ((abs(x - 981) le 10) and (abs(y - 44) le 10)) then begin    ;user clicked enhance overhead mode box
      if (enhanceoverhead eq !true) then begin
       enhanceoverhead = !false
       xyouts, 978, 40, 'X', /device, color = !p.background
     endif else begin
       enhanceoverhead = !true
       xyouts, 978, 40, 'X', /device
     endelse
     autoscaledone = !false     ;need to redo autoscale
  endif ;(x)
  endif ;nrcvrs eq 4


if abs(x - 573) le 10 and abs(y - 44) le 10 then spreadf=0
if abs(x - 623) le 10 and abs(y - 44) le 10 then spreadf=1
if abs(x - 673) le 10 and abs(y - 44) le 10 then spreadf=2
if abs(x - 723) le 10 and abs(y - 44) le 10 then spreadf=3
  if (abs(x - 723) le 10 and abs(y - 44) le 10) OR (abs(x - 673) le 10 and abs(y - 44) le 10) OR (abs(x - 623) le 10 and abs(y - 44) le 10) OR (abs(x - 573) le 10 and abs(y - 44) le 10) then begin
      if not keyword_set(heightfrequencyscale) then datascanned = !false else begin     ;clicked to clear scaling
      ;find out if file for data exists if scaling hight fof2 then click this box to record data
        print, 'directory = ', dir
        ;stop
        if hpheight gt 100 then begin
        datafilename = 'c:\cdata\' + strmid(dir,3,8) + '_F.tfh'
        if cursorfreq le 2.0 then cursorfreq = 'NaN'
        if cursorheight le 150 then cursorheight = 'NaN'
        if cursorheight1 le 150 then cursorheight1 = 'NaN'
        if hpheight le 150 then hpheight = 'NaN'
        print, datafilename
        print, spreadf
        openw, 1, datafilename, /append
      ;put data in file
      ;if keyword_set(spreadf) then
        printf, 1, format = '(a,f6.2,f6.1, f6.1, f6.1, f6.1)',hhmm_str1,cursorfreq,cursorheight,cursorheight1,hpheight, spreadf
      close, 1
      endif else begin
      datafilename = 'c:\cdata\' + strmid(dir,3,8) + '_E.tfh'
        print, datafilename
        openw, 1, datafilename, /append
      ;put data in file
        printf, 1, format = '(a,f6.2,f6.1, f6.2)',hhmm_str1,cursorfreq,cursorheight,cursorfreq1
      close, 1   ;close the file
      end

      boxcolor = 0  ;change click box to black
      goto, forward1
     endelse
  endif

           ;     if abs(x - 573) le 10 and abs(y - 44) le 10 then begin
;      spreadf=0
;        if not keyword_set(heightfrequencyscale) then datascanned = !false else begin     ;clicked to clear scaling
;      ;find out if file for data exists if scaling hight fof2 then click this box to record data
;        print, 'directory = ', dir
;        ;stop
;        if hpheight gt 100 then begin
;        datafilename = 'c:\cdata\' + strmid(dir,3,8) + '_F.tfh'
;        if cursorfreq le 2.0 then cursorfreq = 'NaN'
;        if cursorheight le 150 then cursorheight = 'NaN'
;        if cursorheight1 le 150 then cursorheight1 = 'NaN'
;        if hpheight le 150 then hpheight = 'NaN'
;        print, datafilename
;        print, spreadf
;        openw, 1, datafilename, /append
;      ;put data in file
;      ;if keyword_set(spreadf) then
;        printf, 1, format = '(a,f6.2,f6.1, f6.1, f6.1, f6.1)',hhmm_str1,cursorfreq,cursorheight,cursorheight1,hpheight, spreadf
;      close, 1
;      endif



  if (x gt 49) and (y gt 90) and datascanned eq !false then begin    ;cursor in plot region, start reading data
      dno = 0                    ;array pointer
      dxarr = fltarr(200)             ;frequency array
      dyarr = fltarr(200)             ;height array
      oldx = 0.0
      repeat begin
        cursor, x, y, /nowait, /data       ;read cursor position in data coordinates
        if finc gt 0.2 then fs = 0.2 else fs = 0.1  ;uselargersteps if wide frequency range
        if x ge oldx + fs    then begin     ;if freq incremented by 0.1 save the data
          dxarr(dno) = x              ;in arrays dxarr
          dyarr(dno) = y              ;and dyarr
          dno = dno + 1               ;increment the array pointer
          oldx = x                   ;copy frequency to oldx (could just get value from array)
        endif
      endrep until !mouse.button ne 1    ;mouse button released so datascan finished
                                          ;check if value when mouse button released is useful
      if x gt dxarr(dno-1) then begin     ;if higher frequency it is presumably ok
        dxarr(dno) = x             ;so save the value
        dyarr(dno) = y
        dno = dno + 1
      endif
      pxarr = dxarr(0:dno-1)       ;transfer the non zero freqency value to a plot array
      pyarr = dyarr(0:dno-1)       ;along with the height values; these will be overplotted as a dashed curve
      dxarr(dno) = dxarr(dno-1) +0.05     ;set the highest frequency read as a critical frequencybut need to increment
      ;dyarr(dno -1) = 0.0                                 ;the frequency a little otherwise get an error in polan
      dxarr = dxarr(0:dno +1)     ;trimarray but keep one 0.0 at end
      dyarr = dyarr(0:dno +1)     ;trim array but keep one 0.0 at end
      ;print, 'dxarr = ', dxarr
      ;print, 'dyarr = ', dyarr
      datascanned = !true
      ;stop
   realhtanalysis, dir, site, hhmm_str, hour, minute, dxarr,dyarr,freqs,heights     ;run polan to get real ht analysis
   endif                           ;cursor value: reading plot data

;print,'!false = ', !false, '  autoscaledone = ', autoscaledone
wait, 0.2
goto, startplotpanels               ;do next


endtimex:
 ;endfor ; timex
   end_hourx_loop:
  ;; endfor ; hourx


!p.noerase = 0
!p = psav
!x = xsav
!y = ysav

error_return:
return

end
      program Trajectory
      integer, parameter::itraj=100000
      real*8 V(100000),t(100000),gamma(100000),rho(100000),alt(100000)
     &,S(100000),
     &Q(100000),g,goo,rhoo,B,alti,dt,decel,H,Re,nmax,n(100000),
     &dvdt,LoD,bank,Lf,Po
        real*8 Vt,alt_t,gammat,rhot,Vtold,gammatold,Vnmax,hnmax,alt_told
        real*8 qconv(itraj),qtotal(itraj),qrad(itraj),rnose,dS,acc
        real*8 tatm, p_atms
        integer nc,Tp,tprnt,prntint,iplanet

      real*8 Hrec,ruCh,lam,surfP,patm,Cp,qradbc

      real*8 y1k1,y1k2,y1k3,y1k4,y1k5,y1k6
      real*8 y2k1,y2k2,y2k3,y2k4,y2k5,y2k6
      real*8 y3k1,y3k2,y3k3,y3k4,y3k5,y3k6
      real*8 y4k1,y4k2,y4k3,y4k4,y4k5,y4k6

      real*8 k11,k12,k13,k14,k15,k16
      real*8 k21,k22,k23,k24,k25,k26
      real*8 k31,k32,k33,k34,k35,k36
      real*8 k41,k42,k43,k44,k45,k46

      write(*,*) 'This Program numerically integrates the ballistic'
      write(*,*) 'entry equations of motion'
      write(*,*)
      write(*,*) 'Input Planet: 1) Earth 2) Mars '
        read(*,*) iplanet
        write(*,*) 'input initial flight path angle: '
      read(*,*) gamma(1)
      write(*,*) 'Input entry velocity, km/s: '
      read(*,*) V(1)
      V(1)=V(1)*1000.0 !convert to m/s
      write(*,*) 'Input Ballistic coefficient: '
      read(*,*) B
      write(*,*) 'Input effective nose radius (m)'
      read(*,*) rnose
      write(*,*) 'Input L/D: '
      read(*,*) LoD
      write(*,*) 'Input Bank angle: '
      read(*,*) bank
      write(*,*) 'Input entry interface altitude, km'
      read(*,*) alti
      alti=alti*1000.0 !convert to meters
c     write(*,*) 'Input time step: '
c     read(*,*) dt
        dt=.05
c     write(*,*) 'Input time interval for output: '
c     read(*,*)Tp

      write(*,*)
      write(*,*) 'This program assumes: '
      write(*,*) ' - The scale height is 7.25km for Earth, 11.1km for
     &Mars'
      write(*,*) ' - Exponential atmosphere, rho=rhoi*EXP(-h/H) '
      write(*,*) ' - g varies with altitude, g=goo/(1+h/ro)**2 '
      write(*,*)

      t(1)=0.0
        Cp=1.27
        prntint=1.0
        tprnt=1.0
      alt(1)=alti
      !Earth
      if (iplanet .eq. 1) then
        goo=9.8062
      rhoi=1.226
      !H=6850.0
        H=7250.0
      Re=6378000.0
        Po=101325.0
      endif

        if (iplanet .eq. 2) then
        goo=3.71
      rhoi=0.057
      H=11100.0
      Re=3380000.0
        Po=600.0
        endif

        nmax=100.0
      n(1)=0.0
      qconv(1)=0.0
        qrad(1)=0.0
      qtotal(1)=0.0
          S(1)=0.0
      dS=0.0
        Q(1)=0.0
        Hrec=(V(1)**2./2.)/2325.854324
        ruCh=1.0D-6
        lam=.4

        if (iplanet .eq. 1) then
            tatm = 141.94 + 0.00229 * alt(1)
            patm = 2488 * ((tatm/216.6))**(-11.388)
                    if (tatm .LE. 0) then
                        tatm = 2.73
                    endif
                rhot = ABS((patm/1000) / (.2869 * (tatm)))
        endif

        if (iplanet .eq. 2) then
C        INSERT MARS , PRESSURE, TEMPERATURE - INITIAL HIGH ALTITUDE
            patm = 669*exp(-0.00009*alt(1))
            tatm = 249.75-0.00222*alt(1)
            rhot = ABS((patm/1000)/(.1921*(tatm)))
        endif

        if (tatm .LE. 0) then
                tatm = 2.73
        endif

        if (patm .LE. 0) then
                patm = 0
        endif
c        rhot=rhoi*exp(-alt(1)/H)
c        patm=(goo/(1 + (alt(1)/Re))**2)*rhoi*dexp(-alt(1)/H)*alt(1)
c        patm=po*rhot/rhoi
        surfP=(patm+Q(1)*Cp)/101325
        qradbc=qrad(1)*.881
      p_atms = patm/101325
c     Output trajectory info every 2 seconds
      nc=Tp/dt ! number of calculations per output interval
      open(unit=1,file='output.dat',status='unknown',form='formatted')
      open(unit=11,file='cma_out.dat',status='unknown',form='formatted')
      open(unit=21,file='heat_in.txt',status='unknown',form='formatted')

      i=1
c     write(1,110)'variables="time (s)","Alt (m)","Vel (m/s)"
c     &,"gamma (deg)","
c     &decel (g)","Q conv (W/cm2)","Qrad (W/cm2)","Qtotal (J/cm2)"
c     &,"S (m)","P (N/m2)"'
      write(1,111)'variables="time (s)","Alt (m)","Vel (m/s)"
     &,"Patm (N/m2)","Rho atm (kg/m3)","gamma (deg)","
     &decel (g)","Q conv (W/cm2)","Qrad (W/cm2)","Qtotal (J/cm2)"
     &,"S (m)","Pdyn (N/m2)"'
      write(1,101)'zone T="BC=',B,' V=',V(1),' Gamma=',gamma(1),' LoD='
     &,LoD,' Bnk=',bank,' rnose=',rnose,'"'
      write(1,105)t(i),alt(i),V(i),patm,rhot,gamma(i),n(i),qconv(i),
     &qrad(i),qtotal(i),S(i),Q(i)
      write(21,44) p_atms,tatm,V(i)
      write(11,200) t(i),Hrec,qradbc,ruCh,surfP,lam


      do while (alt(i) .gt. 0.0)

*  Setup starting conditions for the next interval
          alt_t = alt(i)
             Vt = V(i)
        gammat  = gamma(i)


c   do j=1,nc

        g=goo/(1 + (alt_t/Re))**2

        if (iplanet .eq. 1) then

            if (alt_t .LT. 11000) then
                tatm = 288.19-0.00649*alt_t
                patm = 101290 * (tatm/288.08)**5.256
                    if (tatm .LE. 0) then
                        tatm = 2.73
                    endif
                rhot = ABS((patm/1000)/(.2869*(tatm)))
            endif

            if (alt_t .GE. 11000 .AND. alt_t .LE. 25000) then
                tatm = 216.69
                patm = 0.02265 * exp(1.73-((0.000157)* alt_t))
                    if (tatm .LE. 0) then
                        tatm = 2.73
                    endif
                rhot = ABS((patm/1000) / (.2869 * (tatm)))
            endif

            if (alt_t .GT. 25000) then
                tatm = 141.94 + 0.00229 * alt_t
                patm = 2488 * ((tatm/216.6))**(-11.388)

                    if (tatm .LE. 0) then
                        tatm = 2.73
                    endif
                rhot = ABS((patm/1000) / (.2869 * (tatm)))
            endif
c            rhot=rhoi*exp(-alt_t/H)
c           patm=po*rhot/rhoi
        endif

        if (iplanet .eq. 2) then
C        INSERT MARS , PRESSURE, TEMPERATURE - CHECK ALTITUDE

            if (alt_t .LT. 7000) then
                patm = 699 * exp(-0.00009 * alt_t)
                tatm = 242.15-0.00222*alt_t
                    if (tatm .LE. 0) then
                        tatm = 2.73
                    endif
                rhot = ABS((patm/1000)/(.1921*(tatm)))
            endif

            if (alt_t .GE. 7000) then
                patm = 669 * exp(-0.00009 * alt_t)
                tatm = 249.75 - 0.00222 * alt_t
                    if (tatm .LE. 0) then
                        tatm = 2.73
                    endif
                rhot = ABS((patm/1000) / (.1921 * (tatm)))
            endif

        endif

        if (tatm .LE. 0) then
                tatm = 2.73
        endif

        if (patm .LE. 0) then
                patm = 0
        endif

c        rhot=rhoi*exp(-alt_t/H)
        y1k1=Vt
        y3k1=alt_t
        y2k1=gammat

c   if (LoD .gt. 0.0) then

c       Vt=Vt+(-(7702.0**2-Vt**2)/((Re+alt_t)*Lod*cosd(bank))
cc    &     +g*sind(gammat))*dt

c   else
*
        k11=(-rhot*y1k1**2)/(2.0*B) + g*sin(y2k1)
        k21=(-y1k1*cos(y2k1) )/(Re+y3k1)
     & -((rhot*y1k1)/(2.0*B))*(LoD)*cos(bank)
     & + (g/y1k1)*cos(y2k1)
        k31=-y1k1*sin(y2k1)


        y1k2=y1k1+0.25*k11*dt
        y2k2=y2k1+0.25*k21*dt
        y3k2=y3k1+0.25*k31*dt

        k12=(-rhot*y1k2**2)/(2.0*B) + g*sin(y2k2)
        k22=(-y1k2*cos(y2k2) )/(Re+y3k2)
     & -((rhot*y1k2)/(2.0*B))*(LoD)*cos(bank)
     & + (g/y1k2)*cos(y2k2)
        k32=-y1k2*sin(y2k2)

        y1k3=y1k1+0.125*k11*dt+0.125*k12*dt
        y2k3=y2k1+0.125*k21*dt+0.125*k22*dt
        y3k3=y3k1+0.125*k31*dt+0.125*k32*dt

        k13=(-rhot*y1k3**2)/(2.0*B) + g*sin(y2k3)
        k23=(-y1k3*cos(y2k3) )/(Re+y3k3)
     & -((rhot*y1k3)/(2.0*B))*(LoD)*cos(bank)
     & + (g/y1k3)*cos(y2k3)
        k33=-y1k3*sin(y2k3)

        y1k4=y1k1-0.5*k12*dt+k13*dt
        y2k4=y2k1-0.5*k22*dt+k23*dt
        y3k4=y3k1-0.5*k32*dt+k33*dt

        k14=(-rhot*y1k4**2)/(2.0*B) + g*sin(y2k4)
        k24=(-y1k4*cos(y2k4) )/(Re+y3k4)
     & -((rhot*y1k4)/(2.0*B))*(LoD)*cos(bank)
     & + (g/y1k4)*cos(y2k4)
        k34=-y1k4*sin(y2k4)

        y1k5=y1k1+0.1875*k11*dt+0.5625*k14*dt
        y2k5=y2k1+0.1875*k21*dt+0.5625*k24*dt
        y3k5=y3k1+0.1875*k31*dt+0.5625*k34*dt

        k15=(-rhot*y1k5**2)/(2.0*B) + g*sin(y2k5)
        k25=(-y1k5*cos(y2k5) )/(Re+y3k5)
     & -((rhot*y1k5)/(2.0*B))*(LoD)*cos(bank)
     & + (g/y1k5)*cos(y2k5)
        k35=-y1k5*sin(y2k5)

       y1k6=y1k1-0.428571429*k11*dt+0.285714286*k12*dt
     &  +1.71428571429*k13*dt-1.71428571429*k14*dt+1.14285714286*k15*dt
       y2k6=y2k1-0.428571429*k21*dt+0.285714286*k22*dt
     &  +1.71428571429*k23*dt-1.71428571429*k24*dt+1.14285714286*k25*dt
       y3k6=y3k1-0.428571429*k31*dt+0.285714286*k32*dt
     &  +1.71428571429*k33*dt-1.71428571429*k34*dt+1.14285714286*k35*dt

        k16=(-rhot*y1k6**2)/(2.0*B) + g*sin(y2k6)
        k26=(-y1k6*cos(y2k6) )/(Re+y3k6)
     & -((rhot*y1k6)/(2.0*B))*(LoD)*cos(bank)
     & + (g/y1k6)*cos(y2k6)
        k36=-y1k6*sin(y2k6)



        Vt = y1k1 + 0.011111111*(7.0*k11+32.*k13+12.*k14+32.*k15+7.*k16)
     &*dt

        gammat = y2k1
     &   + 0.011111111*(7.0*k21+32.*k23+12.*k24+32.*k25+7.*k26)*dt

        alt_t = y3k1
     &   + 0.011111111*(7.0*k31+32.*k33+12.*k34+32.*k35+7.*k36)*dt

        dS = -1./tan(y2k1)
     &     *( + 0.011111111*(7.0*k31+32.*k33+12.*k34+32.*k35+7.*k36)*dt)

c   endif
        Q(i+1)=.5*rhot*Vt**2
        t=t+dt

c   if(LoD .gt. 0.0) then

c   dvdt=-(7702.0**2-Vtold**2)/((Re+alt_told)*Lod*cosd(bank))
c     &          + g*sind(gammatold)

c   else
        dvdt=-((rhot*Vt**2)/(2.0*B)) + g*sin(gammat)
c   endif
        nmax=min1(nmax,dvdt/goo)

        if(nmax .eq. dvdt/goo) then
            Vnmax=y1k1
            hnmax=y3k1
        endif

c   end do




      S(i+1)=S(i)+dS
      V(i+1)=Vt
      alt(i+1)=alt_t
      gamma(i+1)=gammat
      rho(i+1)=rhot
      if(iplanet .eq. 1) then
c     units in W/cm^2
      qconv(i+1)=(1.74153e-4*(rho(i+1)/rnose)**.5*V(i+1)**3)/10000.0
      endif
      if(iplanet .eq. 2) then
c     units in W/cm^2
      qconv(i+1)=(1.9027e-4*(rho(i+1)/rnose)**.5*V(i+1)**3)/10000.0
      endif


      call RadHeat(iplanet,Vt,rhot,rnose,qrad(i+1))

c   qconv(i+1)=(1.83e-4*sqrt(rho(i+1)/rnose)*V(i+1)**3)/10000.0  !W/cm^2
!     call radheat(
      qtotal(i+1)=qtotal(i)+ ((qconv(i)+qconv(i+1))/2.0)*dt
     & + ((qrad(i)+qrad(i+1))/2.0)*dt  !J/cm^2

      Hrec=V(i+1)**2./2.
      ruCh=(qconv(i+1)*10000)/Hrec

      Hrec=Hrec/2325.854324
      ruCh=ruCh/4.882
c      patm=g*rhot*alt_t
C      patm=Po*rhot/rhoi
      surfP=(patm+Q(1)*Cp)/101325

c   if (LoD .gt. 0.0) then

c   dvdt=-(7702.0**2-V(i+1)**2)/((Re+alt(i+i))*LoD*cosd(bank))
c     &  +g*sind(gamma(i+1))

c   else

c      if(t(i) .ge. 44.5) then
c        write(*,*) "pause"
c      endif

      dvdt=-(rho(i+1)*V(i+1)**2)/(2.0*B) !*LoD*cosd(bank)!*(1+(.8346*cosd(bank)/.77)**2)
     &   +g*sin(gamma(i+1))
c   endif
      n(i+1)=dvdt/goo
      acc=(V(i+1)-V(i))/dt/goo
      if(alt(i+1) .gt. alt(i)+100) then
      write(*,*)'Flight path angle too shallow for
     &the given entry velocity'
      write(*,*)'Vehicle will skip out of the atmosphere,
     &make flight path angle larger'
        stop
      endif

      i=i+1
        p_atms = patm/101325
      if(t(i) >= tprnt) then
      write(1,105) t(i),alt(i),V(i),Patm,rho(i),gamma(i),n(i),qconv(i),
     &qrad(i),qtotal(i),S(i),Q(i)



      write(21,44) p_atms,tatm, V(i)




      write(11,200) t(i),Hrec,qradbc,ruCh,surfP,lam

      tprnt=tprnt+prntint

      endif


      end do

c   write(1,*)
c   write(1,115)'     Nmax in Earth g:   ',nmax
c   write(1,115)'Velocity at Nmax (m/s): ',Vnmax
c   write(1,125)'Altitude at Nmax (m):   ',hnmax
44    format(f15.10,f10.3,f10.3)
100   format(1x,f9.2,3x,f8.1,3x,f8.2,3x,f7.3,4x,f8.2,
     &5x,f10.2,5x,f10.2,5x,f10.2,3x,f10.2,3x,f10.2)
105   format(1x,f9.2,3x,f8.1,3x,f8.2,3x,ES12.5,3x,ES12.5,3x,f7.3,4x,
     &f8.2,5x,f10.2,5x,f10.2,5x,f10.2,3x,f10.2,3x,f10.2)
101   format(a11,f6.2,a3,f8.2,a7,f5.2,a5,f4.2,a5,f5.2,a7,f5.3,a1)
110   format(a140)
111   format(a171)
115   format(1x,a24,f8.2)
125   format(1x,a24,f9.2)
200   format(f9.2,3x,f9.2,3x,f7.2,3x,e11.4,3x,e11.4,3x,f3.1)
      stop
      end

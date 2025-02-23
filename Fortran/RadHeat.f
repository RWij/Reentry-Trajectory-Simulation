*********************************************************************
*	This program calculates the radiative heat transfer to a blunt 
*	body entry for Earth and Mars
*	- Requires that at POST input file has been read and that density, 
*		and velocity are available and stored in Vfi(i), and Rho_P(i)
*********************************************************************

c	subroutine RadHeat
	subroutine RadHeat(iplanet,v,rho,rnose,qrad) 
c
c      Tauber Sutton Radiation Calculation
c
c   "Stagnation Point Radiative Heating Relations for Earth and Mars"
c      From J. Spacecraft, Vol 28. No. 1, pp 40-43
c
	integer size,iplanet
       real*8 v,rho,rnose,qrad

      dimension vedat(19),fevdat(19),vmdat(17),fmvdat(17),rho0(44)
      data vedat/9000,9250,9500,9750,10000,10250,10500,10750,11000,
     & 11500,12000,12500,13000,13500,14000,14500,15000,15500,16000/
      data fevdat/1.5,4.3,9.7,19.5,35.,55.,81.,115.,151.,238.,359.,
     & 495.,660.,850.,1065.,1313.,1550.,1780.,2040./
      data vmdat/6000,6150,6300,6500,6700,6900,7000,7200,7400,7600,
     & 7800,8000,8200,8400,8600,8800,9000/
      data fmvdat/0.2,1.0,1.95,3.42,5.1,7.1,8.1,10.2,12.5,14.8,17.1,
     & 19.2,21.4,24.1,26.0,28.9,32.8/

c	include 'Variables.inc'

c
c    Earth
c
      if(iplanet.eq.1) then
c	do j=1,mTraj

c		v=Vfi(j)
c		rho=Rho_P(j)

       a = 1.072e+06 * v **(-1.88) * rho **(-0.325)
       a = amin1(a,1.0)
       if(rnose.ge.1.0.and.rnose.le.2.0) a = amin1(a,0.6)
       if(rnose.ge.2.0.and.rnose.le.3.0) a = amin1(a,0.5)
       b = 1.22
       c = 4.736e+04
c
c        find Fe(V)
c 
        if(v.le.vedat(1)) then
            fv = fevdat(1)
            go to 100
        endif
        do 10 i = 2,19 
          if(v.le.vedat(i)) then
              dfdv = (fevdat(i)-fevdat(i-1))/(vedat(i)-vedat(i-1))
              fv = fevdat(i-1) + (v-vedat(i-1))*dfdv
              go to 100
          endif
 10     continue
          fv = fevdat(19)
 100    continue 
       if(v.lt.9000) c = 0.

       qrad = c * rnose**a * rho**b * fv
	 
c	 end do
	 
	endif !end Earth
c
c     MARS
c
       if(iplanet.eq.2) then

c	do j=1,mTraj

c		v=Vfi(j)
c		rho=Rho_P(j)

       a = 0.526
       b=1.19
       c=2.35e+04
c
c        find Fm(V)
c 
        if(v.le.vmdat(1)) then
            fv = fmvdat(1)
            go to 200
        endif
        do 20 i = 2,17 
          if(v.le.vmdat(i)) then
              dfdv = (fmvdat(i)-fmvdat(i-1))/(vmdat(i)-vmdat(i-1))
              fv = fmvdat(i-1) + (v-vmdat(i-1))*dfdv
              go to 200
          endif
 20     continue
          fv = fmvdat(17)
 200    continue 
       if(v.lt.5500) c = 0.

	qrad = c * rnose**a * rho**b * fv
c	end do

      endif !End Mars
c
*      qrad(i) = c * rnose**a * rho**b * fv
c      chrad = qrad/(0.5*rho*v**3)
c
c     write(6,*)'a,b,c,fv = ',a,b,c,fv
c     write(6,*)' Radiative Heating = ',qrad,' W/cm^2'
c     write(6,*)' Radiative Heat transfer Coef Chr = ',chrad
c
		



      
	return
      end
c
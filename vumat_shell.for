!************************************************************************
!
! A simple VUMAT for a single fibre bundle material
!
! Sumedh Joshi- Mini thesis
        subroutine vumat(
     +     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
     +     stressNew, stateNew, enerInternNew, enerInelasNew)
      
      include 'vaba_param.inc'

      dimension props(nprops), density(nblock), coordMp(nblock,*),
     +     charLength(nblock), strainInc(nblock,ndir+nshr),
     +     relSpinInc(nblock,nshr), tempOld(nblock),
     +     stretchOld(nblock,ndir+nshr),
     +     defgradOld(nblock,ndir+nshr+nshr),
     +     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     +     stateOld(nblock,nstatev), enerInternOld(nblock),
     +     enerInelasOld(nblock), tempNew(nblock),
     +     stretchNew(nblock,ndir+nshr),
     +     defgradNew(nblock,ndir+nshr+nshr),
     +     fieldNew(nblock,nfieldv),
     +     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     +     enerInternNew(nblock), enerInelasNew(nblock)

      character*80 cmname       
      
      integer km
      
      real*8 Iden(3,3), La(3,3)
      real*8 U(3,3),U_inv(3,3),T_tau(3,3)
      
      real*8 muein,kappaein,wa,a
      
      real*8 zero,one,two,three,half,third,four,Pi,two_third
      real*8 det,det1,det2,det3,detinv,pres
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,half=0.5d0,
     +          third=1.0d0/3.0d0,two_third=2.d0/3.0d0,Pi=3.1415926d0)
      
      
      
      ! Read material properties needed here

      muein =    PROPS(1)
      kappaein=  PROPS(2)
      wa =       PROPS(3)
      a     =    PROPS(4)
	  
     
	  
	   !**************Call the first structural tensor*********************!
          
            
           La(1,1) = cos(Pi*a)**2.d0;
           La(1,2) = sin(Pi*a)*cos(Pi*a);
           La(1,3) = 0;
           La(2,1) = sin(Pi*a)*cos(Pi*a);
           La(2,2) = sin(Pi*a)**2.d0;
           La(2,3) = 0;
           La(3,1) = 0;
           La(3,2) = 0;
           La(3,3) = 0;
           
      !*****************specify the identity tensor*******************!
	      
		   Iden(1,1) = 1.d0
		   Iden(1,2) = 0.d0
		   Iden(1,3) = 0.d0
		   Iden(2,1) = 0.d0
		   Iden(2,2) = 1.d0
		   Iden(2,3) = 0.d0
		   Iden(3,1) = 0.d0
		   Iden(3,2) = 0.d0
		   Iden(3,3) = 1.d0

      
      ! START LOOP OVER MATERIAL POINTS:
      
      do km=1,nblock
          U(1,1) = StretchNew(km,1)
          U(2,2) = StretchNew(km,2)
          U(1,2) = Stretchnew(km,4)
          U(2,1) = U(1,2)
          U(3,3) = 1.d0/((U(2,2)*U(1,1))-(U(1,2)*U(2,1)))
          U(1,3) = zero
          U(2,3) = zero
          U(3,1) = zero
          U(3,2) = zero
           
	 
		  
	  !***********************Compute the Inverse of the Stretch tensor matrix***********!
	  
	  ! Calculate the determinanat inverse
	  
	  
       det1 =   U(1,1)* U(2,2)* U(3,3) - U(1,1)* U(2,3)* U(3,2)
       det2 = - U(1,2)* U(2,1)* U(3,3) + U(1,2)* U(2,3)* U(3,1)
       det3 =   U(1,3)* U(2,1)* U(3,2) - U(1,3)* U(2,2)* U(3,1)
       det  =   det1 + det2 + det3
       detinv = 1.d0/det         
                     
         
      ! Calculate the inverse of the matrix
         U_inv(1,1) =  +detinv * (U(2,2)*U(3,3) - U(2,3)*U(3,2))
         U_inv(2,1) =  -detinv * (U(2,1)*U(3,3) - U(2,3)*U(3,1))
         U_inv(3,1) =  +detinv * (U(2,1)*U(3,2) - U(2,2)*U(3,1))
         U_inv(1,2) =  -detinv * (U(1,2)*U(3,3) - U(1,3)*U(3,2))
         U_inv(2,2) =  +detinv * (U(1,1)*U(3,3) - U(1,3)*U(3,1))
         U_inv(3,2) =  -detinv * (U(1,1)*U(3,2) - U(1,2)*U(3,1))
         U_inv(1,3) =  +detinv * (U(1,2)*U(2,3) - U(1,3)*U(2,2))
         U_inv(2,3) =  -detinv * (U(1,1)*U(2,3) - U(1,3)*U(2,1))
         U_inv(3,3) =  +detinv * (U(1,1)*U(2,2) - U(1,2)*U(2,1))
         
         
         
      !***********************Compute the Cauchy Stress**********************************!
      
         pres  =    (kappaein*2.d0*(1.d0-wa)*third*U(3,3)**2.d0) 
     +             +(muein*2.d0*(wa-1.d0)*third*U_inv(3,3)**2.d0)
         
         T_tau =    kappaein*2.d0*(wa*matmul(U,matmul(La,U)))
     +            + kappaein*2.d0*(1.d0-wa)*third*matmul(U,U)         
     +            - muein*2.d0*wa*matmul(U_inv,matmul(La,U_inv))
     +            + muein*2.d0*(wa-1.d0)*third*matmul(U_inv,U_inv)
     +            - pres*Iden            
        
      !**********************Update the stress vector************************************!
        
		
            stressNew(km,1) = T_tau(1,1)
			stressNew(km,2) = T_tau(2,2)
			stressNew(km,3) = T_tau(3,3)
			stressNew(km,4) = T_tau(1,2)
			stressNew(km,5) = T_tau(2,3)
			stressNew(km,6) = T_tau(1,3)
              
              
              
         
        
         print *, "h" , pres
      !  print *, "T" , T_tau
         print *, "U" , U

        
	
      enddo ! end loop over material points

      end subroutine vumat




      
          
           





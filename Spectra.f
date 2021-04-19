C! ######################################################################
C! THIS PROGRAM SIMULATES THE SCATTERING OF HIGH ENERGY PHOTONS FROM    # 
C! NON-LINEAR THOMSON SCATTERING: THE INTERACTION OF A FREE AND HIGHLY  #
C! RELATIVISTIC ELECTRON WITH A REALISTIC 3-DIMENSONAL LASER FIELD.     #
C! THE LASER FIELD CAN BE CALCULATED PARAXIALLY BY SETTING eps = 0      #
C! ######################################################################

C!#######################################################################
C! THE CODE OUTPUT IS STORED IN TWO FILES "traj.dat" WHICH IS THE       # 
C! ELECTRON TRAJECTORY AND "Efield.dat" WHICH STORES THE RADIATED       # 
C! ELECTRIC FIELD. THE FORMAT OF THE OUTPUT FILES IS AS FOLLOWS         #  
C! retarded time, laser time, X,Y,Z --> FOR "traj.dat"                  # 
C! retarded time, laser time, Ex, Ey, Ez --> for "Efield.dat"           #
C!#######################################################################
      
      PROGRAM NL_THOMSON_SCATTERING
                                                                  
      IMPLICIT          NONE

C----Parameters
      INTEGER           COUNTER,INTPHI,ENDPOINT,NUMF,e_start,e_end,
     &                  e_max,CT,E_N0,TOTAL_ES,J,I,start_freq_iter,
     &                  START_THETA_ITER,MAX_THETA_ITER,MAX_PHI_ITER,
     &                  END_THETA_ITER,END_PHI_ITER,START_PHI_ITER,NF
      DOUBLE PRECISION  FNUM,s0_F,Z_RR,EENERGY,GAMA_I,W_LENGTH,PI_P,
     &                  Z_LASER,Z_ELECTRON,T0_fs,I_Wcm2,BETA_i,df,
     &                  freq_eV,ko(4),kt(4),kth(4),kfr(4),kfv(4),ksx(4),
     &                 GOLD(6),told,dtheta,int_ang,dphi_deg,PHI_e_START,
     &               PHI_e_END,THETA_e_START,THETA_e_END,dphi_rad,dphi,
     &               THETA_e_END_mrad,PHI_e_START_deg,PHI_e_END_deg,
     &               dtheta_mrad,np,THETA_e_START_mrad
      
      INTEGER*4		 POINTS
      
C!######## Setting code Runge-Kutta parameters. #######################      
      
      PARAMETER (POINTS = 180000,df = 368.7750d0,MAX_PHI_ITER = 2000)

C!######## SETTING LASER FIELD PARAMETERS #############################     
      
      PARAMETER (START_THETA_ITER = 1,MAX_THETA_ITER = 2000,
     &          START_PHI_ITER = 1,e_max = 20804,NF = 4000)
	
C!####### LOCAL SCALAR VARIABLES FOR RK-CODE AND MAIN CODE #############
      LOGICAL           CALCULATE, OUTSIDE
      INTEGER           L, JJ ,LL, KJ,KK, K1, KL, MM, NL,TT,LT  
      
      DOUBLE PRECISION  LAMBDA,TEMPSUM,TG, TG_START,TG_INC, TG_WANT, 
     &	                TG_LAST,THETA, S_Z,R0y, R0z,R0,rz,rx,ry,R0x,
     &                  eps,freq_MeV,PERIOD,PI,C,W,E0,T0,Z_R,s0,BcKz,
     &                  PX_e(e_max),PY_e(e_max),PZ_e(e_max),
     &                  X_e(e_max),Y_e(e_max),Z_e(e_max),int_ang_deg,
     &                  X_um,Y_um,Z_um,pz_i,DENSITY,temp_phi1,
     &                  temp_phi2,temp_theta1,temp_theta2,
     &                  theta_l(MAX_PHI_ITER,MAX_THETA_ITER), 
     &                  phi_l(MAX_PHI_ITER,MAX_THETA_ITER),
     &                  theta_e(MAX_PHI_ITER,MAX_THETA_ITER),
     &                  phi_e(MAX_PHI_ITER,MAX_THETA_ITER)
   


C!####### LOCAL ARRAYS FOR RK-CODE AND MAIN PROGRAM ###################
      
      DOUBLE PRECISION  AREA,PHI,PSUM,G(6), DG(6), G_START(6),
     &                  X(POINTS),Y(POINTS),Z(POINTS), PX(POINTS), 
     &	                PY(POINTS), PZ(POINTS),T_R0(POINTS),EW(6,NF),
     &                  EX_RAD(POINTS), EY_RAD(POINTS), EZ_RAD(POINTS),
     &                  T(POINTS),FREQUENCY(NF),POWSPEC(NF),
     &                  EX1,EX2,EY1,EY2,EZ1,EZ2,RMSIN,RMCOS,
     &                  par(16),NTHETA,NPHI,X_f_eframe,Y_f_eframe,
     &                  Z_f_eframe,PX_f_eframe,PY_f_eframe,PZ_f_eframe
    
       CHARACTER     label(16)*250
       REAL etime, elapsed(2),Total                      

C!##### INTRINSIC FUNCTIONS ############################################
      
       INTRINSIC         ATAN, COS, SIN,ACOS

C!##### COMMON PARAMETERS ##############################################

        COMMON          PI, C, W, E0, T0, Z_R, s0,eps, BckZ, 
     &                  DENSITY,THETA,CT 


C!#### Set the initial conditions.######################################  
C #--E in +X and B in +y direction and k along +z direction      !     #
C-#--G1, G2, G3 is PX, PY, PZ                                    !     # 
C-#--The problem is set up in atomic units                       !     #
C!######################################################################

!      Read in command-line arguments
        CHARACTER(LEN=32)  arg,mystart,mystop,myname
        CHARACTER(LEN=100) pow_spec_cohf,pow_spec_incf,
     &  finaltraj_f
!        INTEGER E_START,E_STOP

        call getarg(0,arg)
        read(arg,'(A32)') myname
        if (iargc() /= 2) then
                print*, "Usage is '",trim(myname),"' start_electron",
     &           " stop_electron"
                STOP
        endif

        CALL GETARG(1,arg)
        READ(arg,*) e_start
        READ(arg,*) mystart

        CALL GETARG(2,arg)
        read(arg,*) e_end
        read(arg,*) mystop


        pow_spec_incf = 'powspec_inc.dat.'//trim(mystart)//'.'//
     &  trim(mystop)
        pow_spec_cohf = 'powspec_coh.dat.'//trim(mystart)//'.'//
     &  trim(mystop)
        finaltraj_f = 'e_final_tra.dat.'//trim(mystart)//'.'//
     &  trim(mystop)


        print*,pow_spec_incf
        print*,pow_spec_cohf
        print*,finaltraj_f
        print*,'e_start=',e_start
        print*,'e_end=',e_end

       OPEN (UNIT=65, FILE='parameters.dat',STATUS='OLD',
     &     FORM='FORMATTED')
       DO 100 J = 1,14
           read(65,*) par(J),label(J)
100     CONTINUE
      W_LENGTH = par(1); FNUM = par(2);I_Wcm2 = par(3); T0_fs = par(4);
      INT_ANG_DEG = par(5);THETA_e_START_mrad = par(6);
      theta_e_end_mrad = par(7);NTHETA = NINT(par(8)); 
      phi_e_start_deg = par(9);phi_e_end_deg=par(10);NPHI=NINT(par(11));
      start_freq_iter=NINT(par(12));NUMF=NINT(par(13));
      np=par(14);
      

      PI = 4.0D0*ATAN(1.0D0)     
      PHI_e_START = PHI_e_START_deg*PI/180.0d0
      PHI_e_END = PHI_e_END_deg*PI/180.0d0
      THETA_e_END = THETA_e_END_mrad/1000.0d0
      THETA_e_START = THETA_e_START_mrad/1000.0d0
      IF (NTHETA.gt.1) then
         dtheta = (THETA_e_END_mrad - THETA_e_START_mrad)/(NTHETA - 1)
         dtheta = dtheta/1000.0d0
      ELSE
         dtheta = 0.0d0
      END IF     
      END_THETA_ITER = NINT(NTHETA)
      IF (NPHI.GT.1) then
         dphi = (PHI_e_END - PHI_e_start)/(NPHI - 1)
      ELSE
         dphi = 0.0d0
      END IF
      END_PHI_ITER= NINT(NPHI)
c      int_ang_deg = 180.0d0
      int_ang = PI*int_ang_deg/180.0d0
      RMSIN = dsin(int_ang)
      RMCOS = dcos(int_ang)
      DENSITY = 0.0d0
      BckZ = 8.0d0

C# Read and write electron beam phase space.################### 
C# ebeam.dat is file for e-beam phase space data
       EENERGY = 300.00d0     
       GAMA_I = EENERGY/0.5110d0 + 1.0d0
       BETA_i = dsqrt(1 - 1/GAMA_I**2.0d0)
       DO 47 kk = 1,END_PHI_ITER
          DO 36 J = 1,END_THETA_ITER
          theta_e(kk,j) = theta_e_start + (J-1)*dtheta
          phi_e(kk,j) = phi_e_start + (kk-1)*dphi
          temp_phi1 = dsin(theta_e(kk,j))*dsin(phi_e(kk,j))
          temp_phi2 = dsin(theta_e(kk,j))*dcos(phi_e(kk,j))*
     &                dcos(int_ang) + dcos(theta_e(kk,j))*dsin(int_ang)
          phi_l(kk,j) = ATAN(temp_phi1/temp_phi2)
          temp_theta1 = dcos(theta_e(kk,j))*dcos(int_ang)
      temp_theta2 = dsin(theta_e(kk,j))*dcos(phi_e(kk,j))*dsin(int_ang)
         theta_l(kk,j) = ACOS(temp_theta1 - temp_theta2)
        WRITE(*,*) theta_e(kk,j),phi_l(kk,j)      
36     Continue
47     CONTINUE       
       OPEN (UNIT=25, FILE='ephase_space.dat',STATUS='UNKNOWN',
     &     FORM='FORMATTED')
      OPEN (UNIT=41, FILE=pow_spec_incf,STATUS='UNKNOWN',
     &  FORM='FORMATTED')
      OPEN (UNIT=42, FILE=pow_spec_cohf,STATUS='UNKNOWN',
     &  FORM='FORMATTED')
       OPEN (UNIT=43, FILE=finaltraj_f,STATUS='UNKNOWN',
     &     FORM='FORMATTED')

       
       DO 35 J = 1,e_max
       READ(25,*) X_e(j),Y_e(j),Z_e(j),Px_e(j),Py_e(j),Pz_e(j)
35     CONTINUE 
        

C####### END of reading and writing electron beam phase space.############
      WRITE(*,*) 'f_number used is', FNUM
C!##### Defining Laser Parameters ######################################    
      C = 137.0359895d0           !    SPEED OF LIGHT                 !#    
      LAMBDA = W_LENGTH/0.052917721080d0
      WRITE(*,*) 'Wavelength of laser is', LAMBDA
      W = 2.0d0*PI*C/LAMBDA         !  LASER ANGULAR FREQUENCY        !#    
      PERIOD = 2.0d0*PI/W                                             !#
      E0 = DSQRT((I_Wcm2)/(3.55*1.0d+16)) ! PEAK LASER E FILED        !#
      T0 = (T0_fs)/(0.0241888432650d0*DSQRT(2.0d0*abs(log(2.0d0))))   !#
      s0_F=FNUM*2.0d0*LAMBDA/PI
      s0 = s0_F                                                       !#
      Z_R = W*s0**2.0d0/(c*2.0d0) ! RAYLEIGH RANGE                    !#
      np = 1.0d0
      eps = LAMBDA*np/(2.0d0*PI*s0) ! NP PARAMETER                 !#
      WRITE(*,*) 'For Field = ',E0                                    !#
      WRITE(*,*) 'Temporal duration',T0                         
      WRITE(*,*) 'beam waist radius is',s0
C!#### END OF DEFINING LASER PARAMETERS.################################        
        
C!#### put laser at z_laser when electron is at z_electron ###############

      Z_LASER = -C*2.0d0*T0 
      Z_ELECTRON = C*1.50d0*T0
      WRITE(*,*) 'temporal laser dist. from peak is', Z_LASER/C
      PZ_i = GAMA_I*BETA_i*C


        DO 4000 KK = START_PHI_ITER,END_PHI_ITER
          DO 3000 TT = START_THETA_ITER,END_THETA_ITER
             THETA = theta_l(kk,tt)
             phi = phi_l(kk,tt)
             frequency(start_freq_iter) = (start_freq_iter - 1)*df
             powspec(1) = 0.0d0
             EW(1,1) = 0.0d0
             EW(2,1) = 0.0d0
             EW(3,1) = 0.0d0
             EW(4,1) = 0.0d0
             EW(5,1) = 0.0d0
             EW(6,1) = 0.0d0
      DO 17 J = start_freq_iter + 1,NUMF
         Frequency(J) = Frequency(J-1) + df
         powspec(J) = 0.0d0
         EW(1,J) = 0.0d0
         EW(2,J) = 0.0d0
         EW(3,J) = 0.0d0
         EW(4,J) = 0.0d0
         EW(5,J) = 0.0d0
         EW(6,J) = 0.0d0
17    CONTINUE

C#### Loop across electrons starting from electron number####################
C#### e_start ending at electron number e_end ###############################      
      
      TOTAL_ES = 0
      DO 2000 E_N0 = e_start,e_end 
         TOTAL_ES = TOTAL_ES + 1
         WRITE(*,*) 'total electrons = ', TOTAL_ES
         AREA = 0.0d0


C!########## INITIALIZING ELECTRON INITIAL CONDITIONS.##################
c!########## G_START(1) is Px,G_START(2) is Py, G_START(3) is Pz #######
c!######### G_START(4) X, G_START(5) is Y, and G_START(6) is Z #########
c           Z_e(E_N0) =  Z_LASER 
            G_START(1) =(PX_e(E_N0)*RMCOS + PZ_e(E_N0)*RMSIN)/6.0d0
c           G_START(1) = PX_e(E_N0)
           G_START(2) = PY_e(E_N0)/6.0d0
           G_START(3) = PZ_e(E_N0)*RMCOS - PX_e(e_N0)*RMSIN 
c            G_START(3) = -PZ_e(E_N0)
         G_START(4) =  Z_LASER*RMSIN + X_e(E_N0)*RMCOS + Z_e(E_N0)*RMSIN 
c          G_START(4) = X_e(E_N0)
         G_START(5) = Y_e(E_N0)
         G_START(6) = Z_LASER*RMCOS + Z_e(E_N0)*RMCOS - X_e(E_N0)*RMSIN
c          G_START(6) = -Z_LASER
C!########## END of Electron initial conditions ########################
         
c!####### START LOOP FOR ELECTRON TRAJECTORY CALCULATION.###############        
	  TG_START = Z_LASER/c
          TG_LAST = 2.0d0*abs(TG_START) 
	  TG_INC = (TG_LAST-TG_START)/POINTS
          COUNTER = COUNTER + 1
               WRITE(*,*) 'trajectory calculation begins '
               WRITE(*,*) 'Trajectory count is ', COUNTER
          L=0
C!####### Solve electron dyanmics with a runge-kutta algorithm.############
           OUTSIDE = .false.
           DO WHILE ((L .LT. POINTS).and.(OUTSIDE .eqv. .false.))
            L=L+1
            TG_WANT = TG_LAST + (L-POINTS)*TG_INC

      CALL F(TG_START,told,G_START,DG,ko,kt,kth,kfr,kfv,ksx,GOLD,TG_INC)
            G(1) = GOLD(1) +
     &         (ko(1) + 2.0d0*ko(2) + 2.0d0*ko(3) + ko(4))*TG_INC/6.0d0

             G(2) = GOLD(2) +
     &         (kt(1) + 2.0d0*kt(2) + 2.0d0*kt(3) +kt(4))*TG_INC/6.0d0
            G(3) = GOLD(3) +
     &     (kth(1) + 2.0d0*kth(2) + 2.0d0*kth(3) + kth(4))*TG_INC/6.0d0
         G(4) = GOLD(4) +
     &      (kfr(1) + 2.0d0*kfr(2) + 2.0d0*kfr(3) +kfr(4))*TG_INC/6.0d0
             G(5) = GOLD(5) +
     &      (kfv(1) + 2.0d0*kfv(2) + 2.0d0*kfv(3) +kfv(4))*TG_INC/6.0d0
             G(6) = GOLD(6) +
     &      (ksx(1) + 2.0d0*ksx(2) + 2.0d0*ksx(3) +ksx(4))*TG_INC/6.0d0
            PX(L) = G(1)
            PY(L) = G(2)
            PZ(L) = G(3)
            X(L) = G(4)
            Y(L) = G(5)
            Z(L) = G(6)
            T(L) = TG_WANT
            G_START(1) = G(1)
            G_START(2) = G(2)
            G_START(3) = G(3)
            G_START(4) = G(4)
            G_START(5) = G(5)
            G_START(6) = G(6)
            TG_START = told + TG_INC
          ENDDO
          ENDPOINT = L        
          CALL RAD_ECALC (X,Y,Z, PX,PY,PZ, T,T_R0, EX_RAD,
     &                  EY_RAD, EZ_RAD, THETA,PHI,ENDPOINT)

          CALL fourier(T_R0,EX_RAD,EY_RAD,EZ_RAD,POINTS,NUMF,df,
     &              FREQUENCY,powspec,EW,start_freq_iter)
 
         IF ((KK.EQ.START_PHI_ITER).AND.(TT.EQ.START_THETA_ITER)) THEN
            X_f_eframe = X(ENDPOINT)*RMCOS - Z(ENDPOINT)*RMSIN
            Y_f_eframe = Y(ENDPOINT)
            Z_f_eframe = Z(ENDPOINT)*RMCOS - X(ENDPOINT)*RMSIN
            PX_f_eframe = PX(ENDPOINT)*RMCOS - PZ(ENDPOINT)*RMSIN
            PY_f_eframe = PY(ENDPOINT)
            PZ_f_eframe = PZ(ENDPOINT)*RMCOS - PX(ENDPOINT)*RMSIN  
            WRITE(43,'(6e20.10)') X_f_eframe,Y_f_eframe,Z_f_eframe,
     &                         PX_f_eframe,PY_f_eframe,PZ_f_eframe
          END IF
2000   CONTINUE
      
        DO 670 L = start_freq_iter,NUMF
           powspec(L) = powspec(L)*(c*1.0d+18)/(4*PI)
           freq_eV = frequency(L)*27.21
           freq_MeV = freq_eV/(1.0d+6)
           WRITE(41,'(4e20.10)') phi_e(kk,tt),theta_e(kk,tt),
     &                           freq_MeV,powspec(L) 
           EX1 = EW(1,L)*dsqrt(c/PI)*(1.0d+9)/2.0d0
           EY1 = EW(2,L)*dsqrt(c/PI)*(1.0d+9)/2.0d0
           EZ1 = EW(3,L)*dsqrt(c/PI)*(1.0d+9)/2.0d0
           EX2 = EW(4,L)*dsqrt(c/PI)*(1.0d+9)/2.0d0
           EY2 = EW(5,L)*dsqrt(c/PI)*(1.0d+9)/2.0d0
           EZ2 = EW(6,L)*dsqrt(c/PI)*(1.0d+9)/2.0d0
           WRITE(42,'(9e20.10)') phi_e(kk,tt),theta_e(kk,tt),
     &                           freq_MeV,EX1,EX2,EY1,EY2,EZ1,EZ2
670    CONTINUE
3000   CONTINUE
4000   CONTINUE

      CLOSE(UNIT=21)
      CLOSE(UNIT=20)
      CLOSE(UNIT=41)
      CLOSE (UNIT=42)
      CLOSE(UNIT=43)
      STOP
      END

C---------------------SUBROUTINES


      SUBROUTINE F(t,told,gs,dg,ko,kt,kth,kfr,kfv,ksx,GOLD,dt)

      DOUBLE PRECISION  gs(6), dg(*), bx_field, by_field,bz_field,
     &                  gama_c,r,t,ko(4),kt(4),kth(4),ex_las,
     &                  kfr(4),kfv(4),ksx(4),dt,XOLD,YOLD,ZOLD,
     &                  PXOLD,PYOLD,PZOLD,GOLD(6),e_las(4),ey_las,
     &                  ez_las,bx_las,by_las,bz_las,told

      DOUBLE PRECISION  PI, c, w, E0, T0, X_END, LBfield, s0,
     &                  eps, Qe, Me,PcrsB_x,PcrsB_y,
     &                  PcrsB_z,pdotE,Rx,Ry,Rz,PcrsB_sqr,
     &                  Elas_mag,T1,Rconst,Bckz,density

      INTEGER           i,j,k
      COMMON            PI, c, w, E0, T0, X_END, LBfield, s0,
     &                  eps, Qe, Me,Bckz,density

      PXOLD = gs(1)
      PYOLD = gs(2)
      PZOLD = gs(3)
      XOLD = gs(4)
      YOLD = gs(5)
      ZOLD = gs(6)
      told = t
      GOLD(1) = PXOLD
      GOLD(2) = PYOLD
      GOLD(3) = PZOLD
      GOLD(4) = XOLD
      GOLD(5) = YOLD
      GOLD(6) = ZOLD
      told = t
      DO 10 i = 1,4
           IF ((i.gt.1).and.(i.lt.4)) then
              gs(1) = PXOLD + (dt/2)*ko(i-1)
              gs(2) = PYOLD + (dt/2)*kt(i-1)
              gs(3) = PZOLD + (dt/2)*kth(i-1)
              gs(4) = XOLD + (dt/2)*kfr(i-1)
              gs(5) = YOLD + (dt/2)*kfv(i-1)
              gs(6) = ZOLD + (dt/2)*ksx(i-1)
               t = told + dt/2.0d0
         END IF
         IF (i.eq.4) then
             gs(1) = PXOLD + (dt)*ko(i-1)
             gs(2) = PYOLD + (dt)*kt(i-1)
             gs(3) = PZOLD + (dt)*kth(i-1)
             gs(4) = XOLD + (dt)*kfr(i-1)
              gs(5) = YOLD + (dt)*kfv(i-1)
             gs(6) = ZOLD + (dt)*ksx(i-1)
             t = told + dt
         END IF
       gamma_c = DSQRT(gs(1)**2.0d0 + gs(2)**2.0d0 + gs(3)**2.0d0 +
     &                c**2.0d0 )

      CALL calcEfield(gs(4),gs(5),gs(6),t,e_las)
      ex_las = e_las(1)
      ey_las = e_las(2)
      ez_las = e_las(3)*gs(4)
      bx_las = e_las(2)
      by_las = e_las(4)
      bz_las = e_las(3)*gs(5)

      pdotE = gs(1)*ex_las + GS(2)*ey_las + GS(3)*ez_las
      PcrsB_x = GS(2)*bz_las - GS(3)*by_las
      PcrsB_y = GS(3)*bx_las - GS(1)*bz_las
      PcrsB_z = GS(1)*by_las - GS(2)*bx_las
      PcrsB_sqr = PcrsB_x**2.0d0 + PcrsB_y**2.0d0 + PcrsB_z**2.0d0
      Elas_mag = ex_las**2.0d0 + ey_las**2.0d0 + ez_las**2.0d0
      T1 = 2.0d0*gamma_c*(0.50d0*gamma_c*Elas_mag + ex_las*PcrsB_x +
     &                   ey_las*PcrsB_y + ez_las*PcrsB_z)
      Rconst = (2.0d0/(3.0d0*gamma_c*(c**6.0d0)))*(PcrsB_sqr +T1 -pdotE)
      Rx = -GS(1)*Rconst*0.0d0
      Ry = -GS(2)*Rconst*0.0d0
      Rz = -GS(3)*Rconst*0.0d0
C----differential equations of motion

      DG(1) = -(4.0d0/3.0d0)*PI*BckZ*density*gs(4) -
     &       ex_las -(gs(2)*bz_las-gs(3)*by_las)/gamma_c + Rx
      DG(2) = -(4.0d0/3.0d0)*PI*BckZ*density*gs(5) -
     &       ey_las -(gs(3)*bx_las-gs(1)*bz_las)/gamma_c + Ry
      DG(3) = -(4.0d0/3.0d0)*PI*BckZ*density*gs(6) -
     &       ez_las -(gs(1)*by_las-gs(2)*bx_las)/gamma_c + Rz

      DG(4) = gs(1)*c/gamma_c
      DG(5) = gs(2)*c/gamma_c
      DG(6) = gs(3)*c/gamma_c

C----differential equations of motion
             ko(i) = DG(1)
             kt(i) = DG(2)
             kth(i) = DG(3)
             kfr(i) = DG(4)
             kfv(i) = DG(5)
             ksx(i) = DG(6)
10      CONTINUE
        RETURN
        END
        
         SUBROUTINE calcEfield(X, y, z, t, e_las)

      DOUBLE PRECISION  X,y,z,t,e_las(4),s_z,amplitude,phase,
     &                  alpha, r2, ExDiff, ByDiff
      DOUBLE PRECISION  PI, c, w, E0, T0, Z_R, s0, eps
      COMMON            PI, c, w, E0, T0, Z_R, s0, eps

C fields up to first order in eps = 1/(k0s0) is included
C Z_R is the raleigh length
C e_las(1) is Ex[0,2]
C e_las(2) is Ey[2]=Bx[2]
c e_las(3) is Ez[1,3]/X = Bz[1,3]/y
c e_las(4) is By[0,2]

      s_z = s0*DSQRT( 1.0d0 + ( z/Z_R )**2.0d0 )
      alpha = atan(z/Z_R)
c       alpha = 0.0d0
      r2 = X**2.0d0+y**2.0d0
      amplitude = E0*(s0/s0)*dexp(- r2*1.0d0/s_z**2.0d0 )*
     &            dexp(-((z/c-t)*1.0d0/T0)**2.0d0)
       phase = w*(t-z/c)-z*r2*1.0d0 /
     &        (Z_R*s_z**2.0d0)
c      phase=2.0d0*(w/3)*(1+dexp(-0.028d0*((t/T0)**2.0d0))/2.0d0)*(t-z/c)
c      phase = w*(t-z/c)
      ExDiff = amplitude*eps**2.0d0*(2.0d0*X**2.0d0+r2)*
     &      dcos(phase+3.0d0*alpha)/s_z**2.0d0
      ByDiff = amplitude*eps**2.0d0*(2.0d0*y**2.0d0+r2)*
     &      dcos(phase+3.0d0*alpha)/s_z**2.0d0
      e_las0 = amplitude*( dcos(phase+alpha) - eps**2.0d0*
     &     r2**2.0d0*dcos(phase+4.0d0*alpha)/
     &     (s_z**3.0d0*s0) )
      e_las(1) = e_las0 + ExDiff
      e_las(2) = 2.0d0*amplitude*eps**2.0d0*X*y*dcos(phase+
     &           3.0d0*alpha)/s_z**2.0d0
      e_las(3) = -2.0d0*amplitude*(eps/s_z)*( dsin(phase +
     &    2.0d0*alpha) + (eps/s_z)**2.0d0*r2*
     &    ( 3.0d0*dsin(phase+4.0d0*alpha)-dsin(phase + 5.0d0*
     &    alpha)*r2/(s_z*s0) ) )
        e_las(4) = e_las0 + ByDiff
      RETURN
      END
C-----------------bcw
      SUBROUTINE RAD_ECALC (X,y,z,px,py,pz,t,t_r0,
     &                      ex_rad,ey_rad,ez_rad,theta,PHI,points)
        
      DOUBLE PRECISION  X(*),y(*),z(*),px(*),py(*),pz(*),PHI,
     &                 t(*),t_r0(*), theta, ex_rad(*), ey_rad(*),
     &          ez_rad(*),R0, R0x, R0y, R0z, r, rx, ry, rz, nx, ny, nz,
     &                  nbeta, denom, n_a, numX, numY, numZ
        
      INTEGER 		i, points
      DOUBLE PRECISION 	ax,ay,az,ex_las,ey_las,ez_las,bx_las, by_las,
     &                  bz_las,pdotx, pdoty, pdotz, gamma_c,e_las(4),  
     &                  vdotE
 
      DOUBLE PRECISION  PI, c, w, E0, T0, Z_R, s0,
     &                  eps, BckZ, density
      COMMON            PI, c, w, E0, T0, Z_R, s0,
     &                  eps, BckZ, density
 
C----E_rad observed 10^9 AU length from focua, xz plane, phi=0
      R0  = 1.0d+9
      R0x = R0*dsin(theta)*dcos(PHI)
      R0y = R0*dsin(theta)*dsin(PHI)
      R0z = R0*dcos(theta)

      DO 30, i = 1, points

        rx = R0x - X(i)   
        ry = R0y - y(i)
        rz = R0z - z(i)
        r  = DSQRT(rx*rx + ry*ry + rz*rz)
        nx = rx/r
        ny = ry/r
        nz = rz/r

        gamma_c = DSQRT(px(i)**2.0d0 + py(i)**2.0d0 + 
     &                   pz(i)**2.0d0 + c**2.d0)
      
        CALL calcEfield(X(i), y(i), z(i), t(i), e_las)
        ex_las = e_las(1)
        ey_las = e_las(2)
        ez_las = e_las(3)*X(i)
        bx_las = ey_las 
        by_las = e_las(4)
        bz_las = e_las(3)*y(i)
        vx = px(i)*c/gamma_c
        vy = py(i)*c/gamma_c
        vz = pz(i)*c/gamma_c

        pdotx = -(4.0d0/3.0d0)*PI*BckZ*density*X(i) -
     &           ex_las - ( vy*bz_las-vz*by_las )/c 
        pdoty = -(4.0d0/3.0d0)*PI*BckZ*density*y(i) -
     &           ey_las - ( vz*bx_las-vx*bz_las )/c
        pdotz = -(4.0d0/3.0d0)*PI*BckZ*density*z(i) -
     &          ez_las -  ( vx*by_las-vy*bz_las )/c

        vdotE = vx*ex_las + vy*ey_las + vz*ez_las
        ax = c*( pdotx + vdotE*vx/(c*c) )/gamma_c
        ay = c*( pdoty + vdotE*vy/(c*c) )/gamma_c
        az = c*( pdotz + vdotE*vz/(c*c) )/gamma_c

        nbeta = (nx*vx + ny*vy + nz*vz)/c
        denom = r*c*c*(1.0d0-nbeta)**3.0d0
        n_a = nx*ax + ny*ay + nz*az
        numX = -( (nx - vx/c)*n_a - (1.0d0-nbeta)*ax )
        numY = -( (ny - vy/c)*n_a - (1.0d0-nbeta)*ay )
        numZ = -( (nz - vz/c)*n_a - (1.0d0-nbeta)*az )
        ex_rad(i) = numX/denom
        ey_rad(i) = numY/denom
        ez_rad(i) = numZ/denom   
c        t_r0(i) = -DSQRT( (R0x - rx)**2.0d0 + (R0y - ry)**2.0d0 + 
c     &                  (R0z - rz)**2.0d0)/c + t(i)
         t_r0(i) = (r-R0)/c + t(i)
30    CONTINUE

      RETURN
      END  

C  Subroutine to generate filenames 
C  for storing individual trajectories 
C  this takes an integer 'no' and converts it to a filename 
C  character like Coh00001.dat and returns it to the main program

      SUBROUTINE int_to_str(fname,MC)
      character nmNo(5)
      character*12 fname(*)
        character temp(20)
        integer int0, lennmNo, m, MC, LL, i
        int0 = ichar('0')
     
      DO 88, no = 1, MC
        IF(no .eq. 0)then
          lennmNo=5
          DO LL=1,lennmNo
              nmNo(LL)='0'
          ENDDO
        ELSE
          m = no
          i = 0
          DO WHILE ( m .gt. 0)
             i = i + 1 
             temp(i) = char(int0 + mod(m,10))
             m = m/10
          ENDDO
          DO WHILE ( i .lt. 5 )
              i = i+1
              temp(i) = char(int0)
          ENDDO
          lennmNo = i
          DO i = 1, lennmNo
             nmNo(i) = temp(lennmNo+1-i)
          ENDDO
          fname(no) = 'Rad' // nmNo(1) // nmNo(2) // nmNo(3)
     &                // nmNo(4) // nmNo(5) // '.dat'
        ENDIF 
88      CONTINUE
        
      END
      SUBROUTINE int_to_str_spec(fname,MC)
      character nmNo(5)
      character*12 fname(*)
        character temp(20)
        integer int0, lennmNo, m, MC, LL, i
        int0 = ichar('0')

      DO 88, no = 1, MC
        IF(no .eq. 0)then
          lennmNo=5
          DO LL=1,lennmNo
              nmNo(LL)='0'
          ENDDO
        ELSE
          m = no
          i = 0
          DO WHILE ( m .gt. 0)
             i = i + 1
             temp(i) = char(int0 + mod(m,10))
             m = m/10
          ENDDO
          DO WHILE ( i .lt. 5 )
              i = i+1
              temp(i) = char(int0)
          ENDDO
          lennmNo = i
          DO i = 1, lennmNo
             nmNo(i) = temp(lennmNo+1-i)
          ENDDO
          fname(no) = 'Pol' // nmNo(1) // nmNo(2) // nmNo(3)
     &                // nmNo(4) // nmNo(5) // '.dat'
        ENDIF
   88 CONTINUE

      END
        
       subroutine rad_binning3(bin,bgap,Bmax,Trad,Eradx,
     &                        Erady,Eradz,Imax)
      
      integer          I, Imax, B, Bmax, tally
      
      double precision Trad(*), bgap, Eradx(*),Erady(*),
     &                 Eradz(*),slope(3), TempBin(3),
     &                 bin(Bmax,*)
      logical          Tend, start
CC The subroutine has
CC two counters, B for bin array, I for individual array
CC to be binned
            
CC Bin array has four columns, 1 for bin time, the other
CC three for binned function values, Ex Ey, and Ez
            
      I = 1
      B = 1 
      Tend = .false.
      start = .false.  
      do while (Tend .eqv. .false.)
        if (start .eqv. .false.) then
c          if (Trad(I) .lt. bin(B,1) ) then
          if (Trad(I) .lt. bin(B,1) ) then
            write(*,*) 'field time preceeds bin time'
            I = I + 1
            if (I .gt. Imax) Tend = .true.
          else
            start = .true.
            TempBin(1) = 0.0d0
            TempBin(2) = 0.0d0
            TempBin(3) = 0.0d0
            tally = 0
          endif
        else
c          if ( (Trad(I) .ge. bin(B,1)) .and.  
c     &       (Trad(I) .le. bin(B,1) ) then
           if ( (Trad(I) .ge. bin(B,1)) .and.
     &       (Trad(I) .le. bin(B,1)+bgap) ) then

            tally = tally + 1
            TempBin(1) = TempBin(1) + Eradx(I)
            TempBin(2) = TempBin(2) + Erady(I)
            TempBin(3) = TempBin(3) + Eradz(I)
            I = I + 1
            if (I .gt. Imax) then
              Tend = .true.
              bin(B,2) = bin(B,2)+ TempBin(1)/tally
              bin(B,3) = bin(B,3)+ TempBin(2)/tally
              bin(B,4) = bin(B,4)+ TempBin(3)/tally
            endif   
          else
            if (tally .eq. 0) then 
              if (I .eq. 1) then
                TempBin(1) = 0.0d0
                TempBin(2) = 0.0d0
                TempBin(3) = 0.0d0
              else
                slope(1) = (Eradx(I)-Eradx(I-1))/(Trad(I)-Trad(I-1))
                TempBin(1)=slope(1)*(bin(B,1)-Trad(I-1)) + Eradx(I-1)
                slope(2) = (Erady(I)-Erady(I-1))/(Trad(I)-Trad(I-1))
                TempBin(2)=slope(2)*(bin(B,1)-Trad(I-1)) + Erady(I-1)
                slope(3) = (Eradz(I)-Eradz(I-1))/(Trad(I)-Trad(I-1))
                TempBin(3)=slope(3)*(bin(B,1)-Trad(I-1)) + Eradz(I-1)
              endif
            else
              TempBin(1) = TempBin(1)/tally
              TempBin(2) = TempBin(2)/tally
              TempBin(3) = TempBin(3)/tally
              tally = 0
            endif
            bin(B,2) = bin(B,2)+ TempBin(1)   
            bin(B,3) = bin(B,3)+ TempBin(2)
            bin(B,4) = bin(B,4)+ TempBin(3)
            TempBin(1) = 0.0d0
            TempBin(2) = 0.0d0
            TempBin(3) = 0.0d0
            B = B + 1
            if (B .gt. Bmax) then
              Tend = .true.
c              if (I .le. Imax) write(*,*)
c     &          'field time past bin time'
            endif
          endif
        endif
      enddo
                
      return
      end

       subroutine fourier(Trad,Eradx,Erady,Eradz,Imax,NUMF,df,FREQ,
     &                   powspec,EW,start_freq)

       integer          I, Imax,NUMF,J,L,start_freq,ff

       DOUBLE PRECISION  Eradx(IMAX),Erady(IMAX),Eradz(IMAX),Trad(IMAX),
     &                  FREQ(NUMF),POWSPEC(NUMF),Fcos(NUMF), Fsin(NUMF)
        DOUBLE PRECISION  Icosx(NUMF),Icosy(NUMF),Icosz(NUMF),
     &             Isiny(NUMF),Isinz(NUMF),dt_scalar,Isinx(NUMF),
     &             slopex,slopey,slopez,t_rave(IMAX),temp,
     &             Dt(IMAX),MAX_FREQ,df,FREQ_eV,POWSPEC_eV,
     &             Ex(IMAX),Ey(IMAX),Ez(IMAX),EW(6,NUMF),
     &             alpha,beta

         DOUBLE PRECISION  PI, c
         COMMON            PI, c

       
         DO 300, J = 1, IMAX-1
          Dt(J) = Trad(j+1) - Trad(j)
          slopex = (Eradx(j+1)-Eradx(j))/Dt(J)
          slopey = (Erady(j+1)-Erady(j))/Dt(J)
          slopez = (Eradz(j+1)-Eradz(j))/Dt(J)
          t_rave(J) = (Trad(j+1) + Trad(j))/2.0d0
          Ex(J) = slopex*(t_rave(J) - Trad(j)) + Eradx(j)
          Ey(J) = slopey*(t_rave(J) - Trad(j)) + Erady(j)
          Ez(J) = slopez*(t_rave(J) - Trad(j)) + Eradz(j)
300     Continue
        DO 912, j = start_freq, NUMF
c       df = 13.633889510d0
c        Freq(j) = Freq(j-1) + df
        Icosx(J) = 0.0d0
        Icosy(J) = 0.0d0
        Icosz(J) = 0.0d0
        Isinx(J) = 0.0d0
        Isiny(J) = 0.0d0
        Isinz(J) = 0.0d0
912   CONTINUE
c       alpha = 2.0d0*dsin(df/2.0d0)*dsin(df/2.0d0)
c       beta = dsin(df)
       DO 177, J = 1, IMAX-1
        alpha = 2.0d0*dsin(df*T_rave(J)/2.0d0)*dsin(df*T_rave(J)/2.0d0)          
        beta = dsin(df*T_rave(J))

         DO 188, I = start_freq, NUMF

           IF (I.LE.(start_freq+1)) THEN
                  Fcos(I) = dcos(FREQ(I)*T_rave(J))
                  Fsin(I) = dsin(FREQ(I)*T_rave(J))
           END IF
           IF (I.gt.(start_freq+1)) THEN
c              ff = start_freq + 1
              Fcos(I) = Fcos(I-1) - (alpha*Fcos(I-1) + beta*Fsin(I-1))
              Fsin(I) = Fsin(I-1) - (alpha*Fsin(I-1) - beta*Fcos(I-1))
           END If
           Icosx(I) = Icosx(I) + Ex(J)*Fcos(I)*Dt(J)

           Icosy(I) = Icosy(I) + Ey(J)*Fcos(I)*Dt(J)

           Icosz(I) = Icosz(I) + Ez(J)*Fcos(I)*Dt(J)

           Isinx(I) = Isinx(I) + Ex(J)*Fsin(I)*Dt(J)

           Isiny(I) = Isiny(I) + Ey(J)*Fsin(I)*Dt(J)
           
           Isinz(I) = Isinz(I) + Ez(J)*Fsin(I)*Dt(J)

188      CONTINUE

177   CONTINUE

            DO 211, L = start_freq, NUMF
          powspec(L) = powspec(L) +  Icosx(L)**2.0d0+Icosy(L)**2.0d0+
     &                Icosz(L)**2.0d0+Isinx(L)**2.0d0+Isiny(L)**2.0d0+
     &                Isinz(L)**2.0d0
           EW(1,L) = EW(1,L) + Icosx(L)
           EW(2,L) = EW(2,L) + Icosy(L)
           EW(3,L) = EW(3,L) + Icosz(L)
           EW(4,L) = EW(4,L) + Isinx(L)
           EW(5,L) = EW(5,L) + Isiny(L)
           EW(6,L) = EW(6,L) + Isinz(L)
        
c                  FREQ_ev = FREQ(L)*27.21
cc                  powspec(L) = powspec(L)*1.0d+18*c/(4*PI)
c                  powspec_ev = powspec(L)
211       CONTINUE
          RETURN
          END






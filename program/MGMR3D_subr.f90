!------------------------------
    subroutine LateralInt(idi)
    use BigArrays, only : t_tb, Ex_tb,AxD_tb,Ey_tb,AyD_tb,Ar_tb, Erh_tb,t_ch
    use BigArrays, only : Ex_spld,AxD_spld,Ey_spld,AyD_spld,Ar_spld, Erh_spld
    use BigArrays, only : t_to, Ex_to,Ey_to,Er_to
    use BigArrays, only : ObsDist_dim, tTrace_dim_b, tTrace_dim_o, ObsDist_Step
    use BigArrays, only : Line2Core, CoreDist_Dim
!    use BigArrays, only : Ix,Iy,IQ,Ix_int,Iy_Int
    use constants, only : ci,pi,dp
    implicit none
    integer, intent(in) :: idi
    integer :: Nrs,irs,ith,nth,id,CD_i
    real(dp)  :: dth,theta,drs,rs,y_s,x_s,x_o,ro,Wei,weid,rs_max,gth,did
    real(dp)  :: Int,d0,int2,intx,wrs,wro,dwr,CD_d,tshft ! ,CoreDist,r_s
    real(dp) :: DEx(tTrace_dim_o),DAxD(tTrace_dim_o),DEy(tTrace_dim_o),DAyD(tTrace_dim_o), &
        DAr(tTrace_dim_o),DErh(tTrace_dim_o)
    external W_tc
    real(dp) W_tc
!
! lateral integral
! create grid at a fixed d0 as distance between core and observer
    d0= idi*ObsDist_Step
    int=0.  ; intx=0. ;  int2=0.
    nrs=CoreDist_Dim*2+2
    rs_max= Line2Core(nrs)
 !   Nth=15  ; dth=pi/Nth
    Nth=7  ; dth=pi/Nth
    gth=1.
    Do ith=0,Nth    ! for theta=0 point lies inbetween sh-core and observer
      theta=ith*dth
      do irs = 1, nrs   ! distance from shower core; CoreDist_A(0,CoreDist_Dim)
        rs=Line2Core(irs)
        drs=(Line2Core(irs+1)-Line2Core(irs-1))/2.
        !rs=(irs-0.5)*drs
        x_s=rs*cos(theta)
        if(x_s.gt. d0/2.) exit
        wrs=drs
        !wrs=0.     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        If((rs+drs)*cos(theta) .gt. d0/2.) then
            wrs=d0/(2.*cos(theta))-(rs-drs/2.)  ! remaining distance to the center line betwee observer and core
        endif
        wei=wrs*dth*W_tc(rs,dwr)
        weid=wrs*dth*dwr
        y_s=rs*sin(theta)
        ro=sqrt(y_s*y_s+(x_s-d0)*(x_s-d0))
        !
        id=IDINT(ro/ObsDist_Step + 0.5)  ;  did=ro/ObsDist_Step-(id-0.5)
        if(id.ge.ObsDist_dim) exit
        if(id.lt. 1) then
          id=1  ; did=0.
        endif
        call Find_CD_Interpol(rs,CD_i,CD_d)
        !
        tshft=((1.-did)*t_ch(id)+did*t_ch(id+1))   ! need to apply time-shift
        call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            Ex_tb(1,id,CD_i),   Ex_spld(1,id,CD_i),    Ex_tb(1,id,CD_i+1),   Ex_spld(1,id,CD_i+1), &
            Ex_tb(1,id+1,CD_i), Ex_spld(1,id+1,CD_i),  Ex_tb(1,id+1,CD_i+1), Ex_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DEx )
        call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            AxD_tb(1,id,CD_i),   AxD_spld(1,id,CD_i),    AxD_tb(1,id,CD_i+1),   AxD_spld(1,id,CD_i+1), &
            AxD_tb(1,id+1,CD_i), AxD_spld(1,id+1,CD_i),  AxD_tb(1,id+1,CD_i+1), AxD_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DAxD )
        Ex_to(:,idi)=Ex_to(:,idi) + gth*wei *DEx(:) - gth*DAxD(:)*weid
        call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            Ey_tb(1,id,CD_i),   Ey_spld(1,id,CD_i),    Ey_tb(1,id,CD_i+1),   Ey_spld(1,id,CD_i+1), &
            Ey_tb(1,id+1,CD_i), Ey_spld(1,id+1,CD_i),  Ey_tb(1,id+1,CD_i+1), Ey_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DEy )
        call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            AyD_tb(1,id,CD_i),   AyD_spld(1,id,CD_i),    AyD_tb(1,id,CD_i+1),   AyD_spld(1,id,CD_i+1), &
            AyD_tb(1,id+1,CD_i), AyD_spld(1,id+1,CD_i),  AyD_tb(1,id+1,CD_i+1), AyD_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DAyD )
        Ey_to(:,idi)=Ey_to(:,idi) + gth*wei *DEy(:) - gth*DAyD(:)*weid
        call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            Ar_tb(1,id,CD_i),   Ar_spld(1,id,CD_i),    Ar_tb(1,id,CD_i+1),   Ar_spld(1,id,CD_i+1), &
            Ar_tb(1,id+1,CD_i), Ar_spld(1,id+1,CD_i),  Ar_tb(1,id+1,CD_i+1), Ar_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DAr )
        call  Spline_Average2_Shift(tTrace_dim_b, t_tb(1), &
            Erh_tb(1,id,CD_i),   Erh_spld(1,id,CD_i),    Erh_tb(1,id,CD_i+1),   Erh_spld(1,id,CD_i+1), &
            Erh_tb(1,id+1,CD_i), Erh_spld(1,id+1,CD_i),  Erh_tb(1,id+1,CD_i+1), Erh_spld(1,id+1,CD_i+1), &
            did, CD_d, tshft, tTrace_dim_o, t_to(1), DErh )
        Er_to(:,idi)=Er_to(:,idi) + gth*(DAr(:)*weid + DErh(:)*wei)*cos(theta)
      enddo     ! irs
!
      Do id = 1, ObsDist_dim   ! distance from observer in direction of shower core for theta=0.
        ro=(id-0.5)*ObsDist_Step
        x_o=ro*cos(theta)
        if(x_o.gt. d0/2.) exit
        wro=ObsDist_Step
        !wro=0. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        If((ro+ObsDist_Step)*cos(theta) .gt. d0/2.) then
            wro=d0/(2.*cos(theta))-(ro-ObsDist_Step/2.)  ! remaining distance to the center line betwee observer and core
        endif
        x_s=d0-x_o  ;  y_s=ro*sin(theta)  ; rs=sqrt(x_s*x_s+y_s*y_s)
        if(rs.gt. rs_max) cycle
        wei=wro*dth*W_tc(rs,dwr)*ro/rs  ! replaced wrs by wro
        weid=wro*dth*dwr*ro/rs          ! replaced wrs by wro
!
        call Find_CD_Interpol(rs,CD_i,CD_d)
        tshft=t_ch(id)
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            Ex_tb(1,id,CD_i), Ex_spld(1,id,CD_i), Ex_tb(1,id,CD_i+1), Ex_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DEx )
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            AxD_tb(1,id,CD_i),   AxD_spld(1,id,CD_i),    AxD_tb(1,id,CD_i+1),   AxD_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DAxD )
        Ex_to(:,idi)=Ex_to(:,idi) + gth*wei *DEx(:) - gth*DAxD(:)*weid
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            Ey_tb(1,id,CD_i),   Ey_spld(1,id,CD_i),    Ey_tb(1,id,CD_i+1),   Ey_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DEy )
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            AyD_tb(1,id,CD_i),   AyD_spld(1,id,CD_i),    AyD_tb(1,id,CD_i+1),   AyD_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DAyD )
        Ey_to(:,idi)=Ey_to(:,idi) + gth*wei *DEy(:) - gth*DAyD(:)*weid
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            Ar_tb(1,id,CD_i),   Ar_spld(1,id,CD_i),    Ar_tb(1,id,CD_i+1),   Ar_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DAr )
        call  Spline_Average_Shift(tTrace_dim_b, t_tb(1), &
            Erh_tb(1,id,CD_i),   Erh_spld(1,id,CD_i),    Erh_tb(1,id,CD_i+1),   Erh_spld(1,id,CD_i+1), &
            CD_d, tshft, tTrace_dim_o, t_to(1), DErh )
        Er_to(:,idi)=Er_to(:,idi) + gth*(DAr(:)*weid + DErh(:)*wei)*(x_s/rs)
      enddo     ! id
!      write(2,*) irs,ith,id,did,wei,Ex_nui((inui+inum)/2,idi)
      gth=2. ; if(ith.eq.(Nth-1)) gth=1.
    enddo       ! ith
!    write(2,*) 'Int W(r)=',Ex_nui(inui,idi),Ex_nui(inum,idi)
!    write(2,*) 'd0=',d0,int,int2,int+int2
    ! stop
    return
    end subroutine LateralInt
!-------------------
subroutine Get_ObsPlsDelay(idi)
    use BigArrays, only : t_to, Ex_to,Ey_to,Er_to
    use BigArrays, only : tTrace_dim_o, tTrace_step, ObsPlsTime
    use constants, only : dp
    use RFootPars, only : test
    implicit none
    integer, intent(in) :: idi
    integer :: i,j,k,nth,id,inu,CD_i,m(1)
    real(dp) :: Overlap(100), Norma, Normb
    !
    k=ObsPlsTime(idi-1)/tTrace_step
    If((k+300).ge.tTrace_dim_o) then
        write(2,*) 'ERROR, time-trace too short to contain full peak at distance#',idi
        stop 'ERROR, time-trace too short'
    endif
    Do i=1,100
        Overlap(i)=0.0
        Norma=0.
        Normb=0.
        Do j=1,200
            Overlap(i)=Overlap(i) + Ex_to(k+j,idi-1)*Ex_to(k+j+i,idi) + Ey_to(k+j,idi-1)*Ey_to(k+j+i,idi)
            Norma = Norma + Ex_to(k+j,idi-1)*Ex_to(k+j,idi-1) + Ey_to(k+j,idi-1)*Ey_to(k+j,idi-1)
            Normb = Normb + Ex_to(k+j+i,idi)*Ex_to(k+j+i,idi) + Ey_to(k+j+i,idi)*Ey_to(k+j+i,idi)
        enddo
        Overlap(i)=Overlap(i)/sqrt(Norma*Normb)
    enddo
    m = maxloc(Overlap)
    !if(test) then
    !    write(2,*) idi,m
    !    write(2,"(20f6.3)") Overlap
    !endif
    ObsPlsTime(idi)=ObsPlsTime(idi-1) + m(1)*tTrace_step
end subroutine Get_ObsPlsDelay
!------------------------------
Real*8 Function W_tc(r,dwr)
!  W_tc == W = W_NKG *r
!  dwr == x d(W/x)/dx  !   Wprim= r {d\, w(r)/r \over dr}
!  -----   case  NKG=.true. : NKG for s=2
!    W_tc=Norm_tc*x**(s-1)*(x+1.)**(s-4.5)  ! r*NKG for s
!   for W/x=Norm* x**(s-2)*(x+1.)**(s-4.5)
!        gives dwr= x d(W/x)/dx= x Norm[(s-2)*x**(s-3)*(x+1.)**(s-4.5) + (s-4.5)*x**(s-2)*(x+1.)**(s-5.5)]=
!        thus dwr= x Norm[(s-2)*(x+1.) + (s-4.5)*x]x**(s-3)*(x+1.)**(s-5.5)= [(s-2)/x + (2s-6.5)] w_tc /(x+1.)
!  ----   case  NKG=.false. & tu1=.true. : tu form with t=6, u=1
!    W_tc=Norm_tc*W=x^u/(x+1)^t  ! r* non NKG form t, u
!        Max of W at x=u/(t-u)=a
!        thus set x=a* r/Moli   to have max independent of t and u, normalized for the value for (t=6, u=1) at r=moli
!   for W/x=Norm x^(u-1)/(x+1.)^t
!        gives dwr= x d(W/x)/dx= x Norm[(u-1)/x -t/(x+1)] x^(u-1)/(x+1.)^t= Norm [x(u-1-t)+u-1] w_tc /(x(x+1.))
   real*8 :: r,x,dwr,a
   common /W_TrCurr/ Norm_tc,Moli_tc,s_tc,dalpha,Nrm_alpha(0:41)
   real*8 :: Norm_tc,Moli_tc,s_tc,dalpha,Nrm_alpha ! ,alpha
   !Integer, parameter :: NKG=0     ! NKG, s=2 =(t=2.5,u=1)
   !Integer, parameter :: NKG=1      ! =(t=6,u=1)
   !Integer, parameter :: NKG=2      ! =(t=6,u=1.5)
   Integer, parameter :: NKG=3      ! =(t=6,u=2)
    !---------------------------------------------------
    select case (NKG)
    case (1)
      !    W_tc=Norm_tc*x/((x+1.)**4)  ! r* non NKG form t=4, u=1
      !    dwr= -4 *W_tc/(x+1.) !  for t=4, u=1
      a=1.5/(5.*Moli_tc) !  for t=6, u=1  ! factor 1.5 in a to keep the same value for Moli for the NKG (default) case
      x=r*a     !
      W_tc=Norm_tc*x/((x+1.)**6)  ! r* non NKG form t=6, u=1
      dwr= -6 *a*W_tc/(x+1.) !  for t=6, u=1
    case (2)
      a=1.5/(3.*Moli_tc)    ! for  t=6, u=1.5 ; *1.5 to normalize for NKG
      x=r*a
      W_tc=Norm_tc*x*sqrt(x)/((x+1.)**6)  ! r* non NKG form t=6, u=1.5  ! seems to need to increase J0Q
      dwr= a*(0.5/x - 5.5) *W_tc/(x+1.) !  for t=6, u=1.5
    case (3)
      a=1.5/(2.*Moli_tc) ! for  t=6, u=2 ; *1.5 to normalize for NKG
      x=r*a
      W_tc=Norm_tc*x*x/((x+1.)**6)  ! r* non NKG form t=6, u=2  ! seems to have a too small CE force near core
      dwr= a*(1./x-5) *W_tc/(x+1.) !  for t=6, u=2
    case default
      !    W_tc=Norm_tc*sqrt(x)/((x+1.)**3)  ! r*NKG for s=1.5
      !    dwr= (-0.5/x-3.5) *W_tc/(x+1.) !  for s=1.5
      a=1./Moli_tc !  for t=2.5, u=1 ; NKG, s=2
      x=r*a     !
      W_tc=Norm_tc*x/(x+1.)**2.5  ! r*NKG for s=2
      dwr= -2.5*a *W_tc/(x+1.) !  for s=2
      !    W_tc=Norm_tc*x*sqrt(x)/(x+1.)**2  ! r*NKG for s=2.5
      !    dwr= (0.5/x-1.5) *W_tc/(x+1.) !  for s=2.5
    End Select
   !
End function W_tc
!------------------------------
    subroutine Find_CD_Interpol(CD,CD_i,CD_d)
        use BigArrays, only : CoreDist_Dim, CoreDist_A
        use constants, only : dp
        implicit none
    real(dp), intent(in) :: CD
    integer , intent(out) :: CD_i
    real(dp), intent(out) :: CD_d
    integer :: i
    CD_i=CoreDist_Dim-1
    Do i=1,CoreDist_Dim
      if(CD.lt.CoreDist_A(i)) then
        CD_i=i-1
        exit
      endif
    enddo
    CD_d=(CD-CoreDist_A(CD_i))/(CoreDist_A(CD_i+1)-CoreDist_A(CD_i))
    if(CD_d.lt.0.) CD_d=0.
    ! write(2,*) 'cd',cd,cd_i,cd_d,CoreDist_A(CD_i),CoreDist_A(CD_i+1)
    end subroutine Find_CD_Interpol
!------------------------------
    subroutine fieh(h,CoreDist_i,alpha,fh,dfhl,dfha)
!   Calculate the value of the pancake function, fh, and
!       the derivaties v.s. lambda (=integration variable in line-emission i.e. NOT 'lam'), dfhl, and
!       the derivaties v.s. alpha, dfha.
!   h=distance behind front
    use BigArrays, only : PancakeThi
    use RFootPars, only : SelectFh
    implicit none
    integer, intent(in) :: CoreDist_i
    real*8, intent(in) :: h,alpha
    real*8, intent(out) :: fh,dfhl,dfha
    real*8 :: x,lam,di,nrm,a,b,dfdx,d  ! r,
    real*8, parameter :: c=1.d0, xScale=10.d0 , xScale3=30.d0
    integer :: i
    common /W_TrCurr/ Norm_tc,Moli_tc,s_tc,dalpha,Nrm_alpha(0:41)
      real*8 :: Norm_tc,Moli_tc,s_tc,dalpha,Nrm_alpha
    lam=alpha*PancakeThi(CoreDist_i)
    x=h/lam
    if(x.gt.20.) then
      fh=0.0d0 ; dfhl=0.0d0 ; dfha=0.0d0
      return
    endif
    i=IDINT((alpha-1.)/dalpha)
    di=alpha-1.-i*dalpha
    nrm=((1.-di)*Nrm_alpha(i) + di*Nrm_alpha(i+1))/lam
    !
    select case (SelectFh)
    !---------------------------------------------------
    case (1)
    !------  fh=x^alpha e^(-2x)  ----------------------
        fh =Nrm*x**(alpha-1) *exp(-2.*x) ! =fh/x
        dfhl=fh*(alpha-2*x)/lam   !  =df/dh at fixed alpha
        fh=fh*x
        dfha=(Nrm_alpha(i+1)-Nrm_alpha(i))/dalpha  ! =df/da at fixed h
        dfha=(dfha/Nrm + log(x) - (1-2*x/alpha))*fh
        ! PT=lam/alpha
        ! dfhdPT=fh/PT  - dfdx * x/PT =  alpha/lam *fh (1. - (alpha-2 x))
        ! from line integration:       dfhdr = - (1. + alpha - 2. * h/(alpha*lam) ) * fh*dlambda_hdr/(lam=pancakethicjness)
    !---------------------------------------------------
    case (3)
    !------  fh=x^2 /(exp(sqrt(x))+c*alpha)  -------------
        x=x*xScale3       ! x=10 h / lam
        a=sqrt(x)
        b=exp(a)
        d=c*alpha
        fh =Nrm*x/(b+d)
        dfdx=fh*(2.d0 - 0.5d0*a*b/(b+d))
        dfhl=dfdx*xScale3/lam     ! (dfh/dx) * (dx/dh=10/lam)
        fh=fh*x
        dfha=(Nrm_alpha(i+1)-Nrm_alpha(i))/dalpha  ! =df/da at fixed h
        dfha=(dfha/Nrm - c/(b+d) )*fh - dfdx*x/alpha ! (dx/da=-x/alpha)
    !---------------------------------------------------
    case (4)
    !------  fh=x /(exp(sqrt(x))+1)  -------------
        x=x*xScale
        a=sqrt(x)
        b=exp(a)
        fh =Nrm /(b+c)
        dfdx=fh*(1.d0- 0.5d0*a*b/(b+c))
        dfhl=dfdx*xScale/lam
        fh=fh*x
        dfha=(Nrm_alpha(i+1)-Nrm_alpha(i))/dalpha  ! =df/da at fixed h
        dfha=-fh*xScale/alpha - dfdx*x*xScale/alpha
    !---------------------------------------------------
    case default
    !------  fh=x^alpha /(exp(sqrt(x))+1)  -------------
        x=x*xScale
        a=sqrt(x)
        b=exp(a)
        fh =Nrm*x**(alpha-1) /(b+c)
        dfdx=fh*(alpha- 0.5d0*a*b/(b+c))
        dfhl=dfdx*xScale/lam
        fh=fh*x
        dfha=(Nrm_alpha(i+1)-Nrm_alpha(i))/dalpha  ! =df/da at fixed h
        dfha=(dfha/Nrm + log(x))*fh - dfdx*x*xScale/alpha
    !---------------------------------------------------
    !
    end select
    !
    return
    end subroutine fieh
!------------------------------
    subroutine getZeta(lamb,zetai,zeta)
    ! lamb is a measure for the distance behind the shower front (difference with Cherenkov height)
    ! zeta is the emission height in the atmosphere (including the distance behind the shower front)
    use BigArrays, only : lambda, AtmHei_dim, AtmHei_step
    implicit none
    real*8 :: zetaa,zetab,lamb,da,db,zetai(2),zeta(2),zet
    integer :: i,k
    common / time / zeta_c,hcto_c
      real*8 :: zeta_c,hcto_c
!
!       lower branch
    k=zetai(1)/AtmHei_step
!    write(2,*) 'zetai',zetai,k
    if(k.ge.AtmHei_dim) then
      k=AtmHei_dim
      if(lambda(k).lt.lamb) goto 1
      zeta(1)=AtmHei_dim*AtmHei_step+2
      goto 3
    endif
1   continue
    if(k.lt.0) then
        zeta(1)=-AtmHei_step
        goto 3
    else if(lambda(k).lt.lamb) then
      k=k-1
      goto 1
    else
      zetab=k*AtmHei_step
      db=lambda(k)-lamb
      zetaa=(k+1)*AtmHei_step ! search lower in height
      da=lambda(k+1)-lamb
      if(zetaa .gt. zeta_c) then
        zetaa=zeta_c
        da=-lamb
      endif
!        lambda(i)=hcto + sqrt( (sqrt(d*d+z*z)-(1.+xi(i))*z)**2+((1.+xi(i))**2-1.)*d*d)
      zet=(zetaa*db -zetab*da)/(db-da)
!      write(2,*) k,zet,zetaa,da,zetab,db
      call RefineLam(zetaa,da,zetab,db,zet,lamb)
      zeta(1)=zet
!      write(2,*) 'zeta1',zeta(1),lamb
    endif
!       upper branch
3    k=zetai(2)/AtmHei_step+1
2   continue
    if(k.ge.AtmHei_dim) then
        zeta(2)=(AtmHei_dim+1.)*AtmHei_step
        goto 4
    else if(lambda(k).lt.lamb) then
      k=k+1
      goto 2
    else
      zetab=k*AtmHei_step
      db=lambda(k)-lamb
      zetaa=(k-1)*AtmHei_step ! search larger in height
      da=lambda(k-1)-lamb
      if(zetaa .lt. zeta_c) then
        zetaa=zeta_c
        da=-lamb
      endif
      zet=(zetaa*db -zetab*da)/(db-da)
!      write(2,*) k,zet,zetaa,da,zetab,db
      call RefineLam(zetaa,da,zetab,db,zet,lamb)
      zeta(2)=zet
!      write(2,*) 'zeta2',zeta(2),lamb
    endif
!
4   continue
    return
  end subroutine getZeta
!------------------------------
    subroutine RefineLam(zeta,da,zetb,db,zet,lam)
    ! Solves for \zeta:  n\sqrt{(-\beta t +h)^2 + (1-\beta^2 n^2)d^2} =n(R-n\zeta)
    ! given (-t+h)= hcto_c - lam where lam is input parameter
    ! input: value of zet is used as start value for the search between values zeta and zetb
    ! input: da and db are the differences between lama, lamb values for zeta and zetb with lam
    use BigArrays, only : xi,dxi, AtmHei_dim, AtmHei_step, ObsDist
    implicit none
    real*8 :: zeta,zetb,zet,da,db,dd,lam,NN,R2,di  ! ,dr
    real*8, parameter :: Converg=1.e-5
    integer, parameter :: MaxIter=5
    integer :: i,k
    common / time / zeta_c,hcto_c
      real*8 :: zeta_c,hcto_c
    Real*8 zetx, hxto
    i=zet/AtmHei_step
!
    Do k=1,MaxIter  ! find value of `zet' that zeros `DD' for a given value of `lam', all have units of length
      R2=ObsDist*ObsDist+zet*zet
      if(i.lt.AtmHei_dim) then
        di=zet/AtmHei_step-i
        NN=1.+ di*xi(i+1)+(1.-di)*xi(i)
      else
        zet=AtmHei_dim*AtmHei_step+1.
        return
      endif
      DD=hcto_c + sqrt( (sqrt(ObsDist*ObsDist+zet*zet)-NN*zet)**2+(NN*NN-1.)*ObsDist*ObsDist) - lam
!      write(2,*) 'niter',k,zet,dd,i
      if(Da*DD.le. 0.) then
        zetb=zet
        Db=DD
      else if(Db*DD.le. 0.) then
        zeta=zet
        Da=DD
      else
        write(2,*) 'RefineLam:',zeta,Da,zetb,Db,zet, dd
        exit
      endif
      zet=(zeta*db -zetb*da)/(db-da)
      if(abs(zeta-zetb).lt.converg) exit
    enddo
     !write(2,*) 'RefineLam:',zeta,Da,zetb,Db,zet
     !hxto=hcto_c - lam
     !zetx= (hxto + nn*sqrt(hxto*hxto + (1-nn*nn)*ObsDist*ObsDist))/(1-nn*nn)
     !DD=hxto + sqrt( (sqrt(ObsDist*ObsDist+zetx*zetx)-NN*zetx)**2+(NN*NN-1.)*ObsDist*ObsDist)
     !write(2,*) hxto, zetx, DD
    end subroutine RefineLam
!------------------------------
!------------------------------
    subroutine getZeta1(lamb,zetai,zeta)
    use BigArrays, only : lambda, AtmHei_dim, AtmHei_step
    implicit none
    real*8 :: zetaa,zetab,lamb,da,db,zetai,zeta! ,zet
    integer :: k    ! i,
    common / time / zeta_c,hcto_c
      real*8 :: zeta_c,hcto_c
!
  k=zetai/AtmHei_step
  if(zetai .lt. zeta_c) then
!       lower branch
1   continue
    if(k.gt.0) then
      if(lambda(k).lt.lamb) then
        k=k-1
        goto 1
      endif
    endif
    if(k.ge.AtmHei_dim) then
        zeta=(AtmHei_dim+0.5)*AtmHei_step
        return
    endif
      zetab=k*AtmHei_step
      db=lambda(k)-lamb
      zetaa=(k+1)*AtmHei_step
      da=lambda(k+1)-lamb
      if(zetaa .gt. zeta_c) then
        zetaa=zeta_c
        da=-lamb
      endif
!        lambda(i)=hcto + sqrt( (sqrt(d*d+z*z)-(1.+xi(i))*z)**2+((1.+xi(i))**2-1.)*d*d)
!      zet=(zetaa*db -zetab*da)/(db-da)
!      write(2,*) 'getZeta1-low',k,zet,zetaa,da,zetab,db,zetai,lamb,lambda(k)
      call RefineLam1(zetaa,da,zetab,db,k,zeta,lamb)
!      write(2,*) 'zeta1',zeta(1),lamb
!
  else    !       upper branch
2   continue
!    write(2,*) zetai,lamb
    if(k.ge.AtmHei_dim) then
        zeta=(AtmHei_dim+1.)*AtmHei_step
        return
    else if(lambda(k).lt.lamb) then
      k=k+1
      goto 2
    else
      zetab=k*AtmHei_step
      db=lambda(k)-lamb
      zetaa=(k-1)*AtmHei_step
      da=lambda(k-1)-lamb
      if(zetaa .lt. zeta_c) then
        zetaa=zeta_c
        da=-lamb
      endif
!      zet=(zetaa*db -zetab*da)/(db-da)
!      write(2,*) 'getZeta1-up',zetai,zeta_c,k,zet,zetaa,da,zetab,db
      call RefineLam1(zetaa,da,zetab,db,k-1,zeta,lamb)
!      write(2,*) 'zeta2',zeta(2),lamb
    endif
!
  endif
    return
  end subroutine getZeta1
!------------------------------
    subroutine RefineLam1(zeta,da,zetb,db,i,zet,lam)
    use BigArrays, only : xi,dxi, AtmHei_dim, AtmHei_step, ObsDist
    implicit none
    real*8 :: zeta,zetb,zet,da,db,dd,lam,NN,R2,di,dr
    real*8, parameter :: Converg=1.e-5
    integer, parameter :: MaxIter=5
    integer :: i,k
    common / time / zeta_c,hcto_c
      real*8 :: zeta_c,hcto_c
!
    if(i.ge.AtmHei_dim) then
        zet=AtmHei_dim*AtmHei_step+1.
        return
    endif
    if(Da*Db.gt. 0.) then
        write(2,*) 'RefineLam1:',zeta,Da,zetb,Db,i
        zet=AtmHei_dim*AtmHei_step+2.
        return
    endif
    zet=(zeta*db -zetb*da)/(db-da)
    Do k=1,MaxIter
      R2=ObsDist*ObsDist+zet*zet
      di=zet/AtmHei_step-i
      NN=1.+ di*xi(i+1)+(1.-di)*xi(i)
      DD=hcto_c + sqrt( (sqrt(ObsDist*ObsDist+zet*zet)-NN*zet)**2+(NN*NN-1.)*ObsDist*ObsDist)-lam
!      write(2,*) 'niter',k,zet,dd,i
      if(Da*DD.le. 0.) then
        zetb=zet
        Db=DD
      else
        zeta=zet
        Da=DD
      endif
      zet=(zeta*db -zetb*da)/(db-da)
      if(abs(zeta-zetb).lt.converg) exit
    enddo
    end subroutine RefineLam1
!------------------------------
!------------------------------
  subroutine FindZetaC(k,zeta,hcto)
!
!     In general: zeta= negative retarded time
!     Zeta_C= the zeta (=height from which signal is emitted) for which the retarded distance CalD equals zero
!     Cald is equal to the the denominator of eq 34, and may be different from eq 35 by a factor n.
!     The programmed form for CalD rather follows eq 36 where aslo the change in index of refraction with height is accounted for
!  This subroutine is called at line 124 of      Subroutine MGMR3D
!     where k is index in atmospheric height for which
!            T_obs(i)=-z+sqrt(ObsDist*ObsDist+z*z)*(1.+xi(i))  ! z=i*AtmHei_step  (i=atm height)
!     is minimal (i=k), i.e. close to the cherenkov distance
!
!
    use BigArrays, only : AtmHei_dim, AtmHei_step
    implicit none
    real*8 :: zeta1,zeta2,calDa,calDb,calD,zeta,hcto ! ,N_s
    integer :: i,k
    real*8, parameter :: Converg=1.e-10
!
    if(k.ge.AtmHei_dim) then
! at large height $xi*zeta=constant$ and $xi'=-xi*zeta/zeta^2$
      zeta1=(AtmHei_dim-1)*AtmHei_step
      zeta2=10*AtmHei_dim*AtmHei_step+5.
    else
      zeta1=(k-1)*AtmHei_step
      zeta2=(k+1)*AtmHei_step
    endif
    call GetCalD(zeta1,calDa,hcto)
    call GetCalD(zeta2,calDb,hcto)
    if(k.ge.AtmHei_dim) then
      if(calDa*calDb .gt. 0.) then
        zeta=zeta2
        return
      endif
    endif
    Do i=1,10
      zeta=(zeta1*calDb-zeta2*calDa)/(calDb-calDa)
      call GetCalD(zeta,calD,hcto)
     if(calDa*calD.le. 0.) then
        zeta2=zeta
        calDb=calD
      else if(calDb*calD.le. 0.) then
        zeta1=zeta
        calDa=calD
      else
        write(2,*) 'FindZetaC:',zeta1,calDa,zeta2,calDb,zeta,hcto,calD
        exit
      endif
    enddo
    if(k.ge.AtmHei_dim) write(2,*) 'z',zeta,cald
    return
  end subroutine FindZetaC
!------------------------------
    subroutine GetCalD(zeta,calD,hcto)
    use BigArrays, only : xi,dxi, PenDepth, IndRefOrho, AtmHei_step, AtmHei_dim, TopAtmExpon, ObsDist
    use constants, only : dp
    implicit none
    real(dp) :: zeta,calD,hcto,NN,R2,di,dr,h0  ! ,hcto2
    integer :: i
    i=zeta/AtmHei_step
    di=zeta/AtmHei_step-i
    if(i.ge.AtmHei_dim) then
      if(i.ge. AtmHei_dim*2) then
        hcto=-IndRefOrho*PenDepth(0)  ! (IndRefOrho == Index of refraction over density)(PenDepth == penetration depth in [g/cm^2])
        calD=ObsDist*ObsDist/(2.*zeta)
!  N.B. at d=140 m there is an inflection point (at 12+ km height) in t-ret v.s. t-obs which causes a 3rd branch.
        return
      endif
!      NN=1.+rh0 * (1./zeta)* h0*(1.- exp(-zeta/h0)) ! mean index of refraction along path
      h0= TopAtmExpon
      NN= 1. + IndRefOrho*(PenDepth(0)-PenDepth(AtmHei_dim)*exp(-(zeta-AtmHei_dim*AtmHei_step)/h0))/zeta
      dr=IndRefOrho * (-PenDepth(0)+PenDepth(AtmHei_dim)*(zeta/h0+1)*exp(-(zeta-AtmHei_dim*AtmHei_step)/h0))/(zeta*zeta)   !  ((h0+zeta)*exp(-zeta/h0) - h0)/(zeta*zeta)
    else
      NN=1.+ di*xi(i+1)+(1.-di)*xi(i)   ! Mean index of refraction from height zeta to ground
      dr=   di*dxi(i+1)+(1.-di)*dxi(i)  ! change in refractivity
    endif
    R2=ObsDist*ObsDist+zeta*zeta
    calD=sqrt(R2) - NN*zeta- dr*R2 ! corresponds to eq 35 of the notes where calD=D/n and a term is added to account for varying inderx of refraction
    hcto=-sqrt(dr*R2*dr*R2 + (NN*NN-1.)*ObsDist*ObsDist)  !=-\sqrt{(R-N\zeta)^2+(N^2-1) d^2}
!    hcto2=-sqrt( (sqrt(R2)-NN*zeta)**2+(NN**2-1.)*d*d)
!    write(2,*) 'hcto-hcto2=',hcto-hcto2
    !write(2,*) 'GetCalD:',hcto, 'h0to - hcto =', -ObsDist*sqrt(NN*NN-1.) - hcto
    end subroutine GetCalD
!------------------------------
!------------------------------
    subroutine CurrField(Ex,AxD,Ey,AyD,Ar,Erh, CD_i,T_o)
!   Calculate the field (i.e. integral over pancake thickness and z) for a single ray
    use BigArrays, only : L_I,JL_I,dL_I,zeta_I, PancakeThi, LamInt_dim, LamGrid_Nrt, LamGrid_max
    use BigArrays, only : xi,dxi, AtmHei_dim, AtmHei_step, CoreDist_A, CoreDist_Dim
    use BigArrays, only : Ix,Iy,IQ,alpha_tr,Ix_int,Iy_Int, ObsDist, ObsDist_dim, ObsDist_Step
    use constants, only : dp
    implicit none
!    real*8 :: Ex(0:CoreDist_Dim),AxD(0:CoreDist_Dim),Ey(0:CoreDist_Dim),AyD(0:CoreDist_Dim),Ar(0:CoreDist_Dim),nt_r,t_o
    real(dp), intent(out) :: Ex,AxD,Ey,AyD,Ar,Erh
    real(dp), intent(in) :: t_o
    integer, intent(in) :: CD_i
    real(dp) :: h_c,h,zeta,di,R,NN,dr,calD,Jx,Jx_int,Jy,Jy_int,dJx,dJy,JQ, nt_r
    real(dp) :: Intx,Inty,Intr,B,Lbd,weigh ! ,dh,max_h,fh_ce
    real(dp) :: fh,dfh,dfhl,dfha,alpha,lam, dlambda_hdr, dfhdr, Int_dh
    integer :: i,branch,il
    real(dp), save :: d_prev = -10.
    integer, save :: il_start = 1, CD_i_prev=-1
    common / time / zeta_c,hcto_c
    real(dp) :: zeta_c,hcto_c
    logical :: Frst,subtr
!
    h_c=hcto_c+T_o
    Ex=0.   ; Ey=0. ; Ar=0. ; AxD=0.    ; AyD=0.    ; Erh=0.
    if(ObsDist.ne.d_prev .or. CD_i.ne.CD_i_prev) il_start=1
    d_prev=ObsDist ;  CD_i_prev=CD_i
    !if(int(T_o/0.01)==655 .and. ObsDist.gt.400. .and. CD_i==0) write(2,*) 'a1',Ar
    lam=PancakeThi(CD_i)
    if(CD_i.eq.0) then
        dlambda_hdr=0.5*(PancakeThi(1)-PancakeThi(0))/(CoreDist_A(1)-CoreDist_A(0))
    elseif(CD_i.eq.CoreDist_Dim) then
        dlambda_hdr=(PancakeThi(CD_i)-PancakeThi(CD_i-1))/(CoreDist_A(CD_i)-CoreDist_A(CD_i-1))
    else
        dlambda_hdr=(PancakeThi(CD_i+1)-PancakeThi(CD_i-1))/(CoreDist_A(CD_i+1)-CoreDist_A(CD_i-1))
    endif
  Do branch=1,2
    Intx=0. ;   Inty=0. ; Intr=0. ; B=0.
    Frst=.true.
    subtr=.false.
    Int_dh=0.d0
    Do il=il_start,LamGrid_max-1   ! integration over pancake thickness
        lbd=L_I(il)  ! distance from Cherenkov divergence
        h=h_c-lbd    ! distance behind the shower front
        if(h.le.0.) exit
        if(h.gt.(10.*lam)) then
          il_start=il
          cycle
        endif
        zeta=zeta_I(branch,il)
        if(zeta.gt. AtmHei_dim*AtmHei_step) then
          if(branch.eq.2) exit
          if(branch.eq.1) cycle
        endif
        !        if(zeta.le.h) exit
        if(zeta.le.0.) exit
        weigh=JL_I(il)* dL_I(il)
        if((h_c-L_I(il+1)).le. 0. .and. (zeta_I(branch,il+1).gt. 0.) .and. (zeta_I(branch,il+1).lt. AtmHei_dim*AtmHei_step)) then ! last point on this branch, above ground
          if(il.le. LamGrid_Nrt) then  ! ne 1 then sqrt integral
            fh=(dL_I(il)-JL_I(il))**2/4.   !=ha
            lbd= (h_c+fh)/2.
            weigh=h_c-fh
            !            write(2,*) il,zeta_I(branch,il),zeta_I(branch,il-1),(zeta_c+zeta_I(branch,il))/2.
            !            write(2,*) L_I(il),fh,lbd
          else
            fh=L_I(il)-dL_I(il)/2.    !=ha
            lbd= (h_c + fh)/2.
            weigh=h_c - fh  ! JL_I(il) should be =1
          endif
          fh=zeta_I(branch,il-1)
          call getZeta1(lbd,fh,zeta)
          !          write(2,*) 'end-point',branch,il,lbd,zeta,zeta_c,zeta_I(branch,il),zeta_I(branch,il-1)
          h=h_c-lbd
          !          write(2,*) 'zeta,h',zeta,h,ex
        endif
        i=(zeta/AtmHei_step)
        di=zeta/AtmHei_step-i
        If(i.ge.AtmHei_dim) then
            i=AtmHei_dim-1
            di=1.d0
        endif
        !        write(2,*) 'di,i',di,i
        R=sqrt(ObsDist*ObsDist +zeta*zeta)
        NN=1.+ di*xi(i+1) + (1.-di)*xi(i)
        dr=   di*dxi(i+1) + (1.-di)*dxi(i)
        cald=abs(NN*(R - NN*zeta - dr*R*R)) ! retarded distance, including the needed factor n=NN
        !
        !        nt_r=zeta-h     ! - t_r (=negative t_r) to have a positive number
        !t_retarded(?)=zeta  Reatrded time  = zeta_I(branch,il-1)  ?
        !        nt_r=zeta  !  -h     ! - t_r (=negative t_r) to have a positive number
        !        i=(nt_r/AtmHei_step)
        !        di=nt_r/AtmHei_step-i
        JQ= di*IQ(i+1) + (1.-di)*IQ(i)
        Jx= di*Ix(i+1) + (1.-di)*Ix(i)
        Jx_int= di*Ix_int(i+1) + (1.-di)*Ix_int(i)  ! contribution from moving dipole
        dJx= (Ix(i+1) - Ix(i))/AtmHei_step  !  seems to give a vanishing O(10^{-4}) contribution to E
        Jy= di*Iy(i+1) + (1.-di)*Iy(i)
        !write(2,*) 'jy',jy,i,di
        Jy_int= di*Iy_int(i+1) + (1.-di)*Iy_int(i)
        dJy= (Iy(i+1) - Iy(i))/AtmHei_step  !  seems to give a vanishing contribution to E
        alpha=di*alpha_tr(i+1) + (1.-di)*alpha_tr(i) ! the parameter that determines the width of the pancake function
        !        Do CD_i=CoreDist_Dim,0,-1  ! Distance to the shower core, called rs in radial integral
        !          CoreDist = 0.5*ObsDist_Step + CoreDist_i*CoreDist_Step
        !    write(2,*) 'fieh1',h,CD_i,alpha,fh,dfhl,dfha
        call fieh(h,CD_i,alpha,fh,dfhl,dfha)  ! fh is being defined
        !    write(2,*) 'fieh2',h,CD_i,alpha,fh,dfhl,dfha
        dfh=dfhl + dfha *(alpha_tr(i+1)-alpha_tr(i))/AtmHei_step
        dfhdr = - (fh + dfhl*h ) *dlambda_hdr/lam
        if(Frst) then
            if(h.gt.(9.*lam) .and. branch.eq.1) then
                subtr=.true.
            else
                subtr=.false.
            endif
        endif
        !jx= & jy= terms give by far the dominat contribution (as compared to djx and djy)
        !jx=0.
        !jy=0.
        Ex =Ex  + weigh * (dfh*Jx -fh *dJx)/cald  ! (dfh*Jx/cald - fh *dJx/cald)
        Int_dh=Int_dh + weigh * dfh ! should be =0 when integrated over the complete height, any finite value is due to
        !                               integration-step-size problems which is why it is more stable to include the subtraction
        !  at large observer distances where 'cald' is almost constant over the integration range,
        !    the integral over dfh*Jx should yield zero, but this is from a cancellation of relatively large numbers.
        AxD=AxD + weigh*Jx_int*fh/cald   ! proportional to integrated current, thus like the moving dipole
        Ey =Ey  + (dfh*Jy/cald - fh *dJy/cald) * weigh
        AyD=AyD + weigh*Jy_int*fh/cald
        Ar =Ar  -(fh *JQ/cald) * weigh    ! still needs to be times dw/dr
        Erh =Erh  -(dfhdr *JQ/cald) * weigh    ! still needs to be times w(r)
        !if(ObsDist.lt.10 .and. T_o.lt.1.0 .and. CD_i.eq.0) write(2,*) il, h,zeta,cald,ex,i,weigh, dfh,Jx/cald, fh ,dJx
        !if(ObsDist.gt.410 .and. T_o.gt.44.2 .and. CD_i.eq.0) write(2,*) il, h,zeta,cald,ex,i
        !if(ObsDist.gt.410. .and. T_o.gt.44.2 .and. CD_i.eq.0) write(2,*) weigh , dfh,Jx ,fh ,dJx
        Frst=.false.
    enddo
    if(subtr) then
        Ex=Ex-Int_dh*Jx/cald
        Ey=Ey-Int_dh*Jy/cald
    endif
    !    ex=ex-b+intr
    !    write(2,*) 'int fh, int dfh',ex,CD_i,T_o
    !if(ObsDist.gt.410. .and. T_o.gt.44.2 .and. CD_i.eq.0) write(2,*) T_o,Ex,il
    !if(ObsDist.lt.10 .and. T_o.lt.1.0 .and. CD_i.eq.0) write(2,*) T_o,Ex,il,Int_dh*Jx/cald
    if(il.eq.LamGrid_max) then
        Ex=0.d0 ; AxD=0.d0 ; Ey=0.d0 ; AyD=0.d0 ; Ar=0.d0 ; Erh=0.d0
        return
    endif
  enddo
  !if(ObsDist.gt.410. .and. T_o.gt.45.2 .and. CD_i.eq.0) stop
!    write(3,*) 'cald,Exb=',cald,Ex
!    write(3,*) Jx,djx,cald,Exa,Exb
!    write(2,*) t_o,ex,ey,ar
    return
    end
!------------------------------

module mod_2d_in
    implicit none

    ! 初期変数
    integer,parameter::nt=1000
    real,parameter::dt=1.0e-2
    integer,parameter::nx0=100
    integer,parameter::ny0=100
    real,parameter::dx=1.0e-3
    real,parameter::dy=1.0e-3
    integer,parameter::lbound=1
    
    ! 計算領域
    integer,parameter::nx=nx0+2*lbound
    integer,parameter::ny=ny0+2*lbound

    ! アルミナ99.6%(Al2O3)
    real,parameter::lam=32.0
    real,parameter::rho=3.94e3
    real,parameter::c=0.78e3

    ! 初期温度
    real,parameter::u_init=293.0

    ! norm評価値
    real,parameter::norm_ok=1.0e-4
    integer,parameter::norm_check_num=100

    ! 境界条件
     
    !ディリクレ条件  =(0,温度(K))
    !ノイマン条件    =(1,何か)
    !熱流束          =(2,熱流束(W/m2))
    real,parameter::left_bd(2)=(/0., 573.0/)
    real,parameter::right_bd(2)=(/0., 293.0/)
    real,parameter::upper_bd(2)=(/1., 0.0/)
    real,parameter::lower_bd(2)=(/1., 0.0/)
    
    ! 出力情報
    integer,parameter::every_outnum=100
    character(len=5)::out_file_front="time_"
    character(len=3)::out_file_back="d-3"
    character(len=4)::out_ext=".csv"

    ! 初期発散抑え用
    real,parameter::init_small=1.0e-4

    ! 格子定義
    real::x(1:nx)
    real::y(1:ny)

    ! 温度, rhoc, 境界熱伝導率
    real::u(1:nx,1:ny)
    real::rhoc(1:nx,1:ny)
    real::lam_bdx(1:nx-1,1:ny-1)
    real::lam_bdy(1:nx-1,1:ny-1)

end module mod_2d_in
!-----------------------------------------------
! main program
!-----------------------------------------------
    program main
    use mod_2d_in
    implicit none
    integer::t
    real::t1,t2
    
    call cpu_time(t1)

    call setup()
    call output_T(0)
    do t=1,nt
        call u_cal()
        call boundary_left()
        call boundary_right()
        call boundary_upper()
        call boundary_lower()
        if (mod(t, every_outnum) ==0) then
            write(*,*) t
            call output_T(t)
        end if
    end do

    call cpu_time(t2)
    write(*,*) t2-t1,"sec"
    
    end program main


!-----------------------------------------------
! setup
!-----------------------------------------------
    subroutine setup()
        use mod_2d_in
        implicit none
        integer::i,j
        
        x=0.0
        do i = 1,nx
            x(i)=(i-1)*dx-dx/2
        end do
        y=0.0
        do j = 1,ny
            y(j)=(j-1)*dy-dy/2
        end do

        u=u_init
        u(1+lbound,1+lbound)=u(1+lbound,1+lbound)+init_small
        
        rhoc=rho*c
        lam_bdx=lam
        lam_bdy=lam
    
    end subroutine setup
  

!-----------------------------------------------
! u cal
!-----------------------------------------------
    subroutine u_cal()
        use mod_2d_in
        implicit none
        integer::i,j
        real::delta_u(1:nx,1:ny)
        real::R(1:nx,1:ny)
        real::sum_b_Ax
        real::sum_b
        real::x_numerator,y_numerator,keisu,Rx,Ry,norm2d
        integer::con,roop_num

        delta_u=0.0
        R=0.0
        sum_b=0.0
        
        do i = 1+lbound,nx-lbound
            do j = 1+lbound,ny-lbound
                x_numerator=lam_bdx(i,j)*(u(i+1,j)-u(i,j))-lam_bdx(i-1,j)*(u(i,j)-u(i-1,j))     ! 右辺Rのxの分子
                y_numerator=lam_bdy(i,j)*(u(i,j+1)-u(i,j))-lam_bdy(i,j-1)*(u(i,j)-u(i,j-1))     ! 右辺Rのyの分子
                R(i,j)=x_numerator/(rhoc(i,j)*dx**2)+y_numerator/(rhoc(i,j)*dy**2)                 ! 右辺Rの計算
                sum_b = sum_b + (R(i,j))**2 
            end do
        end do      
        write(*,*) sum_b


        con=0       ! ループ条件
        roop_num=0  ! ループ回数

        do while (con==0)
            do i = 1+lbound,nx-lbound
                do j = 1+lbound,ny-lbound
                    keisu = 1/dt+(lam_bdx(i,j)+lam_bdx(i-1,j))/(rhoc(i,j)*dx**2)+(lam_bdy(i,j)+lam_bdy(i,j-1))/(rhoc(i,j)*dy**2)    ! Δu_i,jの係数
                    Rx=lam_bdx(i,j)*delta_u(i+1,j)+lam_bdx(i-1,j)*delta_u(i-1,j)
                    Ry=lam_bdy(i,j)*delta_u(i,j+1)+lam_bdy(i,j-1)*delta_u(i,j-1)
                    delta_u(i,j)=(Rx/(rhoc(i,j)*dx**2)+Ry/(rhoc(i,j)*dy**2)+R(i,j))/keisu ! Δuの計算
                end do
            end do

            if (mod(roop_num+1 , norm_check_num) == 0) then
                sum_b_Ax=0 ! ノルムチェック用b-Ax
            
                do i = 1+lbound,nx-lbound
                    do j = 1+lbound,ny-lbound
                        ! ノルムチェック用b-Axの計算
                        keisu = 1/dt+(lam_bdx(i,j)+lam_bdx(i-1,j))/(rhoc(i,j)*dx**2)+(lam_bdy(i,j)+lam_bdy(i,j-1))/(rhoc(i,j)*dy**2)    ! Δu_i,jの係数
                        Rx=lam_bdx(i,j)*delta_u(i+1,j)+lam_bdx(i-1,j)*delta_u(i-1,j)
                        Ry=lam_bdy(i,j)*delta_u(i,j+1)+lam_bdy(i,j-1)*delta_u(i,j-1)
                        sum_b_Ax = sum_b_Ax+((Rx/(rhoc(i,j)*dx**2)+Ry/(rhoc(i,j)*dy**2)+R(i,j))/keisu-delta_u(i,j))**2
                    end do
                end do
                
                norm2d=(sum_b_Ax/sum_b)**0.5    ! 2乗ノルムの計算
                write(*,*) norm2d                        ! ノルムのプリント

                if (norm2d < norm_ok ) then                  ! ノルムが規定以下のときループ終了
                    con=1
                end if
            end if

            roop_num = roop_num+1    ! ループカウンター
        end do


        ! 時間の更新
        do i = 1+lbound,nx-lbound
            do j = 1+lbound,ny-lbound
                u(i,j)=u(i,j)+delta_u(i,j)
            end do
        end do

    end subroutine u_cal

!-----------------------------------------------
! boundary cal
!-----------------------------------------------
    subroutine boundary_left()
        use mod_2d_in
        implicit none
        integer::j

        if (left_bd(1)==0) then
            do j = 1+lbound,ny-1
                u(lbound,j)=left_bd(2)
            end do
        else if (left_bd(1)==1) then
            do j = 1+lbound,ny-1
                u(lbound,j)=u(lbound+1,j)
            end do
        else if (left_bd(1)==2) then
            do j = 1+lbound,ny-1
                u(lbound,j)=u(lbound+1,j)+left_bd(2)*dx/lam
            end do 
        else
            write(*,*) "boundary_left miss"
        end if
                
    end subroutine boundary_left

    subroutine boundary_right()
        use mod_2d_in
        implicit none
        integer::j

        if (right_bd(1)==0) then
            do j = 1+lbound,ny-1
                u(nx,j)=right_bd(2)
            end do
        else if (right_bd(1)==1) then
            do j = 1+lbound,ny-1
                u(nx,j)=u(nx-lbound,j)
            end do
        else if (right_bd(1)==2) then
            do j = 1+lbound,ny-1
                u(nx,j)=u(nx-lbound,j)+right_bd(2)*dx/lam
            end do
        else
            write(*,*) "boundary_right miss"
        end if

    end subroutine boundary_right
        
    subroutine boundary_lower()
        use mod_2d_in
        implicit none
        integer::i

        if (lower_bd(1)==0) then
            do i = 1+lbound,nx-1
                u(i,lbound)=lower_bd(2)
            end do
        else if (lower_bd(1)==1) then
            do i = 1+lbound,nx-1
                u(i,lbound)=u(i,lbound+1)
            end do
        else if (lower_bd(1)==2) then
            do i = 1+lbound,nx-1
                u(i,lbound)=u(i,lbound+1)+lower_bd(2)*dx/lam
            end do
        end if
                
    end subroutine boundary_lower

    subroutine boundary_upper()
        use mod_2d_in
        implicit none
        integer::i

        if (upper_bd(1)==0) then
            do i = 1+lbound,nx-1
                u(i,ny)=upper_bd(2)
            end do
        else if (upper_bd(1)==1) then
            do i = 1+lbound,nx-1
                u(i,ny)=u(i,ny-lbound)
            end do
        else if (upper_bd(1)==2) then
            do i = 1+lbound,nx-1
                u(i,ny)=u(i,ny-lbound)+upper_bd(2)*dx/lam
            end do
        end if
                
    end subroutine boundary_upper
  
    
!-----------------------------------------------
! output file
!-----------------------------------------------
    subroutine output_T(t)
    use mod_2d_in
    implicit none
    integer::i,j
    integer::t
    character(len=4)::s
    character(len=20)::sss
    
    write( s, '(I4.4)') t

    sss=out_file_front//s//out_file_back//out_ext

    open (unit=5,file=sss,action='write',status='replace')
    write(5,*) "x(m),y(m),T(K)"
    do i= 1,nx
        do j= 1,ny
             write(5,*) x(i),',',y(j),',',u(i,j)
        end do
    end do

    close(5)

    end subroutine output_T


program exp6
    implicit none
    real,dimension(10) :: w1
    real,dimension(10) :: w2
    real,dimension(10) :: w0
    real,dimension(10) :: w2p
    real,dimension(10) :: deltaw
    real,dimension(10) :: X
    real,dimension(10) :: tao
    real,dimension(10) :: U 
    real,allocatable,dimension(:) :: qv
    real,allocatable,dimension(:) :: qv1
    real,allocatable,dimension(:) :: uu 
    integer,allocatable,dimension(:) :: num2
    integer :: i,t,skip
    integer :: status
    real :: rho
    real :: temp
    real :: p 
    integer,dimension(10) :: num1

    write (*,1)
1   format (/,1x,'本程序用以处理化工实验6:流化床干燥实验')
    write (*,2)
2   format(/,1x,'输入数据时，每个数据用空格或者英文逗号隔开(PS:用英文逗号方便复制到Excel中)')
    write (*,3)
3   format (/,1x,'By 2020级化学基地班步允霆,使用语言-Fortran')
    write (*,4)
4 format (/,'--------------------------------------------------------------------------------',/)
    write (*,*)'输入114514直接进入测定流化曲线,输入1919810保持干燥速率曲线数据处理'
    read (*,*) skip 
    if (skip == 1919810) then 
    write (*,00411)
    00411 format (/,'----------干燥速率曲线----------',/)
    write (*,*) '输入称量瓶质量w0'
    read (*,*) w0
    write (*,*) '输入瓶+湿物料质量w1'
    read (*,*) w1
    write (*,*) '输入瓶+干物料质量w2' 
    read (*,*) w2 
    write (*,*) '输入取样时间(s)' 
    read (*,*) tao 
    do i=1,10
    w2p(i)=w2(i)-w0(i)
    deltaw(i)=(w1(i)-w2(i))
    X(i)=deltaw(i)/(w2p(i))
    if (i==1) then 
    U(i)=0.0
    else
    U(i)=-(X(i)-X(i-1))*3600.0/((tao(i)-tao(i-1)))
    end if 
    num1(i)=i
    end do 
    write (*,1000) 
    1000 format (/,/,2x,'NO.',3x,'  干物料质量  ',3x,'   含水量    ',3x,'      X      ',3x,'      U      ',/,&
                     2x,'===',3x,'=============',3x,'=============',3x,'=============',3x,'=============') 
    write (*,300) num1(1),w2p(1),deltaw(1),X(1)
    300 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5,3x,'       NONE  ')
    write (*,200) (num1(i),w2p(i),deltaw(i),X(i),U(i),i=2,10) 
    200 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5,3x,f13.5)
    else if (skip == 114514) then
    write (*,00011)
    00011 format (/,'----------测定流化曲线----------',/)
    write (*,*) '输入数据数量'
    read (*,*) t
    allocate (qv(t),stat=status)
    allocate (qv1(t),stat=status)
    allocate (uu(t),stat=status)
    allocate (num2(t),stat=status)
    write (*,*) '输入每组数据的流量(m3/h)'
    read (*,*) qv
    write (*,*) '输入当前温度下空气的密度(g/L)(或者输入1我帮你算密度)'
    read (*,*) rho 
    if (rho == 1.0) then
    write (*,*) '输入当时的温度(K)与压强(Pa)'
    read (*,*) temp,p 
    rho=1.293*(p/(10.0)**5)*(273.15/temp)
    write (*,999) rho
    999 format (/,2x,'此时的空气密度为',f8.5,'g/L',/,2x,'公式1.293*(实际压力/标准物理大气压)*(273.15/实际绝对温度)')
    write (*,888)
    888 format (2x,'由20化基孙守康提供')
    end if 
    do i=1,t
    qv1(i)=sqrt(1.205/rho)*qv(i)
    uu(i)=(4.0*qv1(i))/(3.14*((0.110)**2)*3600.0)
    num2(i)=i
    end do 
    write (*,10000) 
    10000 format (/,/,2x,'NO.',3x,'  矫正流量  ',3x,'   u(m/s)  ',/,&
                     2x,'===',3x,'=============',3x,'=============') 
    write (*,2000) (num2(i),qv1(i),uu(i),i=1,t)
    2000 format (2x,i3,3x,f13.5,3x,f13.5)
    end if 
    read (*,*)

end program exp6
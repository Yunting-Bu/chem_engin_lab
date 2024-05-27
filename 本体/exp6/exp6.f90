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
1   format (/,1x,'���������Դ�����ʵ��6:����������ʵ��')
    write (*,2)
2   format(/,1x,'��������ʱ��ÿ�������ÿո����Ӣ�Ķ��Ÿ���(PS:��Ӣ�Ķ��ŷ��㸴�Ƶ�Excel��)')
    write (*,3)
3   format (/,1x,'By 2020����ѧ���ذಽ����,ʹ������-Fortran')
    write (*,4)
4 format (/,'--------------------------------------------------------------------------------',/)
    write (*,*)'����114514ֱ�ӽ���ⶨ��������,����1919810���ָ��������������ݴ���'
    read (*,*) skip 
    if (skip == 1919810) then 
    write (*,00411)
    00411 format (/,'----------������������----------',/)
    write (*,*) '�������ƿ����w0'
    read (*,*) w0
    write (*,*) '����ƿ+ʪ��������w1'
    read (*,*) w1
    write (*,*) '����ƿ+����������w2' 
    read (*,*) w2 
    write (*,*) '����ȡ��ʱ��(s)' 
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
    1000 format (/,/,2x,'NO.',3x,'  ����������  ',3x,'   ��ˮ��    ',3x,'      X      ',3x,'      U      ',/,&
                     2x,'===',3x,'=============',3x,'=============',3x,'=============',3x,'=============') 
    write (*,300) num1(1),w2p(1),deltaw(1),X(1)
    300 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5,3x,'       NONE  ')
    write (*,200) (num1(i),w2p(i),deltaw(i),X(i),U(i),i=2,10) 
    200 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5,3x,f13.5)
    else if (skip == 114514) then
    write (*,00011)
    00011 format (/,'----------�ⶨ��������----------',/)
    write (*,*) '������������'
    read (*,*) t
    allocate (qv(t),stat=status)
    allocate (qv1(t),stat=status)
    allocate (uu(t),stat=status)
    allocate (num2(t),stat=status)
    write (*,*) '����ÿ�����ݵ�����(m3/h)'
    read (*,*) qv
    write (*,*) '���뵱ǰ�¶��¿������ܶ�(g/L)(��������1�Ұ������ܶ�)'
    read (*,*) rho 
    if (rho == 1.0) then
    write (*,*) '���뵱ʱ���¶�(K)��ѹǿ(Pa)'
    read (*,*) temp,p 
    rho=1.293*(p/(10.0)**5)*(273.15/temp)
    write (*,999) rho
    999 format (/,2x,'��ʱ�Ŀ����ܶ�Ϊ',f8.5,'g/L',/,2x,'��ʽ1.293*(ʵ��ѹ��/��׼�������ѹ)*(273.15/ʵ�ʾ����¶�)')
    write (*,888)
    888 format (2x,'��20�������ؿ��ṩ')
    end if 
    do i=1,t
    qv1(i)=sqrt(1.205/rho)*qv(i)
    uu(i)=(4.0*qv1(i))/(3.14*((0.110)**2)*3600.0)
    num2(i)=i
    end do 
    write (*,10000) 
    10000 format (/,/,2x,'NO.',3x,'  ��������  ',3x,'   u(m/s)  ',/,&
                     2x,'===',3x,'=============',3x,'=============') 
    write (*,2000) (num2(i),qv1(i),uu(i),i=1,t)
    2000 format (2x,i3,3x,f13.5,3x,f13.5)
    end if 
    read (*,*)

end program exp6
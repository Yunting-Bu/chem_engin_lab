module main
implicit none
contains

subroutine one(dout,din,Z,eta,d0,dr,rho,mu,H)
    use,intrinsic :: iso_fortran_env, only: wp => real64
    use pyplot_module
    implicit none
    type(pyplot) :: plt
    real,intent(in) :: dout
    real,intent(in) :: din
    real,intent(in) :: Z
    real,intent(in) :: eta 
    real,intent(in) :: d0
    real,intent(in) :: dr
    real,intent(in) :: rho 
    real,intent(in) :: mu 
    real,allocatable,dimension(:) :: Pout
    real,allocatable,dimension(:) :: Pin
    real,allocatable,dimension(:) :: uout
    real,allocatable,dimension(:) :: uin
    real(wp),allocatable,dimension(:) :: qv
    real(wp),allocatable,dimension(:),intent(out) :: H
    real(wp),allocatable,dimension(:) :: etae
    real(wp),allocatable,dimension(:) :: Pe 
    real,allocatable,dimension(:) :: P 
    real,allocatable,dimension(:) :: eread
    real(wp),allocatable,dimension(:) :: C 
    real,allocatable,dimension(:) :: qv0
    real(wp),allocatable,dimension(:) :: Re 
    real,allocatable,dimension(:) :: u0
    integer,allocatable,dimension(:) :: num
    real :: pi=3.14
    real :: A0
    integer :: i,t,status

    write (*,*) '������������'
    read (*,*) t
    allocate (qv(t),stat=status)
    allocate (Pout(t),stat=status)
    allocate (Pin(t),stat=status)
    allocate (uout(t),stat=status)
    allocate (uin(t),stat=status)
    allocate (num(t),stat=status)
    allocate (H(t),stat=status)
    allocate (etae(t),stat=status)
    allocate (eread(t),stat=status)
    allocate (P(t),stat=status)
    allocate (Pe(t),stat=status)
    allocate (Re(t),stat=status)
    allocate (qv0(t),stat=status)
    allocate (C(t),stat=status)
    allocate (u0(t),stat=status)
    write (*,*) '����ÿ�����ݵ�����(m3/h)'
    read (*,*) qv
    write (*,*) '����ÿ���������ѹ��P1(MPa)'
    read (*,*) Pin
    write (*,*) '����ÿ�����ݳ���ѹ��P2(MPa)'
    read (*,*) Pout
    write (*,*) '����ÿ�����ݵ������(kW)'
    read (*,*) eread
    write (*,*) '����ÿ�������������ƶ���(KPa)'
    read (*,*) qv0 
    A0=(pi*(d0**2.0))/4
    do i=1,t 
    uout(i)=qv(i)/(pi*900.0*((dout)**2.0))
    uin(i)=qv(i)/(pi*900.0*((din)**2.0))
    H(i)=Z+((Pout(i)-Pin(i))*1000000.0)/(rho*9.81)+((uout(i)**2.0)-(uin(i)**2.0))/(2.0*9.81)
    P(i)=eread(i)*eta 
    Pe(i)=(H(i)*qv(i)*rho*9.81)/(3600.0*1000.0)
    etae(i)=Pe(i)/P(i)
    C(i)=(qv(i)/(A0*3600.0))*sqrt(rho/(2.0*qv0(i)*1000.0))
    u0(i)=qv(i)/(900.0*pi*dr*dr)
    Re(i)=(dr*u0(i)*rho)/mu 
    num(i)=i 
    end do 
    write (*,100)
    100 format (/,'---------------------���ı����ܲⶨ1---------------------')
    write (*,101) 
    101 format (/,2x,'NO.',3x,'   u��(m/s)  ',3x,'   u��(m/s)  ',3x,'     ���     ',/,&
                  2x,'===',3x,'=============',3x,'=============',3x,'=============') 
    write (*,102) (num(i),uout(i),uin(i),H(i),i=1,t)
    102 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5)
    write (*,103)
    103 format (/,'---------------------���ı����ܲⶨ2---------------------')
    write (*,104) 
    104 format (/,2x,'NO.',3x,'  �Ṧ��(KW)  ',3x,'    Pe(KW)  ',3x,'     Ч��     ',/,&
                  2x,'===',3x,'=============',3x,'=============',3x,'=============') 
    write (*,105) (num(i),P(i),Pe(i),etae(i),i=1,t)
    105 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5)
    write (*,106)
    106 format (/,'---------------------���������ܲ���--------------------')
    write (*,107) 
    107 format (/,2x,'NO.',3x,'   ����ϵ��    ',3x,'  ����(m/s) ',3x,'      Re     ',/,&
                  2x,'===',3x,'=============',3x,'=============',3x,'=============') 
    write (*,108) (num(i),C(i),u0(i),Re(i),i=1,t)
    108 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5)
    call plt%initialize(grid=.true.,xlabel='Re',ylabel='C0',&
                     title='Plot of C0~Re',legend=.true.)
    call plt%add_plot(Re,C,label='C0~Re',linestyle='b-o',markersize=5,linewidth=2)
    call plt%savefig('C0~Re.png', pyfile='C0~Re.py')
    call plt%showfig('C0~Re.png')
    call plt%initialize(grid=.true.,xlabel='qv',&
                     title='total',legend=.true.)
    call plt%add_plot(qv,H/30.0_wp,label='H divide by 30.0~qv',linestyle='b-o',markersize=5,linewidth=2)
    call plt%add_plot(qv,Pe,label='Pe~qv',linestyle='r-o',markersize=5,linewidth=2)
    call plt%add_plot(qv,etae,label='eta~qv',linestyle='g-o',markersize=5,linewidth=2,ylim=(/0.0_wp,0.8_wp/))
    call plt%savefig('total.png', pyfile='total.py') 
    call plt%showfig('total.png')
    write (*,999)
    999 format(/,2x,'ͼ��total.png��C0~Re.png',/)

    end subroutine one 

    subroutine two(dout,din,Z,rho,mu,Hin)
    use,intrinsic :: iso_fortran_env, only: wp => real64
    use pyplot_module
    implicit none
    type(pyplot) :: plt
    real,intent(in) :: dout
    real,intent(in) :: din
    real,intent(in) :: Z
    real,intent(in) :: rho 
    real,intent(in) :: mu 
    real,allocatable,dimension(:) :: Pout
    real,allocatable,dimension(:) :: Pin
    real,allocatable,dimension(:) :: uout
    real,allocatable,dimension(:) :: uin
    real(wp),allocatable,dimension(:) :: qv
    real(wp),allocatable,dimension(:) :: H
    real(wp),allocatable,dimension(:) :: H1
    real(wp),allocatable,dimension(:),intent(in) :: Hin
    integer,allocatable,dimension(:) :: num
    real :: pi=3.14
    integer :: i,t,status

    write (*,*) '������������'
    read (*,*) t
    allocate (qv(t),stat=status)
    allocate (Pout(t),stat=status)
    allocate (Pin(t),stat=status)
    allocate (uout(t),stat=status)
    allocate (uin(t),stat=status)
    allocate (num(t),stat=status)
    allocate (H(t),stat=status)
    allocate (H1(t),stat=status)
    write (*,*) '����ÿ�����ݵ�����(m3/h)'
    read (*,*) qv
    write (*,*) '����ÿ���������ѹ��P1(MPa)'
    read (*,*) Pin
    write (*,*) '����ÿ�����ݳ���ѹ��P2(MPa)'
    read (*,*) Pout
    do i=1,t 
    uout(i)=qv(i)/(pi*900.0*((dout)**2.0))
    uin(i)=qv(i)/(pi*900.0*((din)**2.0))
    H(i)=Z+((Pout(i)-Pin(i))*1000000.0)/(rho*9.81)+((uout(i)**2.0)-(uin(i)**2.0))/(2.0*9.81)
    num(i)=i 
    H1(i)=Hin(i)
    end do 
    write (*,200)
    200 format (/,'---------------------��·�������߲ⶨ---------------------')
    write (*,201) 
    201 format (/,2x,'NO.',3x,'   u��(m/s)  ',3x,'   u��(m/s)  ',3x,'     ���     ',/,&
                  2x,'===',3x,'=============',3x,'=============',3x,'=============') 
    write (*,202) (num(i),uout(i),uin(i),H(i),i=1,t)
    202 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5)
    call plt%initialize(grid=.true.,xlabel='qv',ylabel='H',&
                     title='pipeline',legend=.true.)
    call plt%add_plot(qv,H,label='need',linestyle='b-o',markersize=5,linewidth=2)
    call plt%add_plot(qv,H1,label='provide',linestyle='r-o',markersize=5,linewidth=2)
    call plt%savefig('pipeline.png', pyfile='pipeline.py') 
    call plt%showfig('pipeline.png')
    write (*,998)
    998 format(/,2x,'ͼ��pipeline.png',/)
    end subroutine two 

    subroutine three(d,mu,rho,L)
    use,intrinsic :: iso_fortran_env, only: wp => real64
    use pyplot_module

    implicit none
    type(pyplot) :: plt
    real,intent(in) :: d
    real,intent(in) :: rho 
    real,intent(in) :: mu 
    real,intent(in) :: L
    real,allocatable,dimension(:) :: P1 
    real,allocatable,dimension(:) :: P2 
    real,allocatable,dimension(:) :: u
    real,allocatable,dimension(:) :: qv
    real(wp),allocatable,dimension(:) :: lamda
    real(wp),allocatable,dimension(:) :: Re
    integer,allocatable,dimension(:) :: num
    real :: pi=3.14
    integer :: i,t1,t2,status

    write (*,*) '������mmH2OΪ��λ����������'
    read (*,*) t1
    write (*,*) '������kPaΪ��λ����������'
    read (*,*) t2
    allocate (P1(t1),stat=status)
    allocate (P2(t2),stat=status)
    allocate (u(t1+t2),stat=status)
    allocate (lamda(t1+t2),stat=status)
    allocate (num(t1+t2),stat=status)
    allocate (qv(t1+t2),stat=status)
    allocate (Re(t1+t2),stat=status)
    write (*,*) '����ÿ�����ݵ�����(L/h)'
    read (*,*) qv
    write (*,*) '����ÿ������ѹ��P(mmH2O)'
    read (*,*) P1
    write (*,*) '����ÿ������ѹ��P(kPa)'
    read (*,*) P2 
    do i=1,t1
    u(i)=(qv(i)*0.001)/(900.0*pi*d*d)
    lamda(i)=((2.0*d)/(rho*L))*((P1(i)*0.001*101325)/(10.33*u(i)*u(i)))
    Re(i)=d*u(i)*rho/mu 
    num(i)=i
    end do
    do i=1,t2 
    u(i+t1)=(qv(i+t1)*0.001)/(900.0*pi*d*d)
    lamda(i+t1)=((2.0*d)/(rho*L))*((P2(i)*1000.0)/(u(i+t1)*u(i+t1)))
    Re(i+t1)=d*u(i+t1)*rho/mu 
    num(i+t1)=i+t1 
    end do 
    write (*,300)
    300 format (/,'---------------------ֱ������ʵ��---------------------')
    write (*,301) 
    301 format (/,2x,'NO.',3x,'    u(m/s)   ',3x,'   Ħ��ϵ��  ',3x,'      Re      ',/,&
                  2x,'===',3x,'=============',3x,'=============',3x,'=============') 
    write (*,302) (num(i),u(i),lamda(i),Re(i),i=1,t1+t2)
    302 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5)
    call plt%initialize(grid=.true.,xlabel='Re',ylabel='Lambda',&
                     title='Plot of Lambda~Re',legend=.true.)
    call plt%add_plot(Re,lamda,label='Lambda~Re',linestyle='b-o',markersize=5,linewidth=2,xscale='log',yscale='log')
    call plt%savefig('Lambda~Re1.png', pyfile='Lambda~Re1.py')
    call plt%showfig('Lambda~Re1.png')
    write (*,997)
    997 format(/,2x,'ͼ��Lambda~Re1.png',/)

    end subroutine three 

    subroutine four(d,mu,rho,L)
    use,intrinsic :: iso_fortran_env, only: wp => real64
    use pyplot_module

    implicit none
    type(pyplot) :: plt
    real,intent(in) :: d
    real,intent(in) :: rho 
    real,intent(in) :: mu 
    real,intent(in) :: L
    real,allocatable,dimension(:) :: P1 
    real,allocatable,dimension(:) :: P2 
    real,allocatable,dimension(:) :: u
    real,allocatable,dimension(:) :: qv
    real(wp),allocatable,dimension(:) :: lamda
    real(wp),allocatable,dimension(:) :: Re
    integer,allocatable,dimension(:) :: num
    real :: pi=3.14
    integer :: i,t1,t2,status

    write (*,*) '������mmH2OΪ��λ����������'
    read (*,*) t1
    write (*,*) '������kPaΪ��λ����������'
    read (*,*) t2
    allocate (P1(t1),stat=status)
    allocate (P2(t2),stat=status)
    allocate (u(t1+t2),stat=status)
    allocate (lamda(t1+t2),stat=status)
    allocate (num(t1+t2),stat=status)
    allocate (qv(t1+t2),stat=status)
    allocate (Re(t1+t2),stat=status)
    write (*,*) '����ÿ�����ݵ�����(L/h)'
    read (*,*) qv
    write (*,*) '����ÿ������ѹ��P(mmH2O)'
    read (*,*) P1
    write (*,*) '����ÿ������ѹ��P(kPa)'
    read (*,*) P2 
    do i=1,t1
    u(i)=(qv(i)*0.001)/(900.0*pi*d*d)
    lamda(i)=((2.0*d)/(rho*L))*((P1(i)*0.001*101325)/(10.33*u(i)*u(i)))
    Re(i)=d*u(i)*rho/mu 
    num(i)=i
    end do
    do i=1,t2 
    u(i+t1)=(qv(i+t1)*0.001)/(900.0*pi*d*d)
    lamda(i+t1)=((2.0*d)/(rho*L))*((P2(i)*1000.0)/(u(i+t1)*u(i+t1)))
    Re(i+t1)=d*u(i+t1)*rho/mu 
    num(i+t1)=i+t1 
    end do 
    write (*,322)
    322 format (/,'---------------------ֱ������ʵ��---------------------')
    write (*,321) 
    321 format (/,2x,'NO.',3x,'    u(m/s)   ',3x,'   Ħ��ϵ��  ',3x,'      Re      ',/,&
                  2x,'===',3x,'=============',3x,'=============',3x,'=============') 
    write (*,332) (num(i),u(i),lamda(i),Re(i),i=1,t1+t2)
    332 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5)
    call plt%initialize(grid=.true.,xlabel='Re',ylabel='Lambda',&
                     title='Plot of Lambda~Re',legend=.true.)
    call plt%add_plot(Re,lamda,label='Lambda~Re',linestyle='b-o',markersize=5,linewidth=2,xscale='log',yscale='log')
    call plt%savefig('Lambda~Re2.png', pyfile='Lambda~Re2.py')
    call plt%showfig('Lambda~Re2.png')
    write (*,996)
    996 format(/,2x,'ͼ��Lambda~Re2.png',/)

    end subroutine four

    subroutine five(df,rho)
    implicit none
    real,intent(in) :: df
    real,intent(in) :: rho
    real,allocatable,dimension(:) :: qv 
    real,allocatable,dimension(:) :: Pb
    real,allocatable,dimension(:) :: Pa 
    real,allocatable,dimension(:) :: u 
    real,allocatable,dimension(:) :: zeta
    real,allocatable,dimension(:) :: Pf
    integer,allocatable,dimension(:) :: num
    real :: pi=3.14
    integer :: i,t,status

    write (*,*) '������������'
    read (*,*) t
    allocate (qv(t),stat=status)
    allocate (Pb(t),stat=status)
    allocate (Pa(t),stat=status)
    allocate (u(t),stat=status)
    allocate (Pf(t),stat=status)
    allocate (zeta(t),stat=status)
    allocate (num(t),stat=status)
    write (*,*) '����ÿ�����ݵ�����(L/h)'
    read (*,*) qv
    write (*,*) '����ÿ�����ݵĽ���ѹ��(kPa)'
    read (*,*) Pb 
    write (*,*) '����ÿ�����ݵ�Զ��ѹ��(kPa)'
    read (*,*) Pa 
    do i=1,t 
    Pf(i)=2.0*Pb(i)-Pa(i)
    u(i)=(qv(i)*0.001)/(900.0*pi*df*df)
    zeta(i)=(2.0*Pf(i)*1000.0)/(rho*u(i)*u(i))
    num(i)=i
    end do 
    write (*,400)
    400 format (/,'---------------------�ֲ�����ʵ��--------------------')
    write (*,401) 
    401 format (/,2x,'NO.',3x,'    u(m/s)   ',3x,' �ֲ�����ϵ�� ',3x,' ѹǿ��(kPa) ',/,&
                  2x,'===',3x,'=============',3x,'=============',3x,'=============') 
    write (*,402) (num(i),u(i),zeta(i),Pf(i),i=1,t)
    402 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5)
    end subroutine five 
end module main


program exp1
    use main
    use pyplot_module
    use,intrinsic :: iso_fortran_env, only: wp => real64
    implicit none

    real(wp),allocatable,dimension(:) :: H    
    real :: dout,din
    real :: d0,dr
    real :: dt,ds
    real :: df
    real :: Z,eta
    real :: L
    real :: mu,rho
    integer :: choose

    write (*,1)
    1   format (/,1x,'���������Դ�����ʵ��1:���������ۺ�ʵ������')
    write (*,2)
    2   format(/,1x,'��������ʱ��ÿ�������ÿո����Ӣ�Ķ��Ÿ���')
    write (*,3)
    3   format (/,1x,'By 2020����ѧ���ذಽ����,ʹ������-Fortran')
    write (*,4)
    4 format (/,'--------------------------------------------------------------',/)
    open (10,file='parameters.txt',status='old',action='read')
    read (10,*) dout,din,d0,dr,ds,dt,df,Z,L,eta
    write (*,99) 
    99 format (/,2x,'!!!!!�˳������ò�����¼��parameters.txt��!!!!!',/)
    write (*,*) '���������ܶ�(kg/m3)��ճ��(Pa��s)'
    read (*,*) rho,mu
    write (*,5) 
    5 format (/,2x,'------------��ѡ�����¹���------------',/)
    write (*,6)
    6 format (2x,'1.���ı����ܲⶨ�����������ܲ���ʵ��',/,&
                2x,'2.��·�������߲ⶨʵ��',/,&
                2x,'3.ֱ������ʵ��(�⻬��)',/,&
                2x,'4.ֱ������ʵ��(�ֲڹ�)',/,&
                2x,'5.�ֲ�����ʵ��(ȫ��)',/,&
                2x,'6.�ֲ�����ʵ��(�뿪)',/,&
                2x,'0.�˳�')
    read (*,*) choose 
    do 
    if (choose == 0) exit
    if (choose == 1) then
    call one(dout,din,Z,eta,d0,dr,rho,mu,H)
    else if(choose == 2) then 
    call two(dout,din,Z,rho,mu,H)
    else if(choose == 3) then 
    call three(ds,mu,rho,L)
    else if(choose == 4) then 
    call four(dt,mu,rho,L)
    else if(choose == 5) then 
    call five(df,rho)
    else if(choose == 6) then 
    call five(df,rho)
    end if 
    write (*,7) 
    7 format (/,2x,'------------��ѡ�����¹���------------',/)
    write (*,8)
    8 format (2x,'1.���ı����ܲⶨ�����������ܲ���ʵ��',/,&
                2x,'2.��·�������߲ⶨʵ��',/,&
                2x,'3.ֱ������ʵ��(�⻬��)',/,&
                2x,'4.ֱ������ʵ��(�ֲڹ�)',/,&
                2x,'5.�ֲ�����ʵ��(ȫ��)',/,&
                2x,'6.�ֲ�����ʵ��(�뿪)',/,&
                2x,'0.�˳�')
    read (*,*) choose 
    end do 
    read (*,*)


end program exp1













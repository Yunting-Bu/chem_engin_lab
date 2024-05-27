program exp2
    implicit none

    real :: d17,d2,L,cp,c0,d16
    integer :: choose,t,n1,n2

    write (*,1)
    1   format (/,1x,'���������Դ�����ʵ��2:��������ʵ������')
    write (*,2)
    2   format(/,1x,'��������ʱ��ÿ�������ÿո����Ӣ�Ķ��Ÿ���')
    write (*,3)
    3   format (/,1x,'By 2020����ѧ���ذಽ����,ʹ������-Fortran')
    write (*,4)
    4 format (/,'--------------------------------------------------------------',/)
    open (10,file='parameters.txt',status='old',action='read')
    read (10,*) d2,d17,L,cp,c0,d16,n1,n2
    write (*,99) 
    99 format (/,2x,'!!!!!�˳������ò�����¼��parameters.txt��!!!!!',/)
    write (*,98) 
    98 format (/,2x,'���Բ������㹫ʽ����',/,&
                2x,'�ܶ�:1.128+(1.093-1.128)*(tm-40)/10',/,&
                2x,'ճ��(10^5):1.91+(1.96-1.91)*(tm-40)/10',/,&
                2x,'��(10^2):2.756+(2.826-2.756)*(tm-40)/10')

    write (*,5) 
    5 format (/,2x,'------------��ѡ�����¹���------------',/)
    write (*,6)
    6 format (2x,'1.�⻬�ܻ�����',/,&
                2x,'2.ǿ���ܻ�����',/,&
                2x,'3.�йܻ�����ȫ��ͨ',/,&
                2x,'4.�йܻ���������ͨ',/,&
                2x,'0.�˳�')
    read (*,*) choose 
    do 
    if (choose == 0) exit
    if (choose == 1) then
    call one(d17,d2,L,cp,c0,choose)
    else if (choose == 2) then
    call one(d17,d2,L,cp,c0,choose)
    else if (choose == 3) then
    call two(d17,d2,L,cp,c0,d16,n1)
    else if (choose == 4) then
    call two(d17,d2,L,cp,c0,d16,n2)
    end if 
    write (*,7) 
    7 format (/,2x,'------------��ѡ�����¹���------------',/)
    write (*,8)
    8 format (2x,'1.�⻬�ܻ�����',/,&
                2x,'2.ǿ���ܻ�����',/,&
                2x,'3.�йܻ�����ȫ��ͨ',/,&
                2x,'4.�йܻ���������ͨ',/,&
                2x,'0.�˳�')
    read (*,*) choose 
    end do 
    read (*,*)



end program exp2


subroutine one(d17,d2,L,cp,c0,choose)
    implicit none
    real,intent(in) :: d17,d2,L,cp,c0
    real,allocatable,dimension(:) :: p 
    real,allocatable,dimension(:) :: t1
    real,allocatable,dimension(:) :: t2
    real,allocatable,dimension(:) :: tw
    real,allocatable,dimension(:) :: u 
    real,allocatable,dimension(:) :: W 
    real,allocatable,dimension(:) :: Q 
    real,allocatable,dimension(:) :: Nu 
    real,allocatable,dimension(:) :: Nu0
    real,allocatable,dimension(:) :: Re 
    real,allocatable,dimension(:) :: Pr 
    real,allocatable,dimension(:) :: tm 
    real,allocatable,dimension(:) :: deltatm
    real,allocatable,dimension(:) :: rhotm
    real,allocatable,dimension(:) :: lamdatm
    real,allocatable,dimension(:) :: mutm
    real,allocatable,dimension(:) :: alpha
    real,allocatable,dimension(:) :: lnRe
    real,allocatable,dimension(:) :: lnPr
    real,allocatable,dimension(:) :: lnNu
    real,allocatable,dimension(:) :: div
    real,allocatable,dimension(:) :: mins
    integer,allocatable,dimension(:) :: num
    real :: pi=3.14
    integer :: i,status,t
    integer,intent(in) :: choose

    write (*,*) '������������'
    read (*,*) t
    allocate (p(t),stat=status)
    allocate (t1(t),stat=status)
    allocate (t2(t),stat=status)
    allocate (tw(t),stat=status)
    allocate (u(t),stat=status)
    allocate (W(t),stat=status)
    allocate (Q(t),stat=status)
    allocate (Nu(t),stat=status)
    allocate (Re(t),stat=status)
    allocate (Pr(t),stat=status)
    allocate (tm(t),stat=status)
    allocate (deltatm(t),stat=status)
    allocate (rhotm(t),stat=status)
    allocate (lamdatm(t),stat=status)
    allocate (mutm(t),stat=status)
    allocate (num(t),stat=status)
    allocate (alpha(t),stat=status)
    allocate (lnRe(t),stat=status)
    allocate (lnPr(t),stat=status)
    allocate (lnNu(t),stat=status)
    allocate (Nu0(t),stat=status)
    allocate (div(t),stat=status)
    allocate (mins(t),stat=status)
    write (*,*) '����ÿ�����ݵĿ�������ѹ��(kPa)'
    read (*,*) p 
    write (*,*) '����ÿ�����ݵ�����¶�(��)'
    read (*,*) t1
    write (*,*) '����ÿ�����ݵĳ����¶�(��)'
    read (*,*) t2
    write (*,*) '����ÿ�����ݵ�tw(��)'
    read (*,*) tw 
    do i=1,t
    tm(i)=(t2(i)+t1(i))/2
    rhotm(i)=1.128+(1.093-1.128)*(tm(i)-40)/10
    lamdatm(i)=(2.756+(2.826-2.756)*(tm(i)-40)/10)
    mutm(i)=(1.91+(1.96-1.91)*(tm(i)-40)/10)
    deltatm(i)=tw(i)-tm(i)
    u(i)=c0*((2000*p(i)/rhotm(i))**0.5)*d17*d17/(d2*d2)
    W(i)=0.25*pi*d2*d2*u(i)*rhotm(i)
    Q(i)=W(i)*cp*(t2(i)-t1(i))
    alpha(i)=Q(i)/(deltatm(i)*pi*d2*L)
    Re(i)=(d2*u(i)*rhotm(i))/(mutm(i)*0.00001)
    Nu(i)=(alpha(i)*0.02)/(lamdatm(i)*0.01)
    Pr(i)=(cp*mutm(i)*0.00001)/(lamdatm(i)*0.01)
    lnRe(i)=log(Re(i))
    lnPr(i)=log((Pr(i))**0.4)
    lnNu(i)=log(Nu(i))
    mins(i)=lnNu(i)-lnPr(i)
    num(i)=i
    end do 
    write (*,100)
    100 format (/,'---------------------���Բ���---------------------')
    write (*,101) 
    101 format (/,2x,'NO.',3x,'     �ܶ�    ',3x,'      ��      ',3x,'     ��     ',/,&
                  2x,'===',3x,'=============',3x,'=============',3x,'=============') 
    write (*,102) (num(i),rhotm(i),lamdatm(i),mutm(i),i=1,t)
    102 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5)
    write (*,103)
    103 format (/,'---------------------��ͨ�׹ܻ�����1---------------------')
    write (*,104) 
    104 format (/,2x,'NO.',3x,'    u(m/s)   ',3x,'��������(kg/s)',3x,'   ����(J)   ',/,&
                  2x,'===',3x,'=============',3x,'=============',3x,'=============') 
    write (*,105) (num(i),u(i),W(i),Q(i),i=1,t)
    105 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5)
    write (*,106)
    106 format (/,'---------------------��ͨ�׹ܻ�����2---------------------')
    write (*,107) 
    107 format (/,2x,'NO.',3x,'      ��      ',3x,'      Nu     ',3x,'      Re     ',3x,'      Pr     ',/,&
                  2x,'===',3x,'=============',3x,'=============',3x,'=============',3x,'=============') 
    write (*,108) (num(i),alpha(i),Nu(i),Re(i),Pr(i),i=1,t)
    108 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5,3x,f13.5)
     write (*,109)
    109 format (/,'---------------------��ͨ�׹ܻ�����3---------------------')
    write (*,110) 
    110 format (/,2x,'NO.',3x,'     lnNu    ',3x,'     lnRe    ',3x,' ln(Pr**0.4) ',3x,'lnNu/(Pr**0.4)'/,&
                  2x,'===',3x,'=============',3x,'=============',3x,'=============',3x,'=============') 
    write (*,111) (num(i),lnNu(i),lnRe(i),lnPr(i),mins(i),i=1,t)
    111 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5,3x,f13.5)
    if (choose == 1) then
    open (8,file='Nu.txt',status='old',action='readwrite')
    write (8,*) Nu
    close (8,status='keep')
    else if (choose == 2) then
    open (8,file='Nu.txt',status='old',action='readwrite')
    read (8,*) Nu0
    do i=1,t 
    div(i)=Nu(i)/Nu0(i)
    end do 
    write (*,112)
    112 format (/,2x,'ǿ���ȷֱ�Ϊ',/)
    write (*,113) (div(i),i=1,t)
    113 format (2x,f13.5)
    end if
    end subroutine one

subroutine two(d17,d2,L,cp,c0,d16,n)
    implicit none
    real,intent(in) :: d17,d2,L,cp,c0,d16
    real,allocatable,dimension(:) :: p 
    real,allocatable,dimension(:) :: t1
    real,allocatable,dimension(:) :: t2
    real,allocatable,dimension(:) :: TT1 
    real,allocatable,dimension(:) :: TT2 
    real,allocatable,dimension(:) :: W 
    real,allocatable,dimension(:) :: Q 
    real,allocatable,dimension(:) :: tm 
    real,allocatable,dimension(:) :: deltatm
    real,allocatable,dimension(:) :: rhotm
    real,allocatable,dimension(:) :: lamdatm
    real,allocatable,dimension(:) :: mutm
    real,allocatable,dimension(:) :: K 
    real,allocatable,dimension(:) :: u      
    integer,allocatable,dimension(:) :: num
    real :: pi=3.14
    integer :: i,status,t,n

    write (*,*) '������������'
    read (*,*) t
    allocate (p(t),stat=status)
    allocate (t1(t),stat=status)
    allocate (t2(t),stat=status)
    allocate (u(t),stat=status)
    allocate (W(t),stat=status)
    allocate (Q(t),stat=status)
    allocate (K(t),stat=status)
    allocate (tm(t),stat=status)
    allocate (deltatm(t),stat=status)
    allocate (rhotm(t),stat=status)
    allocate (lamdatm(t),stat=status)
    allocate (mutm(t),stat=status)
    allocate (TT1(t),stat=status)
    allocate (TT2(t),stat=status)
    allocate (num(t),stat=status)   
    write (*,*) '����ÿ�����ݵĿ�������ѹ��(kPa)'
    read (*,*) p 
    write (*,*) '����ÿ�����ݵĿ�������¶�(��)'
    read (*,*) t1
    write (*,*) '����ÿ�����ݵĿ��������¶�(��)'
    read (*,*) t2
    write (*,*) '����ÿ�����ݵ���������¶�(��)'
    read (*,*) TT1
    write (*,*) '����ÿ�����ݵ����������¶�(��)'
    read (*,*) TT2 
    do i=1,t 
    tm(i)=(t2(i)+t1(i))/2
    rhotm(i)=1.128+(1.093-1.128)*(tm(i)-40)/10
    lamdatm(i)=(2.756+(2.826-2.756)*(tm(i)-40)/10)
    mutm(i)=(1.91+(1.96-1.91)*(tm(i)-40)/10)
    u(i)=(c0*d17*d17*((2000*p(i)/rhotm(i))**0.5))/(d16*d16*n)
    W(i)=0.25*pi*d16*d16*n*u(i)*rhotm(i)
    Q(i)=W(i)*cp*(t2(i)-t1(i))    
    deltatm(i)=((TT2(i)-t2(i))-(TT1(i)-t1(i)))/log((TT2(i)-t2(i))/(TT1(i)-t1(i)))
    K(i)=(cp*W(i)*(t2(i)-t1(i)))/(deltatm(i)*pi*n*L*d16)
    num(i)=i
    end do 
    write (*,200)
    200 format (/,'---------------------���Բ���---------------------')
    write (*,201) 
    201 format (/,2x,'NO.',3x,'     �ܶ�    ',3x,'      ��      ',3x,'     ��     ',/,&
                  2x,'===',3x,'=============',3x,'=============',3x,'=============') 
    write (*,202) (num(i),rhotm(i),lamdatm(i),mutm(i),i=1,t)
    202 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5)
    write (*,203)
    203 format (/,'---------------------�йܻ�����1---------------------')
    write (*,204) 
    204 format (/,2x,'NO.',3x,'    u(m/s)   ',3x,'��������(kg/s)',3x,'   ����(J)   ',/,&
                  2x,'===',3x,'=============',3x,'=============',3x,'=============') 
    write (*,205) (num(i),u(i),W(i),Q(i),i=1,t)
    205 format (2x,i3,3x,f13.5,3x,f13.5,3x,f13.5)
    write (*,206)
    206 format (/,'---------------------�йܻ�����2---------------------')
    write (*,207) 
    207 format (/,2x,'NO.',3x,'       K     ',3x,'   deltatm   ',/,&
                  2x,'===',3x,'=============',3x,'=============') 
    write (*,208) (num(i),K(i),deltatm(i),i=1,t)
    208 format (2x,i3,3x,f13.5,3x,f13.5)
    end subroutine two



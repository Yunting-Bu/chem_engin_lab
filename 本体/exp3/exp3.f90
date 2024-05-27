program exp3

    implicit none
    real :: d
    real :: K,R
    real :: rou
    real :: mu
    integer :: t,j,ifc
    integer :: status
    real,allocatable,dimension(:) :: V
    real,allocatable,dimension(:) :: I
    real,allocatable,dimension(:) :: n
    real,allocatable,dimension(:) :: Np
    real,allocatable,dimension(:) :: NN
    real,allocatable,dimension(:) :: Re
    integer,allocatable,dimension(:) :: num
    write (*,1)
1   format (/,1x,'本程序用以处理化工实验3:搅拌器性能测定实验数据')
    write (*,2)
2   format(/,1x,'输入数据时，每个数据用空格或者英文逗号隔开')
    write (*,3)
3   format (/,1x,'By 2020级化学基地班步允霆,使用语言-Fortran')
    write (*,4)
4 format (/,'--------------------------------------------------------------',/)
    write (*,*) '输入数据数量'
    read (*,*) t
    allocate (n(t),stat=status)
    allocate (I(t),stat=status)
    allocate (V(t),stat=status)
    allocate (NN(t),stat=status)
    allocate (Np(t),stat=status)
    allocate (Re(t),stat=status)
    allocate (num(t),stat=status)
    open (10,file='parameters.txt',status='old',action='read')
    read (10,*) R,K 
    write (*,*) '依次输入密度(kg/m3)、粘度(Pa·s)、直径(m)'
    read (*,*) rou,mu,d
    write (*,*) '输入每组数据的转速(rpm)'
    read (*,*) n
    write (*,*) '输入每组数据的电流(A)'
    read (*,*) I
    write (*,*) '输入每组数据的电压(V)'
    read (*,*) V
    R=30.0
    K=0.177
    do j=1,t
        NN(j)=I(j)*V(j)-((I(j)**2.0)*R+K*(n(j)/60.0)**1.2)
        Np(j)=NN(j)/((rou*((n(j)/60.0))**3.0)*(d**5.0))
        Re(j)=((d**2.0)*(n(j)/60.0)*rou)/mu
        num(j)=j
    end do
    write (*,1000) 
    1000 format (/,/,2x,'NO.',3x,'    N/P     ',3x,'     Re     ',3x,'     Np     ',/,&
                     2x,'===',3x,'============',3x,'============',3x,'============') 
    write (*,200) (num(j),NN(j),Re(j),Np(j),j=1,t)
    200 format (2x,i3,3x,f12.4,3x,f12.4,3x,f12.4)
    write (*,222)
    222 format (/,1x,'是否计算下一组数据?是输入1,否输入0')
    read (*,*) ifc
    do 
    if (ifc == 0) exit
    write (*,*) '输入数据数量'
    read (*,*) t
    allocate (n(t),stat=status)
    allocate (I(t),stat=status)
    allocate (V(t),stat=status)
    allocate (NN(t),stat=status)
    allocate (Np(t),stat=status)
    allocate (Re(t),stat=status)
    allocate (num(t),stat=status)
    write (*,*) '输入每组数据的转速(rpm)'
    read (*,*) n
    write (*,*) '输入每组数据的电流(A)'
    read (*,*) I
    write (*,*) '输入每组数据的电压(V)'
    read (*,*) V
    do j=1,t
        NN(j)=I(j)*V(j)-((I(j)**2.0)*R+K*(n(j)/60.0)**1.2)
        Np(j)=NN(j)/((rou*((n(j)/60.0))**3.0)*(d**5.0))
        Re(j)=((d**2.0)*(n(j)/60.0)*rou)/mu
        num(j)=j
    end do
    write (*,1900) 
    1900 format (/,/,2x,'NO.',3x,'    N/P     ',3x,'     Re     ',3x,'     Np     ',/,&
                     2x,'===',3x,'============',3x,'============',3x,'============') 
    write (*,2100) (num(j),NN(j),Re(j),Np(j),j=1,t)
    2100 format (2x,i3,3x,f12.4,3x,f12.4,3x,f12.4)
    write (*,2122)
    2122 format (/,1x,'是否计算下一组数据?是输入1,否输入0')
    read (*,*) ifc
    end do
                 read (*,*)
    end program exp3


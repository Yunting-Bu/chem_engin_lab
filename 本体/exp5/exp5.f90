program exp5

    implicit none
    real :: Dd,Dw,inD,outD
    real :: Xd,Xw,inX,outX
    real :: aa,bb,cc,xn
    real :: Eml
    real :: N,Ne,E
    real :: wd,ww,inw,outw
    integer :: ifc,ifnext,i
    real :: x,y,temp

    write (*,1)
1   format (/,1x,'���������Դ�����ʵ��5: ��ʽ���������ܲⶨ����')
    write (*,2)
2   format(/,1x,'��������ʱ��ÿ�������ÿո����Ӣ�Ķ��Ÿ���')
    write (*,3)
3   format (/,1x,'By 2020����ѧ���ذಽ����,ʹ������-Fortran')
    write (*,4)
4 format (/,'--------------------------------------------------------------------------------',/)
    write(*,*)'�������������۹��ʡ���Һ�۹��ʡ������n����Һ���۹��ʡ��뿪��n����Һ���۹���'
    read (*,*) Dd,Dw,inD,outD
    write (*,*)'����ʵ��������'
    read (*,*) Ne
    wd=58.8441-42.6132*Dd
    ww=58.8441-42.6132*Dw
    inw=58.8441-42.6132*inD
    outw=58.8441-42.6132*outD
    Xd=(30.0*wd)/(23.0+7.0*wd)
    Xw=(30.0*ww)/(23.0+7.0*ww)
    inX=(30.0*inw)/(23.0+7.0*inw)
    outX=(30.0*outw)/(23.0+7.0*outw)
    aa=-0.72001
    bb=1.66931
    cc=0.02466-inX
    call equ(aa,bb,cc,xn)
    Eml=(inX-outX)/(inX-xn)
    y=Xd 
    i=1
    write (*,110) 
110 format (/,'-----------����������-----------',/)
    do 
    cc=0.02466-y 
    call equ(aa,bb,cc,x)
    write (*,123) i,x,y
    if (x < Xw) exit 
123 format (/,2x,'NO.',I2,'   x=',f10.5,'   y=',f10.5)
    temp=x  
    y=temp 
    i=i+1 
    end do 
    N=i 
    E=N/Ne
    write (*,124) i 
124 format (/,2x,'��',i2,'����������,���xy����������',/)
    write(*,5)
5 format (/,2x,'!!!!!ע��!!!!!',/)
    write(*,6)
6 format (2x,'���������y=-0.72001x**2+1.66931x+0.02466����Ҵ�-��������ƽ���ϵ',/)
    write(*,7)
7 format (/,'-----------�ܽ�-----------',/)
    write (*,8) wd,ww,inw,outw,Xd,Xw,inX,outX,xn,N,E,Eml
8 format (2x,'��������Һ���Wd:',f12.6,/,&
        2x,'��Һ���Ww:',f12.6,/,&
        2x,'�����n����Һ�����Wn-1:',f12.6,/,&
        2x,'�뿪��n����Һ�����Wn:',f12.6,/,&
        2x,'��������Һ���Xd:',f12.6,/,&
        2x,'��Һ���Xw:',f12.6,/,&
        2x,'�����n����Һ�����Xn-1:',f12.6,/,&
        2x,'�뿪��n����Һ�����Xn:',f12.6,/,&
        2x,'��yn(yn=Xn-1)����ƽ���Һ�����xn*:',f12.6,/,&
        2x,'����������N:',f12.6,/,&
        2x,'ȫ��Ч��E:',f12.6,/,&
        2x,'����Ч��Eml:',f12.6,/,&
        2x,'д������,W��x�Ļ��㹫ʽΪX=(30.0*W)/(23.0+7.0*W),��20�������ؿ��ṩ',/)
    write (*,10)
10 format (2x,'ͼ�ⷨ����:��origin�ж��Ҵ�-��������ƽ���ϵ��x��y���ж���ʽ���,',&
                '����ʹ��y=-0.72001x**2+1.66931x+0.02466;',&
                '����y=x,��XdΪ��㻭����,��xС��Xw��ֹ,������������������')
    read (*,*)

end program

subroutine equ(a,b,c,x)

    implicit none
    real,intent(in) :: a
    real,intent(in) :: b
    real,intent(in) :: c
    real,intent(out) :: x
    real :: x1,x2,delta
    delta=sqrt(b**2.0-4.0*a*c)
    x1=(-b-delta)/(2.0*a)
    x2=(-b+delta)/(2.0*a)
    if (x1 < x2) then
    x=x1
    else
    x=x2
    end if
    end subroutine equ

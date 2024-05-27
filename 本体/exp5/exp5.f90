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
1   format (/,1x,'本程序用以处理化工实验5: 板式精馏塔性能测定数据')
    write (*,2)
2   format(/,1x,'输入数据时，每个数据用空格或者英文逗号隔开')
    write (*,3)
3   format (/,1x,'By 2020级化学基地班步允霆,使用语言-Fortran')
    write (*,4)
4 format (/,'--------------------------------------------------------------------------------',/)
    write(*,*)'输入塔顶回流折光率、釜液折光率、进入第n块板的液体折光率、离开第n块板的液体折光率'
    read (*,*) Dd,Dw,inD,outD
    write (*,*)'输入实际塔板数'
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
110 format (/,'-----------理论塔板数-----------',/)
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
124 format (/,2x,'共',i2,'块理论塔板,输出xy即曲线坐标',/)
    write(*,5)
5 format (/,2x,'!!!!!注意!!!!!',/)
    write(*,6)
6 format (2x,'本程序采用y=-0.72001x**2+1.66931x+0.02466拟合乙醇-正丙醇相平衡关系',/)
    write(*,7)
7 format (/,'-----------总结-----------',/)
    write (*,8) wd,ww,inw,outw,Xd,Xw,inX,outX,xn,N,E,Eml
8 format (2x,'塔顶回流液组成Wd:',f12.6,/,&
        2x,'釜液组成Ww:',f12.6,/,&
        2x,'进入第n块板的液体组成Wn-1:',f12.6,/,&
        2x,'离开第n块板的液体组成Wn:',f12.6,/,&
        2x,'塔顶回流液组成Xd:',f12.6,/,&
        2x,'釜液组成Xw:',f12.6,/,&
        2x,'进入第n块板的液体组成Xn-1:',f12.6,/,&
        2x,'离开第n块板的液体组成Xn:',f12.6,/,&
        2x,'与yn(yn=Xn-1)成相平衡的液体组成xn*:',f12.6,/,&
        2x,'理论塔板数N:',f12.6,/,&
        2x,'全塔效率E:',f12.6,/,&
        2x,'单板效率Eml:',f12.6,/,&
        2x,'写报告用,W与x的换算公式为X=(30.0*W)/(23.0+7.0*W),由20化基孙守康提供',/)
    write (*,10)
10 format (2x,'图解法步骤:在origin中对乙醇-正丙醇相平衡关系中x、y进行多项式拟合,',&
                '个人使用y=-0.72001x**2+1.66931x+0.02466;',&
                '后做y=x,以Xd为起点画阶梯,当x小于Xw终止,阶梯数即理论塔板数')
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

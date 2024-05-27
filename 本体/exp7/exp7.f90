program exp7
    implicit none

    integer :: temp
    real :: M,row,roo,Vw,Vo
    real :: aa,bb,cc
    real :: XRb,XRt,YEb
    real :: VRb,VRt,VEb  
    real :: slope,intercept
    real :: cnaoh
    real,dimension(20) :: YE
    real,dimension(20) :: YE1
    real,dimension(20) :: XR
    real,dimension(20) :: mins
    integer,dimension(20) :: num
    integer :: i,ifc

    write (*,1)
1   format (/,1x,'本程序用以处理化工实验7:萃取塔性能测定')
    write (*,2)
2   format(/,1x,'输入数据时，每个数据用空格或者英文逗号隔开')
    write (*,3)
3   format (/,1x,'By 2020级化学基地班步允霆,使用语言-Fortran')
    write (*,4)
4 format (/,'--------------------------------------------------------------------------------',/)
    open (10,file='parameters.txt',status='old',action='read')
    read (10,*) M,row,roo,Vw,Vo
    write (*,*) '输入1采用15℃的平衡曲线,输入2采用25℃的平衡曲线'
    read (*,*) temp 
    if (temp == 1) then
    write (*,5) 
5 format (/,2x,'使用y=-245.54113x^2+1.08508x+2.22121E-5进行拟合,单位无10^4',/)
    aa=-245.54113
    bb=1.08508
    cc=2.22121/100000.0
    else if (temp == 2) then
    write (*,6) 
6 format (/,2x,'使用y=-261.125x^2+1.07758x+2.51212E-5进行拟合,单位无10^4',/)
    aa=-261.125
    bb=1.07758
    cc=2.51212/100000.0
end if
    write (*,*) '输入氢氧化钠浓度'
    read (*,*) cnaoh 
    write (*,*) '输入塔底轻相、塔顶轻相、塔底重相的氢氧化钠用量(mL)'
    read (*,*) VRb,VRt,VEb
    XRb=(VRb*cnaoh*M)/(Vo*roo)
    XRt=(VRt*cnaoh*M)/(Vo*roo)
    YEb=(VEb*cnaoh*M)/(Vw*row)
    slope=YEb/(XRb-XRt)
    intercept=(XRt*YEb)/(XRt-XRb)
    write (*,7) slope,intercept
7 format (/,2x,'操作线方程:y=',f8.5,'x',f8.5,/)
    write (*,8)
8 format (2x,'接下来输出绘制1/(YE*-YE)~YE关系所用的点,采用20组YE,数值记录在YE.txt中(关掉程序后可以手动修改数值,但总数要是20个)')
    open (10,file='YE.txt',status='old',action='read')
    read (10,*) YE 
do i=1,20
    XR(i)=(YE(i)-intercept)/slope
    YE1(i)=aa*((XR(i))**2.0)+bb*XR(i)+cc 
    mins(i)=1.0/(YE1(i)-YE(i))
    num(i)=i 
end do 
    write (*,10000) 
10000 format (/,/,2x,'NO.',3x,'      YE     ',3x,'  1/(YE*-YE) ',/,&
                  2x,'===',3x,'=============',3x,'=============') 
    write (*,2000) (num(i),YE(i),mins(i),i=1,20)
2000 format (2x,i3,3x,f13.5,3x,f13.5)
    write (*,11)
11 format (/,2x,'将以上数据拖入origin中作图,拟合后,选择Analysis→Mathematics→Integrate→Open Dialog进行积分,面积即NOE')
    write (*,12)
12 format (/,2x,'HOE与K特别好求,自己敲敲计算器吧(拟合建议最优拟合)')   
    write (*,1911)
1911 format (/,'----------------总结----------------',/)
    write (*,1144) XRb,XRt,YEb
    1144 format (2x,'塔底轻相入口浓度XRb',f12.8,/,&
                   2x,'塔顶轻相出口浓度XRt',f12.8,/,&
                   2x,'塔底重相出口浓度YEb',f12.8,/)
    close(10)
write (*,2223)
    2223 format (/,1x,'是否计算下一组数据?是输入1,否输入0')
    read (*,*) ifc 
    do
    if (ifc == 0) exit
    write (*,*) '输入塔底轻相、塔顶轻相、塔底重相的氢氧化钠用量(mL)'
    read (*,*) VRb,VRt,VEb
    XRb=(VRb*cnaoh*M)/(Vo*roo)
    XRt=(VRt*cnaoh*M)/(Vo*roo)
    YEb=(VEb*cnaoh*M)/(Vw*row)
    slope=YEb/(XRb-XRt)
    intercept=(XRt*YEb)/(XRt-XRb)
    write (*,70) slope,intercept
70 format (/,2x,'操作线方程:y=',f8.5,'x',f8.5,/)
    write (*,80)
80 format (2x,'接下来输出绘制1/(YE*-YE)~YE关系所用的点,采用20组YE,数值记录在YE.txt中(关掉程序后可以手动修改数值,但总数要是20个)')
    open (11,file='YE.txt',status='old',action='read')
    read (11,*) YE 
do i=1,20
    XR(i)=(YE(i)-intercept)/slope
    YE1(i)=aa*((XR(i))**2.0)+bb*XR(i)+cc 
    mins(i)=1.0/(YE1(i)-YE(i))
    num(i)=i 
end do 
    write (*,10001) 
10001 format (/,/,2x,'NO.',3x,'      YE     ',3x,'  1/(YE*-YE) ',/,&
                  2x,'===',3x,'=============',3x,'=============') 
    write (*,2001) (num(i),YE(i),mins(i),i=1,20)
2001 format (2x,i3,3x,f13.5,3x,f13.5)
    write (*,1912)
1912 format (/,'----------------总结----------------',/)
    write (*,1148) XRb,XRt,YEb
    1148 format (2x,'塔底轻相入口浓度XRb',f12.8,/,&
                   2x,'塔顶轻相出口浓度XRt',f12.8,/,&
                   2x,'塔底重相出口浓度YEb',f12.8,/)
    close(11)               
    write (*,2233)
    2233 format (/,1x,'是否计算下一组数据?是输入1,否输入0')
    read (*,*) ifc 
end do                   
read (*,*)
end program exp7

 
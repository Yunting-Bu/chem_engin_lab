program exp4

    implicit none
    real :: D,V
    real :: rof 
    real :: Z
    real :: ro1,ro2,roc,roa
    real :: pi = 3.14159
    real :: qvc,qva,qvcr,qvw
    real :: cbaoh2,chcl
    real :: Vd,Vt
    real :: ca1,ca2,cam !a1塔底,a2塔顶
    real :: H,E
    real :: row
    real :: Mw
    real :: y1,y2
    real :: VsL,LsL
    real :: atm
    real :: kla
    real :: per
    integer :: ifc
    real,allocatable,dimension(:) :: qv1
    real,allocatable,dimension(:) :: qv2
    real,allocatable,dimension(:) :: u
    real,allocatable,dimension(:) :: P0
    real,allocatable,dimension(:) :: PdZ
    integer :: t,i
    integer :: status
    integer,allocatable,dimension(:) :: num

    open (10,file='parameters.txt',status='old',action='read')
    read (10,*) row,atm,Mw,rof,roa,V
    write (*,1)
1   format (/,1x,'本程序用以处理化工实验4:二氧化碳吸收实验数据')
    write (*,2)
2   format(/,1x,'输入数据时，每个数据用空格或者英文逗号隔开(PS:用英文逗号方便复制到Excel中)')
    write (*,3)
3   format (/,1x,'By 2020级化学基地班步允霆,使用语言-Fortran')
    write (*,4)
4 format (/,'--------------------------------------------------------------------------------',/)
    write (*,*) '输入数据数量'
    read (*,*) t
    allocate (qv1(t),stat=status)
    allocate (qv2(t),stat=status)
    allocate (u(t),stat=status)
    allocate (P0(t),stat=status)
    allocate (PdZ(t),stat=status)
    allocate (num(t),stat=status)
     write (*,*) '输入每组数据的空气转子流量计读数(m3/h)'
     read (*,*) qv1
     write (*,*) '输入每组数据的填料层压强降(mmH2O)'
     read (*,*) P0
     write (*,*) '依次输入塔高Z(m)、塔径D(m)、出厂时标定转子流量计所用的流体密度(kg/m3)、实际测定流量时流体的密度(kg/m3)'
     read (*,*) Z,D,ro1,ro2
do i=1,t
    qv2(i)=sqrt((ro1*(rof-ro2))/(ro2*(rof-ro1)))*qv1(i)
    u(i)=(4.0*qv2(i))/(pi*(D**2.0)*3600.0)
    PdZ(i)=P0(i)/Z
    num(i)=i
end do
    write (*,1000) 
    1000 format (/,/,2x,'NO.',3x,'     P/Z     ',3x,'   u(m/s)  ',/,&
                     2x,'===',3x,'=============',3x,'=============') 
    write (*,200) (num(i),PdZ(i),u(i),i=1,t)
    200 format (2x,i3,3x,f13.5,3x,f13.5)
    write (*,222)
    222 format (/,1x,'是否计算下一组数据?是输入1,否输入0,并进入二氧化碳吸收传质系数计算')
    read (*,*) ifc
do
    if (ifc == 0) exit
     write (*,*) '输入数据数量'
    read (*,*) t
    allocate (qv1(t),stat=status)
    allocate (qv2(t),stat=status)
    allocate (u(t),stat=status)
    allocate (P0(t),stat=status)
    allocate (PdZ(t),stat=status)
    allocate (num(t),stat=status)
    write (*,*) '输入每组数据的空气转子流量计读数(m3/h)'
     read (*,*) qv1
     write (*,*) '输入每组数据的填料层压强降(mmH2O)'
     read (*,*) P0
     write (*,*) '依次输入塔高Z(m)、塔径D(m)、出厂时标定转子流量计所用的流体密度(kg/m3)、实际测定流量时流体的密度(kg/m3)'
     read (*,*) Z,D,ro1,ro2
do i=1,t
    qv2(i)=sqrt((ro1*(rof-ro2))/(ro2*(rof-ro1)))*qv1(i)
    u(i)=(4.0*qv2(i))/(pi*(D**2.0)*3600.0)
    PdZ(i)=P0(i)/Z
    num(i)=i
end do
    write (*,10000) 
    10000 format (/,/,2x,'NO.',3x,'     P/Z     ',3x,'   u(m/s)  ',/,&
                     2x,'===',3x,'=============',3x,'=============') 
    write (*,2000) (num(i),PdZ(i),u(i),i=1,t)
    2000 format (2x,i3,3x,f13.5,3x,f13.5)
    write (*,2222)
    2222 format (/,1x,'是否计算下一组数据?是输入1,否输入0,并进入二氧化碳吸收传质系数计算')
    read (*,*) ifc
end do
    write (*,00011)
    00011 format (/,'----------二氧化碳吸收传质系数测定----------',/)
    write (*,*) '依次输入CO2转子流量计读数(m3/h)、空气转子流量计读数(m3/h)、二氧化碳密度(g/L)'
    read (*,*) qvc,qva,roc
    write (*,*) '依次输入Ba(OH)2浓度(mol/L)、盐酸浓度(mol/L)、滴定塔底吸收液用盐酸的体积(mL)、滴定塔顶液用盐酸的体积(mL)'
    read (*,*) cbaoh2,chcl,Vd,Vt
    qvcr=sqrt((roa*(rof-roc))/(roc*(rof-roa)))*qvc
    ca1=(2.0*cbaoh2*V-chcl*Vd)/(2.0*V)
    ca2=(2.0*cbaoh2*V-chcl*Vt)/(2.0*V)
    write (*,*) '输入该温度下CO2在水中的亨利系数(忽略数量级,如输入1.66)'
    read (*,*) E 
    write (*,*) '输入水转子流量计读数(l/h)(如40.0,后面的.0不能省!)'
    read (*,*) qvw
    H=row/(Mw*(E*100000000.0) )
    y1=qvcr/(qva+qvcr)
    VsL=(qva*1000.0)/22.4
    LsL=(qvw*row)/Mw
    y2=y1-(LsL*(ca1-ca2)*Mw)/(VsL*row)
    cam=((H*y1*atm-ca1)-(H*y2*atm-ca2))/(log((H*y1*atm-ca1)/(H*y2*atm-ca2)))
    kla=((LsL/(Z*0.25*pi*(D**2.0)))*((ca1-ca2)/cam))*0.001
    per=((ca1-ca2)/ca1)*100.0
    write (*,1911)
1911 format (/,'----------总结----------',/)
    write (*,1144) qvcr,ca1,ca2,H,y1,y2,VsL,LsL,cam,kla,per
    1144 format (2x,'CO2实际流量=',f8.3,'m3/h',/,&
                   2x,'Ca1=',e12.3,'mol/L',/,&
                   2x,'Ca2=',e12.3,'mol/L',/,&
                   2x,'H=',e12.4,'kmol·m-3·Pa-1',/,&
                   2x,'y1=',f12.6,/,&
                   2x,'y2=',f12.6,/,&
                   2x,'VsL=',f8.4,'mol/h',/,&
                   2x,'LsL=',f12.4,'mol/h',/,&
                   2x,'deltaCam=',f12.5,'kmol/m3',/,&
                   2x,'KLa=',e18.5,'mol/m3·h',/,&
                   2x,'吸收率=',f12.2,'%')
    write (*,2223)
    2223 format (/,1x,'是否计算下一组数据?是输入1,否输入0')
    read (*,*) ifc 
    do
    if (ifc == 0) exit
    write (*,*) '依次输入CO2转子流量计读数(m3/h)、空气转子流量计读数(m3/h)、二氧化碳密度(g/L)'
    read (*,*) qvc,qva,roc
    write (*,*) '依次输入Ba(OH)2浓度(mol/L)、盐酸浓度(mol/L)、滴定塔底吸收液用盐酸的体积(mL)、滴定塔顶液用盐酸的体积(mL)'
    read (*,*) cbaoh2,chcl,Vd,Vt
    qvcr=sqrt((roa*(rof-roc))/(roc*(rof-roa)))*qvc
    ca1=(2.0*cbaoh2*V-chcl*Vd)/(2.0*V)
    ca2=(2.0*cbaoh2*V-chcl*Vt)/(2.0*V)
    write (*,*) '输入该温度下CO2在水中的亨利系数(忽略数量级,如输入1.66)'
    read (*,*) E 
    write (*,*) '输入水转子流量计读数(l/h)(如40.0,后面的.0不能省!)'
    read (*,*) qvw
    H=row/(Mw*(E*100000000.0) )
    y1=qvcr/(qva+qvcr)
    VsL=(qva*1000.0)/22.4
    LsL=(qvw*row)/Mw
    y2=y1-(LsL*(ca1-ca2)*Mw)/(VsL*row)
    cam=((H*y1*atm-ca1)-(H*y2*atm-ca2))/(log((H*y1*atm-ca1)/(H*y2*atm-ca2)))
    kla=((LsL/(Z*0.25*pi*(D**2.0)))*((ca1-ca2)/cam))*0.001
    per=((ca1-ca2)/ca1)*100.0
    write (*,1811)
1811 format (/,'----------总结----------',/)
    write (*,1444) qvcr,ca1,ca2,H,y1,y2,VsL,LsL,cam,kla,per
    1444 format (2x,'CO2实际流量=',f8.3,'m3/h',/,&
                   2x,'Ca1=',e12.3,'mol/L',/,&
                   2x,'Ca2=',e12.3,'mol/L',/,&
                   2x,'H=',e12.4,'kmol·m-3·Pa-1',/,&
                   2x,'y1=',f12.6,/,&
                   2x,'y2=',f12.6,/,&
                   2x,'VsL=',f8.4,'mol/h',/,&
                   2x,'LsL=',f12.4,'mol/h',/,&
                   2x,'deltaCam=',e12.5,'kmol/m3',/,&
                   2x,'KLa=',f18.5,'mol/m3·h',/,&
                   2x,'吸收率=',f8.2,'%')
    write (*,2233)
    2233 format (/,1x,'是否计算下一组数据?是输入1,否输入0')
    read (*,*) ifc 
    end do
    read (*,*)
end program exp4
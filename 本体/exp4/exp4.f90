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
    real :: ca1,ca2,cam !a1����,a2����
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
1   format (/,1x,'���������Դ�����ʵ��4:������̼����ʵ������')
    write (*,2)
2   format(/,1x,'��������ʱ��ÿ�������ÿո����Ӣ�Ķ��Ÿ���(PS:��Ӣ�Ķ��ŷ��㸴�Ƶ�Excel��)')
    write (*,3)
3   format (/,1x,'By 2020����ѧ���ذಽ����,ʹ������-Fortran')
    write (*,4)
4 format (/,'--------------------------------------------------------------------------------',/)
    write (*,*) '������������'
    read (*,*) t
    allocate (qv1(t),stat=status)
    allocate (qv2(t),stat=status)
    allocate (u(t),stat=status)
    allocate (P0(t),stat=status)
    allocate (PdZ(t),stat=status)
    allocate (num(t),stat=status)
     write (*,*) '����ÿ�����ݵĿ���ת�������ƶ���(m3/h)'
     read (*,*) qv1
     write (*,*) '����ÿ�����ݵ����ϲ�ѹǿ��(mmH2O)'
     read (*,*) P0
     write (*,*) '������������Z(m)������D(m)������ʱ�궨ת�����������õ������ܶ�(kg/m3)��ʵ�ʲⶨ����ʱ������ܶ�(kg/m3)'
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
    222 format (/,1x,'�Ƿ������һ������?������1,������0,�����������̼���մ���ϵ������')
    read (*,*) ifc
do
    if (ifc == 0) exit
     write (*,*) '������������'
    read (*,*) t
    allocate (qv1(t),stat=status)
    allocate (qv2(t),stat=status)
    allocate (u(t),stat=status)
    allocate (P0(t),stat=status)
    allocate (PdZ(t),stat=status)
    allocate (num(t),stat=status)
    write (*,*) '����ÿ�����ݵĿ���ת�������ƶ���(m3/h)'
     read (*,*) qv1
     write (*,*) '����ÿ�����ݵ����ϲ�ѹǿ��(mmH2O)'
     read (*,*) P0
     write (*,*) '������������Z(m)������D(m)������ʱ�궨ת�����������õ������ܶ�(kg/m3)��ʵ�ʲⶨ����ʱ������ܶ�(kg/m3)'
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
    2222 format (/,1x,'�Ƿ������һ������?������1,������0,�����������̼���մ���ϵ������')
    read (*,*) ifc
end do
    write (*,00011)
    00011 format (/,'----------������̼���մ���ϵ���ⶨ----------',/)
    write (*,*) '��������CO2ת�������ƶ���(m3/h)������ת�������ƶ���(m3/h)��������̼�ܶ�(g/L)'
    read (*,*) qvc,qva,roc
    write (*,*) '��������Ba(OH)2Ũ��(mol/L)������Ũ��(mol/L)���ζ���������Һ����������(mL)���ζ�����Һ����������(mL)'
    read (*,*) cbaoh2,chcl,Vd,Vt
    qvcr=sqrt((roa*(rof-roc))/(roc*(rof-roa)))*qvc
    ca1=(2.0*cbaoh2*V-chcl*Vd)/(2.0*V)
    ca2=(2.0*cbaoh2*V-chcl*Vt)/(2.0*V)
    write (*,*) '������¶���CO2��ˮ�еĺ���ϵ��(����������,������1.66)'
    read (*,*) E 
    write (*,*) '����ˮת�������ƶ���(l/h)(��40.0,�����.0����ʡ!)'
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
1911 format (/,'----------�ܽ�----------',/)
    write (*,1144) qvcr,ca1,ca2,H,y1,y2,VsL,LsL,cam,kla,per
    1144 format (2x,'CO2ʵ������=',f8.3,'m3/h',/,&
                   2x,'Ca1=',e12.3,'mol/L',/,&
                   2x,'Ca2=',e12.3,'mol/L',/,&
                   2x,'H=',e12.4,'kmol��m-3��Pa-1',/,&
                   2x,'y1=',f12.6,/,&
                   2x,'y2=',f12.6,/,&
                   2x,'VsL=',f8.4,'mol/h',/,&
                   2x,'LsL=',f12.4,'mol/h',/,&
                   2x,'deltaCam=',f12.5,'kmol/m3',/,&
                   2x,'KLa=',e18.5,'mol/m3��h',/,&
                   2x,'������=',f12.2,'%')
    write (*,2223)
    2223 format (/,1x,'�Ƿ������һ������?������1,������0')
    read (*,*) ifc 
    do
    if (ifc == 0) exit
    write (*,*) '��������CO2ת�������ƶ���(m3/h)������ת�������ƶ���(m3/h)��������̼�ܶ�(g/L)'
    read (*,*) qvc,qva,roc
    write (*,*) '��������Ba(OH)2Ũ��(mol/L)������Ũ��(mol/L)���ζ���������Һ����������(mL)���ζ�����Һ����������(mL)'
    read (*,*) cbaoh2,chcl,Vd,Vt
    qvcr=sqrt((roa*(rof-roc))/(roc*(rof-roa)))*qvc
    ca1=(2.0*cbaoh2*V-chcl*Vd)/(2.0*V)
    ca2=(2.0*cbaoh2*V-chcl*Vt)/(2.0*V)
    write (*,*) '������¶���CO2��ˮ�еĺ���ϵ��(����������,������1.66)'
    read (*,*) E 
    write (*,*) '����ˮת�������ƶ���(l/h)(��40.0,�����.0����ʡ!)'
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
1811 format (/,'----------�ܽ�----------',/)
    write (*,1444) qvcr,ca1,ca2,H,y1,y2,VsL,LsL,cam,kla,per
    1444 format (2x,'CO2ʵ������=',f8.3,'m3/h',/,&
                   2x,'Ca1=',e12.3,'mol/L',/,&
                   2x,'Ca2=',e12.3,'mol/L',/,&
                   2x,'H=',e12.4,'kmol��m-3��Pa-1',/,&
                   2x,'y1=',f12.6,/,&
                   2x,'y2=',f12.6,/,&
                   2x,'VsL=',f8.4,'mol/h',/,&
                   2x,'LsL=',f12.4,'mol/h',/,&
                   2x,'deltaCam=',e12.5,'kmol/m3',/,&
                   2x,'KLa=',f18.5,'mol/m3��h',/,&
                   2x,'������=',f8.2,'%')
    write (*,2233)
    2233 format (/,1x,'�Ƿ������һ������?������1,������0')
    read (*,*) ifc 
    end do
    read (*,*)
end program exp4
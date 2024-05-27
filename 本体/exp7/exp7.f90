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
1   format (/,1x,'���������Դ�����ʵ��7:��ȡ�����ܲⶨ')
    write (*,2)
2   format(/,1x,'��������ʱ��ÿ�������ÿո����Ӣ�Ķ��Ÿ���')
    write (*,3)
3   format (/,1x,'By 2020����ѧ���ذಽ����,ʹ������-Fortran')
    write (*,4)
4 format (/,'--------------------------------------------------------------------------------',/)
    open (10,file='parameters.txt',status='old',action='read')
    read (10,*) M,row,roo,Vw,Vo
    write (*,*) '����1����15���ƽ������,����2����25���ƽ������'
    read (*,*) temp 
    if (temp == 1) then
    write (*,5) 
5 format (/,2x,'ʹ��y=-245.54113x^2+1.08508x+2.22121E-5�������,��λ��10^4',/)
    aa=-245.54113
    bb=1.08508
    cc=2.22121/100000.0
    else if (temp == 2) then
    write (*,6) 
6 format (/,2x,'ʹ��y=-261.125x^2+1.07758x+2.51212E-5�������,��λ��10^4',/)
    aa=-261.125
    bb=1.07758
    cc=2.51212/100000.0
end if
    write (*,*) '������������Ũ��'
    read (*,*) cnaoh 
    write (*,*) '�����������ࡢ�������ࡢ���������������������(mL)'
    read (*,*) VRb,VRt,VEb
    XRb=(VRb*cnaoh*M)/(Vo*roo)
    XRt=(VRt*cnaoh*M)/(Vo*roo)
    YEb=(VEb*cnaoh*M)/(Vw*row)
    slope=YEb/(XRb-XRt)
    intercept=(XRt*YEb)/(XRt-XRb)
    write (*,7) slope,intercept
7 format (/,2x,'�����߷���:y=',f8.5,'x',f8.5,/)
    write (*,8)
8 format (2x,'�������������1/(YE*-YE)~YE��ϵ���õĵ�,����20��YE,��ֵ��¼��YE.txt��(�ص����������ֶ��޸���ֵ,������Ҫ��20��)')
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
11 format (/,2x,'��������������origin����ͼ,��Ϻ�,ѡ��Analysis��Mathematics��Integrate��Open Dialog���л���,�����NOE')
    write (*,12)
12 format (/,2x,'HOE��K�ر����,�Լ����ü�������(��Ͻ����������)')   
    write (*,1911)
1911 format (/,'----------------�ܽ�----------------',/)
    write (*,1144) XRb,XRt,YEb
    1144 format (2x,'�����������Ũ��XRb',f12.8,/,&
                   2x,'�����������Ũ��XRt',f12.8,/,&
                   2x,'�����������Ũ��YEb',f12.8,/)
    close(10)
write (*,2223)
    2223 format (/,1x,'�Ƿ������һ������?������1,������0')
    read (*,*) ifc 
    do
    if (ifc == 0) exit
    write (*,*) '�����������ࡢ�������ࡢ���������������������(mL)'
    read (*,*) VRb,VRt,VEb
    XRb=(VRb*cnaoh*M)/(Vo*roo)
    XRt=(VRt*cnaoh*M)/(Vo*roo)
    YEb=(VEb*cnaoh*M)/(Vw*row)
    slope=YEb/(XRb-XRt)
    intercept=(XRt*YEb)/(XRt-XRb)
    write (*,70) slope,intercept
70 format (/,2x,'�����߷���:y=',f8.5,'x',f8.5,/)
    write (*,80)
80 format (2x,'�������������1/(YE*-YE)~YE��ϵ���õĵ�,����20��YE,��ֵ��¼��YE.txt��(�ص����������ֶ��޸���ֵ,������Ҫ��20��)')
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
1912 format (/,'----------------�ܽ�----------------',/)
    write (*,1148) XRb,XRt,YEb
    1148 format (2x,'�����������Ũ��XRb',f12.8,/,&
                   2x,'�����������Ũ��XRt',f12.8,/,&
                   2x,'�����������Ũ��YEb',f12.8,/)
    close(11)               
    write (*,2233)
    2233 format (/,1x,'�Ƿ������һ������?������1,������0')
    read (*,*) ifc 
end do                   
read (*,*)
end program exp7

 
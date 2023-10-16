!     Last change:  NH   14 Apr 2015    8:58 pm
!***********************************************************************
!Development of a 2D Compressible flow solver using PVU-M Scheme
!Under the Guidance of Prof. Nadeem Hasan
!Group Members: Parvez Ahmad, Ayush Gupta, Aqib Khan
!This Program is yet to write formatted output to a file
!***********************************************************************
program StartPVUM
implicit none
real(kind=8),parameter::timeStep=2.0d-05,Gamma=1.4d0,Mach=0.7d0,Pr=0.70d0,kay=0.1d0,Re=100000.0d0
real(kind=8),allocatable,dimension(:)::      X,Y
real(kind=8),dimension(:,:),allocatable::    Rho,totEnrgy,Enrgy,Temp,Press,Mu,kappa,hT
real(kind=8),dimension(:,:,:),allocatable::  U,Vel,FcIntCel,GcIntCel,Fnc,Gnc,Q,U1
!Temporary Variables
integer::                                    i,j,k,m,nPoints,iter,stats,counter,nX,nY,wrfl,nYb
real(kind=8),dimension(10000)::              Xtemp,Ytemp
real(kind=8)::                               tmp1,tmp2,xMid_F,xMid_B,yMid_F,yMid_B,dx,sum_shear,sum_heat,Drag,Heat_Transfer
real(kind=8)::                               percent,time,Yb,p1,p2,p3,uu1,uu2,uu3,dudy1,dTdy1,shear,heat_flux
integer(kind=4)::                            kstart 
!***********************************************************************
call noOfPoints()
call makeGrid()
allocate(X(nX),Y(nY))
X=Xtemp(1:nX)
Y=Ytemp(1:nY)
allocate(U(nX,nY,4),U1(nX,nY,4),Fnc(nX,nY,4),Gnc(nX,nY,4))
allocate(FcIntCel(nX-1,2:nY-1,4),GcIntCel(2:nX-1,nY-1,4),Q(nX,nY,4))
allocate(Rho(nX,nY),Vel(nX,nY,2),totEnrgy(nX,nY))
allocate(Enrgy(nX,nY),Temp(nX,nY),Press(nX,nY),Mu(nX,nY),kappa(nX,nY),hT(nX,nY))

iter=201000
Yb=1.0
do j=1,nY
if(Y(j)>Yb)then
nYb=j
go to 10
endif
end do
10 write(*,*)nYb
!pause
write(*,*)"Start from initial data (0) or latest saved data (1)"
read*,kstart
if(kstart==0)then
time=0.0
counter=0
call initializeAll()
else
open(16,file="Restart.dat",status="unknown", form="unformatted") !RESULT

read(16)counter,time,Rho,Vel,Temp,press,enrgy,totenrgy,hT,mu,kappa

close(16)

endif   ! start from t=0 or previous saved data

!Main time iteration loop starts
do
if(iter==0) exit
do i=2,nX-1
do j=2,nYb-1
U(i,j,1)=Rho(i,j)
U(i,j,2)=Rho(i,j)*Vel(i,j,1)
U(i,j,3)=Rho(i,j)*Vel(i,j,2)
U(i,j,4)=Rho(i,j)*totEnrgy(i,j)
end do
end do
call findFlux()
!*****************Predictor Step starts*********************************
do i=2,nX-1
xMid_F=(X(i)+X(i+1))/2.0
xMid_B=(X(i)+X(i-1))/2.0
do j=2,nYb-1
yMid_F=(Y(j)+Y(j+1))/2.0
yMid_B=(Y(j)+Y(j-1))/2.0
do k=1,4
tmp1=(FcIntCel(i,j,k)-FcIntCel(i-1,j,k))/(xMid_F-xMid_B)+(Fnc(i+1,j,k)-Fnc(i,j,k))/(X(i+1)-X(i))
tmp2=(GcIntCel(i,j,k)-GcIntCel(i,j-1,k))/(yMid_F-yMid_B)+(Gnc(i,j+1,k)-Gnc(i,j,k))/(Y(j+1)-Y(j))
U1(i,j,k)=U(i,j,k)
U(i,j,k)=U(i,j,k)-timeStep*(tmp1+tmp2)
end do!loop k ends
end do!loop j ends
end do!loop i ends
!*****************Predictor Step ends***********************************
call findPrim()
call findFlux()
!*****************Corrector Step starts*********************************
do i=2,nX-1
xMid_F=(X(i)+X(i+1))/2.0
xMid_B=(X(i)+X(i-1))/2.0
do j=2,nYb-1
yMid_F=(Y(j)+Y(j+1))/2.0
yMid_B=(Y(j)+Y(j-1))/2.0
do k=1,4
tmp1=(FcIntCel(i,j,k)-FcIntCel(i-1,j,k))/(xMid_F-xMid_B)+(Fnc(i,j,k)-Fnc(i-1,j,k))/(X(i)-X(i-1))
tmp2=(GcIntCel(i,j,k)-GcIntCel(i,j-1,k))/(yMid_F-yMid_B)+(Gnc(i,j,k)-Gnc(i,j-1,k))/(Y(j)-Y(j-1))
U(i,j,k)=(U1(i,j,k)+U(i,j,k)-timeStep*(tmp1+tmp2))/2.0
end do!loop k ends
end do!loop j ends
end do!loop i ends
call findPrim()
call boundryCond()
!*****************Corrector Step ends***********************************
counter=counter+1
time=time+timeStep

if(mod(counter,100)==0) then        !PERCENTAGE

open(15,file="Result.txt",status="unknown") !RESULT
write(15,*)'zone'
write(15,*)'i=',nX,'j=',nYb
do j=1,nYb                          
do i=1,nX
write(15,*)x(i),y(j),Rho(i,j),Vel(i,j,1),Vel(i,j,2),Temp(i,j),Press(i,j)
end do
write(15,*) 
end do
close(15)

open(15,file="Restart.dat",status="unknown", form="unformatted") !RESULT
write(15)counter,time,Rho,Vel,Temp,press,enrgy,totenrgy,hT,mu,kappa
close(15)

! calculating  shear and wall heat flux

open(15,file="wallshear_heatflux.txt",status="unknown")
sum_shear=0.0
sum_heat=0.0
do i=1,nX

p1=Y(1)
p2=Y(2)
p3=Y(3)
uu1=vel(i,1,1)
uu2=vel(i,2,1)
uu3=vel(i,3,1)

dudy1=uu1*(1./(p1-p3)+1./(p1-p2))+uu2*(p1-p3)/((p2-p1)*(p2-p3))+uu3*(p1-p2)/((p3-p1)*(p3-p2))

uu1=temp(i,1)
uu2=temp(i,2)
uu3=temp(i,3)

dTdy1=uu1*(1./(p1-p3)+1./(p1-p2))+uu2*(p1-p3)/((p2-p1)*(p2-p3))+uu3*(p1-p2)/((p3-p1)*(p3-p2))

shear=mu(i,1)*dudy1
heat_flux=kappa(i,1)*dTdy1

dx=0.5*(x(i+1)-x(i-1))
if(i==1)dx=0.5*(x(2)-x(1))
if(i==nX)dx=0.5*(x(nX)-x(nX-1))

sum_shear=sum_shear+shear*dx
sum_heat=sum_heat+heat_flux*dx

write(15,*)x(i),shear,heat_flux

end do

close(15)

open(17,file="Drag_Heat_time.txt",status="unknown",access="append")

Drag=2.*sum_shear/Re
Heat_Transfer=2.*sum_heat/(Re*Pr*(gamma-1)*Mach*Mach)
write(17,*)time,Drag,Heat_Transfer
close(17)

endif

if(counter==iter) exit

WRITE(*,*)counter
end do
print*,'Time Loop exited'
!Main time iteration loop ends

contains
subroutine findFlux()
implicit none
real(kind=8)::vAvg,vC,vU,psi,phiC,phiU,eta,vIntCel,vIntCel1,vIntCel2,phiL,phiIntCel,Z(2)
real(kind=8)::Wf,vWf,phiWf,usn
real(kind=8)::p1,p2,p3,p4,u1,u2,u3,u4,t1,t2,t3,t4,s1,s2,s3,s4,s5,s6
real(kind=8)::dudx,dvdx,dTdx,dudy,dvdy,divVel,dTdy,Df,Dg,uu1,uu2,uu3,uu4,umax,umin
logical::logic1,logic2
!***********************************************************************
Q(:,:,1)=Rho(:,:)
Q(:,:,2)=Rho(:,:)*Vel(:,:,1)    
Q(:,:,3)=Rho(:,:)*Vel(:,:,2)
Q(:,:,4)=Rho(:,:)*hT(:,:)
!***************Calculation of FcIntCel starts**************************
do i=1,nX-1 
if(.not.(i==1 .or. i==nX-1)) then !weights and spacings
p1=X(i-1)
p2=X(i)
p3=X(i+1)
p4=X(i+2)
!weights for higher order
t1=(((p3-p2)**2)*(2*p4-p3-p2))/((p4-p1)*(p3-p1)*(p2-p1))
t2=((p2+p3-2*p1)*(2*p4-p3-p2))/((p4-p2)*(p2-p1))
t3=((p2+p3-2*p1)*(2*p4-p3-p2))/((p4-p3)*(p3-p1))
t4=((p2+p3-2*p1)*(p3-p2)**2)/((p4-p1)*(p4-p2)*(p4-p3))
!weights for lower order
s1=((p3-p2)**2)/((p2-p1)*(p3-p1))
s2=1+(p3-p1)/(p2-p1)
s3=1+(p2-p1)/(p3-p1)
s4=1+(p4-p3)/(p4-p2)
s5=1+(p4-p2)/(p4-p3)
s6=((p3-p2)**2)/((p4-p2)*(p4-p3))
endif
do j=2,nYb-1  
u2=Vel(i,j,1)
u3=Vel(i+1,j,1)  
vAvg=(u2+u3)/2.0
!Intercel velocity calculation
if(i==1 .or. i==nX-1) then  
vIntcel=(u2+u3)/2.0+Sig(vAvg)*(u2-u3)/2.0
else
u1=Vel(i-1,j,1)
u4=Vel(i+2,j,1)
vC=(-u1*t1+u2*t2+u3*t3-u4*t4)/8                   !cubic central interpolation
vU=(-u1*s1+u2*(s2+s4)+u3*(s3+s5)-u4*s6)/8+Sig(vAvg)*(-u1*s1+u2*(s2-s4)+u3*(s3-s5)+u4*s6)/8  !second order forward upwinding
usn=max(abs(u2),abs(u3))
Wf=usn/(usn+kay)  
vWf=vC+Wf*(vU-vC) 
if((u3-u2)/=0.0 .or. (u2-u1)/=0.0) then
Z(1)=(abs(u3-u2)-abs(u2-u1))/(abs(u3-u2)+abs(u2-u1)) 
else 
Z(1)=0.0
end if
if((u4-u3)/=0.0 .or. (u3-u2)/=0.0) then
Z(2)=(abs(u4-u3)-abs(u3-u2))/(abs(u4-u3)+abs(u3-u2))
else 
Z(2)=0.0
end if
psi=max(Z(1),Z(2));

vIntCel1=vWf+psi*(vU-vWf) 
vIntCel2=vU+psi*(vWf-vU)
!Range Boundedness Criteria
umax=u2
umin=u3
if(u3.gt.umax)then
umax=u3
umin=u2
endif

!if(umax.gt.0)umax=0.7*umax
!if(umax.lt.0)umax=1.3*umax
!if(umin.gt.0)umin=1.3*umin
!if(umin.lt.0)umin=0.7*umin

logic1=vIntCel1>=umin.and.vIntCel1<=umax
logic2=vIntCel2>=umin.and.vIntCel2<=umax

 ! if(logic1 .or. logic2) then
 !  if(logic1 .and. logic2) then
 !  vIntCel=(vIntCel1+vIntCel2)/2.0;
 !  elseif(logic1) then
 !  vIntCel=vIntCel1
 !  else
 !  vIntCel=vIntCel2
 !  end if
 ! else
 ! vIntCel=vAvg
 ! end if
  IF(logic1.or.logic2)then
     IF(logic1)vIntCel=vIntCel1
     IF(logic2)vIntCel=vIntCel2
     if(logic1.and.logic2)vIntCel=(vIntCel1+vIntCel2)/2.0
  else
     vIntCel=vAvg
  endif

end if
!Intercel property calculation
do m=1,4
u2=Q(i,j,m)
u3=Q(i+1,j,m)

uu2=U(i,j,m)
uu3=U(i+1,j,m)

phiL=(u2+u3)/2.0+Sig(vAvg)*(u2-u3)/2.0
if(i==1 .or. i==nX-1) then
phiIntcel=phiL
else 
u1=Q(i-1,j,m)
u4=Q(i+2,j,m)

uu1=U(i-1,j,m)
uu4=U(i+2,j,m)

phiC=(-u1*t1+u2*t2+u3*t3-u4*t4)/8.0  !cubic central interpolation
phiU=(-u1*s1+u2*(s2+s4)+u3*(s3+s5)-u4*s6)/8.0+Sig(vAvg)*(-u1*s1+u2*(s2-s4)+u3*(s3-s5)+u4*s6)/8.0  !second order forward upwinding  
phiWf=phiC+Wf*(phiU-phiC)
  if((uu3-uu2)/=0 .or. (uu2-uu1)/=0) then
  Z(1)=(abs(uu3-uu2)-abs(uu2-uu1))/(abs(uu3-uu2)+abs(uu2-uu1));
  else 
  Z(1)=0.0
  end if
  if((uu4-uu3)/=0 .or. (uu3-uu2)/=0) then
  Z(2)=(abs(uu4-uu3)-abs(uu3-uu2))/(abs(uu4-uu3)+abs(uu3-uu2));
  else 
  Z(2)=0.0
  end if
eta=max(Z(1),Z(2))

  
phiIntCel=phiWf+eta*(phiL-phiWf)
!Range Boundedness Criteria
logic1=phiIntCel>=min(u2,u3).and.phiIntCel<=max(u2,u3)
  if(.not.logic1) then
  phiIntCel=(u2+u3)/2.0
  endif
endif
FcIntCel(i,j,m)=vIntCel*phiIntCel
!FcIntCel(i,j,m)=vAvg*phiIntCel

end do!loop m ends
end do!loop j ends
end do!loop i ends
!***************Calculation of FcIntCel ends*************************
!***************Calculation of GcIntCel starts***********************
do j=1,nYb-1
if(.not.(j==1 .or. j==nYb-1)) then !weights and spacings
p1=Y(j-1)
p2=Y(j)
p3=Y(j+1)
p4=Y(j+2)
!weights for higher order
t1=(((p3-p2)**2)*(2*p4-p3-p2))/((p4-p1)*(p3-p1)*(p2-p1))  
t2=((p2+p3-2*p1)*(2*p4-p3-p2))/((p4-p2)*(p2-p1))
t3=((p2+p3-2*p1)*(2*p4-p3-p2))/((p4-p3)*(p3-p1))
t4=((p2+p3-2*p1)*(p3-p2)**2)/((p4-p1)*(p4-p2)*(p4-p3))
!weights for lower order
s1=((p3-p2)**2)/((p2-p1)*(p3-p1))
s2=1+(p3-p1)/(p2-p1)
s3=1+(p2-p1)/(p3-p1)
s4=1+(p4-p3)/(p4-p2)
s5=1+(p4-p2)/(p4-p3)
s6=((p3-p2)**2)/((p4-p2)*(p4-p3))
endif
do i=2,nX-1 
u2=Vel(i,j,2)
u3=Vel(i,j+1,2) 
vAvg=(u2+u3)/2.0 
!Intercell velocity calculation
if(j==1 .or. j==nYb-1) then  
vIntCel=(u2+u3)/2.0+Sig(vAvg)*(u2-u3)/2.0
else
u1=Vel(i,j-1,2)
u4=Vel(i,j+2,2)
vC=(-u1*t1+u2*t2+u3*t3-u4*t4)/8.0                     !cubic central interpolation
vU=(-u1*s1+u2*(s2+s4)+u3*(s3+s5)-u4*s6)/8.0+Sig(vAvg)*(-u1*s1+u2*(s2-s4)+u3*(s3-s5)+u4*s6)/8.0   !second order forward upwinding
usn=max(abs(u2),abs(u3))
Wf=usn/(usn+kay)

!if(Wf.lt.0.8)wf=0.8  

vWf=vC+Wf*(vU-vC) 
if((u3-u2)/=0 .or. (u2-u1)/=0) then
Z(1)=(abs(u3-u2)-abs(u2-u1))/(abs(u3-u2)+abs(u2-u1)) 
else 
Z(1)=0.0
end if
if((u4-u3)/=0 .or. (u3-u2)/=0) then
Z(2)=(abs(u4-u3)-abs(u3-u2))/(abs(u4-u3)+abs(u3-u2))
else 
Z(2)=0.0
end if
psi=max(Z(1),Z(2))

vIntCel1=vWf+psi*(vU-vWf) 
vIntCel2=vU+psi*(vWf-vU)

!Range Boundedness Criterion
umax=u2
umin=u3
if(u3.gt.umax)then
umax=u3
umin=u2
endif

if(umax.gt.0)umax=0.7*umax
if(umax.lt.0)umax=1.3*umax
if(umin.gt.0)umin=1.3*umin
if(umin.lt.0)umin=0.7*umin

logic1=vIntCel1>umin.and.vIntCel1<umax
logic2=vIntCel2>umin.and.vIntCel2<umax

!  if(logic1 .or. logic2) then
!   if(logic1 .and. logic2) then
!   vIntCel=(vIntCel1+vIntCel2)/2.0
!   elseif(logic1) then
!   vIntCel=vIntCel1
!   else
!   vIntCel=vIntCel2
!   end if
!  else
!  vIntCel=vAvg
!  endif

  IF(logic1.or.logic2)then
     IF(logic1)vIntCel=vIntCel1
     IF(logic2)vIntCel=vIntCel2
     if(logic1.and.logic2)vIntCel=(vIntCel1+vIntCel2)/2.0
  else
     vIntCel=vAvg
  endif


endif
!Intercell property calculation
do m=1,4
u2=Q(i,j,m)
u3=Q(i,j+1,m)

uu2=U(i,j,m)
uu3=U(i,j+1,m)

phiL=(u2+u3)/2.0+Sig(vAvg)*(u2-u3)/2.0
if(j==1 .or. j==nYb-1) then
phiIntCel=phiL
else 
u1=Q(i,j-1,m)
u4=Q(i,j+2,m)

uu1=U(i,j-1,m)
uu4=U(i,j+2,m)

phiC=(-u1*t1+u2*t2+u3*t3-u4*t4)/8                                                               !cubic central interpolation
phiU=(-u1*s1+u2*(s2+s4)+u3*(s3+s5)-u4*s6)/8+Sig(vAvg)*(-u1*s1+u2*(s2-s4)+u3*(s3-s5)+u4*s6)/8  !second order forward upwinding
phiWf=phiC+Wf*(phiU-phiC)
  if((uu3-uu2)/=0 .or. (uu2-uu1)/=0) then
  Z(1)=(abs(uu3-uu2)-abs(uu2-uu1))/(abs(uu3-uu2)+abs(uu2-uu1));
  else 
  Z(1)=0;
  end if
  if((uu4-uu3)/=0 .or. (uu3-uu2)/=0) then
  Z(2)=(abs(uu4-uu3)-abs(uu3-uu2))/(abs(uu4-uu3)+abs(uu3-uu2));
  else 
  Z(2)=0;
  end if
eta=max(Z(1),Z(2))
  
phiIntCel=phiWf+eta*(phiL-phiWf)
!Range Boundedness Criteria
logic1=phiIntCel>=min(u2,u3).and.phiIntCel<=max(u2,u3)
  if(.not.logic1) then
  phiIntCel=(u2+u3)/2.0
  endif
endif
!GcIntCel(i,j,m)=vIntCel*phiIntCel
GcIntCel(i,j,m)=vAvg*phiIntCel

end do!loop m ends
end do!loop i ends
end do!loop j ends
!****************Calculation of GcIntCel ends************************
!************Calculation of non convective flux starts***************
do i=1,nX
do j=1,nYb

if(i==1) then
dudx=(Vel(i+1,j,1)-Vel(i,j,1))/(X(i+1)-X(i))
dvdx=(Vel(i+1,j,2)-Vel(i,j,2))/(X(i+1)-X(i))
dTdx=(Temp(i+1,j)-Temp(i,j))/(X(i+1)-X(i))
endif

if(i==nX) then
dudx=(Vel(i,j,1)-Vel(i-1,j,1))/(X(i)-X(i-1))
dvdx=(Vel(i,j,2)-Vel(i-1,j,2))/(X(i)-X(i-1))
dTdx=(Temp(i,j)-Temp(i-1,j))/(X(i)-X(i-1))
endif

IF(i.gt.1.and.i.lt.nx)then
p1=X(i-1)
p2=X(i)
p3=X(i+1)
dudx=(p2-p1)/(p3-p1)*(Vel(i+1,j,1)-Vel(i,j,1))/(p3-p2)+(p3-p2)/(p3-p1)*(Vel(i,j,1)-Vel(i-1,j,1))/(p2-p1)
dvdx=(p2-p1)/(p3-p1)*(Vel(i+1,j,2)-Vel(i,j,2))/(p3-p2)+(p3-p2)/(p3-p1)*(Vel(i,j,2)-Vel(i-1,j,2))/(p2-p1)
dTdx=(p2-p1)/(p3-p1)*(Temp(i+1,j)-Temp(i,j))/(p3-p2)+(p3-p2)/(p3-p1)*(Temp(i,j)-Temp(i-1,j))/(p2-p1)
endif

if(j==1) then
dudy=(Vel(i,j+1,1)-Vel(i,j,1))/(Y(j+1)-Y(j))
dvdy=(Vel(i,j+1,2)-Vel(i,j,2))/(Y(j+1)-Y(j))
dTdy=(Temp(i,j+1)-Temp(i,j))/(Y(j+1)-Y(j))
endif

if(j==nYb) then
dudy=(Vel(i,j,1)-Vel(i,j-1,1))/(Y(j)-Y(j-1))
dvdy=(Vel(i,j,2)-Vel(i,j-1,2))/(Y(j)-Y(j-1))
dTdy=(Temp(i,j)-Temp(i,j-1))/(Y(j)-Y(j-1))
endif

IF(j.gt.1.and.j.lt.nyb)then
p1=Y(j-1)
p2=Y(j)
p3=Y(j+1)
dudy=(p2-p1)/(p3-p1)*(Vel(i,j+1,1)-Vel(i,j,1))/(p3-p2)+(p3-p2)/(p3-p1)*(Vel(i,j,1)-Vel(i,j-1,1))/(p2-p1)
dvdy=(p2-p1)/(p3-p1)*(Vel(i,j+1,2)-Vel(i,j,2))/(p3-p2)+(p3-p2)/(p3-p1)*(Vel(i,j,2)-Vel(i,j-1,2))/(p2-p1)
dTdy=(p2-p1)/(p3-p1)*(Temp(i,j+1)-Temp(i,j))/(p3-p2)+(p3-p2)/(p3-p1)*(Temp(i,j)-Temp(i,j-1))/(p2-p1)
endif

divVel=dudx+dvdy
Df=(2.0*Vel(i,j,1)/3.0)*(dvdy-2.0*dudx)-Vel(i,j,2)*(dvdx+dudy)
Dg=(2.0*Vel(i,j,2)/3.0)*(dudx-2.0*dvdy)-Vel(i,j,1)*(dvdx+dudy)

Fnc(i,j,1)=0.0
Fnc(i,j,2)=Press(i,j)/(Gamma*Mach*Mach)-2.0*Mu(i,j)*(dudx-divVel/3.0)/Re
Fnc(i,j,3)=-Mu(i,j)/Re*(dvdx+dudy)
Fnc(i,j,4)=-Gamma*kappa(i,j)*dTdx/(Re*Pr)+Gamma*(Gamma-1)*Mu(i,j)*Mach*Mach*Df/Re

Gnc(i,j,1)=0.0
Gnc(i,j,2)=-Mu(i,j)*(dvdx+dudy)/Re
Gnc(i,j,3)=Press(i,j)/(Gamma*Mach*Mach)-2.0*Mu(i,j)*(dvdy-divVel/3.0)/Re
Gnc(i,j,4)=-Gamma*kappa(i,j)*dTdy/(Re*Pr)+Gamma*(Gamma-1)*Mu(i,j)*Mach*Mach*Dg/Re
end do
end do
!***************Calculation of non convective flux ends*****************

end subroutine findFlux

subroutine initializeAll()
implicit none
!Non dimensional initialization
do i=1,nX
do j=1,nYb
if(j==1)then
Vel(i,j,1)=0.0d0
else
Vel(i,j,1)=1.0d0
endif
Vel(i,j,2)=0.0d0
Rho(i,j)=1.0d0
Temp(i,j)=1.0d0
Mu(i,j)=(1.0/1.0761)*(Temp(i,j)/273.0*300.0)**1.5*(273.0+110.0)/(Temp(i,j)*300+110.0)
kappa(i,j)=Mu(i,j)
Press(i,j)=Rho(i,j)*Temp(i,j)
Enrgy(i,j)=Temp(i,j)
totEnrgy(i,j)=Temp(i,j)+0.5*Gamma*(Gamma-1)*Mach**2*(Vel(i,j,1)**2+Vel(i,j,2)**2)
hT(i,j)=totEnrgy(i,j)+(Gamma-1)*Press(i,j)/Rho(i,j)
end do
end do
end subroutine initializeAll

subroutine boundryCond()
implicit none
real(kind=8)::dpdy1,dpdy2,dpdy3,R_pos,R_neg,t2,t3,Entpy,p1,p2,p3,dtxy,rru,rrv,qq,divvel,dtyy,dpdy
real(kind=8),dimension(2)::dvdy,dudx,tyy
real(kind=8),dimension(3)::dudy,txy
integer::i,j
!Wall BC starts
do i=2,nX-1
dvdy(1)=(Vel(i,2,2)-Vel(i,1,2))/(Y(2)-Y(1))
p1=Y(1)
p2=Y(2)
p3=Y(3)
dvdy(2)=(p2-p1)/(p3-p1)*(Vel(i,3,2)-Vel(i,2,2))/(p3-p2)+(p3-p2)/(p3-p1)*(Vel(i,2,2)-Vel(i,1,2))/(p2-p1)

!tyy(1)=(4.0/3.0)*Mu(i,1)*dvdy(1)
!tyy(2)=(4.0/3.0)*Mu(i,2)*dvdy(2)

dudy(1)=(Vel(i-1,2,1)-Vel(i-1,1,1))/(p2-p1)
dudy(2)=(Vel(i,2,1)-Vel(i,1,1))/(p2-p1)
dudy(3)=(Vel(i+1,2,1)-Vel(i+1,1,1))/(p2-p1)

txy(1)=dudy(1)*Mu(i-1,1)
txy(2)=dudy(2)*Mu(i,1)
txy(3)=dudy(3)*Mu(i+1,1)


p1=X(i-1)
p2=X(i)
p3=X(i+1)
dudx(2)=(p2-p1)/(p3-p1)*(Vel(i+1,2,1)-Vel(i,2,1))/(p3-p2)+(p3-p2)/(p3-p1)*(Vel(i,2,1)-Vel(i-1,2,1))/(p2-p1)
dudx(1)=0.0  !dudx at wall
divvel=dvdy(1)+dudx(1)
tyy(1)=2*Mu(i,1)*(dvdy(1)-divvel/3.0)
divvel=dvdy(2)+dudx(2)
tyy(2)=2*Mu(i,2)*(dvdy(2)-divvel/3.0)
dtyy=(tyy(2)-tyy(1))/(Y(2)-Y(1))

dtxy=(p2-p1)/(p3-p1)*(txy(3)-txy(2))/(p3-p2)+(p3-p2)/(p3-p1)*(txy(2)-txy(1))/(p2-p1)

dpdy=((dtxy+dtyy)/Re-Rho(i,2)*(vel(i,2,2)**2)/(Y(2)-Y(1)))*Gamma*Mach*Mach

!dpdy=(press(i,3)-press(i,2))/(Y(3)-Y(2))
p1=Y(2)
p2=Y(3)
p3=Y(4)
dpdy3=(p2-p1)/(p3-p1)*(Press(i,4)-Press(i,3))/(p3-p2)+(p3-p2)/(p3-p1)*(Press(i,3)-Press(i,2))/(p2-p1)
dpdy2=press(i,2)*(1./(p1-p3)+1./(p1-p2))+press(i,3)*(p1-p3)/((p2-p1)*(p2-p3))+press(i,4)*(p1-p2)/((p3-p1)*(p3-p2))
dpdy1=dpdy2-(dpdy3-dpdy2)*p1/(p2-p1)

Press(i,1)=Press(i,2)-dpdy*(Y(2)-Y(1))

!Press(i,1)=press(i,2)
Rho(i,1)=Press(i,1)/Temp(i,1)
end do

!Wall BC ends
Press(1,1)=Press(2,1)
Press(nX,1)=Press(nX-1,1)
Rho(1,1)=Rho(2,1)
Rho(nX,1)=Rho(nX-1,1)
do i=1,nX
Mu(i,1)=(1.0/1.0761)*((Temp(i,1)/273.0*300.0)**1.5)*(273.0+110.0)/(Temp(i,1)*300+110.0)
kappa(i,1)=Mu(i,1)
Enrgy(i,1)=Temp(i,1)
totEnrgy(i,1)=Temp(i,1)+0.5*Gamma*(Gamma-1)*Mach*Mach*(Vel(i,1,1)**2+Vel(i,1,2)**2)
hT(i,1)=totEnrgy(i,1)+(Gamma-1)*Press(i,1)/Rho(i,1)
end do


!Left Boundary BC starts
do j=2,nYb-1 !Left Boundary
if(Vel(1,j,1)>0.0) then !inflow
  R_neg=-1-2/(Mach*(Gamma-1))
  if(abs(Vel(1,j,1)*Mach/sqrt(Temp(1,j)))>1.0) then
  R_pos=-1+2/(Mach*(gamma-1))
  else
  t2=-Vel(2,j,1)+2.0*sqrt(Temp(2,j))/(Mach*(gamma-1.0))
  t3=-Vel(3,j,1)+2.0*sqrt(Temp(3,j))/(Mach*(gamma-1.0))  
  R_pos=t2+(t2-t3)/(X(2)-X(3))*(X(1)-X(2))
  end if
  Vel(1,j,2)=0.0
  Vel(1,j,1)=-(R_pos+R_neg)/2.0
  !Vel(1,j,1)=1.0d0
  Temp(1,j)=((Gamma-1)*Mach*(R_pos-R_neg)/4.0)**2
  Rho(1,j)=Temp(1,j)**(1.0/(Gamma-1.0))
!  Rho(1,j)=1.0
else 
  t2=-Vel(2,j,1)+2.0*sqrt(Temp(2,j))/(Mach*(gamma-1.0))
  t3=-Vel(3,j,1)+2.0*sqrt(Temp(3,j))/(Mach*(gamma-1.0))
  R_pos=t2+(t2-t3)/(X(2)-X(3))*(X(1)-X(2))
  if(abs(Vel(1,j,1)*Mach/sqrt(Temp(1,j)))>1.0) then
  t2=-Vel(2,j,1)-2.0*sqrt(Temp(2,j))/(Mach*(gamma-1.0))
  t3=-Vel(3,j,1)-2.0*sqrt(Temp(3,j))/(Mach*(gamma-1.0))
  R_neg=t3-(t3-t2)/(X(3)-X(2))*(X(3)-X(1))
  else
  R_neg=-1-2/(Mach*(gamma-1))
  end if
  Vel(1,j,2)=Vel(3,j,2)-(Vel(3,j,2)-Vel(2,j,2))/(X(3)-X(2))*(X(3)-X(1))
  !Vel(1,j,2)=0.0d0
  Vel(1,j,1)=-(R_pos+R_neg)/2.0
  !Vel(1,j,1)=1.0d0
  Temp(1,j)=((Gamma-1)*Mach*(R_pos-R_neg)/4.0)**2
  t2=log(Temp(2,j))/(Gamma-1)-log(Rho(2,j))
  t3=log(Temp(3,j))/(Gamma-1)-log(Rho(3,j))
  Entpy=t3-(t3-t2)/(X(3)-X(2))*(X(3)-X(1))
  Rho(1,j)=Temp(1,j)**(1/(Gamma-1))/exp(Entpy)
  !Rho(1,j)=1.0d0
end if
end do
!patch for 1,nY
Vel(1,nYb,1)=Vel(1,nYb-1,1)
Vel(1,nYb,2)=Vel(1,nYb-1,2)
Temp(1,nYb)=Temp(1,nYb-1)
Rho(1,nYb)=Rho(1,nYb-1)

do j=1,nYb
Press(1,j)=Temp(1,j)*Rho(1,j)
Mu(1,j)=(1.0/1.0761)*((Temp(1,j)/273.0*300.0)**1.5)*(273.0+110.0)/(Temp(1,j)*300+110.0)
kappa(1,j)=Mu(1,j)
Enrgy(1,j)=Temp(1,j)
totEnrgy(1,j)=Temp(1,j)+0.5*Gamma*(Gamma-1)*Mach*Mach*(Vel(1,j,1)**2+Vel(1,j,2)**2)
hT(1,j)=totEnrgy(1,j)+(Gamma-1)*Press(1,j)/Rho(1,j)
end do
!Left Boundary BC ends

!Upper Boundary BC starts
do i=2,nX-1  !Upper Boundary
if(Vel(i,nYb,2)<=0.0) then !inflow
  R_neg=0.0-2/(Mach*(gamma-1.0)) !Vn-c
  if(abs(Vel(i,nYb,2)*Mach/sqrt(Temp(i,nYb)))>1.0) then!changed nh
  R_pos=0.0+2.0/(Mach*(gamma-1.0))
  else
  t2=Vel(i,nYb-1,2)+2.0*sqrt(Temp(i,nYb-1))/(Mach*(gamma-1.0))
  t3=Vel(i,nYb-2,2)+2.0*sqrt(Temp(i,nYb-2))/(Mach*(gamma-1.0))
  R_pos=t3-(t3-t2)/(Y(nYb-2)-Y(nYb-1))*(Y(nYb-2)-Y(nYb))
  end if
  Vel(i,nYb,1)=1.0d0 
  Vel(i,nYb,2)=(R_pos+R_neg)/2.0
  Temp(i,nYb)=((Gamma-1.0)*Mach*(R_pos-R_neg)/4.0)**2.0
  Rho(i,nYb)=1.0d0
else
  t2=Vel(i,nYb-1,2)+2.0*sqrt(Temp(i,nYb-1))/(Mach*(gamma-1.0))
  t3=Vel(i,nYb-2,2)+2.0*sqrt(Temp(i,nYb-2))/(Mach*(gamma-1.0))
  R_pos=t3-(t3-t2)/(Y(nYb-2)-Y(nYb-1))*(Y(nYb-2)-Y(nYb))

  if(abs(Vel(1,nYb,2)*Mach/sqrt(Temp(1,nYb)))>1.0) then!changed NH
  t2=Vel(i,nYb-1,2)-2.0*sqrt(Temp(i,nYb-1))/(Mach*(gamma-1.0))
  t3=Vel(i,nYb-2,2)-2.0*sqrt(Temp(i,nYb-2))/(Mach*(gamma-1.0))
  R_neg=t3-(t3-t2)/(Y(nYb-2)-Y(nYb-1))*(Y(nYb-2)-Y(nYb))
  else
  R_neg=0.0-2.0/(Mach*(gamma-1.0))
  end if
  Vel(i,nYb,1)=Vel(i,nYb-2,1)-(Vel(i,nYb-2,1)-Vel(i,nYb-1,1))/(Y(nYb-2)-Y(nYb-1))*(Y(nYb-2)-Y(nYb))  
  Vel(i,nYb,2)=(R_pos+R_neg)/2.0
  Temp(i,nYb)=((Gamma-1.0)*Mach*(R_pos-R_neg)/4.0)**2.0
  t2=log(Temp(i,nYb-1))/(Gamma-1)-log(Rho(i,nYb-1))
  t3=log(Temp(i,nYb-2))/(Gamma-1)-log(Rho(i,nYb-2))
  Entpy=t3-(t3-t2)/(Y(nYb-2)-Y(nYb-1))*(Y(nYb-2)-Y(nYb)) 
  Rho(i,nYb)=Temp(i,nYb)**(1.0/(Gamma-1.0))/exp(Entpy)
end if
end do
!patch for nX,nY
Vel(nX,nYb,1)=Vel(nX-1,nYb,1)
Vel(nX,nYb,2)=Vel(nX-1,nYb,2)
Temp(nX,nYb)=Temp(nX-1,nYb)
Rho(nX,nYb)=Rho(nX-1,nYb)
do i=1,nX
Press(i,nYb)=Temp(i,nYb)*Rho(i,nYb)
Mu(i,nYb)=(1.0/1.0761)*((Temp(i,nYb)/273.0*300.0)**1.5)*(273.0+110.0)/(Temp(i,nYb)*300+110.0)
kappa(i,nYb)=Mu(i,nYb)
Enrgy(i,nYb)=Temp(i,nYb)
totEnrgy(i,nYb)=Temp(i,nYb)+0.5*Gamma*(Gamma-1)*Mach*Mach*(Vel(i,nYb,1)**2+Vel(i,nYb,2)**2)
hT(i,nYb)=totEnrgy(i,nYb)+(Gamma-1)*Press(i,nYb)/Rho(i,nYb)
end do
!Upper Boundary BC ends

!Right Boundary BC starts
do j=2,nYb-1 !Right Boundary
if(Vel(nX,j,1)<0.0) then !inflow
  R_neg=1.0-2.0/(Mach*(gamma-1.0)) !Vn-c
  if(abs(Vel(nX,j,1)*Mach/sqrt(Temp(nX,j)))>1.0) then
  R_pos=1.0+2.0/(Mach*(gamma-1.0))
  else
  t2=Vel(nX-1,j,1)+2.0*sqrt(Temp(nX-1,j))/(Mach*(gamma-1.0))
  t3=Vel(nX-2,j,1)+2.0*sqrt(Temp(nX-2,j))/(Mach*(gamma-1.0))
  R_pos=t3-(t3-t2)/(X(nX-2)-X(nX-1))*(X(nX-2)-X(nX))
  end if
  !Vel(nX,j,2)=0.0
  Vel(nX,j,2)=Vel(nX-1,j,2)
  !Vel(nX,j,1)=(R_pos+R_neg)/2.0
  rru=U(nX-2,j,2)+((U(nX-1,j,2)-U(nX-2,j,2))/(X(nX-1)-X(nX-2)))*(X(nX)-X(nX-2))
  
  Temp(nX,j)=((Gamma-1.0)*Mach*(R_pos-R_neg)/4.0)**2.0
  !Rho(nX,j)=1.0d0
  Rho(nX,j)=Rho(nX-1,j)
  Vel(nX,j,1)=rru/Rho(nX,j)
else      !outflow
  t2=Vel(nX-1,j,1)+2.0*sqrt(Temp(nX-1,j))/(Mach*(gamma-1.0))
  t3=Vel(nX-2,j,1)+2.0*sqrt(Temp(nX-2,j))/(Mach*(gamma-1.0))
  R_pos=t3-(t3-t2)/(X(nX-2)-X(nX-1))*(X(nX-2)-X(nX))

  if(abs(Vel(nX,j,1)*Mach/sqrt(Temp(nX,j)))>1.0) then
  t2=Vel(nX-1,j,1)-2.0*sqrt(Temp(nX-1,j))/(Mach*(gamma-1.0))
  t3=Vel(nX-2,j,1)-2.0*sqrt(Temp(nX-2,j))/(Mach*(gamma-1.0))
  R_neg=t3-(t3-t2)/(X(nX-2)-X(nX-1))*(X(nX-2)-X(nX))
  else
  R_neg=1-2.0/(Mach*(gamma-1.0))
  end if

! Vel(nX,j,2)=Vel(nX-2,j,2)+((Vel(nX-1,j,2)-Vel(nX-2,j,2))/(X(nX-1)-X(nX-2)))*(X(nX)-X(nX-2))
 rrv=U(nX-2,j,3)+((U(nX-1,j,3)-U(nX-2,j,3))/(X(nX-1)-X(nX-2)))*(X(nX)-X(nX-2))
 !Vel(nX,j,2)=Vel(nX-1,j,2)
 !Vel(nX,j,1)=(R_pos+R_neg)/2.0
 rru=U(nX-2,j,2)+((U(nX-1,j,2)-U(nX-2,j,2))/(X(nX-1)-X(nX-2)))*(X(nX)-X(nX-2))
 !Temp(nX,j)=((Gamma-1.0)*Mach*(R_pos-R_neg)/4.0)**2.0
 !qq=U(nX-2,j,4)+((U(nX-1,j,4)-U(nX-2,j,4))/(X(nX-1)-X(nX-2)))*(X(nX)-X(nX-2))  
 Temp(nX,j)=Temp(nX-2,j)+((Temp(nX-1,j)-Temp(nX-2,j))/(X(nX-1)-X(nX-2)))*(X(nX)-X(nX-2))

 t2=log(Temp(nX-1,j))/(Gamma-1.0)-log(Rho(nX-1,j))
 t3=log(Temp(nX-2,j))/(Gamma-1.0)-log(Rho(nX-2,j))
 Entpy=t3-(t3-t2)/(X(nX-2)-X(nX-1))*(X(nX-2)-X(nX))
 Rho(nX,j)=Temp(nX,j)**(1/(Gamma-1.0))/exp(Entpy)
 !Rho(nX,j)=Rho(nX-1,j)
 Vel(nX,j,1)=rru/Rho(nX,j)
 Vel(nX,j,2)=rrv/Rho(nX,j)
  
end if

end do
do j=1,nYb
Press(nX,j)=Temp(nX,j)*Rho(nX,j)
Mu(nX,j)=(1.0/1.0761)*((Temp(nX,j)/273.0*300.0)**1.5)*(273.0+110.0)/(Temp(nX,j)*300+110.0)
kappa(nX,j)=Mu(nX,j)
Enrgy(nX,j)=Temp(nX,j)
totEnrgy(nX,j)=Temp(nX,j)+0.5*Gamma*(Gamma-1)*Mach*Mach*(Vel(nX,j,1)**2+Vel(nX,j,2)**2)
hT(nX,j)=totEnrgy(nX,j)+(Gamma-1)*Press(nX,j)/Rho(nX,j)
end do
!Right Boundary BC ends
end subroutine boundryCond

subroutine noOfPoints()
implicit none
open (10,file="fileGrid.txt",status="old")
counter=0
stats=0
do
read(10,*,iostat=stats)
if(stats>0) then
!print*,'Problem reading file'
exit
elseif(stats<0) then
!print*,'Successful End of file reached'
exit
else
counter=counter+1
end if
end do
close(10)
nPoints=counter+1
!print*,nPoints,'grid points successfully read'
end subroutine noOfPoints

subroutine makeGrid()
implicit none
real(kind=8)::temp(0:nPoints,2)
temp(0,1)=1234
open (10,file="fileGrid.txt",status="old")
nX=0;
nY=0;
do i=1,nPoints
read(10,*,iostat=stats) temp(i,1),temp(i,2)
if(stats/=0) exit

if(temp(1,1)==temp(i,1)) then
nY=nY+1
Ytemp(nY)=temp(i,2)
end if
if(temp(i,1)/=temp(i-1,1)) then
nX=nX+1
Xtemp(nX)=temp(i,1)
end if
end do

print*,'The file contains ',nX,' points in X direction'
print*,'The file contains ',nY,' points in Y direction'
close(10)
end subroutine makeGrid

subroutine findPrim()
implicit none
do i=2,nX-1
do j=2,nYb-1
Rho(i,j)=U(i,j,1)
Vel(i,j,1)=U(i,j,2)/Rho(i,j)
Vel(i,j,2)=U(i,j,3)/Rho(i,j)
totEnrgy(i,j)=U(i,j,4)/Rho(i,j)
Enrgy(i,j)=totEnrgy(i,j)-0.5*Gamma*(Gamma-1)*Mach*Mach*(Vel(i,j,1)**2+Vel(i,j,2)**2) 
Temp(i,j)=Enrgy(i,j)
Press(i,j)=Rho(i,j)*Temp(i,j)
Mu(i,j)=(1.0/1.0761)*((Temp(i,j)/273.0*300.0)**1.5)*(273.0+110.0)/(Temp(i,j)*300+110.0)
kappa(i,j)=Mu(i,j)
hT(i,j)=totEnrgy(i,j)+(Gamma-1)*Press(i,j)/Rho(i,j)
end do
end do
end subroutine findPrim

function Sig(input1)
Implicit none
real(kind=8)::input1
integer::Sig
if(input1>0) then
   Sig=1
elseif(input1<0) then
   Sig=-1
 else
   Sig=0
endif
end function Sig

end program

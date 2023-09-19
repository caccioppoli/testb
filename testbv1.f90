parameter (nsa=100,iter=500000,ndi=15000)
!parameter (iter=500000,ndi=15000)
dimension q(700000),qsa(500),nb(ndi),nsb(ndi),nra(ndi),b(iter),sb(iter),ra(iter)
character*20 nome1, nome2, nome3

open(60,file='set.dat')
open(70,file='sc_conv.dat')





qc=2.2

kk=0
do i=1,700000 
  read(70,*,end=99)t,x,y,z,qq 
  if(qq<qc)cycle
  kk=kk+1
  q(kk)=qq
end do
99 continue
print*,kk

rath=.5
!nsa=0
do ira=1,4
  do i=1,ndi
    nb(i)=0
    nsb(i)=0
    nra(i)=0
  end do
  rath=rath+.5
!  nsa=nsa+50
  read(60,*)nome1
  read(60,*)nome2
  read(60,*)nome3
  open(75,file=nome1)
  open(80,file=nome2)
  open(85,file=nome3)
  bsum=0.
  sdsum=0.
  rasum=0.
  ne=0
  do nrand=1,iter
    qmin=10.
    qmax=-10.
    do i=1,nsa
      l=int(rand()*kk)
      qsa(i)=q(l)
      if(qsa(i)<qmin)qmin=qsa(i)
      if(qsa(i)>qmax)qmax=qsa(i)  
    end do
    ra(nrand)=qmax-qmin
    if(ra(nrand)<rath)then
      ne=ne+1
      cycle
    end if
    rasum=rasum+ra(nrand)
    call bvalue(qsa,nsa,bs,sbs)
    bsum=bsum+bs
    sdsum=sdsum+sbs
    mb=int(bs*100.)
    msd=int(sbs*300)
    mra=int(ra(nrand)*3.)
    nra(mra)=nra(mra)+1
    nb(mb)=nb(mb)+1
    nsb(msd)=nsb(msd)+1
    b(nrand)=bs
    sb(nrand)=sbs
!    write(90,*)b(nrand),sb(nrand),ra(nrand)
  end do
  bm=bsum/float(iter-ne)
  sbm=sdsum/float(iter-ne)
  ram=rasum/float(iter-ne)

!!$   do i=1,iter
!!$    write(75,*)b(i)
!!$    write(80,*)sb(i)
!!$    write(85,*)ra(i)
!!$  end do  
!!$  
  do i=1,ndi
    if(nb(i)>0)write(75,*)i/100.,float(nb(i))/float(iter-ne*0)
    if(nsb(i)>0)write(80,*)i/300.,float(nsb(i))/float(iter-ne)
    if(nra(i)>0)write(85,*)i/3.,float(nra(i))/float(iter-ne)
  end do  
  close(75)
  close(80)
  close(85)
  
  sdbsum=0.
  sdsdsum=0.
  sdrasum=0.
  do i=1,iter
    if(ra(i)<rath)cycle
    sdbsum=sdbsum+(bm-b(i))**2.
    sdsdsum=sdsdsum+(sbm-sb(i))**2.
    sdrasum=sdrasum+(ram-ra(i))**2.   
  end do
  sdb=sqrt(sdbsum/float(iter-1-ne))
  sdsd=sqrt(sdsdsum/float(iter-1-ne))
  sdra=sqrt(sdrasum/float(iter-1-ne))

    ! Calcola la skewness
  skebsum=0. !sk
  skesbsum=0. !sk
  skerasum=0. !sk
  do i=1,iter
     if(ra(i)<rath)cycle
     skebsum=skebsum+((b(i)-bm)/sdb)**3.
     skesbsum=skesbsum+((sb(i)-sbm)/sdsd)**3.
     skerasum=skerasum+((ra(i)-ram)/sdra)**3.
     !write(21,*),skebsum
  end do
  
  nsk=float(iter-1-ne)
  skbsum = skebsum*1./nsk
  sksbsum = skesbsum*1./nsk
  skrasum = skerasum*1./nsk
  
  print*,rath,bm,sdb,skbsum,sbm,sdsd,sksbsum,ram,sdra,skrasum,ne
!  print*,'bvalue',bm,sdb
!  print*,'berr',sbm,sdsd
!  print*,'range',ram,sdra
!  print*,'esclusi',ne
end do

end

subroutine bvalue(q,nsa,bs,sbs)
dimension q(nsa)
qc=2.2
qs=0.
k=0
do i=1,nsa
  qs=qs+q(i)
end do
qm=qs/float(nsa)
bs=(1./log(10.))*(1./(qm-(qc-.05)))
!print*,qs,qm,bs,nsa
qqs=0.
k=0
do i=1,nsa
  qqs=qqs+(q(i)-qm)**2.
  k=k+1
end do
de=float(k)*float(k-1)
sbs=2.3*bs*bs*sqrt(qqs/de)

return    
end

subroutine cv(q,nsa,qc)
dimension q(nsa),q1(nsa)

qc=0.
do iqc=1,35
  qc=qc+.1
  if(qc>4.5)exit
  qs=0.
  k=0
  do i=1,nsa
    if(q(i)<qc)cycle
    k=k+1
    q1(k)=q(i)-qc
    qs=qs+q(i)-qc
  end do
  qm=qs/float(k)
  qqs=0.
  do i=1,k
    qqs=qqs+(q1(i)-qm)**2.
  end do
  sd1=sqrt(qqs/float(k-1))
  cvs=sd1/qm
!  print*,qc,cv
  if(cvs>=.93)exit
end do

return
end



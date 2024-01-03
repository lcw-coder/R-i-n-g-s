program ring_analysis
implicit none
real, allocatable, dimension(:,:) :: x
integer, allocatable, dimension(:,:) :: edgelist,ring3list,ring4list,ring5list,ring6list,ring7list,ring8list
integer, allocatable, dimension(:) :: temp_array,arr
integer :: i,j,k,ierr,cnt3,cnt4,cnt5,cnt6,cnt7,nat,l,m,n,o,p,curr1,curr2,curr3,curr4,curr5,curr6,curr7,cnt_1,flag,curr8,q,cnt8
real :: rij,a,b,c,alpha,beta,gama
character(len=2),allocatable,dimension(:) :: atype
character(len=100) :: fname
open(unit=10,file='process2C.xyz',status='old',action='read')
a=38.3
b=38.3
c=38.3
alpha=90
beta=90
gama=90
read(10,*,iostat=ierr) nat
read(10,*,iostat=ierr)
allocate(x(nat,3),edgelist(nat,7),atype(nat))
allocate(ring3list(5000,3),ring4list(5000,4),ring5list(5000,5),ring6list(5000,6),ring7list(5000,7),ring8list(5000,8))
edgelist=0
ring3list=0
ring4list=0
ring5list=0
ring6list=0
ring7list=0
ring8list=0
cnt3=0
cnt4=0
cnt5=0
cnt6=0
cnt7=0
cnt8=0
open(unit=20,file='ring_count.txt',status='unknown',access='append')
write(20,'("#3mem #4mem  #5mem  #6mem  #7mem  #8mem")')
open(unit=40,file='3ringlist.txt',status='unknown',access='append')
open(unit=50,file='4ringlist.txt',status='unknown',access='append')
open(unit=60,file='5ringlist.txt',status='unknown',access='append')
open(unit=70,file='6ringlist.txt',status='unknown',access='append')
open(unit=80,file='7ringlist.txt',status='unknown',access='append')
open(unit=90,file='8ringlist.txt',status='unknown',access='append')
do while (ierr == 0)
    do i=1,nat
       read(10,*) atype(i),x(i,1),x(i,2),x(i,3)
    enddo
    do i=1,nat
       cnt_1=1
       edgelist(i,cnt_1) = i
       do j=1,nat
          if((i .ne. j).and.(atype(i)=='C').and.(atype(j)=='C')) then
             call calc_rij(x(i,1:3),x(j,1:3),rij,a,b,c,alpha,beta,gama)
             if (rij .le. 1.85) then
                cnt_1=cnt_1+1
                edgelist(i,cnt_1) = j
             endif
          endif
       enddo
    enddo
    do i=1,nat
       curr1 = i
       do j=2,7
          curr2 = edgelist(i,j)
          do k=2,7
             if((edgelist(curr2,k) .ne. curr1) .and. (curr2 .ne. 0)) then
                  curr3 = edgelist(curr2,k)
                  do l=2,7
                     if((edgelist(curr3,l) == curr1) .and. (curr3 .ne. 0)) then !start 3 memeber ring count
                          allocate(temp_array(3),arr(3))
                          arr(1) = curr1
                          arr(2) = curr2
                          arr(3) = curr3
                          temp_array = arr
                          call sort(3,arr)
                          call ring_analyze(3,cnt3,temp_array,arr,ring3list,ring4list,ring5list,ring6list,ring7list,ring8list,x,nat,a,b,c)
                          deallocate(temp_array,arr)
                     endif   !end 3 member ring count
                     if ((edgelist(curr3,l) .ne. curr2) .and. (edgelist(curr3,l) .ne. curr1) .and. (curr3 .ne. 0)) then
                          curr4 = edgelist(curr3,l)
                          do m=2,7
                             if((edgelist(curr4,m) == curr1) .and. (curr4 .ne. 0)) then !start 4 memeber ring count
                                allocate(temp_array(4),arr(4))
                                arr(1) = curr1
                                arr(2) = curr2
                                arr(3) = curr3
                                arr(4) = curr4
                                temp_array = arr
                                call sort(4,arr)
                                call ring_analyze(4,cnt4,temp_array,arr,ring3list,ring4list,ring5list,ring6list,ring7list,ring8list,x,nat,a,b,c)
                                deallocate(temp_array,arr)
                             endif   !end 4 member ring count
                             if((edgelist(curr4,m) .ne. curr3) .and. (edgelist(curr4,m) .ne. curr2) .and.(edgelist(curr4,m) .ne. curr1) .and. (curr4 .ne. 0)) then
                                curr5=edgelist(curr4,m)
                                do n=2,7
                                   if((edgelist(curr5,n) == curr1) .and. (curr5 .ne. 0)) then !start 5 memeber ring count
                                       allocate(temp_array(5),arr(5))
                                       arr(1) = curr1
                                       arr(2) = curr2
                                       arr(3) = curr3
                                       arr(4) = curr4
                                       arr(5) = curr5
                                       temp_array = arr
                                       call sort(5,arr)
                                       call ring_analyze(5,cnt5,temp_array,arr,ring3list,ring4list,ring5list,ring6list,ring7list,ring8list,x,nat,a,b,c)
                                       deallocate(temp_array,arr)
                                   endif   !end 5 member ring count
                                   if((edgelist(curr5,n).ne.curr4).and.(edgelist(curr5,n).ne.curr3).and.(edgelist(curr5,n).ne.curr2).and.(edgelist(curr5,n).ne.curr1).and.(curr5 .ne. 0)) then
                                      curr6=edgelist(curr5,n)
                                      do o=2,7
                                         if((edgelist(curr6,o) == curr1) .and. (curr6 .ne. 0)) then !start 6 memeber ring count
                                            allocate(temp_array(6),arr(6))
                                            arr(1) = curr1
                                            arr(2) = curr2
                                            arr(3) = curr3
                                            arr(4) = curr4
                                            arr(5) = curr5
                                            arr(6) = curr6
                                            temp_array = arr
                                            call sort(6,arr)
                                            call ring_analyze(6,cnt6,temp_array,arr,ring3list,ring4list,ring5list,ring6list,ring7list,ring8list,x,nat,a,b,c)
                                            deallocate(temp_array,arr)
                                         endif ! end 6 member ring count
                                         if((edgelist(curr6,o).ne.curr5).and.(edgelist(curr6,o).ne.curr4).and.(edgelist(curr6,o).ne.curr3).and.(edgelist(curr6,o)&
                                         &.ne.curr2).and.(edgelist(curr6,o).ne.curr1).and.(curr6.ne. 0)) then
                                             curr7=edgelist(curr6,o)
                                             do p=2,7
                                                if((edgelist(curr7,p) == curr1) .and. (curr7 .ne. 0)) then !start 7 memeber ring count
                                                     allocate(temp_array(7),arr(7))
                                                     arr(1) = curr1
                                                     arr(2) = curr2
                                                     arr(3) = curr3
                                                     arr(4) = curr4
                                                     arr(5) = curr5
                                                     arr(6) = curr6
                                                     arr(7) = curr7
                                                     temp_array = arr
                                                     call sort(7,arr)
                                                     call ring_analyze(7,cnt7,temp_array,arr,ring3list,ring4list,ring5list,ring6list,ring7list,ring8list,x,nat,a,b,c)
                                                     deallocate(temp_array,arr)
                                                endif ! end 7 member ring count
                                                if((edgelist(curr7,p).ne.curr6).and.(edgelist(curr7,p).ne.curr5).and.(edgelist(curr7,p).ne.curr4).and.(edgelist(curr7,p).ne.curr3).and.&
                                                &(edgelist(curr7,p).ne.curr2).and.(edgelist(curr7,p).ne.curr1).and.(curr7.ne. 0)) then       
                                                   curr8=edgelist(curr7,p)
                                                   do q=2,7
                                                     if((edgelist(curr8,q) == curr1) .and. (curr8 .ne. 0)) then   !start 8 member ring count
                                                        allocate(temp_array(8),arr(8))
                                                        arr(1) = curr1 
                                                        arr(2) = curr2 
                                                        arr(3) = curr3 
                                                        arr(4) = curr4 
                                                        arr(5) = curr5 
                                                        arr(6) = curr6 
                                                        arr(7) = curr7 
                                                        arr(8) = curr8 
                                                        temp_array = arr
                                                        call sort(8,arr)
                                                        call ring_analyze(8,cnt8,temp_array,arr,ring3list,ring4list,ring5list,ring6list,ring7list,ring8list,x,nat,a,b,c)
                                                        deallocate(temp_array,arr)
                                                      endif
                                                    enddo
                                                 endif
                                             enddo
                                         endif
                                      enddo
                                   endif
                                enddo
                             endif
                           enddo
                       endif  
                   enddo
               endif
           enddo
        enddo
    enddo
    write(20,*) cnt3,cnt4,cnt5,cnt6,cnt7,cnt8
    cnt3=0
    cnt4=0
    cnt5=0
    cnt6=0 
    cnt7=0
    cnt8=0
    ring3list=0
    ring4list=0
    ring5list=0
    ring6list=0 
    ring7list=0
    ring8list=0
    write(40,*)
    write(50,*)
    write(60,*)
    write(70,*)
    write(80,*)
    write(90,*)
    read(10,*,iostat=ierr) nat
    read(10,*,iostat=ierr)
    if (ierr == 0) then
       deallocate(x,edgelist,atype)
       allocate(x(nat,3),atype(nat),edgelist(nat,7))
       edgelist=0
    endif
enddo
!do i=1,nat
    !write(*,*) edgelist(i,1),edgelist(i,2),edgelist(i,3),edgelist(i,4),edgelist(i,5),edgelist(i,6)
!enddo
end program ring_analysis

subroutine calc_rij(x10,x20,r,a,b,c,alpha,beta,gama)
implicit none
real,dimension(3),intent(in) :: x10,x20
real,intent(out)::r
real, intent(in) :: a,b,c,alpha,beta,gama
real,dimension(3) :: x1

x1(:)=x10(:) - x20(:)
x1(1) = x1(1) - a*nint(x1(1)/a)
x1(2) = x1(2) - b*nint(x1(2)/b)
x1(3) = x1(3) - c*nint(x1(3)/c)
   
r = sqrt(x1(1)**2 + x1(2)**2 + x1(3)**2)
end subroutine calc_rij


subroutine sort(n,arr)
  integer,intent(inout) :: arr(n),n
  integer :: i,j,a
  do j=2, n
    a=arr(j)
    do i=j-1,1,-1
      if (arr(i)<=a) goto 10
      arr(i+1)=arr(i)
    end do
    i=0
10  arr(i+1)=a
  end do
  return
end subroutine sort

subroutine ring_analyze(n1,cnt,temp_array,arr,ring3list,ring4list,ring5list,ring6list,ring7list,ring8list,x,nat,a,b,c)
implicit none
integer,intent(in) :: n1,nat
integer,intent(inout)::cnt
real,intent(in)::a,b,c
real,dimension(nat,3),intent(in) :: x
integer,dimension(5000,3),intent(inout) :: ring3list
integer,dimension(5000,4),intent(inout) :: ring4list
integer,dimension(5000,5),intent(inout) :: ring5list
integer,dimension(5000,6),intent(inout) :: ring6list
integer,dimension(5000,7),intent(inout) :: ring7list
integer,dimension(5000,8),intent(inout) :: ring8list
integer, dimension(n1),intent(in) :: temp_array,arr
integer :: v,y,z,flag
real :: r,alpha=90.0,beta=90.0,gama=90.0

if (n1==3) then  !3 member ring start
  if (cnt == 0) then
     call calc_rij(x(temp_array(3),1:3),x(temp_array(1),1:3),r,a,b,c,alpha,beta,gama)
     if(r<=1.85) then
        ring3list(1,1:3) = arr(1:3)
        write(40,'(3(3X,I4))') temp_array
        cnt = cnt + 1
     endif
  else 
     flag = 0
     do z=1,cnt
        if ((ring3list(z,1)==arr(1)) .and. (ring3list(z,2)==arr(2)) .and. (ring3list(z,3)==arr(3))) flag=1
     enddo
     call calc_rij(x(temp_array(3),1:3),x(temp_array(1),1:3),r,a,b,c,alpha,beta,gama)
     if(r>1.85) flag=1
     if(flag == 0) then
        cnt = cnt+1
        ring3list(cnt,1:3) = arr(1:3)
        write(40,'(3(3X,I4))') temp_array
     endif
  endif
endif   ! 3 member ring end

if (n1 == 4) then   !4 member ring start
  if (cnt == 0) then
     flag = 0
     do y = 1,(n1-2)
        do z = y+2,n1
           if ((y.ne.1) .or. (z.ne.n1)) then
              call calc_rij(x(temp_array(y),1:3),x(temp_array(z),1:3),r,a,b,c,alpha,beta,gama)
              if(r<=1.85) flag=1
           endif
        enddo
     enddo
     call calc_rij(x(temp_array(1),1:3),x(temp_array(4),1:3),r,a,b,c,alpha,beta,gama)
     if(r > 1.85) flag=1
     if (flag==0) then
         ring4list(1,1:4) = arr(1:4)
        write(50,'(4(3X,I4))') temp_array
         cnt = cnt + 1
     endif
  else
     flag = 0
     do z=1,cnt
       if ((ring4list(z,1)==arr(1)) .and. (ring4list(z,2)==arr(2)) .and. (ring4list(z,3)==arr(3)) .and. (ring4list(z,4)==arr(4))) flag=1
     enddo
     do y = 1,(n1-2)
        do z = y+2,n1
           if ((y.ne.1) .or. (z.ne.n1)) then
              call calc_rij(x(temp_array(y),1:3),x(temp_array(z),1:3),r,a,b,c,alpha,beta,gama)
              if(r<=1.85) flag=1
           endif
        enddo
     enddo
     call calc_rij(x(temp_array(1),1:3),x(temp_array(4),1:3),r,a,b,c,alpha,beta,gama)
     if(r > 1.85) flag=1
     if(flag == 0) then
       cnt = cnt+1
       ring4list(cnt,1:4) = arr(1:4)
       write(50,'(4(3X,I4))') temp_array
     endif
  endif
endif        ! 4 member ring end

if (n1 == 5) then   !5 member ring start
  if (cnt == 0) then
     flag = 0
     do y = 1,(n1-2)
        do z = y+2,n1
           if ((y.ne.1) .or. (z.ne.n1)) then
              call calc_rij(x(temp_array(y),1:3),x(temp_array(z),1:3),r,a,b,c,alpha,beta,gama)
              if(r<=1.85) flag=1
           endif
        enddo
     enddo
     call calc_rij(x(temp_array(1),1:3),x(temp_array(5),1:3),r,a,b,c,alpha,beta,gama)
     if(r > 1.85) flag=1
     if (flag==0) then
         ring5list(1,1:5) = arr(1:5)
        write(60,'(5(3X,I4))') temp_array
         cnt = cnt + 1
     endif
  else
     flag = 0
     do z=1,cnt
       if ((ring5list(z,1)==arr(1)) .and. (ring5list(z,2)==arr(2)) .and. (ring5list(z,3)==arr(3)) .and. (ring5list(z,4)==arr(4)).and.(ring5list(z,5)==arr(5))) flag=1
     enddo
     do y = 1,(n1-2)
        do z = y+2,n1
           if ((y.ne.1) .or. (z.ne.n1)) then
              call calc_rij(x(temp_array(y),1:3),x(temp_array(z),1:3),r,a,b,c,alpha,beta,gama)
              if(r<=1.85) flag=1
           endif
        enddo
     enddo
     call calc_rij(x(temp_array(1),1:3),x(temp_array(5),1:3),r,a,b,c,alpha,beta,gama)
     if(r > 1.85) flag=1
     if(flag == 0) then
       cnt = cnt+1
       ring5list(cnt,1:5) = arr(1:5)
        write(60,'(5(3X,I4))') temp_array
     endif
  endif
endif        ! 5 member ring end

if (n1 == 6) then   !6 member ring start
  if (cnt == 0) then
     flag = 0
     do y = 1,(n1-2)
        do z = y+2,n1
           if ((y.ne.1) .or. (z.ne.n1)) then
              call calc_rij(x(temp_array(y),1:3),x(temp_array(z),1:3),r,a,b,c,alpha,beta,gama)
              if(r<=1.85) flag=1
           endif
        enddo
     enddo
     call calc_rij(x(temp_array(1),1:3),x(temp_array(6),1:3),r,a,b,c,alpha,beta,gama)
     if(r > 1.85) flag=1
     if (flag==0) then
         ring6list(1,1:6) = arr(1:6)
        write(70,'(6(3X,I4))') temp_array
         cnt = cnt + 1
     endif
  else
     flag = 0
     do z=1,cnt
       if ((ring6list(z,1)==arr(1)).and.(ring6list(z,2)==arr(2)).and.(ring6list(z,3)==arr(3)).and.(ring6list(z,4)==arr(4)).and.(ring6list(z,5)==arr(5)).and.(ring6list(z,6)==arr(6))) flag=1
     enddo
     do y = 1,(n1-2)
        do z = y+2,n1
           if ((y.ne.1) .or. (z.ne.n1)) then
              call calc_rij(x(temp_array(y),1:3),x(temp_array(z),1:3),r,a,b,c,alpha,beta,gama)
              if(r<=1.85) flag=1
           endif
        enddo
     enddo
     call calc_rij(x(temp_array(1),1:3),x(temp_array(6),1:3),r,a,b,c,alpha,beta,gama)
     if(r > 1.85) flag=1
     if(flag == 0) then
       cnt = cnt+1
       ring6list(cnt,1:6) = arr(1:6)
        write(70,'(6(3X,I4))') temp_array
     endif
  endif
endif        ! 6 member ring end

if (n1 == 7) then   !7 member ring start
  if (cnt == 0) then
     flag = 0
     do y = 1,(n1-2)
        do z = y+2,n1
           if ((y.ne.1) .or. (z.ne.n1)) then
              call calc_rij(x(temp_array(y),1:3),x(temp_array(z),1:3),r,a,b,c,alpha,beta,gama)
              if(r<=1.85) flag=1
           endif
        enddo
     enddo
     call calc_rij(x(temp_array(1),1:3),x(temp_array(7),1:3),r,a,b,c,alpha,beta,gama)
     if(r > 1.85) flag=1
     if (flag==0) then
         ring7list(1,1:7) = arr(1:7)
        write(80,'(7(3X,I4))') temp_array
         cnt = cnt + 1
     endif
  else
     flag = 0
     do z=1,cnt
       if ((ring7list(z,1)==arr(1)).and.(ring7list(z,2)==arr(2)).and.(ring7list(z,3)==arr(3)).and.(ring7list(z,4)==arr(4)).and.(ring7list(z,5)==arr(5)).and.(ring7list(z,6)==arr(6)).and.(ring7list(z,7)==arr(7))) flag=1
     enddo
     do y = 1,(n1-2)
        do z = y+2,n1
           if ((y.ne.1) .or. (z.ne.n1)) then
              call calc_rij(x(temp_array(y),1:3),x(temp_array(z),1:3),r,a,b,c,alpha,beta,gama)
              if(r<=1.85) flag=1
           endif
        enddo
     enddo
     call calc_rij(x(temp_array(1),1:3),x(temp_array(7),1:3),r,a,b,c,alpha,beta,gama)
     if(r > 1.85) flag=1
     if(flag == 0) then
       cnt = cnt+1
       ring7list(cnt,1:7) = arr(1:7)
       write(80,'(7(3X,I4))') temp_array
     endif
  endif
endif        ! 7 member ring end


if (n1 == 8) then   !8 member ring start
  if (cnt == 0) then
     flag = 0
     do y = 1,(n1-2)
        do z = y+2,n1
           if ((y.ne.1) .or. (z.ne.n1)) then
              call calc_rij(x(temp_array(y),1:3),x(temp_array(z),1:3),r,a,b,c,alpha,beta,gama)
              if(r<=1.85) flag=1
           endif
        enddo
     enddo
     call calc_rij(x(temp_array(1),1:3),x(temp_array(8),1:3),r,a,b,c,alpha,beta,gama)
     if(r > 1.85) flag=1
     if (flag==0) then
         ring8list(1,1:8) = arr(1:8)
        write(90,'(8(3X,I4))') temp_array
         cnt = cnt + 1
     endif
  else
     flag = 0
     do z=1,cnt
       if ((ring8list(z,1)==arr(1)).and.(ring8list(z,2)==arr(2)).and.(ring8list(z,3)==arr(3)).and.(ring8list(z,4)==arr(4)).and.(ring8list(z,5)==arr(5)).and.(ring8list(z,6)==arr(6))&
       &.and.(ring8list(z,7)==arr(7)).and.(ring8list(z,8)==arr(8))) flag=1
     enddo
     do y = 1,(n1-2)
        do z = y+2,n1
           if ((y.ne.1) .or. (z.ne.n1)) then
              call calc_rij(x(temp_array(y),1:3),x(temp_array(z),1:3),r,a,b,c,alpha,beta,gama)
              if(r<=1.85) flag=1
           endif
        enddo
     enddo
     call calc_rij(x(temp_array(1),1:3),x(temp_array(8),1:3),r,a,b,c,alpha,beta,gama)
     if(r > 1.85) flag=1
     if(flag == 0) then
       cnt = cnt+1
       ring8list(cnt,1:8) = arr(1:8)
        write(90,'(8(3X,I4))') temp_array
     endif
  endif
endif        ! 8 member ring end
end subroutine ring_analyze

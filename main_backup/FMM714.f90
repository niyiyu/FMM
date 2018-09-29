        program sundx
     
      integer n,m,k,si,sj,imin,jmin,x,y,i,j,size,maxsize,isroot,kg,rcsize,tj
!imin jmin
!x y
!i j
!size number of points in the narrow band
!isroot distinguish the source point
!kg switch
!rcsize
!tj
      real h,remin,tmin,isbandwt,dsh
      parameter(n=700,m=1000,maxsize=10000,si=500,sj=500,h=10)
!m n discretized grid size
!maxsize max size of the band
!si sj source
!h grid side length
      real,dimension(m)::z0
!z0 position of the surface
      real,dimension(n,m)::v,ttime
!v speed parameter
!ttime
        real,dimension(8)::time_temp,temp
        integer,dimension(n,m)::kinds,band
      real,dimension(maxsize)::bandwt
      integer,dimension(maxsize)::bandwx,bandwz 
!

       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%1 Input the slowness parameter !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! input the slowness parameter
      open(10,file='z0.txt',status='old')
      read(10,*)(z0(j),j=1,m)
      close(10)
!define the surface position
      z0=h*z0

      open(11,file='v.txt',status='old')
      read(11,*)((v(i,j),j=1,m),i=1,n)
      close(11) 
   
!      kg=0
!      if(kg.eq.1)then
!      open(15,file='v0.txt',status='new')
!      write(15,13)((v(i,j),j=1,m),i=1,n)
!13    format(1201f10.3)
!      close(15)            
!      endif
!!! input the initize ray location 

    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%2Initialize Alive,Narrow,Far point !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!initialize all of the points with a big enough number  
      time_temp=100.0
      ttime=100.0
      band=0
!
      bandwz=0
      bandwx=0
      bandwt=100.0
!
      size=0
      isroot=1
!initialize the travel time of the source point :0
      ttime(si,sj)=0.0  
!initialize the characteristic matrix
        do i=1,n
        do j=1,m
            if ( v(i,j).eq.0.0 )then
                kinds(i,j)=4
            else
                if(ttime(i,j).eq.100)then
                kinds(i,j)=1
                else
                    if(band(i,j).eq.1)then
                    kinds(i,j)=2
                    else
                    kinds(i,j)=3
                    endif
                endif
            endif
        end do
        end do
        kinds(si,sj)=3
!initialize the eight points around the source point
!left up
      if(((si-1).ge.1).and.((sj-1).ge.1).and.(v(si-1,sj-1).ne.0.0))then
          ttime(si-1,sj-1)=sqrt(2.0)*h/v(si-1,sj-1)
          call bandw_insert(si-1,sj-1,ttime(si-1,sj-1),bandwx,bandwz,bandwt,isroot,size,maxsize)  
          kinds(si-1,sj-1)=2 
          band(si-1,sj-1)=1
      endif
!up
      if(((si-1).ge.1).and.(v(si-1,sj).ne.0.0))then
          ttime(si-1,sj)=h/v(si-1,sj)
          call bandw_insert(si-1,sj,ttime(si-1,sj),bandwx,bandwz,bandwt,isroot,size,maxsize) 
        kinds(si-1,sj)=2  
        band(si-1,sj)=1 
      endif  
!right up
      if(((si-1).ge.1).and.((sj+1).le.m).and.(v(si-1,sj+1).ne.0.0))then
          ttime(si-1,sj+1)=sqrt(2.0)*h/v(si-1,sj+1)
          call bandw_insert(si-1,sj+1,ttime(si-1,sj+1),bandwx,bandwz,bandwt,isroot,size,maxsize)  
          kinds(si-1,sj+1)=2 
          band(si-1,sj+1)=1
      endif 
!left down
      if(((si+1).le.n).and.((sj-1).ge.1).and.(v(si+1,sj-1).ne.0.0))then
          ttime(si+1,sj-1)=sqrt(2.0)*h/v(si+1,sj-1)
          call bandw_insert(si+1,sj-1,ttime(si+1,sj-1),bandwx,bandwz,bandwt,isroot,size,maxsize) 
          kinds(si+1,sj-1)=2  
          band(si+1,sj-1)=1
      endif
!down
      if(((si+1).le.n).and.(v(si+1,sj).ne.0.0))then
          ttime(si+1,sj)=h/v(si+1,sj) 
          call bandw_insert(si+1,sj,ttime(si+1,sj),bandwx,bandwz,bandwt,isroot,size,maxsize)
          kinds(si+1,sj)=2 
          band(si+1,sj)=1
      endif 
!right down
      if(((si+1).le.n).and.((sj+1).le.m).and.(v(si+1,sj+1).ne.0.0))then
          ttime(si+1,sj+1)=sqrt(2.0)*h/v(si+1,sj+1)
          call bandw_insert(si+1,sj+1,ttime(si+1,sj+1),bandwx,bandwz,bandwt,isroot,size,maxsize) 
          kinds(si+1,sj+1)=2 
          band(si+1,sj+1)=1 
      endif
!left
      if(((sj-1).ge.1).and.(v(si,sj-1).ne.0.0))then
          ttime(si,sj-1)=h/v(si,sj-1)        
          call bandw_insert(si,sj-1,ttime(si,sj-1),bandwx,bandwz,bandwt,isroot,size,maxsize) 
          kinds(si,sj-1)=2 
          band(si,sj-1)=1
      endif
!right
      if(((sj+1).le.m).and.(v(si,sj+1).ne.0.0))then
          ttime(si,sj+1)=h/v(si,sj+1)      
          call bandw_insert(si,sj+1,ttime(si,sj+1),bandwx,bandwz,bandwt,isroot,size,maxsize) 

          band(si,sj+1)=1
      endif     


100     call bandw_remove(imin,jmin,tmin,bandwx,bandwz,bandwt,isroot,size,maxsize)
            kinds(imin,jmin)=3
            band(imin,jmin)=0
!up
        if((imin-1).ge.1)then
            call computate(ttime,n,m,(imin-1),jmin,v,h,time_temp)
            
            call sort_time(time_temp)
            temp=time_temp(1)
            select case(kinds((imin-1),jmin))
            case (4)
            case (3)
            case (2)
                if (temp(1).le.ttime((imin-1),j))then
                    ttime((imin-1),j)=temp(1)
                endif
            case (1)
                ttime((imin-1),jmin)=temp(1)
                call bandw_insert((imin-1),jmin,ttime((imin-1),jmin),bandwx,bandwz,bandwt,isroot,size,maxsize)
                kinds((imin-1),jmin)=2
                band((imin-1),jmin)=1
            end select
            time_temp=100.0
        endif

!down
        if((imin+1).le.n)then
            call computate(ttime,n,m,(imin+1),jmin,v,h,time_temp)
            call sort_time(time_temp)
            temp=time_temp(1)
            select case(kinds((imin+1),jmin))
            case (4)
            case (3)
            case (2)
                if (temp(1).le.ttime((imin+1),j))then
                    ttime((imin+1),j)=temp(1)
                endif
            case (1)
                ttime((imin+1),jmin)=temp(1)
                call bandw_insert((imin+1),jmin,ttime((imin+1),jmin),bandwx,bandwz,bandwt,isroot,size,maxsize)
                kinds((imin+1),jmin)=2
                band((imin+1),jmin)=1
            end select
            time_temp=100.0
        endif

!left
        if((jmin-1).ge.1)then
            call computate(ttime,n,m,imin,(jmin-1),v,h,time_temp)
            call sort_time(time_temp)
            temp=time_temp(1)
            select case(kinds(imin,(jmin-1)))
            case (4)
            case (3)
            case (2)
                if (temp(1).le.ttime(imin,(jmin-1)))then
                    ttime(imin,(jmin-1))=temp(1)
                endif
            case (1)
                ttime(imin,(jmin-1))=temp(1)
                call bandw_insert(imin,(jmin-1),ttime(imin,(jmin-1)),bandwx,bandwz,bandwt,isroot,size,maxsize)
                kinds(imin,(jmin-1))=2
                band(imin,(jmin-1))=1
            end select
            time_temp=100.0
        endif

 
!right
        if((jmin+1).le.m)then
            call computate(ttime,n,m,imin,(jmin+1),v,h,time_temp)
            call sort_time(time_temp)
            temp=time_temp(1)
            select case(kinds(imin,(jmin+1)))
            case (4)
            case (3)
            case (2)
                if (temp(1).le.ttime(imin,(jmin+1)))then
                    ttime(imin,(jmin+1))=temp(1)
                endif
            case (1)
                ttime(imin,(jmin+1))=temp(1)
                call bandw_insert(imin,(jmin+1),ttime(imin,(jmin+1)),bandwx,bandwz,bandwt,isroot,size,maxsize)
                kinds(imin,(jmin+1))=2
                band(imin,(jmin+1))=1
            end select
            time_temp=100.0
        endif
                  !!!!! output czht.txt
      
      if(size.gt.0)goto 100
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%3 Output to the data of Recpt/Narrow/Alive !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      kg=1  
      if(kg.eq.1)then
      open(17,file='ttime.txt',status='replace')
      write(17,35)((ttime(i,j),j=1,m),i=1,n)
      !write(17,35)(time_temp(i),i=1,8)
      !write(17,35)temp(1)
      !write(17,35)(time_temp(i),i=1,8)
      !write(17,35)(kinds(imin,(jmin+1)))
      !write(17,35)((kinds(i,j),j=1,m),i=1,n)
      close(17)
35    format(1000F10.3)  
!35    format(500I3)  
      endif 






      end




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine bandw_insert(insertz,insertx,insertt,bandwx,bandwz,bandwt,isroot,size,maxsize)    
     
      integer insertz,insertx,size,maxsize,isroot
      real insertt
      integer,dimension(maxsize)::bandwx,bandwz
      real,dimension(maxsize)::bandwt 
     
      if(isroot.eq.1)then       
        bandwz(1)=insertx
        bandwx(1)=insertz
        bandwt(1)=insertt    
        size=size+1
        isroot=0
        call filterdown(bandwx,bandwz,bandwt,size,maxsize)
      else
        bandwz(size+1)=insertx
        bandwx(size+1)=insertz               
        bandwt(size+1)=insertt
        size=size+1
        call filterup(bandwx,bandwz,bandwt,size,maxsize)
      endif
      
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!remove
      subroutine bandw_remove(imin,jmin,tmin,bandwx,bandwz,bandwt,isroot,csize,msize)

      integer imin,jmin,isroot,csize,msize
      real tmin
      integer,dimension(msize)::bandwx,bandwz
      real,dimension(msize)::bandwt 

      if(isroot.eq.1)then      
        bandwx(1)=bandwx(csize+1) 
        bandwz(1)=bandwz(csize+1)      
        bandwt(1)=bandwt(csize+1)
        bandwz(csize+1)=0
        bandwx(csize+1)=0              
        bandwt(csize+1)=100.0
        call filterdown(bandwx,bandwz,bandwt,csize,msize)
      endif
!!!!!!output point
        imin=bandwx(1)
        jmin=bandwz(1)
        tmin=bandwt(1)
!number of point in the band
        csize=csize-1
        isroot=1
     
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!升序排列
      subroutine filterdown(bandwx,bandwz,bandwt,size,maxsize)

      integer i,j,size,maxsize,tempx,tempz
      real tempt
      integer,dimension(maxsize)::bandwx,bandwz
      real,dimension(maxsize)::bandwt  
     
      i=1
      j=2*i
      tempx=bandwx(i)
      tempz=bandwz(i)
      tempt=bandwt(i)
101   if(j.le.size)then       
        if(((j+1).le.size).and.(bandwt(j).ge.bandwt(j+1)))j=j+1
          if(tempt.ge.bandwt(j))then
              bandwx(i)=bandwx(j)
              bandwz(i)=bandwz(j)
              bandwt(i)=bandwt(j)
              i=j
              j=2*i
           if(j.le.size)goto 101
           endif
       endif        
       bandwx(i)=tempx
       bandwz(i)=tempz
       bandwt(i)=tempt 

      end  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!降序排列
      subroutine filterup(bandwx,bandwz,bandwt,size,maxsize)

      integer i,j,size,maxsize,tempx,tempz
      real tempt
      integer,dimension(maxsize)::bandwx,bandwz
      real,dimension(maxsize)::bandwt


      j=size
      i=j/2
      tempx=bandwx(j)
      tempz=bandwz(j)
      tempt=bandwt(j)
102   if(j.ge.1)then
        if(bandwt(i).ge.tempt)then
           bandwx(j)=bandwx(i)
           bandwz(j)=bandwz(i)
           bandwt(j)=bandwt(i)
           j=i
           i=j/2
           if(j.ge.1)goto 102
        endif 
      endif
      bandwx(j)=tempx
      bandwz(j)=tempz
      bandwt(j)=tempt 

      end



    subroutine computate(ttime,n,m,i,j,v,h,time_temp)
    integer m,n,i,j
    real h
    real,dimension(n,m)::v,ttime
    real,dimension(8)::time_temp

    if((i-1).ge.1)then
        time_temp(1)=ttime((i-1),j)+h/v((i-1),j)
    endif
    if((i+1).le.n)then
        time_temp(2)=ttime((i+1),j)+h/v((i+1),j)
    endif
    if((j-1).ge.1)then
        time_temp(3)=ttime(i,(j-1))+h/v(i,(j-1))
    endif
    if((j+1).le.m)then
        time_temp(4)=ttime(i,(j+1))+h/v(i,(j+1))
    endif
    if(((i-1).ge.1).and.((j-1).ge.1))then
        time_temp(5)=(ttime((i-1),j)+ttime(i,(j-1))+sqrt(abs((2*(h/v(i,j))**2-(ttime((i-1),j)-ttime(i,(j-1)))**2))))/2
    endif
    if(((i+1).le.n).and.((j-1).ge.1))then
        time_temp(6)=(ttime((i+1),j)+ttime(i,(j-1))+sqrt(abs((2*(h/v(i,j))**2-(ttime((i+1),j)-ttime(i,(j-1)))**2))))/2
    endif
    if(((i-1).ge.1).and.((j+1).le.m))then
        time_temp(7)=(ttime((i-1),j)+ttime(i,(j+1))+sqrt(abs((2*(h/v(i,j))**2-(ttime((i-1),j)-ttime(i,(j+1)))**2))))/2
    endif
    if(((i+1).le.n).and.((j+1).le.m))then
        time_temp(8)=(ttime((i+1),j)+ttime(i,(j+1))+sqrt(abs((2*(h/v(i,j))**2-(ttime((i+1),j)-ttime(i,(j+1)))**2))))/2
    endif

    end 


     
    subroutine sort_time(time_temp)
      real time_temp(8),tmp
      integer length
      integer min,i,j


      length=8
      do i=1,length-1
        min=i
        do j=i+1,length
            if (time_temp(j).lt.time_temp(min)) then
                min=j

            end if

        end do
        tmp = time_temp(min)
        time_temp(min) = time_temp(i)
        time_temp(i) = tmp 
      end do
        return
    end subroutine
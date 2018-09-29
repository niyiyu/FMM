        program FMM_niyiyu
        integer n,m,si,sj,imin,jmin,i,j,size,maxsize,isroot,breaker
!imin jmin
!i j
!size number of points in the narrow band
!isroot distinguish the source point
!breaker switch
        real(8) h,tmin
        parameter(n=350,m=500,si=125,sj=250,h=20,maxsize=10000)
!m n discretized grid size
!maxsize max size of the band
!si sj source
!h grid side length
        real(8),dimension(n,m)::v,travel_time
!v speed parameter
!travel_time
        real(8),dimension(12)::time_temp,temp
        integer,dimension(n,m)::kinds,band
        real(8),dimension(maxsize)::bandwt
        integer,dimension(maxsize)::bandwx,bandwz
        
        call inverse(n,m,si,sj,imin,jmin,i,j,size,maxsize,isroot,breaker,h,tmin,&
        v,travel_time,time_temp,temp,kinds,band,bandwt,bandwx,bandwz)

        end program


        subroutine inverse(n,m,si,sj,imin,jmin,i,j,size,maxsize,isroot,breaker,h,tmin,&
        v,travel_time,time_temp,temp,kinds,band,bandwt,bandwx,bandwz)
        integer si,sj,imin,jmin,i,j,size,maxsize,isroot,breaker
        real(8) h,tmin
        real(8),dimension(n,m)::v,travel_time
        real(8),dimension(12)::time_temp,temp
        integer,dimension(n,m)::kinds,band
        real(8),dimension(maxsize)::bandwt
        integer,dimension(maxsize)::bandwx,bandwz 
        integer radius

        print*, 'Input the radius:'
        read(*,*) radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%1 Input the slowness parameter !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        open(11,file='v.txt',status='old')
        read(11,*)((v(i,j),j=1,m),i=1,n)
        close(11) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%2Initialize the matrix!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        time_temp=50.0
        bandwz=0
        bandwx=0
        bandwt=50.0
        size=0
        isroot=1
        
        breaker=1

        if(breaker.eq.0)then
            travel_time=50.0
            travel_time(si,sj)=0.0 
        endif
        if(breaker.eq.1)then
            open(13,file='travel_time.txt',status='old')
            read(13,*)((travel_time(i,j),j=1,m),i=1,n)
            close(13)
        endif

        if(breaker.eq.1)then
            open(14,file='band.txt',status='old')
            read(14,*)((band(i,j),j=1,m),i=1,n)
            close(14)
        endif
        if(breaker.eq.0)then
            band=0
        endif
!initialize the characteristic matrix
        do i=1,n
        do j=1,m
            if ( v(i,j).eq.0.0 )then
                kinds(i,j)=4
            else
                if(travel_time(i,j).eq.50)then
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

        if(breaker.eq.0)then
!initialize the eight points around the source point
!left up
        if(((si-1).ge.1).and.((sj-1).ge.1).and.(v(si-1,sj-1).ne.0.0))then
           travel_time(si-1,sj-1)=sqrt(2.0)*h/v(si-1,sj-1) 
           kinds(si-1,sj-1)=2 
           band(si-1,sj-1)=1
        endif
!up
        if(((si-1).ge.1).and.(v(si-1,sj).ne.0.0))then
           travel_time(si-1,sj)=h/v(si-1,sj) 
           kinds(si-1,sj)=2  
           band(si-1,sj)=1 
        endif  
!right up
        if(((si-1).ge.1).and.((sj+1).le.m).and.(v(si-1,sj+1).ne.0.0))then
            travel_time(si-1,sj+1)=sqrt(2.0)*h/v(si-1,sj+1)
            kinds(si-1,sj+1)=2 
            band(si-1,sj+1)=1
        endif 
!left down
        if(((si+1).le.n).and.((sj-1).ge.1).and.(v(si+1,sj-1).ne.0.0))then
            travel_time(si+1,sj-1)=sqrt(2.0)*h/v(si+1,sj-1)
            kinds(si+1,sj-1)=2  
            band(si+1,sj-1)=1
        endif
!down
        if(((si+1).le.n).and.(v(si+1,sj).ne.0.0))then
            travel_time(si+1,sj)=h/v(si+1,sj) 
            kinds(si+1,sj)=2 
            band(si+1,sj)=1
        endif 
!right down
        if(((si+1).le.n).and.((sj+1).le.m).and.(v(si+1,sj+1).ne.0.0))then
            travel_time(si+1,sj+1)=sqrt(2.0)*h/v(si+1,sj+1)    
            kinds(si+1,sj+1)=2 
            band(si+1,sj+1)=1 
        endif
!left
        if(((sj-1).ge.1).and.(v(si,sj-1).ne.0.0))then
            travel_time(si,sj-1)=h/v(si,sj-1)        
            kinds(si,sj-1)=2 
            band(si,sj-1)=1
        endif
!right
        if(((sj+1).le.m).and.(v(si,sj+1).ne.0.0))then
            travel_time(si,sj+1)=h/v(si,sj+1)      
            kinds(si,(sj+1))=2
            band(si,sj+1)=1
        endif     

        endif

        do i=1,n
        do j=1,m
            if ( band(i,j).eq.1)then
                call bandw_insert(i,j,travel_time(i,j),bandwx,bandwz,bandwt,isroot,size,maxsize) 
            endif
        end do
        end do



100     call bandw_remove(imin,jmin,tmin,bandwx,bandwz,bandwt,isroot,size,maxsize)
            kinds(imin,jmin)=3
            band(imin,jmin)=0   
!up
        if((imin-1).ge.1)then
            call computate(travel_time,n,m,(imin-1),jmin,v,h,time_temp)
            
            call sort_time(time_temp)
            temp=time_temp(1)
            select case(kinds((imin-1),jmin))
            case (4)
            case (3)
            case (2)
                if(temp(1).ge.travel_time((imin-1),jmin))then
                    travel_time((imin-1),jmin)=temp(1)
                endif
            case (1)
                travel_time((imin-1),jmin)=temp(1)
                if(((imin-1-si)**2+(jmin-sj)**2).le.radius**2)then
                call bandw_insert((imin-1),jmin,travel_time((imin-1),jmin),bandwx,bandwz,bandwt,isroot,size,maxsize)
                endif
                kinds((imin-1),jmin)=2
                band((imin-1),jmin)=1
            end select
            time_temp=50.0
        endif
!down
        if((imin+1).le.n)then
            call computate(travel_time,n,m,(imin+1),jmin,v,h,time_temp)
            call sort_time(time_temp)
            temp=time_temp(1)
            select case(kinds((imin+1),jmin))
            case (4)
            case (3)
            case (2)
                if(temp(1).ge.travel_time((imin+1),jmin))then
                    travel_time((imin+1),jmin)=temp(1)
                endif
            case (1)
                travel_time((imin+1),jmin)=temp(1)
                if(((imin+1-si)**2+(jmin-sj)**2).le.radius**2)then
                call bandw_insert((imin+1),jmin,travel_time((imin+1),jmin),bandwx,bandwz,bandwt,isroot,size,maxsize)
                endif
                kinds((imin+1),jmin)=2
                band((imin+1),jmin)=1
            end select
            time_temp=50.0
        endif
!left
        if((jmin-1).ge.1)then
            call computate(travel_time,n,m,imin,(jmin-1),v,h,time_temp)
            call sort_time(time_temp)
            temp=time_temp(1)
            select case(kinds(imin,(jmin-1)))
            case (4)
            case (3)
            case (2)
                if(temp(1).ge.travel_time(imin,(jmin-1)))then
                    travel_time(imin,(jmin-1))=temp(1)
                endif
            case (1)
                travel_time(imin,(jmin-1))=temp(1)
                if(((imin-si)**2+(jmin-1-sj)**2).le.radius**2)then
                call bandw_insert(imin,(jmin-1),travel_time(imin,(jmin-1)),bandwx,bandwz,bandwt,isroot,size,maxsize)
                endif
                kinds(imin,(jmin-1))=2
                band(imin,(jmin-1))=1
            end select
            time_temp=50.0
        endif
!right
        if((jmin+1).le.m)then
            call computate(travel_time,n,m,imin,(jmin+1),v,h,time_temp)
            call sort_time(time_temp)
            temp=time_temp(1)
            select case(kinds(imin,(jmin+1)))
            case (4)
            case (3)
            case (2)
                if(temp(1).ge.travel_time(imin,(jmin+1)))then
                    travel_time(imin,(jmin+1))=temp(1)
                endif
            case (1)
                travel_time(imin,(jmin+1))=temp(1)
                if(((imin-si)**2+(jmin+1-sj)**2).le.radius**2)then
                call bandw_insert(imin,(jmin+1),travel_time(imin,(jmin+1)),bandwx,bandwz,bandwt,isroot,size,maxsize)
                endif
                kinds(imin,(jmin+1))=2
                band(imin,(jmin+1))=1
            end select
            time_temp=50.0
        endif
        
        if(size.ge.0)goto 100
        !call band_write(bandwx,bandwz,bandwt,maxsize)
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%3 Output!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        breaker=1
        if(breaker.eq.1)then
            open(17,file='travel_time.txt',status='replace')
            write(17,35)((travel_time(i,j),j=1,m),i=1,n)
            close(17)
35          format(500F12.8)  
        endif 

        end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine bandw_insert(insertz,insertx,insertt,bandwx,bandwz,bandwt,isroot,size,maxsize)
        integer insertz,insertx,size,maxsize,isroot
        real(8) insertt
        integer,dimension(maxsize)::bandwx,bandwz
        real(8),dimension(maxsize)::bandwt 
     
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
      
        end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine bandw_remove(imin,jmin,tmin,bandwx,bandwz,bandwt,isroot,csize,msize)
        integer imin,jmin,isroot,csize,msize
        real(8) tmin
        integer,dimension(msize)::bandwx,bandwz
        real(8),dimension(msize)::bandwt 

        if(isroot.eq.1)then      
            bandwx(1)=bandwx(csize+1) 
            bandwz(1)=bandwz(csize+1)      
            bandwt(1)=bandwt(csize+1)
            bandwz(csize+1)=0
            bandwx(csize+1)=0              
            bandwt(csize+1)=50.0
            call filterdown(bandwx,bandwz,bandwt,csize,msize)
        endif
        imin=bandwx(1)
        jmin=bandwz(1)
        tmin=bandwt(1)
        csize=csize-1
        isroot=1
     
        end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine filterdown(bandwx,bandwz,bandwt,size,maxsize)
        integer i,j,size,maxsize,tempx,tempz
        real(8) tempt
        integer,dimension(maxsize)::bandwx,bandwz
        real(8),dimension(maxsize)::bandwt  
     
        i=1
        j=2*i
        tempx=bandwx(i)
        tempz=bandwz(i)
        tempt=bandwt(i)
101     if(j.le.size)then       
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

        end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine filterup(bandwx,bandwz,bandwt,size,maxsize)
        integer i,j,size,maxsize,tempx,tempz
        real(8) tempt
        integer,dimension(maxsize)::bandwx,bandwz
        real(8),dimension(maxsize)::bandwt

        j=size
        i=j/2
        tempx=bandwx(j)
        tempz=bandwz(j)
        tempt=bandwt(j)
102     if(j.ge.1)then
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

        end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine computate(travel_time,n,m,i,j,v,h,time_temp)
        integer m,n,i,j
        real(8) h
        real(8),dimension(n,m)::v,travel_time
        real(8),dimension(12)::time_temp

        if((i-1).ge.1)then
            time_temp(1)=travel_time((i-1),j)+h/v((i-1),j)
        endif
        if((i+1).le.n)then
            time_temp(2)=travel_time((i+1),j)+h/v((i+1),j)
        endif
        if((j-1).ge.1)then
            time_temp(3)=travel_time(i,(j-1))+h/v(i,(j-1))
        endif
        if((j+1).le.m)then
            time_temp(4)=travel_time(i,(j+1))+h/v(i,(j+1))
        endif
        if(((i-1).ge.1).and.((j-1).ge.1))then
            time_temp(5)=(travel_time((i-1),j)+travel_time(i,(j-1))+sqrt(abs((2*(h/v(i,j))**2-(travel_time((i-1),j)-&
            travel_time(i,(j-1)))**2))))/2
            time_temp(9)=travel_time((i-1),(j-1))+sqrt(2.0)*h/v(i,j)
        endif
        if(((i+1).le.n).and.((j-1).ge.1))then
            time_temp(6)=(travel_time((i+1),j)+travel_time(i,(j-1))+sqrt(abs((2*(h/v(i,j))**2-(travel_time((i+1),j)-&
            travel_time(i,(j-1)))**2))))/2
            time_temp(10)=travel_time((i+1),(j-1))+sqrt(2.0)*h/v(i,j)
        endif
        if(((i-1).ge.1).and.((j+1).le.m))then
            time_temp(7)=(travel_time((i-1),j)+travel_time(i,(j+1))+sqrt(abs((2*(h/v(i,j))**2-(travel_time((i-1),j)-&
            travel_time(i,(j+1)))**2))))/2
            time_temp(11)=travel_time((i-1),(j+1))+sqrt(2.0)*h/v(i,j)
        endif
        if(((i+1).le.n).and.((j+1).le.m))then
            time_temp(8)=(travel_time((i+1),j)+travel_time(i,(j+1))+sqrt(abs((2*(h/v(i,j))**2-(travel_time((i+1),j)-&
            travel_time(i,(j+1)))**2))))/2
            time_temp(12)=travel_time((i+1),(j+1))+sqrt(2.0)*h/v(i,j)
        endif

        end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        subroutine sort_time(time_temp)
        real(8) time_temp(8),tmp
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine band_write(bandwx,bandwz,bandwt,maxsize)
        integer maxsize
        real(8),dimension(maxsize)::bandwt
        integer,dimension(maxsize)::bandwx,bandwz
        breaker=1  
        if(breaker.eq.1)then
            open(17,file='band_output.txt',status='replace')
            write(17,37)(bandwx(i),i=1,maxsize)
            write(17,37)(bandwz(i),i=1,maxsize)
            write(17,38)(bandwt(i),i=1,maxsize)
            close(17)
37          format(10000I5)
38          format(10000F12.8)    
        endif 
        end subroutine
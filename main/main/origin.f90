      program sundx
     
      integer n,m,si,sj,imin,jmin,x,y,i,j,csize,msize,isroot,kg,rcsize,tj
!imin jmin
!x y
!i j
!csize number of points in the narrow band
!isroot distinguish the source point
!kg switch
!rcsize
!tj
      real h,remin,tmin,isbandwt,dsh
      parameter(n=1000,m=1000,msize=10000,si=500,sj=500,h=10)
!m n discretized grid size
!msize max size of the band
!si sj source
!h grid side length
      real,dimension(m)::z0
!z0 position of the surface
      real,dimension(n,m)::v,recpt,alive
!v speed parameter
!recpt
!alive

      real,dimension(msize)::bandwt
      integer,dimension(msize)::bandwx,bandwz 
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
      recpt=100.0
      alive=100.0
!
      bandwz=0
      bandwx=0
      bandwt=100.0
!
      csize=0
      isroot=1
!initialize the travel time of the source point :0
      recpt(si,sj)=0.0
      alive(si,sj)=0.0   
!initialize the eight points around the source point
!left up
      if(((si-1).ge.1).and.((sj-1).ge.1).and.(v(si-1,sj-1).ne.0.0))then
          recpt(si-1,sj-1)=sqrt(2.0)*h/v(si-1,sj-1)
          alive(si-1,sj-1)=sqrt(2.0)*h/v(si-1,sj-1)
          call bandw_insert(si-1,sj-1,recpt(si-1,sj-1),bandwx,bandwz,bandwt,isroot,csize,msize)   
      endif
!up
      if(((si-1).ge.1).and.(v(si-1,sj).ne.0.0))then
          recpt(si-1,sj)=h/v(si-1,sj)
          alive(si-1,sj)=h/v(si-1,sj)
          call bandw_insert(si-1,sj,recpt(si-1,sj),bandwx,bandwz,bandwt,isroot,csize,msize)   
      endif  
!right up
      if(((si-1).ge.1).and.((sj+1).le.m).and.(v(si-1,sj+1).ne.0.0))then
          recpt(si-1,sj+1)=sqrt(2.0)*h/v(si-1,sj+1)
          alive(si-1,sj+1)=sqrt(2.0)*h/v(si-1,sj+1)
          call bandw_insert(si-1,sj+1,recpt(si-1,sj+1),bandwx,bandwz,bandwt,isroot,csize,msize)   
      endif 
!left down
      if(((si+1).le.n).and.((sj-1).ge.1).and.(v(si+1,sj-1).ne.0.0))then
          recpt(si+1,sj-1)=sqrt(2.0)*h/v(si+1,sj-1)
          alive(si+1,sj-1)=sqrt(2.0)*h/v(si+1,sj-1)
          call bandw_insert(si+1,sj-1,recpt(si+1,sj-1),bandwx,bandwz,bandwt,isroot,csize,msize)   
      endif
!down
      if(((si+1).le.n).and.(v(si+1,sj).ne.0.0))then
          recpt(si+1,sj)=h/v(si+1,sj) 
          alive(si+1,sj)=h/v(si+1,sj)
          call bandw_insert(si+1,sj,recpt(si+1,sj),bandwx,bandwz,bandwt,isroot,csize,msize) 
      endif 
!right down
      if(((si+1).le.n).and.((sj+1).le.m).and.(v(si+1,sj+1).ne.0.0))then
          recpt(si+1,sj+1)=sqrt(2.0)*h/v(si+1,sj+1)
          alive(si+1,sj+1)=sqrt(2.0)*h/v(si+1,sj+1)
          call bandw_insert(si+1,sj+1,recpt(si+1,sj+1),bandwx,bandwz,bandwt,isroot,csize,msize)   
      endif
!left
      if(((sj-1).ge.1).and.(v(si,sj-1).ne.0.0))then
          recpt(si,sj-1)=h/v(si,sj-1)
          alive(si,sj-1)=h/v(si,sj-1)          
          call bandw_insert(si,sj-1,recpt(si,sj-1),bandwx,bandwz,bandwt,isroot,csize,msize)  
      endif
!right
      if(((sj+1).le.m).and.(v(si,sj+1).ne.0.0))then
          recpt(si,sj+1)=h/v(si,sj+1)    
          alive(si,sj+1)=h/v(si,sj+1)     
          call bandw_insert(si,sj+1,recpt(si,sj+1),bandwx,bandwz,bandwt,isroot,csize,msize) 
      endif      
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Choose the smallest travel time in Narrow Band !!!!!!!!!!!!!!!!!!!!!!!!!! 

      tj=0  
!choose and then remove the point(csize-1)    
100   call bandw_remove(imin,jmin,tmin,bandwx,bandwz,bandwt,isroot,csize,msize)
        
!!!!!!!!!!!!! Add the point(imin,jmin) to Alive and remove it from Narrow Band !!!!!!!!!!!!!!!!!!!!!

      alive(imin,jmin)=tmin      
      
!!!!!!! Rcompute the valuses of T at neighbor(imin-1,jmin) if it is in Narrow Band or Far Away !!!!!
!!!!!!!!! If the neighbor piont is in Far Away,add it to the set of Narrow Band !!!!!!!!!!!!!!!!!!!!

!!!boundry condition      
      if((imin-1).ge.1)then
!!!hasn't been calculated and isn't above the ground
      if((alive(imin-1,jmin).eq.100.0).and.(v(imin-1,jmin).ne.0.0))then
       
        isbandwt=recpt(imin-1,jmin)

!!!down-left
!!!if this point has been calculated but in the calculating area
        if(((jmin-1).ge.1).and.(alive(imin,jmin-1).ne.100.0))then
!!!re-calculate and choose the smaller one(call the recompute2 function)
          recpt(imin-1,jmin)=recompute2(alive(imin,jmin),alive(imin,jmin-1),1.0/v(imin-1,jmin))
        else
!!!else just calculate!
          recpt(imin-1,jmin)=alive(imin,jmin)+h/v(imin-1,jmin)
        endif
        remin=recpt(imin-1,jmin)

!!!down-right
        if(((jmin+1).le.m).and.(alive(imin,jmin+1).ne.100.0))then
          recpt(imin-1,jmin)=recompute2(alive(imin,jmin),alive(imin,jmin+1),1.0/v(imin-1,jmin))
          if(recpt(imin-1,jmin).lt.remin)then
            remin=recpt(imin-1,jmin)
          endif
        endif      
!!!right-down+up
        if(((jmin+1).le.m).and.(alive(imin-1,jmin+1).ne.100.0))then

          if(alive(imin,jmin+1).ne.100.0)then
            recpt(imin-1,jmin)=recompute2(alive(imin-1,jmin+1),alive(imin,jmin+1),1.0/v(imin-1,jmin))
          else
            recpt(imin-1,jmin)=alive(imin-1,jmin+1)+h/v(imin-1,jmin)
          endif
          if(recpt(imin-1,jmin).lt.remin)then
            remin=recpt(imin-1,jmin)
          endif

          if(((imin-2).ge.1).and.(alive(imin-2,jmin+1).ne.100.0))then
            recpt(imin-1,jmin)=recompute2(alive(imin-1,jmin+1),alive(imin-2,jmin+1),1.0/v(imin-1,jmin))
            if(recpt(imin-1,jmin).lt.remin)then
              remin=recpt(imin-1,jmin)
            endif
          endif
        endif
!!!left-down+up
        if(((jmin-1).ge.1).and.(alive(imin-1,jmin-1).ne.100.0))then

          if(alive(imin,jmin-1).ne.100.0)then
            recpt(imin-1,jmin)=recompute2(alive(imin-1,jmin-1),alive(imin,jmin-1),1.0/v(imin-1,jmin))
          else
            recpt(imin-1,jmin)=alive(imin-1,jmin-1)+h/v(imin-1,jmin)
          endif
          if(recpt(imin-1,jmin).lt.remin)then
            remin=recpt(imin-1,jmin)
          endif

          if(((imin-2).ge.1).and.(alive(imin-2,jmin-1).ne.100.0))then
            recpt(imin-1,jmin)=recompute2(alive(imin-1,jmin-1),alive(imin-2,jmin-1),1.0/v(imin-1,jmin))
            if(recpt(imin-1,jmin).lt.remin)then
              remin=recpt(imin-1,jmin)
            endif
          endif

        endif
!!!up-right+left
        if(((imin-2).ge.1).and.(alive(imin-2,jmin).ne.100.0))then

          if(((jmin+1).le.m).and.(alive(imin-2,jmin+1).ne.100.0))then
            recpt(imin-1,jmin)=recompute2(alive(imin-2,jmin),alive(imin-2,jmin+1),1.0/v(imin-1,jmin))			
          else
            recpt(imin-1,jmin)=alive(imin-2,jmin)+h/v(imin-1,jmin) 
          endif
          if(recpt(imin-1,jmin).lt.remin)then
            remin=recpt(imin-1,jmin)
          endif

          if(((jmin-1).ge.1).and.(alive(imin-2,jmin-1).ne.100.0))then
            recpt(imin-1,jmin)=recompute2(alive(imin-2,jmin),alive(imin-2,jmin-1),1.0/v(imin-1,jmin))
            if(recpt(imin-1,jmin).lt.remin)then
              remin=recpt(imin-1,jmin)
            endif
          endif

        endif
!!!8-ok
        recpt(imin-1,jmin)=remin
                                   
        if(isbandwt.eq.100.0)then 
          call bandw_insert(imin-1,jmin,recpt(imin-1,jmin),bandwx,bandwz,bandwt,isroot,csize,msize) 
        else
          if(recpt(imin-1,jmin).lt.isbandwt)then
!           tj=tj+1
            rcsize=csize+1
            do 21 i=1,rcsize
              if((bandwx(i).eq.(imin-1)).and.(bandwz(i).eq.jmin))then
                bandwt(i)=recpt(imin-1,jmin)
                call filterup(bandwx,bandwz,bandwt,i,msize)
              endif 
21         continue
          else
            recpt(imin-1,jmin)=isbandwt
          endif 
        endif 

      endif
      endif

!!!!!!! Rcompute the valuses of T at neighbor(imin+1,jmin) if it is in Narrow Band or Far Away !!!!!
!!!!!!!!! If the neighbor piont is in Far Away,add it to the set of Narrow Band !!!!!!!!!!!!!!!!!!!!

      if((imin+1).le.n)then
      if((alive(imin+1,jmin).eq.100.0).and.(v(imin+1,jmin).ne.0.0))then

        isbandwt=recpt(imin+1,jmin)
!!!up-left
        if(((jmin-1).ge.1).and.(alive(imin,jmin-1).ne.100.0))then
          recpt(imin+1,jmin)=recompute2(alive(imin,jmin),alive(imin,jmin-1),1.0/v(imin+1,jmin))
        else
          recpt(imin+1,jmin)=alive(imin,jmin)+h/v(imin+1,jmin)
        endif
        remin=recpt(imin+1,jmin)
!!!up-right
        if(((jmin+1).le.m).and.(alive(imin,jmin+1).ne.100.0))then
          recpt(imin+1,jmin)=recompute2(alive(imin,jmin),alive(imin,jmin+1),1.0/v(imin+1,jmin))
          if(recpt(imin+1,jmin).lt.remin)then
            remin=recpt(imin+1,jmin)
          endif
        endif
!!!right-up+down
        if(((jmin+1).le.m).and.(alive(imin+1,jmin+1).ne.100.0))then

          if(alive(imin,jmin+1).ne.100.0)then
            recpt(imin+1,jmin)=recompute2(alive(imin+1,jmin+1),alive(imin,jmin+1),1.0/v(imin+1,jmin))
          else
            recpt(imin+1,jmin)=alive(imin+1,jmin+1)+h/v(imin+1,jmin)
          endif
          if(recpt(imin+1,jmin).lt.remin)then
            remin=recpt(imin+1,jmin)
          endif

          if(((imin+2).le.n).and.(alive(imin+2,jmin+1).ne.100.0))then
            recpt(imin+1,jmin)=recompute2(alive(imin+1,jmin+1),alive(imin+2,jmin+1),1.0/v(imin+1,jmin))
            if(recpt(imin+1,jmin).lt.remin)then
              remin=recpt(imin+1,jmin)
            endif
          endif

        endif
!!!left-up+down
        if(((jmin-1).ge.1).and.(alive(imin+1,jmin-1).ne.100.0))then

          if(alive(imin,jmin-1).ne.100.0)then
            recpt(imin+1,jmin)=recompute2(alive(imin+1,jmin-1),alive(imin,jmin-1),1.0/v(imin+1,jmin))
          else
            recpt(imin+1,jmin)=alive(imin+1,jmin-1)+h/v(imin+1,jmin)
          endif 
          if(recpt(imin+1,jmin).lt.remin)then
            remin=recpt(imin+1,jmin)
          endif  
      
          if(((imin+2).le.n).and.(alive(imin+2,jmin-1).ne.100.0))then
            recpt(imin+1,jmin)=recompute2(alive(imin+1,jmin-1),alive(imin+2,jmin-1),1.0/v(imin+1,jmin))
            if(recpt(imin+1,jmin).lt.remin)then
              remin=recpt(imin+1,jmin)
            endif 
          endif

        endif
!!!down-right+left
        if(((imin+2).le.n).and.(alive(imin+2,jmin).ne.100.0))then

          if(((jmin+1).le.m).and.(alive(imin+2,jmin+1).ne.100.0))then
            recpt(imin+1,jmin)=recompute2(alive(imin+2,jmin),alive(imin+2,jmin+1),1.0/v(imin+1,jmin))	    
          else
            recpt(imin+1,jmin)=alive(imin+2,jmin)+h/v(imin+1,jmin)
          endif
          if(recpt(imin+1,jmin).lt.remin)then
            remin=recpt(imin+1,jmin)
          endif        
  
          if(((jmin-1).ge.1).and.(alive(imin+2,jmin-1).ne.100.0))then
            recpt(imin+1,jmin)=recompute2(alive(imin+2,jmin),alive(imin+2,jmin-1),1.0/v(imin+1,jmin))
            if(recpt(imin+1,jmin).lt.remin)then
              remin=recpt(imin+1,jmin)
            endif 
          endif 

        endif
!!!8-ok 
        recpt(imin+1,jmin)=remin
      
        if(isbandwt.eq.100.0)then 
          call bandw_insert(imin+1,jmin,recpt(imin+1,jmin),bandwx,bandwz,bandwt,isroot,csize,msize) 
        else
          if(recpt(imin+1,jmin).lt.isbandwt)then
!           tj=tj+1
            rcsize=csize+1
            do 22 i=1,rcsize
              if((bandwx(i).eq.(imin+1)).and.(bandwz(i).eq.jmin))then
                bandwt(i)=recpt(imin+1,jmin)
                call filterup(bandwx,bandwz,bandwt,i,msize)
              endif 
22         continue
          else
            recpt(imin+1,jmin)=isbandwt
          endif
        endif
 
      endif
      endif               
         
!!!!!!! Rcompute the valuses of T at neighbor(imin,jmin+1) if it is in Narrow Band or Far Away !!!!!
!!!!!!!!! If the neighbor piont is in Far Away,add it to the set of Narrow Band !!!!!!!!!!!!!!!!!!!!

      if((jmin+1).le.m)then
      if((alive(imin,jmin+1).eq.100.0).and.(v(imin,jmin+1).ne.0.0))then

        isbandwt=recpt(imin,jmin+1)

!!!left-up
        if(((imin-1).ge.1).and.(alive(imin-1,jmin).ne.100.0))then
          recpt(imin,jmin+1)=recompute2(alive(imin,jmin),alive(imin-1,jmin),1.0/v(imin,jmin+1))
        else
          recpt(imin,jmin+1)=alive(imin,jmin)+h/v(imin,jmin+1)
        endif
        remin=recpt(imin,jmin+1)
!!!left-down
        if(((imin+1).le.n).and.(alive(imin+1,jmin).ne.100.0))then
          recpt(imin,jmin+1)=recompute2(alive(imin,jmin),alive(imin+1,jmin),1.0/v(imin,jmin+1))
          if(recpt(imin,jmin+1).lt.remin)then
            remin=recpt(imin,jmin+1)
          endif
        endif
!!!down-left+right
        if(((imin+1).le.n).and.(alive(imin+1,jmin+1).ne.100.0))then

          if(alive(imin+1,jmin).ne.100.0)then
            recpt(imin,jmin+1)=recompute2(alive(imin+1,jmin+1),alive(imin+1,jmin),1.0/v(imin,jmin+1))
          else
            recpt(imin,jmin+1)=alive(imin+1,jmin+1)+h/v(imin,jmin+1) 
          endif
          if(recpt(imin,jmin+1).lt.remin)then
            remin=recpt(imin,jmin+1)
          endif

          if(((jmin+2).le.m).and.(alive(imin+1,jmin+2).ne.100.0))then
            recpt(imin,jmin+1)=recompute2(alive(imin+1,jmin+1),alive(imin+1,jmin+2),1.0/v(imin,jmin+1))
            if(recpt(imin,jmin+1).lt.remin)then
              remin=recpt(imin,jmin+1)
            endif
          endif

        endif
!!!up-left+right
        if(((imin-1).ge.1).and.(alive(imin-1,jmin+1).ne.100.0))then
         
          if(alive(imin-1,jmin).ne.100.0)then
            recpt(imin,jmin+1)=recompute2(alive(imin-1,jmin+1),alive(imin-1,jmin),1.0/v(imin,jmin+1))
          else
            recpt(imin,jmin+1)=alive(imin-1,jmin+1)+h/v(imin,jmin+1)
          endif
          if(recpt(imin,jmin+1).lt.remin)then
            remin=recpt(imin,jmin+1)
          endif

          if(((jmin+2).le.m).and.(alive(imin-1,jmin+2).ne.100.0))then
            recpt(imin,jmin+1)=recompute2(alive(imin-1,jmin+1),alive(imin-1,jmin+2),1.0/v(imin,jmin+1))
            if(recpt(imin,jmin+1).lt.remin)then
              remin=recpt(imin,jmin+1)
            endif
          endif

        endif
!!!right-down+up
        if(((jmin+2).le.m).and.(alive(imin,jmin+2).ne.100.0))then

          if(((imin+1).le.n).and.(alive(imin+1,jmin+2).ne.100.0))then
            recpt(imin,jmin+1)=recompute2(alive(imin,jmin+2),alive(imin+1,jmin+2),1.0/v(imin,jmin+1))
          else
            recpt(imin,jmin+1)=alive(imin,jmin+2)+h/v(imin,jmin+1)
          endif
          if(recpt(imin,jmin+1).lt.remin)then
            remin=recpt(imin,jmin+1)
          endif

          if(((imin-1).ge.1).and.(alive(imin-1,jmin+2).ne.100.0))then
            recpt(imin,jmin+1)=recompute2(alive(imin,jmin+2),alive(imin-1,jmin+2),1.0/v(imin,jmin+1))
            if(recpt(imin,jmin+1).lt.remin)then
              remin=recpt(imin,jmin+1)
            endif
          endif

        endif
!!!8-ok 
        recpt(imin,jmin+1)=remin       
        
        if(isbandwt.eq.100.0)then 
          call bandw_insert(imin,jmin+1,recpt(imin,jmin+1),bandwx,bandwz,bandwt,isroot,csize,msize)
        else
          if(recpt(imin,jmin+1).lt.isbandwt)then
!           tj=tj+1
            rcsize=csize+1
            do 23 i=1,rcsize
              if((bandwx(i).eq.(imin)).and.(bandwz(i).eq.(jmin+1)))then
                bandwt(i)=recpt(imin,jmin+1)
                call filterup(bandwx,bandwz,bandwt,i,msize)
              endif 
23         continue
          else
            recpt(imin,jmin+1)=isbandwt
          endif 
        endif

      endif
      endif
      
!!!!!!! Rcompute the valuses of T at neighbor(imin,jmin-1) if it is in Narrow Band or Far Away !!!!!
!!!!!!!!! If the neighbor piont is in Far Away,add it to the set of Narrow Band !!!!!!!!!!!!!!!!!!!!
  
      if((jmin-1).ge.1)then
      if((alive(imin,jmin-1).eq.100.0).and.(v(imin,jmin-1).ne.0.0))then

        dsh=sqrt(2.0)*h/(2.0*v(imin,jmin-1))

        isbandwt=recpt(imin,jmin-1)

!!!right-up
        if(((imin-1).ge.1).and.(alive(imin-1,jmin).ne.100.0))then
          recpt(imin,jmin-1)=recompute2(alive(imin,jmin),alive(imin-1,jmin),1.0/v(imin,jmin-1))
        else
          recpt(imin,jmin-1)=alive(imin,jmin)+h/v(imin,jmin-1)
        endif
        remin=recpt(imin,jmin-1)
!!!right-down
        if(((imin+1).le.n).and.(alive(imin+1,jmin).ne.100.0))then
          recpt(imin,jmin-1)=recompute2(alive(imin,jmin),alive(imin+1,jmin),1.0/v(imin,jmin-1))
          if(recpt(imin,jmin-1).lt.remin)then
            remin=recpt(imin,jmin-1)
          endif
        endif
!!!down-right+left
        if(((imin+1).le.n).and.(alive(imin+1,jmin-1).ne.100.0))then

          if(alive(imin+1,jmin).ne.100.0)then
            recpt(imin,jmin-1)=recompute2(alive(imin+1,jmin-1),alive(imin+1,jmin),1.0/v(imin,jmin-1))
          else
            recpt(imin,jmin-1)=alive(imin+1,jmin-1)+h/v(imin,jmin-1)
          endif
          if(recpt(imin,jmin-1).lt.remin)then
              remin=recpt(imin,jmin-1)
          endif
  
          if(((jmin-2).ge.1).and.(alive(imin+1,jmin-2).ne.100.0))then  
            recpt(imin,jmin-1)=recompute2(alive(imin+1,jmin-1),alive(imin+1,jmin-2),1.0/v(imin,jmin-1))
            if(recpt(imin,jmin-1).lt.remin)then
              remin=recpt(imin,jmin-1)
            endif
          endif

        endif
!!!up-left+right
        if(((imin-1).ge.1).and.(alive(imin-1,jmin-1).ne.100.0))then

          if(alive(imin-1,jmin).ne.100.0)then
            recpt(imin,jmin-1)=recompute2(alive(imin-1,jmin-1),alive(imin-1,jmin),1.0/v(imin,jmin-1))
          else
            recpt(imin,jmin-1)=alive(imin-1,jmin-1)+h/v(imin,jmin-1)
          endif 
          if(recpt(imin,jmin-1).lt.remin)then
            remin=recpt(imin,jmin-1)
          endif 
 
          if(((jmin-2).ge.1).and.(alive(imin-1,jmin-2).ne.100.0))then
            recpt(imin,jmin-1)=recompute2(alive(imin-1,jmin-1),alive(imin-1,jmin-2),1.0/v(imin,jmin-1))
            if(recpt(imin,jmin-1).lt.remin)then
              remin=recpt(imin,jmin-1)
            endif
          endif

        endif
!!!left-down+up
        if(((jmin-2).ge.1).and.(alive(imin,jmin-2).ne.100.0))then

          if(((imin+1).le.n).and.(alive(imin+1,jmin-2).ne.100.0))then
            recpt(imin,jmin-1)=recompute2(alive(imin,jmin-2),alive(imin+1,jmin-2),1.0/v(imin,jmin-1))
          else
            recpt(imin,jmin-1)=alive(imin,jmin-2)+h/v(imin,jmin-1)
          endif 
          if(recpt(imin,jmin-1).lt.remin)then
            remin=recpt(imin,jmin-1)
          endif

          if(((imin-1).ge.1).and.(alive(imin-1,jmin-2).ne.100.0))then          
            recpt(imin,jmin-1)=recompute2(alive(imin,jmin-2),alive(imin-1,jmin-2),1.0/v(imin,jmin-1))
            if(recpt(imin,jmin-1).lt.remin)then
              remin=recpt(imin,jmin-1)
            endif
          endif

        endif
!!!8-ok 
        recpt(imin,jmin-1)=remin  

        if(isbandwt.eq.100.0)then 
          call bandw_insert(imin,jmin-1,recpt(imin,jmin-1),bandwx,bandwz,bandwt,isroot,csize,msize) 
        else
          if(recpt(imin,jmin-1).lt.isbandwt)then
!           tj=tj+1
            rcsize=csize+1
            do 24 i=1,rcsize
              if((bandwx(i).eq.(imin)).and.(bandwz(i).eq.(jmin-1)))then
                bandwt(i)=recpt(imin,jmin-1)
                call filterup(bandwx,bandwz,bandwt,i,msize)
              endif 
24         continue
          else
            recpt(imin,jmin-1)=isbandwt 
          endif 
        endif
        
      endif
      endif
      
!!!!!!!!!!!!!! Return to the top of Loop (If the narrow is not empty)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
      if(csize.gt.0)goto 100
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%3 Output to the data of Recpt/Narrow/Alive !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      kg=1              !!!!! output czht.txt
      if(kg.eq.1)then
      open(17,file='ttime.txt',status='replace')
      write(17,35)((alive(i,j),j=1,m),i=1,n)
      close(17)
!!!!!important!!!!!  
35    format(1000f10.5) 
      endif
      end

     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!remove
      subroutine bandw_insert(insertz,insertx,insertt,bandwx,bandwz,bandwt,isroot,csize,msize)    
     
      integer insertz,insertx,csize,msize,isroot
      real insertt
      integer,dimension(msize)::bandwx,bandwz
      real,dimension(msize)::bandwt 
     
      if(isroot.eq.1)then       
        bandwz(1)=insertx
        bandwx(1)=insertz
        bandwt(1)=insertt    
        csize=csize+1
        isroot=0
        call filterdown(bandwx,bandwz,bandwt,csize,msize)
      else
        bandwz(csize+1)=insertx
        bandwx(csize+1)=insertz               
        bandwt(csize+1)=insertt
        csize=csize+1
        call filterup(bandwx,bandwz,bandwt,csize,msize)
      endif
      
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!output
        imin=bandwx(1)
        jmin=bandwz(1)
        tmin=bandwt(1)
!number of point in the band
        csize=csize-1
        isroot=1
     
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!升序排列
      subroutine filterdown(bandwx,bandwz,bandwt,csize,msize)

      integer i,j,csize,msize,tempx,tempz
      real tempt
      integer,dimension(msize)::bandwx,bandwz
      real,dimension(msize)::bandwt  
     
      i=1
      j=2*i
      tempx=bandwx(i)
      tempz=bandwz(i)
      tempt=bandwt(i)
101   if(j.le.csize)then       
        if(((j+1).le.csize).and.(bandwt(j).gt.bandwt(j+1)))j=j+1
          if(tempt.gt.bandwt(j))then
              bandwx(i)=bandwx(j)
              bandwz(i)=bandwz(j)
              bandwt(i)=bandwt(j)
              i=j
              j=2*i
           if(j.le.csize)goto 101
           endif
       endif        
       bandwx(i)=tempx
       bandwz(i)=tempz
       bandwt(i)=tempt 

      end  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!降序排列
      subroutine filterup(bandwx,bandwz,bandwt,csize,msize)

      integer i,j,csize,msize,tempx,tempz
      real tempt
      integer,dimension(msize)::bandwx,bandwz
      real,dimension(msize)::bandwt


      j=csize
      i=j/2
      tempx=bandwx(j)
      tempz=bandwz(j)
      tempt=bandwt(j)
102   if(j.gt.1)then
        if(bandwt(i).gt.tempt)then
           bandwx(j)=bandwx(i)
           bandwz(j)=bandwz(i)
           bandwt(j)=bandwt(i)
           j=i
           i=j/2
           if(j.gt.1)goto 102
        endif 
      endif
      bandwx(j)=tempx
      bandwz(j)=tempz
      bandwt(j)=tempt 

      end
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function recompute1(x,z)
      real recompute1,x,z,h
      h=10
      recompute1=x+z*h
      end

      function recompute2(x,y,z)
      real recompute2,x,y,z,h
      h=10
      if((x.gt.y).and.(x.lt.(y+sqrt(2.0)*z*h/2.0)))then
        recompute2=x+sqrt((z*h)**2-(x-y)**2)
      else
        if(x.le.y)then
          recompute2=x+z*h
        else
          recompute2=y+sqrt(2.0)*z*h
        endif      
      endif
      end

      


      

    

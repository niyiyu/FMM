        program niyiyu
        
        integer i,j
        parameter(n=1,m=1000)
        real,dimension(1,m)::v

        v=1
        write(*,*) "generation of z0(height)"
        open(15,file='z0.txt',status='replace')
        write(15,13)((v(i,j),j=1,m),i=1,n)
13      format(1000f10.3)
        stop
        end

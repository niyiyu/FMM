        program niyiyu
        
        integer i,j
        parameter(n=1000,m=1000)
        real,dimension(n,m)::v

        v=1.0e+03
        write(*,*) "generation of v"
        open(15,file='v.txt',status='replace')
        write(15,13)((v(i,j),j=1,m),i=1,n)
13      format(1000f10.3)
        stop
        end

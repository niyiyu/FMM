load /Users/niyiyu/documents/fortran/fmm/main/travel_time.txt;
[n,m]=size(travel_time);
travel_time=50*ones(700,1000);
isband=zeros(n,m);
for i=2:n-1
    for j=2:m-1
        if((travel_time(i,j)~=50)&&((travel_time((i-1),(j-1))==50)||(travel_time((i+1),(j-1))==50)||(travel_time((i-1),(j+1))==50)||(travel_time((i+1),(j+1))==50)))
            isband(i,j)=1;
        end
    end
end
dlmwrite('/Users/niyiyu/documents/fortran/fmm/main/band.txt',isband,' ');


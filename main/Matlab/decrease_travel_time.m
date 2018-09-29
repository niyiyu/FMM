load /Users/niyiyu/documents/fortran/fmm/main/travel_time.txt;
[n,m]=size(travel_time);
for i=1:n
    for j=1:m
        if travel_time(i,j)==50
            travel_time(i,j)=NaN;
        end
    end
end
travel_time=imresize(travel_time,0.5,'bilinear');
[n,m]=size(travel_time);
for i=1:n
    for j=1:m
        if isnan(travel_time(i,j))
            travel_time(i,j)=50;
        end
    end
end
isband=zeros(n,m);
for i=2:n-1
    for j=2:m-1
        if((travel_time(i,j)~=50)&&((travel_time((i-1),(j-1))==50)||(travel_time((i+1),(j-1))==50)||(travel_time((i-1),(j+1))==50)||(travel_time((i+1),(j+1))==50)))
            isband(i,j)=1;
        end
    end
end
dlmwrite('/Users/niyiyu/documents/fortran/fmm/main/band.txt',isband,' ');
dlmwrite('/Users/niyiyu/documents/fortran/fmm/main/travel_time.txt',travel_time,' ');
clear all
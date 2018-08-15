load ttime2.txt;
[n,m]=size(ttime2);
n=ceil(n/2);
m=ceil(m/2);
ttime=50*ones(700,1000);
band=zeros(700,1000);
for i=1:n
    for j=1:m
        ttime(i+24,j+374)=ttime2(2*i-1,2*j-1);
    end
end
for i=2:699
    for j=2:999
        if((ttime(i,j)~=50)&&((ttime((i-1),(j-1))==50)||(ttime((i+1),(j-1))==50)||(ttime((i-1),(j+1))==50)||(ttime((i+1),(j+1))==50)))
            band(i,j)=1;
        end
    end
end
dlmwrite('ttime.txt',ttime,' ');
dlmwrite('band.txt',band,' ');


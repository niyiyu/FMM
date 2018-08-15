si=n/2;sj=m/2;
A=zeros(n,m);B=zeros(n,m);
for i=1:n
    for j=1:m
        A(i,j)=abs(h*sqrt((i-si)^2+(j-sj)^2)/1000-ttime(i,j))/ttime(i,j);%/ttime(i,j)
        B(i,j)=h*sqrt((i-si)^2+(j-sj)^2)/1000;
    end
    i
end
A(si,sj)=NaN;

%image(A,'CDataMapping','scaled');colorbar;

%A(si,sj)=0;
%delta=sum(sum(A));
%digits(5);
%delta=vpa(delta/(m*n))
mesh(A);colorbar
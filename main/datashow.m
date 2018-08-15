%clear all;
%clc;

h=10;
load v.txt;
[n,m]=size(v);
load z0.txt;

load ttime.txt;

for i=1:n
    for j=1:m
        if ttime(i,j)==50
            ttime(i,j)=NaN;
        end
    end
end

imagesc(v);
figure(gcf);
hold on
plot(z0,'DisplayName','z0','YDataSource','z0');figure(gcf);hold on 
contour(ttime,20);hold on 

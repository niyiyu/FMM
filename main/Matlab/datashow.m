load /Users/niyiyu/documents/fortran/fmm/main/v.txt;
load /Users/niyiyu/documents/fortran/fmm/main/travel_time.txt;
[n,m]=size(v);
load z0.txt;



for i=1:n
    for j=1:m
        if travel_time(i,j)==50
            travel_time(i,j)=NaN;
        end
    end
end

imagesc(v);
figure(gcf);
hold on
plot(z0,'DisplayName','z0','YDataSource','z0');figure(gcf);hold on 
contour(travel_time,20);hold on 
clear i j m n z0
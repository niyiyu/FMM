load v.txt;
[n,m]=size(v);
si=112;
sj=500;
is=si-floor(n/8);
ie=si+floor(n/8);
js=sj-floor(m/8);
je=sj+floor(m/8);

n=2*(ie-is)+1;
m=2*(je-js)+1;
v2=zeros(n,m);
for i=0:(ie-is)
    for j=0:(je-js)
        v2((2*i+1),(2*j+1))=v((is+i),(js+j));
    end
end
for i=1:(ie-is)
    for j=0:(je-js)
        v2((2*i),(2*j+1))=(v2((2*i-1),(2*j+1))+v2((2*i+1),(2*j+1)))/2;
    end
end 
for i=0:(ie-is)
    for j=1:(je-js)
        v2((2*i+1),(2*j))=(v2((2*i+1),(2*j-1))+v2((2*i+1),(2*j+1)))/2;
    end
end
for i=1:(ie-is)
    for j=1:(je-js)
        v2((2*i),(2*j))=(v2((2*i-1),(2*j-1))+v2((2*i-1),(2*j+1))+v2((2*i+1),(2*j-1))+v2((2*i+1),(2*j+1)))/4;
    end
end
dlmwrite('v2.txt',v2,' ');
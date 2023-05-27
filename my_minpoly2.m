function [Flag,v1] = my_minpoly2(A,n)
%Gaussian elimination
m=n+1;
s=[1;randi([0,1],n-2,1);1];
a=zeros(n,m);
p=2;
Flag=true;
a(:,1)=s;
min_poly=zeros(1,m);
min_poly(m)=1;
for i=1:n
    s=A*s;
    s=mod(s,p);
    a(:,i+1)=s;   
end

for i=1:n-1
    k=find(a(i+1:n,i));
    k=k+i;
    if a(i,i)==1
        a(k,i:m)=xor(a(k,i:m),a(i,i:m));
    else
        if isempty(k)
            Flag=false;
            break
        end
        a(i,i:m)=xor(a(k(1),i:m),a(i,i:m));
        a(k,i:m)=xor(a(k,i:m),a(i,i:m));
    end
end 
if a(n,n)==0
    Flag=false;
end


if Flag==true
    min_poly(m-1)=a(n,m);
    for i=1:n-1
        min_poly(m-1-i)=mod(a(n-i,m)-a(n-i,n-i+1:n)*min_poly(m-i:n).',2);
    end
    
end
  v1=flip(min_poly);
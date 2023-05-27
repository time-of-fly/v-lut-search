function [result] = new_primCheck_u50(f)


p=2;
n=length(f)-1;

r=(p^n-1)/(p-1);

result=0;

if(f(n+1)==0)
     result=4;
     return
end

if(mod(sum(f),2)==0)
    result=6;
    return
end
    z_res=zeros(1,n+1);
    z_sq=zeros(1,p*n);
    z_res(n+1)=1;
    z_sq(p*n-p)=1;
    row=zeros(n,n);
    row(1,1)=1;
for k=1:n-1
    [~,z_res]=deconv(conv(z_res,z_sq),f);
    z_res=mod(z_res,p);
    z_res=z_res(end-n+1:end);
    row(k+1,:)=flip(z_res);
end

m_nullity=n-rank(row.'-eye(n));
if m_nullity>=2
    result=8;
    return
end

factorization=factor(r);
uni_factor=unique(factorization);
s_max=length(uni_factor);

power_table=zeros(n,n+1);
z_pt=zeros(1,n+1);
z_pt(n)=1;

for k=1:n-1
    power_table(k,1:n+1)=z_pt;
    [~,z_pt]=deconv(conv(z_pt,z_pt),f);
    z_pt=mod(z_pt,p);
    z_pt=z_pt(end-n:end);
end
power_table(n,1:n+1)=z_pt;

x_r=polypow(r,f,p,n,power_table);

if ~isequal(x_r,[zeros(1,n) 1])
    result=10;
    return
end

for s=1:s_max
    i=uint64(uni_factor(s));
    tmp=zeros(n+1);
    if mod(p-1,i)>0
        tmp=polypow(r/i,f,p,n,power_table);
        tmp(end)=0;
    end
    if tmp==0
        result=i; 
        return
    end
end


end


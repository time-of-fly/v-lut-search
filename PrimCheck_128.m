function [result] = PrimCheck_128(f,t2_128)
%PRIMPCHECK_U50 此处显示有关此函数的摘要
%   此处显示详细说明
%p2=hexToBinaryVector('100000005B');

p=2;
n=length(f)-1;
% n_per=vpi(n);
% p_per=vpi(p);
% r= (p_per^n_per-1)/(p_per-1);
r=(p^vpi(n)-1)/(p-1);
%str_r=string(r);
result=0;
%result=[0 0 0];%result=0 if primitive
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
    z_res=mod(z_res,p);%n+1 row
    z_res=z_res(end-n+1:end);
    row(k+1,:)=flip(z_res);
end
%m_mullity=dim(null(row-I));
m_nullity=n-rank(row.'-eye(n));
if m_nullity>=2
    result=8;
    return
end


uni_factor=[3, 5, 17, 257, 641, 65537, 274177, 6700417, 67280421310721];
s_max=length(uni_factor);

power_table=zeros(n,n+1);%only works if p==2
z_pt=zeros(1,n+1);
z_pt(n)=1;

for k=1:n-1
    power_table(k,1:n+1)=z_pt;
    [~,z_pt]=deconv(conv(z_pt,z_pt),f);
    z_pt=mod(z_pt,p);%n+1 row
    z_pt=z_pt(end-n:end);
end
power_table(n,1:n+1)=z_pt;

x_r=polypow128(r,f,p,n,power_table,t2_128);
%x_r(end)=0;
if ~isequal(x_r,[zeros(1,n) 1])
    result=10;
    return
end

for s=1:s_max
    i=vpi(uni_factor(s));
    tmp=zeros(n+1);
    if mod(p-1,i)>0
        tmp=polypow128(r/i,f,p,n,power_table,t2_128);
        tmp(end)=0;
    end
    if tmp==0
        result=i; 
        return
    end
end


end


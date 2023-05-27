num_l=12;

fullrank=zeros(1,20);
prim_num=zeros(1,20);
for l_count=1:20
    lfsr_lim=l_count+23;
    for batch=1:5000
        A=zeros(lfsr_lim);
        for i=1:lfsr_lim-num_l
            A(i+num_l,i)=1;
        end
        r_tmp2=randperm(num_l);
        for i=1:num_l
            

            r_tmp=randperm(lfsr_lim,5);

            A(i,r_tmp)=1;
            A(i,r_tmp2(i)+lfsr_lim-num_l)=1;
        end
        
        [poly_Flag,m_poly]=my_minpoly2(A,lfsr_lim);
        
        if poly_Flag
            
            fullrank(l_count)=fullrank(l_count)+1;
            prt_r=new_primCheck_u50(m_poly);

                
            
            if prt_r==0
                prim_num(l_count)=prim_num(l_count)+1;
            end
            
        end
    end
end


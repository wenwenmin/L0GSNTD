% Parameter update function
function [core,var,LC,L,btc,bt]=L0GSNTDparaupdate(core,var,coreze,varze,num,ngmar,coreaa,aa,r,coreK,varK,wk,L,LC,LK,LCK,alpha,Lapk,btmax)
    btc=min(wk,btmax);
    core=core+btc*(core-coreK);
    [V,LC]=gradcore(core,var,ngmar,r,num);
    core=PROX(coreze,V,coreaa);

    for j=1:num
    bt(j)=min(wk,btmax);
    if(j<num)
    var{j}=var{j}+bt(j)*(var{j}-varK{j});
    [V,L(j)]=gradIBPL(core,var,ngmar,r,j,num,Lapk,alpha(j));
    var{j}=PROX(varze{j},V,aa(j));
    else
    var{j}=var{j}+bt(j)*(var{j}-varK{j});
    [V,L(j)]=gradIBPL(core,var,ngmar,r,j,num,Lapk,alpha(j));
    var{j}=PROXn1(V);
    end
    end
end


function [V,L]=gradIBPL(core,var,ngmar,r,n,num,Lapk,alpha)
core=tensor(core);
index=1:num;
index(n)=[];
coreg=ttm(core,var,index);
tempB=double(tenmat(coreg,n));
temp=tempB*tempB';
Xn=double(tenmat(ngmar,n));
U=var{n}*temp-Xn*tempB'+alpha*Lapk{n}*var{n};
L=norm(temp,'fro')+alpha*norm(Lapk{n},'fro');
V=var{n}-1/(r*L)*U;
end



function x=PROXn1(x)
x=double(x);
x(x<0)=0;
end

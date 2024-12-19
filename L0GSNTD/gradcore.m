%% Core's gradient calculation function
function [V,LC,U]=gradcore(core,var,ngmar,r,num)
ngmar=tensor(ngmar);
core=tensor(core);
Vc=core;
Vx=ngmar;
temp=1;

for i=1:num
    vart=var{i}'*var{i};
    Vc=ttm(Vc,vart,i);
    Vx=ttm(Vx,var{i}',i);
    temp=temp*norm(vart,'fro');
end
U=Vc-Vx;
LC=temp;
V=core-1/(r*temp)*U;

end
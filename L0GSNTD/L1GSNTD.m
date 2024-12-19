%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions

function [core,var,loss,timerun]=L1GSNTD(core,var,ngmar,maxiteropt,stopindex,r,lamda,btmax)
%% initialization algorithm
loss=[];
timerun=[0];
num=length(size(ngmar));
LK=zeros(1,num);
LCK=0;
L=ones(1,num);
LC=1;
tk=1;

Lapk=LLaplace(ngmar);
varK=var;
coreK=core;
wk=(tk-1)/(tk);

returnloss=norm(tensor(ngmar));
for i=1:num
    alpha(i)=0;
end
alpha(num)=returnloss^2/(2*norm(Lapk{num},'fro')^2);
% alpha(num)=1;
loss(1)=computeloss(ngmar,core,var,lamda,alpha,Lapk);


rate=[0];
t1=clock;


for i=1:maxiteropt
%% update parameters
fprintf("%d\n",i);

LtempK=LK;
LCtempK=LCK;
vartempK=varK;
coretempK=coreK;    
varK=var;
coreK=core;
LK=L;
LCK=LC;    
[core,var,LC,L,btc,bt]=L1GSNTDupdate(core,var,num,ngmar,r,coretempK,vartempK,wk,L,LC,LtempK,LCtempK,lamda,alpha,Lapk,btmax);
loss(i+1)=computeloss(ngmar,core,var,lamda,alpha,Lapk);

check=0;
for j=1:num
    check=check+norm(var{j}-varK{j},'fro')^2;
end
check=check+norm(tensor(core-coreK))^2;
rho=min([L,LC])/(1e+10);
%% Judging whether to extrapolate
if(loss(i+1)>=loss(i)-rho*check)
    var=varK;
    core=coreK;
    L=LK;
    LC=LCK;
    [core,var,LC,L,btc,bt]=L1GSNTDupdate(core,var,num,ngmar,r,coretempK,vartempK,0,L,LC,LtempK,LCtempK,lamda,alpha,Lapk,btmax);
    loss(i+1)=computeloss(ngmar,core,var,lamda,alpha,Lapk);
end


%% Check if termination condition is met

fprintf("L1GSNTD\n");

fprintf("nonzero:%d\n",nnz(core));  

check1=norm(tensor(core));
check2=norm(tensor(coreK));
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}));  
    check1=check1+norm(var{j}-varK{j},'fro');
    check2=check2+norm(varK{j},'fro');
end

bts{i}=[btc,bt];
t2=clock;
timerun(i+1)=etime(t2,t1);
% Res=abs(loss(i+1)-loss(i));
Res=check1/check2;
fprintf("cri：%d\n",Res);
stop=stopcheck(Res,timerun,stopindex);
if(stop==1)
    fprintf("Number of terminations：%d\n",i);
    pause(4);
    break;
end

tk=(1+sqrt(1+4*tk^2))/2;
wk=(tk-1)/(tk);



end
end


function loss=computeloss(ngmar,core,var,lamda,alpha,Lapk)
    loss=compute(core,var,ngmar)+lamda/2*norm(tensor(core));
    for i=1:length(size(ngmar))
    loss=loss+alpha(i)/2*trace(var{i}'*Lapk{i}*var{i});
    end
end

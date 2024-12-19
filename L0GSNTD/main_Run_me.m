clearvars -except 
clc
warning('off');
%% Parameter.
%   index     : The dataset to be used, when index=1, use ConcreteCrackImages dataset
%               when index=2, use orlraws10P dataset. 
%   r         : Step factor
%   maxiteropt: Maximum iteration alloted to the method
%   trigger   : Whether to enable the indicator array of each method, where
%               when 1∈trigger, enable the 𝓁0-SNTD method
%               when 2∈trigger, enable the 𝓁1-GSNTD method 
%               when 3∈trigger, enable the 𝓁0-GSNTD method 
%   percent   ：The proportion of non-zero elements allowed in each decomposition matrix
%   percore   ：The proportion of non-zero elements allowed in the core tensor
%   alphat    : The graph regularization parameter of 𝓁0-GSNTD and 𝓁1-GSNTD
%   stopindex : The indicator of the stop condition.  
%               To set the specific termination condition, see the 'stopcheck' function for details.  
%               The default termination condition is: ϵ<1e-6 or maxiteropt>12000 
%% Display
%   nonzero   ：The number of non-zero elements in each component.
%   Rel       : The difference in the variable value between two iterations.

%% Parameter settings
rng('shuffle')
index=2;
r=1.01;
maxiteropt=12000;
trigger=[1,2,3];
percent=[0.4,0.6,0.5,0.7];
percore=0.4;
alphat=1;
stopindex=4;


%% Select dataset
[ngmar,R,Rdims,y]=readfile(index);
num=length(size(ngmar));
N=R;



for j=1:5
%% Init
var=[];
core=[];
for i=1:num
    var{i}=rand(size(ngmar,i),Rdims(i));
end
core=tenrand(Rdims);
core=tensor(core);

for i=1:length(size(ngmar))
    a(i)=round(Rdims(i)*size(ngmar,i)*percent(i));
end
aa=a;
coreaa=ceil(prod(size(core))*percore);




%% Solving
for i=1:length(trigger)
[datas{i},varss{i}]=ALGOchoose(core,var,ngmar,coreaa,aa,maxiteropt,Rdims,trigger(i),stopindex,r,alphat);
end
vart=var;
vart{num+1}=core;
datas{length(trigger)+1}=varss;
datas{length(trigger)+2}=vart;
for i=1:length(trigger)
    vars1=datas{length(trigger)+1};
    vars11=vars1{i};
    vartemp=vars11{end};
    [acc(j,i),rdx(j,i),NMIs(j,i)]=clustermeans(vartemp{end},N,y);
end
    datass{j}=datas;

end




% Drawing
accmean=mean(acc)
rdxmean=mean(rdx)
nmimean=mean(NMIs)
plt0=plotplt(trigger,acc,rdx,NMIs); plt0=plotplt(trigger,acc,rdx,NMIs);



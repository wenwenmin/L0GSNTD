
function x=PROXL1(U,lamda,tao)
U=double(U);
indbig=find(U>tao*lamda);
indsmall=find(U<-tao*lamda);
x=zeros(size(U));
x(indbig)=U(indbig)-tao*lamda;
x(indsmall)=U(indsmall)+tao*lamda;
x(x<0)=0;
end

function [X,IX]=PROX(X,U,a)
S=U(:);
[~,IX] = sort(S(:),'descend');
if(a>length(IX))
    a=length(IX);
end
X(IX(1:a))=U(IX(1:a));
X(X<0)=0;
end


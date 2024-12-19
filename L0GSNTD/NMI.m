function nmi = NMI(labels1, labels2)
% if(size(labels2,2)>1)
%     [~,labels2]=max(labels2,[],2);
% end
    C = confusionmat(labels1, labels2);
    N = sum(C(:));
    pi = sum(C, 2) / N;
    uj = sum(C, 1) / N;
    piuj = C / N;
    % Add a small constant to avoid log2(0)
    piuj = piuj +1e-10;
    I = sum(sum(piuj .* log2(piuj ./ (pi * uj+1e-10))));
    Hpi = -sum(pi .* log2(pi + 1e-10));
    Huj = -sum(uj .* log2(uj + 1e-10));
    nmi = 2 * I / (Hpi + Huj);
end
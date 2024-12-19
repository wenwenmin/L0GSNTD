%   R         : Rank of Tensor Tucker Decomposition

function [ngmar,R,Rdims,tlabel]=readfile(i)
  if(i==1) 
   E=load('.\Data\ConcreteCrackImages.mat');
   ngmar= double(E.tensor_var);
   ngmar=squeeze(ngmar);
   label= double(E.label);
    M1=tenmat(ngmar,length(size(ngmar)));
    M11=double(M1)';
    M11=normalize(M11,'range');
    M11=M11';
    rdims =M1.rdims; 
    cdims =M1.cdims; 
    tsize = size(ngmar);  
    M10=tenmat(M11,rdims,cdims,tsize);
    ngmar=tensor(M10);
   tlabel=double(label);
    R=length(unique(tlabel));
     if(find(tlabel==0)~=0)
            tlabel=tlabel(:,1)+1;
     end
    len=100;
    labelu=unique(tlabel);  
    label=[];
    id=[];
     for j=1:R
         temp=find(tlabel==labelu(j));
         id=[id;temp(1:len)];
         label=[label,ones(1,len)*j];
     end
     tlabel=tlabel(id);
     ngmart=ngmar(:,:,:,id);
     ngmar=ngmart;
     tsize=size(ngmar);
     Rdims=[ceil(tsize(1:end-1)/2),R];




   elseif(i==2)    
   E=load('.\Data\orlraws10P.mat');
   ngmar= double(E.X);
   ngmar=reshape(ngmar',92,112,size(ngmar,1));
   label= double(E.Y);
    M1=tenmat(ngmar,3);
    M11=double(M1)';
    M11=normalize(M11,'range');
    M11=M11';
    rdims =M1.rdims; 
    cdims =M1.cdims; 
    tsize = size(ngmar);  
    M10=tenmat(M11,rdims,cdims,tsize);
    ngmar=tensor(M10);
   tlabel=double(label);
    R=length(unique(tlabel));
     tsize=size(ngmar);
     Rdims=[ceil(tsize(1:end-1)/2),R];

  end

end



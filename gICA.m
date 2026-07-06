function [Zs,Zt,SSt,SSs] = gICA(Y,nS,rdim1,rdim2)

FY = [];
for j=1:nS
    [F{j},~, ~]= svds(Y{j},rdim1);
    Xs= F{j}'*Y{j};
    FY= [FY;Xs];
end
N = size(Y{1},1);

[G,~,~] = svds(FY,rdim2);
Ss = G'*FY;
[SSs,A,~] = fastica(Ss, 'numOfIC', rdim2,'approach','symm', 'g', 'tanh','verbose', 'off'); 

for j=1:nS
    Zs{j} = inv(A)*G(rdim1*(j-1)+1:rdim1*j,:)'*F{j}'*Y{j};
    Zt{j} = F{j}*G(rdim1*(j-1)+1:rdim1*j,:)*A;
    for i=1:rdim2
        tmpZt(:,i,j)= Zt{j}(:,i);
    end
end

for i=1:rdim2
    [U(:,i),~,~] = svds(reshape(tmpZt(:,i,:),N,nS),1);
end
SSt = U;








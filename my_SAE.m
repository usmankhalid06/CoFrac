function [Wdec,Z] = my_SAE(Y,K,kk,epochs,lr,batch)
    [d,N]=size(Y);
    Wenc=0.01*randn(K,d); Wdec=Wenc'; benc=zeros(K,1); bdec=mean(Y,2);
    b1=0.9;b2=0.999;e0=1e-8;t=0; mWe=0;vWe=0;mWd=0;vWd=0;mbe=0;vbe=0;mbd=0;vbd=0;
    for ep=1:epochs
        idx=randperm(N);
        for s=1:batch:N
            bi=idx(s:min(s+batch-1,N)); xb=Y(:,bi); B=numel(bi);
            h=Wenc*xb+benc; hr=max(h,0);
            sorted=sort(hr,1,'descend'); thr=sorted(min(kk,K),:);
            mask=(hr>=thr)&(hr>0); z=hr.*mask;
            xhat=Wdec*z+bdec; dxhat=(2/B)*(xhat-xb);
            dWdec=dxhat*z'; dbdec=sum(dxhat,2);
            dz=Wdec'*dxhat; dh=dz.*mask.*(h>0);
            dWenc=dh*xb'; dbenc=sum(dh,2); t=t+1;
            [Wenc,mWe,vWe]=adam(Wenc,dWenc,mWe,vWe,lr,b1,b2,e0,t);
            [Wdec,mWd,vWd]=adam(Wdec,dWdec,mWd,vWd,lr,b1,b2,e0,t);
            [benc,mbe,vbe]=adam(benc,dbenc,mbe,vbe,lr,b1,b2,e0,t);
            [bdec,mbd,vbd]=adam(bdec,dbdec,mbd,vbd,lr,b1,b2,e0,t);
        end
    end
    H=max(Wenc*Y+benc,0); sorted=sort(H,1,'descend'); thr=sorted(min(kk,K),:);
    Z=H.*((H>=thr)&(H>0));
end

function [W,m,v]=adam(W,g,m,v,lr,b1,b2,e,t)
    m=b1*m+(1-b1)*g; v=b2*v+(1-b2)*(g.^2);
    W=W-lr*(m/(1-b1^t))./(sqrt(v/(1-b2^t))+e);
end
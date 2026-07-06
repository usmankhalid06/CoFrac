function [E,D,X] = sgBACES(Y,nIter, Kc, Ks, M,lambda1,lambda2)
% Common + individual dictionary learning across M layers (NON-robust).
% Model per layer j:  Y{j} ~= Dc*Xc (shared) + Ds{j}*Xs{j} (layer-specific).
% FIX vs your paste: the two my_sBACES calls dropped the stray nIter arg so the
% call (5 args) matches the function signature (Y,D,X,Kc,spa). Now it runs.
    N = size(Y{1},1);                        % = feature dim d (unused, kept as-is)
    for j=1:M
        d = size(Y{j},1);
        if j == 1                            % init common dictionary ONCE
            Dc = normc(randn(d,Kc)); 
            Dc = Dc*diag(1./sqrt(sum(Dc.*Dc)));
            Xc = zeros(Kc,size(Y{j},2));
        end
        Ds{j} = normc(randn(d,Ks)); 
        Ds{j} = Ds{j}*diag(1./sqrt(sum(Ds{j}.*Ds{j})));
        Xs{j} = zeros(Ks,size(Y{j},2));
    end
    fprintf('Iteration:     ');
    for iter = 1:nIter
        fprintf('\b\b\b\b\b%5i',iter);
        Dcl = Dc;
        %% individual level (per layer)
        for j=1:M
            R{j} = Y{j}-Dc*Xc;
            [Ds{j},Xs{j}] = my_sBACES(R{j},Ds{j},Xs{j},Ks,lambda2);   % <-- 5 args
            % E(iter,j+1)  = sqrt(trace((Ds{j}-Dsp{j})'*(Ds{j}-Dsp{j})))/sqrt(trace(Dsp{j}'*Dsp{j}));
        end
        %% common level (aggregate across layers via median)
        for j =1:M
            tmpEc(:,:,j) = Y{j}-Ds{j}*Xs{j};
        end
        Ec = sum(tmpEc,3)./M;
        % Ec = median(tmpEc,3);
        [Dc,Xc]= my_sBACES(Ec,Dc,Xc,Kc,lambda1);                      % <-- 5 args
        E(iter,1) = sqrt(trace((Dc-Dcl)'*(Dc-Dcl)))/sqrt(trace(Dcl'*Dcl)+eps);
    end
    Dss=[]; Xss=[];
    for i =1:M
        Dss = [Dss Ds{i}];
        Xss = [Xss;Xs{i}];
    end
    D = [Dc Dss];
    X = [Xc;Xss];
end

function [D,X]= my_sBACES(Y,D,X,Kc, spa)
    for j =1:Kc
        X(j,:) = 0;
        E = Y-D*X;
        xk = D(:,j)'*E;
        thr = spa./abs(xk);
        X(j,:) = sign(xk).*max(0, bsxfun(@minus,abs(xk),thr/2));
        rInd = find(X(j,:));
        if (length(rInd)<1)
            [~,ind] = max(sum((Y-D*X).^2,1));
            D(:,j) = Y(:,ind)/(norm(Y(:,ind))+eps);
            continue
        end
        D(:,j) = E(:,rInd)*X(j,rInd)'./norm(E(:,rInd)*X(j,rInd)'+eps);
    end
end
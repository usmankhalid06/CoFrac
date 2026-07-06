clear; close all; clc; rng(0);

% ============================ CONFIG =======================================
FILENAME='esm2_35M_multilayer.mat'; METAFILE='protein_metadata_12.tsv';
N_SUB=60000; PREP='center';
Kc=100; Ks=30;                     % common atoms, per-layer individual atoms
LAMBDA1=3.8; LAMBDA2=7.6;
NITER=15;
TOPN=200; ENR_THR=3; NLAB=5;       % ENR_THR = cysteine gate (unified to 3x)
PUR_MIN=0.25;                      % purity floor for headline fractionation count
ENR_THR4=3;                        % identity-atom threshold (matches SAE 3x reporting)
NULL_NPERM=500;                    % permutations for the chance-partition null
% ===========================================================================
AA='ACDEFGHIKLMNPQRSTVWY'; clean=@(s) s(isstrprop(s,'alphanum'));

%% ---- load ALL layers + preprocess each separately ----
S=load(FILENAME);
assert(ndims(S.X)==3,'Need the multilayer X = N x d x L.');
Nfull=size(S.X,1); d=size(S.X,2); L=size(S.X,3);
if isfield(S,'layers_used'), layers=double(S.layers_used(:)); else, layers=(1:L)'; end
pidx=double(S.protein_idx(:)); rpos=double(S.residue_pos(:));
accAll=cellstr(string(S.accessions));
if numel(accAll)==Nfull, accRes=accAll;
else
    uids=unique(pidx); amap=containers.Map('KeyType','double','ValueType','char');
    for i=1:numel(uids), if i<=numel(accAll), amap(uids(i))=accAll{i}; end, end
    accRes=cell(Nfull,1);
    for i=1:Nfull, if isKey(amap,pidx(i)), accRes{i}=amap(pidx(i)); else, accRes{i}=''; end, end
end
if N_SUB<Nfull, sel=randperm(Nfull,N_SUB); else, sel=1:Nfull; end
pidxS=pidx(sel); rposS=rpos(sel); accResS=accRes(sel); N=numel(sel);
Y=cell(1,L);
for l=1:L
    Xl=double(S.X(:,:,l)); Xl=Xl-mean(Xl,1);
    Yl=Xl(sel,:)'; Yl=zscore(Yl);
    Y{l}=Yl;
end
fprintf('Loaded %d layers (%s), d=%d, N=%d\n', L, num2str(layers'), d, N);

%% ---- COMMON + INDIVIDUAL dictionary learning (sgBACES, non-robust) ----
[~, Dall, Xall] = sgBACES(Y, NITER, Kc, Ks, L, LAMBDA1, LAMBDA2);
Dc = Dall(:,1:Kc); Xc = Xall(1:Kc,:);
num=0; den=0;
for l=1:L
    Dsj=Dall(:, Kc+(l-1)*Ks+(1:Ks)); Xsj=Xall(Kc+(l-1)*Ks+(1:Ks),:);
    num=num+norm(Y{l}-Dc*Xc-Dsj*Xsj,'fro')^2; den=den+norm(Y{l},'fro')^2;
end
fprintf('CID done. relErr=%.3f | %d common atoms, %d x %d layer-specific\n', sqrt(num/den), Kc, Ks, L);

%% ---- metadata + background + per-residue amino acid ----
T=readtable(METAFILE,'FileType','text','Delimiter','\t'); vn=T.Properties.VariableNames;
fc=@(p) find(~cellfun('isempty',regexpi(vn,p,'once')),1);
ci_acc=fc('Entry'); if isempty(ci_acc), ci_acc=1; end; ci_seq=fc('Sequence'); ci_nam=fc('Protein.*ames');
accs=string(T{:,ci_acc}); seqs=string(T{:,ci_seq}); nams=string(T{:,ci_nam}); [ua,ia]=unique(accs,'stable');
m_seq=containers.Map('KeyType','char','ValueType','char'); m_nam=containers.Map('KeyType','char','ValueType','char');
bgc=zeros(1,20);
for i=1:numel(ua), key=clean(char(ua(i)));
    if ~isempty(key), s=char(seqs(ia(i))); m_seq(key)=s; m_nam(key)=char(nams(ia(i)));
        for j=1:20, bgc(j)=bgc(j)+sum(s==AA(j)); end, end
end
bgfreq=bgc/max(sum(bgc),1);
aaLetter=repmat(' ',1,N);
for jj=1:N
    key=clean(accResS{jj});
    if isKey(m_seq,key), s=m_seq(key); p=rposS(jj);
        if p>=1 && p<=numel(s), aaLetter(jj)=s(p); end
    end
end
isCys=(aaLetter=='C');

%% ---- Swiss-Prot label setup (up top so Part 2b can use it) ----
assert(isfield(S,'labels') && isfield(S,'label_names'), 'Need S.labels and S.label_names.');
LAB=double(S.labels);
if size(LAB,1)~=Nfull && size(LAB,2)==Nfull, LAB=LAB.'; end
labS=LAB(sel,:); Cn=size(labS,2);
lnames=cellstr(string(S.label_names(:)));
if numel(lnames)~=Cn, lnames=arrayfun(@(c)sprintf('label%d',c),1:Cn,'uni',0); end
labBg=mean(labS>0,1)+1e-12;
fprintf('Detected labels: %d residues x %d concepts: %s\n', size(LAB,1), Cn, strjoin(lnames,', '));
dcol = find(~cellfun('isempty',regexpi(lnames,'disulf','once')),1);

%% ---- cysteine atoms per block: COMMON vs each LAYER ----
fprintf('\nCysteine atoms (enr>=%g x) per block:\n', ENR_THR);
[~, cysCommon] = block_cys(Xc, isCys, bgfreq, TOPN, ENR_THR);
fprintf('  COMMON           : %d cysteine atoms  %s\n', numel(cysCommon), mat2str(cysCommon(:)'));
for l=1:L
    Xsj=Xall(Kc+(l-1)*Ks+(1:Ks),:);
    [~, cl]=block_cys(Xsj, isCys, bgfreq, TOPN, ENR_THR);
    fprintf('  layer %2d specific : %d cysteine atoms\n', layers(l), numel(cl));
end

[~,ord]=sort(pidxS); pS=pidxS(ord); accS=accResS(ord);

%% ---- fractionation analysis on the COMMON cysteine atoms ----
if numel(cysCommon)>=2
    ATOMS=cysCommon(:).'; Zs=Xc; Wdec=Dc; nA=numel(ATOMS);
    fprintf('\nCysteine fractionates in the COMMON (cross-layer) dictionary -> %d atoms.\n', nA);

    TOPset=cell(nA,1); Pset=cell(nA,1);
    for sN=1:nA, a=abs(Zs(ATOMS(sN),:)); [~,o]=sort(a,'descend'); t=o(1:min(TOPN,N)); TOPset{sN}=t; Pset{sN}=unique(pidxS(t)); end
    fprintf('Pairwise overlap (residue J / protein J):\n');
    fprintf('%8s',''); for j=1:nA, fprintf('  atom%-3d',ATOMS(j)); end; fprintf('\n');
    for i=1:nA
        fprintf('atom%-4d',ATOMS(i));
        for j=1:nA
            if i==j, fprintf('     --    '); continue; end
            rj=numel(intersect(TOPset{i},TOPset{j}))/numel(union(TOPset{i},TOPset{j}));
            pj=numel(intersect(Pset{i},Pset{j}))/numel(union(Pset{i},Pset{j}));
            fprintf('  %.2f/%.2f',rj,pj);
        end
        fprintf('\n');
    end

    %% ---- Part 2b: COMMON cysteine atoms vs DISULFIDE annotation [added] ----
    if ~isempty(dcol)
        isDisulf=labS(:,dcol)>0; nDisCys=sum(isDisulf(isCys)); cysBg=mean(isDisulf(isCys))+1e-9;
        fprintf('\n=== COMMON CYSTEINE ATOMS vs DISULFIDE (label col %d) ===\n', dcol);
        fprintf('Background: %.1f%% of cysteines disulfide-annotated (n=%d)\n', 100*cysBg, nDisCys);
        if nDisCys<10, fprintf('** WARNING: only %d disulfide cysteines -- UNDERPOWERED. **\n', nDisCys); end
        fprintf('%8s  %10s  %10s  %12s\n','atom','top-cys','%disulf','enr_vs_cysbg');
        for sN=1:nA
            a=abs(Zs(ATOMS(sN),:)); [~,o]=sort(a,'descend'); top=o(1:min(TOPN,N));
            topCys=top(isCys(top)); if isempty(topCys), continue; end
            pD=mean(isDisulf(topCys));
            fprintf('%8d  %10d  %9.1f%%  %11.2fx\n', ATOMS(sN), numel(topCys), 100*pD, pD/cysBg);
        end
        fprintf('(>1x disulfide-biased; <1x disulfide-depleted; a split in the COMMON dict = CROSS-LAYER functional partition)\n');

        %% === what ELSE do the cysteine atoms concentrate? (run after main script) ===
        otherCols = setdiff(1:numel(lnames), dcol);
        fprintf('\n--- cysteine atoms vs ALL concepts (top cysteines only) ---\n');
        fprintf('%6s', 'atom');
        for c = otherCols, fprintf('  %16s', lnames{c}); end
        fprintf('   %12s\n', 'unannot%');
        for sN = 1:numel(ATOMS)
            a = abs(Zs(ATOMS(sN),:)); [~,o] = sort(a,'descend'); top = o(1:min(TOPN,N));
            topCys = top(isCys(top));
            if isempty(topCys), continue; end
            fprintf('%6d', ATOMS(sN));
            anyLab = labS(topCys, dcol) > 0;               % start with disulfide
            for c = otherCols
                frac = mean(labS(topCys,c) > 0);
                bg   = mean(labS(isCys,c) > 0) + 1e-9;     % that concept's rate among ALL cysteines
                fprintf('  %5.1f%%(%4.1fx)', 100*frac, frac/bg);
                anyLab = anyLab | (labS(topCys,c) > 0);
            end
            fprintf('   %11.1f%%\n', 100*mean(~anyLab));
        end
        fprintf(['(each cell = %% of the atom''s top cysteines with that label, and enrichment vs\n' ...
                 ' the cysteine-wide rate; unannot%% = carry NONE of the five concepts)\n']);
    end

    Dn=Wdec(:,ATOMS); Dn=Dn./sqrt(sum(Dn.^2,1)+eps); Ccos=abs(Dn'*Dn);
    fprintf('\nDecoder cosines (reference only):\n');
    for i=1:nA
        fprintf('atom%-4d',ATOMS(i));
        for j=1:nA, if i==j, fprintf('    --  '); else, fprintf('  %.2f',Ccos(i,j)); end, end
        fprintf('\n');
    end

    figure('Color','w','Position',[80 50 1200 220*nA]);
    for sN=1:nA
        at=ATOMS(sN); a=abs(Zs(at,:)); a=a/max(a+eps); a=a(ord);
        subplot(nA,1,sN); stem(1:N,a,'Marker','none','Color',[.2 .5 .35]); hold on;
        ylim([0 1.15]); xlim([1 N]); ylabel('|act|'); title(sprintf('COMMON cysteine atom %d',at));
        up=unique(pS,'stable'); mass=zeros(numel(up),1);
        for i=1:numel(up), mass(i)=sum(a(pS==up(i))); end
        [~,oi]=sort(mass,'descend');
        for r=1:min(NLAB,numel(up))
            p=up(oi(r)); idx=find(pS==p); xc=mean(idx); [~,im]=max(a(idx)); yc=a(idx(im));
            key=clean(accS{idx(1)}); nm=key;
            if isKey(m_nam,key), nm=m_nam(key); if numel(nm)>22, nm=nm(1:22); end, end
            text(xc,min(1.1,yc+0.06),nm,'Rotation',40,'FontSize',7,'Color',[.6 0 .1],'Interpreter','none');
        end
    end
    xlabel('residue (voxel), grouped by protein  \rightarrow');
else
    fprintf('\n<2 common cysteine atoms -> cysteine is LAYER-SPECIFIC, not shared across layers.\n');
end

%% ============ Part 3: layering test on COMMON cysteine atoms (mean) =========
if numel(cysCommon)>=2
    ATOMS=cysCommon(:).'; Zs=Xc;
    focal=ATOMS(1);
    aF=abs(Zs(focal,:)); [~,oF]=sort(aF,'descend'); topF=oF(1:min(TOPN,N));
    structCys=topF(isCys(topF));
    cysOther =setdiff(find(isCys),structCys);
    nonCys   =find(~isCys);
    others=ATOMS(2:end);
    fprintf('\nLayering test (COMMON): do atom-%d structural cysteines (n=%d) load on the others?\n', focal, numel(structCys));
    figure('Color','w','Position',[80 80 max(560*numel(others),560) 360]);
    for sN=1:numel(others)
        at=others(sN); a=abs(Zs(at,:));
        mS=mean(a(structCys)); mC=mean(a(cysOther)); mN=mean(a(nonCys));   % mean
        fprintf('  atom %d:  mean|act| struct-cys=%.4f | other-cys=%.4f | non-cys=%.4f\n', at,mS,mC,mN);
        if mS>2*mN, fprintf('     -> ABOVE non-cys baseline: identity+structure LAYERED on the same residue.\n');
        else,       fprintf('     -> not above non-cys baseline: separate pools, not layered.\n'); end
        subplot(1,numel(others),sN); hold on;
        grp={{structCys,[.85 .1 .1]},{cysOther,[.1 .4 .85]},{nonCys,[.5 .5 .5]}};
        for g=1:3, v=sort(a(grp{g}{1})); y=(1:numel(v))/max(numel(v),1); stairs(v,y,'Color',grp{g}{2},'LineWidth',1.3); end
        xlabel(sprintf('|act| on atom %d',at)); ylabel('cum. frac'); ylim([0 1]); title(sprintf('atom %d',at));
        legend({sprintf('atom-%d struct-cys',focal),'other cys','non-cys'},'Location','southeast','FontSize',7);
    end
end

%% ============ Part 4: ALL-ATOM CENSUS over the COMMON dictionary ============
Zs=Xc; K=Kc;
domAA=repmat(' ',K,1); domEnr=zeros(K,1); domPur=zeros(K,1);
for kk=1:K
    a=abs(Zs(kk,:)); [~,o]=sort(a,'descend'); top=o(1:min(TOPN,N));
    lets=aaLetter(top); lets=lets(lets~=' ');
    if isempty(lets), domAA(kk)='-'; continue; end
    cnt=zeros(1,20); for j=1:20, cnt(j)=sum(lets==AA(j)); end
    frac=cnt/sum(cnt); enr=frac./(bgfreq+1e-9);
    [domEnr(kk),ix]=max(enr); domAA(kk)=AA(ix); domPur(kk)=frac(ix);
end
isIdentity = domEnr>=ENR_THR4 & domAA~='-';
ord_atoms=[];
for j=1:20, idx=find(isIdentity & domAA==AA(j)); [~,si]=sort(domEnr(idx),'descend'); ord_atoms=[ord_atoms; idx(si)]; end
df=find(~isIdentity); [~,sd]=sort(domEnr(df),'descend'); ord_atoms=[ord_atoms; df(sd)];
Hmap=abs(Zs(ord_atoms,ord)); Hmap=Hmap./(max(Hmap,[],2)+eps);
figure('Color','w','Position',[60 50 1300 820]);
imagesc(Hmap.^0.5); colormap(flipud(gray)); cb=colorbar; cb.Label.String='|act| (row-norm, sqrt)';
ylab=cell(K,1);
for r=1:K, at=ord_atoms(r);
    if isIdentity(at), ylab{r}=sprintf('%s  a%d  %0.0fx',domAA(at),at,domEnr(at)); else, ylab{r}=sprintf('-   a%d',at); end
end
set(gca,'YTick',1:K,'YTickLabel',ylab,'FontSize',6,'TickLength',[0 0]);
xlabel('residue (voxel), grouped by protein  \rightarrow'); ylabel('common atom (dominant amino acid)');
nIdent=sum(isIdentity);
title(sprintf('COMMON: all %d atoms by dominant amino acid  [%d identity / %d diffuse]',K,nIdent,K-nIdent));
yline(nIdent+0.5,'r-','LineWidth',1.2);
fprintf('\n=== COMMON ATOM CENSUS (Kc=%d, enr>=%g x) ===\n',K,ENR_THR4);
fprintf('Identity atoms: %d  |  Diffuse: %d\n', nIdent, K-nIdent);
fprintf('\n AA  #atoms  meanEnr  bestEnr  meanPurity\n'); covered=0;
for j=1:20
    idx=find(isIdentity & domAA==AA(j)); if isempty(idx), continue; end; covered=covered+1;
    fprintf('  %s   %3d    %5.1fx   %5.1fx     %.2f\n', AA(j),numel(idx),mean(domEnr(idx)),max(domEnr(idx)),mean(domPur(idx)));
end
fprintf('\nAmino acids with a clean identity atom: %d of 20\n', covered);

%% ============ Part 5: STRUCTURE CENSUS vs Swiss-Prot (label sgBACES) ========
STR_THR=3; MINPOS=5; PUR_THR=0.25;
domCon=zeros(K,1); conEnr=zeros(K,1); conPur=zeros(K,1); conCnt=zeros(K,1);
for kk=1:K
    a=abs(Zs(kk,:)); [~,o]=sort(a,'descend'); top=o(1:min(TOPN,N));
    cnt=sum(labS(top,:)>0,1); fr=cnt/numel(top); en=fr./labBg;
    [conEnr(kk),ic]=max(en); domCon(kk)=ic; conPur(kk)=fr(ic); conCnt(kk)=cnt(ic);
end
isStruct = conEnr>=STR_THR & conCnt>=MINPOS & conPur>=PUR_THR;
fprintf('\n=== COMMON STRUCTURE CENSUS (sgBACES, Kc=%d, enr>=%gx, >=%d hits, purity>=%.2f) ===\n', K, STR_THR, MINPOS, PUR_THR);
fprintf('Structure atoms: %d\n\n concept            #atoms  meanEnr  bestEnr  meanPurity\n', sum(isStruct));
for c=1:Cn
    idx=find(isStruct & domCon==c); if isempty(idx), continue; end
    fprintf('  %-18s  %3d    %5.1fx   %5.1fx     %.3f\n', ...
        lnames{c}, numel(idx), mean(conEnr(idx)), max(conEnr(idx)), mean(conPur(idx)));
end
hidden=find(isStruct & ~isIdentity); [~,sh]=sort(conEnr(hidden),'descend');
fprintf('\nStructure atoms hidden in the diffuse pile (no dominant AA, but functional): %d\n', numel(hidden));
fprintf(' atom  concept            enr    purity\n');
for r=hidden(sh).'
    fprintf('  a%-4d %-18s %5.1fx  %.3f\n', r, lnames{domCon(r)}, conEnr(r), conPur(r));
end
pick=[]; for c=1:Cn, idx=find(isStruct&domCon==c); if ~isempty(idx),[~,b]=max(conEnr(idx)); pick=[pick;idx(b)]; end, end
nP=numel(pick);
if nP>0
    figure('Color','w','Position',[80 50 1200 200*nP]);
    for sN=1:nP
        at=pick(sN); a=abs(Zs(at,:)); a=a/max(a+eps); a=a(ord);
        subplot(nP,1,sN); stem(1:N,a,'Marker','none','Color',[.2 .5 .35]); hold on;
        ylim([0 1.15]); xlim([1 N]); ylabel('|act|');
        title(sprintf('sgBACES common atom %d  ->  %s   (enr %.0fx, purity %.2f)', at, lnames{domCon(at)}, conEnr(at), conPur(at)));
        up=unique(pS,'stable'); mass=zeros(numel(up),1);
        for i=1:numel(up), mass(i)=sum(a(pS==up(i))); end
        [~,oi]=sort(mass,'descend');
        for r=1:min(NLAB,numel(up))
            p=up(oi(r)); ix=find(pS==p); xc=mean(ix); [~,im]=max(a(ix)); yc=a(ix(im));
            key=clean(accS{ix(1)}); nm=key;
            if isKey(m_nam,key), nm=m_nam(key); if numel(nm)>22, nm=nm(1:22); end, end
            text(xc,min(1.1,yc+0.06),nm,'Rotation',40,'FontSize',7,'Color',[.6 0 .1],'Interpreter','none');
        end
    end
    xlabel('residue (voxel), grouped by protein  \rightarrow');
end

%% ===== Part 5b: layer-specific structure-atom counts =====
fprintf('\n=== LAYER-SPECIFIC STRUCTURE ATOMS (Ks=%d, enr>=%gx, >=%d hits, purity>=%.2f) ===\n', ...
        Ks, STR_THR, MINPOS, PUR_THR);
for l=1:L
    Xsj = Xall(Kc+(l-1)*Ks+(1:Ks),:);
    ce=zeros(Ks,1); cc=zeros(Ks,1); cp=zeros(Ks,1);
    for kk=1:Ks
        a=abs(Xsj(kk,:)); [~,o]=sort(a,'descend'); top=o(1:min(TOPN,N));
        cnt=sum(labS(top,:)>0,1); fr=cnt/numel(top); en=fr./labBg;
        [ce(kk),ic]=max(en); cc(kk)=cnt(ic); cp(kk)=fr(ic);
    end
    nStr = sum(ce>=STR_THR & cc>=MINPOS & cp>=PUR_THR);
    fprintf('  layer %2d specific : %d structure atoms (of %d)\n', layers(l), nStr, Ks);
end

%% ============ Part 6+7: FRACTIONATION SWEEP + NULL on COMMON dict [added] ===
fprintf('\n=== COMMON FRACTIONATION SWEEP + CHANCE-PARTITION NULL (Kc=%d) ===\n', Kc);
fprintf('%3s %7s %9s %10s %7s %8s   %s\n','AA','#atoms','obsResJ','nullResJ','z','meanPur','verdict');
fracTable = nan(20,6);
for j = 1:20
    idx = find(isIdentity & domAA==AA(j)); nAj=numel(idx);
    if nAj<1, continue; end
    meanPur=mean(domPur(idx));
    if nAj==1
        fprintf('%3s %7d %9s %10s %7s %8.2f   single atom (no partition)\n', AA(j),nAj,'--','--','--',meanPur);
        fracTable(j,:)=[nAj NaN NaN NaN meanPur 0]; continue;
    end
    TS=cell(nAj,1); pool=[];
    for u=1:nAj, a=abs(Xc(idx(u),:)); [~,o]=sort(a,'descend'); TS{u}=o(1:min(TOPN,N)); pool=[pool TS{u}]; end
    obs=[]; for p=1:nAj, for q=p+1:nAj, obs(end+1)=numel(intersect(TS{p},TS{q}))/numel(union(TS{p},TS{q})); end, end %#ok
    obsJ=mean(obs);
    upool=unique(pool); Tn=min(TOPN,N); nullJ=zeros(NULL_NPERM,1);
    for b=1:NULL_NPERM
        sets=cell(nAj,1);
        for u=1:nAj, r=randperm(numel(upool),min(Tn,numel(upool))); sets{u}=upool(r); end
        nj=[]; for p=1:nAj, for q=p+1:nAj, nj(end+1)=numel(intersect(sets{p},sets{q}))/numel(union(sets{p},sets{q})); end, end %#ok
        nullJ(b)=mean(nj);
    end
    z=(obsJ-mean(nullJ))/(std(nullJ)+eps);
    passNull = obsJ < mean(nullJ)-2*std(nullJ);
    if passNull && meanPur>=PUR_MIN, verd='FRACTIONATES (below-chance disjoint, pure)';
    elseif passNull, verd=sprintf('below-chance but LOW PURITY (%.2f) -> noise',meanPur);
    else, verd='not below chance -> redundant'; end
    fprintf('%3s %7d %9.3f %10.3f %7.1f %8.2f   %s\n', AA(j),nAj,obsJ,mean(nullJ),z,meanPur,verd);
    fracTable(j,:)=[nAj obsJ mean(nullJ) z meanPur passNull];
end
isFrac = fracTable(:,1)>=2 & fracTable(:,6)==1 & fracTable(:,5)>=PUR_MIN;
fprintf('\n>>> COMMON (cross-layer shared): %d of 20 amino acids fractionate (>=2 atoms, below-null, purity>=%.2f)\n', sum(isFrac), PUR_MIN);
fprintf('>>> Amino acids: %s\n', AA(isFrac));

%% ============ export focal COMMON cysteine atom ============
if numel(cysCommon)>=1
    ATOM=cysCommon(1); a=abs(Xc(ATOM,:));
    accClean=cellfun(@(s) clean(s), accResS, 'UniformOutput', false);
    Tcsv=table(string(accClean(:)), rposS(:), a(:), 'VariableNames', {'accession','resnum','activation'});
    writetable(Tcsv, sprintf('common_atom%d_all.csv', ATOM));
    fprintf('\nwrote common_atom%d_all.csv\n', ATOM);
end

%% ---- sparsity (avg nonzeros per residue) ----
tol = 1e-7;
nnzC = mean(sum(abs(Xc)>tol,1));
fprintf('avg nnz/col  COMMON (Kc=%d): %.2f  (density %.3f)\n', Kc, nnzC, nnzC/Kc);

specPerLayer = zeros(1,L);
totPerLayer  = zeros(1,L);
for l=1:L
    Xsj = Xall(Kc+(l-1)*Ks+(1:Ks),:);
    specPerLayer(l) = mean(sum(abs(Xsj)>tol,1));      % layer-specific block only
    totPerLayer(l)  = mean(sum(abs([Xc;Xsj])>tol,1)); % common + layer-specific
end
nnzS = mean(specPerLayer);
nnzT = mean(totPerLayer);
fprintf('avg nnz/col  SPECIFIC (Ks=%d): %.2f  (density %.3f)\n', Ks, nnzS, nnzS/Ks);
fprintf('avg nnz/col  COMMON+SPECIFIC (SAE-comparable L0): %.2f\n', nnzT);





% %% ---- local function ----
% function [enr, atoms] = block_cys(Z, isCys, bgfreq, TOPN, THR)
%     K=size(Z,1); enr=zeros(K,1);
%     for kk=1:K
%         a=abs(Z(kk,:)); [~,o]=sort(a,'descend'); top=o(1:min(TOPN,numel(a)));
%         enr(kk)=mean(isCys(top))/(bgfreq(2)+1e-9);
%     end
%     [sv,si]=sort(enr,'descend'); atoms=si(sv>=THR);
% end
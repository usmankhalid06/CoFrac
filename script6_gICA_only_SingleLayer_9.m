clear; close all; clc; rng(0);

% ============================ CONFIG =======================================
FILENAME='esm2_35M_multilayer.mat'; METAFILE='protein_metadata_12.tsv';
N_SUB=60000; PREP='center';
RDIM1=200; RDIM2=200;      % per-layer reduction (rdim1), group component count (rdim2)
TOPN=200; NLAB=5;
PUR_MIN=0.25;              % purity floor for headline fractionation count
ENR_THR=3;                 % identity-atom threshold (matches SAE 3x reporting)
ENR_CYS=3;                 % cysteine auto-find threshold (unified to 3x)
NULL_NPERM=500;            % permutations for the chance-partition null
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
    Xl=double(S.X(:,:,l)); Xl=Xl-mean(Xl,1);     % center features
    Yl=Xl(sel,:)'; Yl=zscore(Yl);                % d x N, per-residue zscore
    Y{l}=Yl;
end
fprintf('Loaded %d layers (%s), d=%d, N=%d\n', L, num2str(layers'), d, N);

%% ---- TRAIN GROUP ICA across layers ----
% gICA(Y,nS,rdim1,rdim2):
%   SSs (rdim2 x N) : group-common components (shared across layers)  -> Zs (group)
%   SSt (N x rdim2) : consensus spatial map across layers
%   Zs{j}/Zt{j}     : per-layer back-reconstructions (layer-specific projections)
[Zs_layer, Zt_layer, SSt, SSs] = gICA(Y, L, RDIM1, RDIM2);
Kg = size(SSs,1);
fprintf('Group ICA done. %d common components across %d layers.\n', Kg, L);

% group decoder (d x Kg): least-squares map from components to a representative
% layer's data so decoder cosines are defined. Use the mean over layers.
Ybar = zeros(d,N); for l=1:L, Ybar=Ybar+Y{l}; end; Ybar=Ybar/L;
Wg = Ybar*pinv(SSs);                            % d x Kg group mixing
Wg = Wg./(sqrt(sum(Wg.^2,1))+eps);

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

%% ---- Swiss-Prot label setup ----
assert(isfield(S,'labels') && isfield(S,'label_names'), 'Need S.labels and S.label_names.');
LAB=double(S.labels);
if size(LAB,1)~=Nfull && size(LAB,2)==Nfull, LAB=LAB.'; end
labS=LAB(sel,:); Cn=size(labS,2);
lnames=cellstr(string(S.label_names(:)));
if numel(lnames)~=Cn, lnames=arrayfun(@(c)sprintf('label%d',c),1:Cn,'uni',0); end
labBg=mean(labS>0,1)+1e-12;
fprintf('Detected labels: %d residues x %d concepts: %s\n', size(LAB,1), Cn, strjoin(lnames,', '));
dcol = find(~cellfun('isempty',regexpi(lnames,'disulf','once')),1);

[~,ord]=sort(pidxS); pS=pidxS(ord); accS=accResS(ord);

%% =====================================================================
%  (A) GROUP / CONSENSUS ANALYSIS  -- runs on SSs (shared-across-layers)
%  =====================================================================
Zs = SSs; Wdec = Wg; K = Kg;
fprintf('\n################ (A) GROUP-COMMON COMPONENTS (SSs) ################\n');

%% ---- auto-find cysteine atoms in the group-common components ----
enrC=zeros(K,1);
for kk=1:K
    a=abs(Zs(kk,:)); nz=find(a>0); if isempty(nz), continue; end
    [~,o]=sort(a(nz),'descend'); top=nz(o(1:min(TOPN,numel(o))));
    enrC(kk)=mean(isCys(top))/(bgfreq(2)+1e-9);
end
[sv,si]=sort(enrC,'descend');
cysAtoms = si(sv>=ENR_CYS);
fprintf('Group cysteine components (enr>=%gx): %s\n', ENR_CYS, mat2str(cysAtoms(:)'));

if numel(cysAtoms)>=2
    ATOMS=cysAtoms(:).'; nA=numel(ATOMS);

    %% Part 2: pairwise overlap
    TOPset=cell(nA,1); Pset=cell(nA,1);
    for sN=1:nA, a=abs(Zs(ATOMS(sN),:)); [~,o]=sort(a,'descend'); t=o(1:min(TOPN,N)); TOPset{sN}=t; Pset{sN}=unique(pidxS(t)); end
    fprintf('\nPairwise overlap (residue J / protein J):\n');
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

    %% Part 2b: DISULFIDE split on group-common cysteine components
    if ~isempty(dcol)
        isDisulf=labS(:,dcol)>0; nDisCys=sum(isDisulf(isCys)); cysBg=mean(isDisulf(isCys))+1e-9;
        fprintf('\n=== GROUP CYSTEINE COMPONENTS vs DISULFIDE (label col %d) ===\n', dcol);
        fprintf('Background: %.1f%% of cysteines disulfide-annotated (n=%d)\n', 100*cysBg, nDisCys);
        if nDisCys<10, fprintf('** WARNING: only %d disulfide cysteines -- UNDERPOWERED. **\n', nDisCys); end
        fprintf('%8s  %10s  %10s  %12s\n','atom','top-cys','%disulf','enr_vs_cysbg');
        for sN=1:nA
            a=abs(Zs(ATOMS(sN),:)); [~,o]=sort(a,'descend'); top=o(1:min(TOPN,N));
            topCys=top(isCys(top)); if isempty(topCys), continue; end
            pD=mean(isDisulf(topCys));
            fprintf('%8d  %10d  %9.1f%%  %11.2fx\n', ATOMS(sN), numel(topCys), 100*pD, pD/cysBg);
        end
        fprintf('(>1x disulfide-biased; <1x disulfide-depleted; a split here = partition is a CROSS-LAYER consensus property)\n');

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

    %% stacked maps
    figure('Color','w','Position',[80 50 1200 220*nA]);
    for sN=1:nA
        at=ATOMS(sN); a=abs(Zs(at,:)); a=a/max(a+eps); a=a(ord);
        subplot(nA,1,sN); stem(1:N,a,'Marker','none','Color',[.2 .5 .35]); hold on;
        ylim([0 1.15]); xlim([1 N]); ylabel('|act|'); title(sprintf('GROUP-COMMON cysteine atom %d',at));
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
    fprintf('<2 group cysteine components -> cysteine is not a clean cross-layer consensus.\n');
end

%% Part 4: all-atom census on group-common components
domAA=repmat(' ',K,1); domEnr=zeros(K,1); domPur=zeros(K,1);
for kk=1:K
    a=abs(Zs(kk,:)); [~,o]=sort(a,'descend'); top=o(1:min(TOPN,N));
    lets=aaLetter(top); lets=lets(lets~=' ');
    if isempty(lets), domAA(kk)='-'; continue; end
    cnt=zeros(1,20); for j=1:20, cnt(j)=sum(lets==AA(j)); end
    frac=cnt/sum(cnt); enr=frac./(bgfreq+1e-9);
    [domEnr(kk),ix]=max(enr); domAA(kk)=AA(ix); domPur(kk)=frac(ix);
end
isIdentity = domEnr>=ENR_THR & domAA~='-';
fprintf('\n=== GROUP ATOM CENSUS (K=%d, enr>=%g x) ===\n',K,ENR_THR);
fprintf('Identity atoms: %d  |  Diffuse: %d\n', sum(isIdentity), K-sum(isIdentity));
fprintf('\n AA  #atoms  meanEnr  bestEnr  meanPurity\n'); covered=0;
for j=1:20
    idx=find(isIdentity & domAA==AA(j)); if isempty(idx), continue; end; covered=covered+1;
    fprintf('  %s   %3d    %5.1fx   %5.1fx     %.2f\n', AA(j),numel(idx),mean(domEnr(idx)),max(domEnr(idx)),mean(domPur(idx)));
end
fprintf('\nAmino acids with a clean identity atom: %d of 20\n', covered);

%% Part 5: structure census on group-common components
STR_THR=3; MINPOS=5; PUR_THR=0.25;
domCon=zeros(K,1); conEnr=zeros(K,1); conPur=zeros(K,1); conCnt=zeros(K,1);
for kk=1:K
    a=abs(Zs(kk,:)); [~,o]=sort(a,'descend'); top=o(1:min(TOPN,N));
    cnt=sum(labS(top,:)>0,1); fr=cnt/numel(top); en=fr./labBg;
    [conEnr(kk),ic]=max(en); domCon(kk)=ic; conPur(kk)=fr(ic); conCnt(kk)=cnt(ic);
end
isStruct = conEnr>=STR_THR & conCnt>=MINPOS & conPur>=PUR_THR;
fprintf('\n=== GROUP STRUCTURE CENSUS (K=%d, enr>=%gx, >=%d hits, purity>=%.2f) ===\n', K, STR_THR, MINPOS, PUR_THR);
fprintf('Structure atoms: %d\n\n concept            #atoms  meanEnr  bestEnr  meanPurity\n', sum(isStruct));
for c=1:Cn
    idx=find(isStruct & domCon==c); if isempty(idx), continue; end
    fprintf('  %-18s  %3d    %5.1fx   %5.1fx     %.3f\n', ...
        lnames{c}, numel(idx), mean(conEnr(idx)), max(conEnr(idx)), mean(conPur(idx)));
end

%% Part 6+7: fractionation sweep + null on group-common components
fprintf('\n=== GROUP FRACTIONATION SWEEP + NULL (K=%d) ===\n', K);
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
    for u=1:nAj, a=abs(Zs(idx(u),:)); [~,o]=sort(a,'descend'); TS{u}=o(1:min(TOPN,N)); pool=[pool TS{u}]; end
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
fprintf('\n>>> GROUP (consensus): %d of 20 amino acids fractionate (>=2 atoms, below-null, purity>=%.2f)\n', sum(isFrac), PUR_MIN);
fprintf('>>> Amino acids: %s\n', AA(isFrac));

%% Consensus spatial map SSt: report how peaked each cysteine component's
%  across-layer consensus is (singular-value concentration = depth-agreement).
if numel(cysAtoms)>=1 && exist('SSt','var')
    fprintf('\n=== CONSENSUS-ACROSS-LAYERS (SSt) for group cysteine components ===\n');
    fprintf('(higher = component''s spatial pattern is more consistent across the 12 layers)\n');
    for sN=1:min(numel(cysAtoms),nA)
        at=cysAtoms(sN);
        v=SSt(:,at); v=v/ (norm(v)+eps);
        % concentration of the consensus map: fraction of energy in top-TOPN residues
        [~,o]=sort(abs(v),'descend'); conc=sum(v(o(1:min(TOPN,N))).^2);
        fprintf('  component %3d : consensus top-%d energy = %.3f\n', at, TOPN, conc);
    end
end

%% =====================================================================
%  (B) PER-LAYER ANALYSIS -- runs on each Zs{j} (layer-specific projections)
%  NOTE: these are the COMMON components re-projected per layer, NOT
%  independently-learned layer-specific atoms. Expected to be WEAKER than
%  sgBACES, which learns genuine per-layer dictionaries.
%  =====================================================================
fprintf('\n################ (B) PER-LAYER PROJECTIONS (Zs{j}) ################\n');
fprintf('(reminder: gICA per-layer = common comps re-projected; sgBACES learns true layer-specific atoms)\n');
fprintf('\n%6s  %12s  %14s  %18s\n','layer','#cys comps','#frac AAs/20','disulfide top-atom enr');
for l=1:L
    Zl = Zs_layer{l};                      % Kg x N, this layer's projection
    Kl = size(Zl,1);
    % cysteine components in this layer
    eC=zeros(Kl,1);
    for kk=1:Kl
        a=abs(Zl(kk,:)); [~,o]=sort(a,'descend'); top=o(1:min(TOPN,N));
        eC(kk)=mean(isCys(top))/(bgfreq(2)+1e-9);
    end
    cysL=find(eC>=ENR_CYS); nCysL=numel(cysL);
    % disulfide top-atom enrichment in this layer
    dEnr=NaN;
    if ~isempty(dcol) && nCysL>=1
        isDisulf=labS(:,dcol)>0; cysBg=mean(isDisulf(isCys))+1e-9; best=0;
        for cc=cysL(:)'
            a=abs(Zl(cc,:)); [~,o]=sort(a,'descend'); top=o(1:min(TOPN,N));
            topCys=top(isCys(top)); if isempty(topCys), continue; end
            best=max(best, mean(isDisulf(topCys))/cysBg);
        end
        dEnr=best;
    end
    % per-layer fractionation count (identity atoms, below-null, pure) -- quick version
    dA=repmat(' ',Kl,1); dE=zeros(Kl,1); dP=zeros(Kl,1);
    for kk=1:Kl
        a=abs(Zl(kk,:)); [~,o]=sort(a,'descend'); top=o(1:min(TOPN,N));
        lets=aaLetter(top); lets=lets(lets~=' ');
        if isempty(lets), dA(kk)='-'; continue; end
        cnt=zeros(1,20); for j=1:20, cnt(j)=sum(lets==AA(j)); end
        fr=cnt/sum(cnt); en=fr./(bgfreq+1e-9); [dE(kk),ix]=max(en); dA(kk)=AA(ix); dP(kk)=fr(ix);
    end
    isId=dE>=ENR_THR & dA~='-'; nfrac=0;
    for j=1:20
        idx=find(isId & dA==AA(j)); if numel(idx)<2, continue; end
        if mean(dP(idx))<PUR_MIN, continue; end
        TS=cell(numel(idx),1); pool=[];
        for u=1:numel(idx), a=abs(Zl(idx(u),:)); [~,o]=sort(a,'descend'); TS{u}=o(1:min(TOPN,N)); pool=[pool TS{u}]; end
        obs=[]; for p=1:numel(idx), for q=p+1:numel(idx), obs(end+1)=numel(intersect(TS{p},TS{q}))/numel(union(TS{p},TS{q})); end, end %#ok
        upool=unique(pool); nullv=zeros(100,1);
        for b=1:100
            st=cell(numel(idx),1); for u=1:numel(idx), r=randperm(numel(upool),min(TOPN,numel(upool))); st{u}=upool(r); end
            nj=[]; for p=1:numel(idx), for q=p+1:numel(idx), nj(end+1)=numel(intersect(st{p},st{q}))/numel(union(st{p},st{q})); end, end %#ok
            nullv(b)=mean(nj);
        end
        if mean(obs) < mean(nullv)-2*std(nullv), nfrac=nfrac+1; end
    end
    if isnan(dEnr)
        fprintf('%6d  %12d  %14d  %18s\n', layers(l), nCysL, nfrac, 'n/a');
    else
        fprintf('%6d  %12d  %14d  %16.2fx\n', layers(l), nCysL, nfrac, dEnr);
    end
end
fprintf('\n(If per-layer disulfide enr / frac counts are noisier or lower than the GROUP consensus\n');
fprintf(' and than sgBACES layer-specific atoms, that is the expected gICA limitation: its per-layer\n');
fprintf(' view is the shared components re-projected, not genuinely layer-specific structure.)\n');

%% ---- CSV export: focal group cysteine component ----
if numel(cysAtoms)>=1
    ATOM=cysAtoms(1); a=abs(SSs(ATOM,:));
    accClean=cellfun(@(s) clean(s), accResS, 'UniformOutput', false);
    Tcsv=table(string(accClean(:)), rposS(:), a(:), 'VariableNames', {'accession','resnum','activation'});
    writetable(Tcsv, sprintf('group_atom%d_gICA.csv', ATOM));
    fprintf('\nwrote group_atom%d_gICA.csv\n', ATOM);
end
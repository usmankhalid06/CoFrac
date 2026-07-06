clear; close all; clc; rng(0);

% ============================ CONFIG =======================================
FILENAME='esm2_35M_multilayer.mat'; METAFILE='protein_metadata_12.tsv';
LAYER=9; K=200; N_SUB=60000; PREP='center';
TOPK=5; EPOCHS=100; LR=1e-2; BATCH=512;
TOPN=200; NLAB=5;
PUR_MIN=0.25;          % purity floor: below this, "dominant AA" is noise
NULL_NPERM=500;        % permutations for the chance-partition null
% ===========================================================================
AA='ACDEFGHIKLMNPQRSTVWY'; clean=@(s) s(isstrprop(s,'alphanum'));


%% ---- load + preprocess + TRAIN SAE ----
S=load(FILENAME);
if ndims(S.X)==3
    if isfield(S,'layers_used'), k=find(S.layers_used(:)==LAYER,1); else, k=LAYER; end
    Xl=double(S.X(:,:,k));
else, Xl=double(S.X); end
Nfull=size(Xl,1); pidx=double(S.protein_idx(:)); rpos=double(S.residue_pos(:));
accAll=cellstr(string(S.accessions));
if numel(accAll)==Nfull, accRes=accAll;
else
    uids=unique(pidx); amap=containers.Map('KeyType','double','ValueType','char');
    for i=1:numel(uids), if i<=numel(accAll), amap(uids(i))=accAll{i}; end, end
    accRes=cell(Nfull,1);
    for i=1:Nfull, if isKey(amap,pidx(i)), accRes{i}=amap(pidx(i)); else, accRes{i}=''; end, end
end
mu=mean(Xl,1);
switch PREP
    case 'center', Xl=Xl-mu;
    case 'zscore', sg=std(Xl,0,1); sg(sg==0)=1; Xl=(Xl-mu)./sg;
end
if N_SUB<Nfull, sel=randperm(Nfull,N_SUB); else, sel=1:Nfull; end
Y=Xl(sel,:)'; pidxS=pidx(sel); rposS=rpos(sel); accResS=accRes(sel); [d,N]=size(Y); Y=zscore(Y);
tic; [Wdec,Zs]=my_SAE(Y,K,TOPK,EPOCHS,LR,BATCH); toc
fprintf('SAE trained. relErr=%.3f\n', norm(Y-Wdec*Zs,'fro')/norm(Y,'fro'));

%% ---- metadata sequences + background ----
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

%% ---- FIX (1): Swiss-Prot LABEL SETUP moved up here (before Part 1/2b) ----
assert(isfield(S,'labels') && isfield(S,'label_names'), ...
    'Need S.labels and S.label_names from the .mat.');
LAB=double(S.labels);
if size(LAB,1)~=Nfull && size(LAB,2)==Nfull, LAB=LAB.'; end   % rows = residues
labS=LAB(sel,:); Cn=size(labS,2);                             % align to subsample
lnames=cellstr(string(S.label_names(:)));
if numel(lnames)~=Cn, lnames=arrayfun(@(c)sprintf('label%d',c),1:Cn,'uni',0); end
labBg=mean(labS>0,1)+1e-12;
fprintf('Detected labels: %d residues x %d concepts: %s\n', size(LAB,1), Cn, strjoin(lnames,', '));

%% ---- per-atom cysteine enrichment -> AUTO-FIND cysteine atoms ----
enrC=zeros(K,1);
for kk=1:K
    a=abs(Zs(kk,:)); nz=find(a>0); if isempty(nz), continue; end
    [~,o]=sort(a(nz),'descend'); top=nz(o(1:min(TOPN,numel(o))));
    nC=0; tot=0;
    for jj=top(:).'
        key=clean(accResS{jj});
        if isKey(m_seq,key), s=m_seq(key); p=rposS(jj);
            if p>=1 && p<=numel(s) && any(s(p)==AA), tot=tot+1; if s(p)=='C', nC=nC+1; end, end
        end
    end
    if tot>0, enrC(kk)=(nC/tot)/(bgfreq(2)+1e-9); end
end
[sv,si]=sort(enrC,'descend');
cysAtoms = si(sv>=3);
if numel(cysAtoms)<2
    fprintf('Found only %d strong cysteine atom(s); nothing to compare.\n', numel(cysAtoms)); return;
end
ATOMS = cysAtoms(:).'; nA=numel(ATOMS);
fprintf('Cysteine atoms found: '); fprintf('%d ',ATOMS); fprintf('\n');

% amino-acid letter for every residue (voxel)
aaLetter=repmat(' ',1,N);
for jj=1:N
    key=clean(accResS{jj});
    if isKey(m_seq,key), s=m_seq(key); p=rposS(jj);
        if p>=1 && p<=numel(s), aaLetter(jj)=s(p); end
    end
end
isCys=(aaLetter=='C');

%% ============ Part 1: stacked maps, one panel per cysteine atom ============
[~,ord]=sort(pidxS); pS=pidxS(ord); accS=accResS(ord);
figure('Color','w','Position',[80 50 1200 220*nA]);
for sN=1:nA
    at=ATOMS(sN); a=abs(Zs(at,:)); a=a/max(a+eps); a=a(ord);
    subplot(nA,1,sN); stem(1:N,a,'Marker','none','Color',[.25 .45 .8]); hold on;
    ylim([0 1.15]); xlim([1 N]); ylabel('|act|');
    title(sprintf('atom %d (cysteine)   enrC = %.0fx', at, enrC(at)));
    up=unique(pS,'stable'); mass=zeros(numel(up),1);
    for i=1:numel(up), mass(i)=sum(a(pS==up(i))); end
    [~,oi]=sort(mass,'descend');
    for r=1:min(NLAB,numel(up))
        p=up(oi(r)); idx=find(pS==p); xc=mean(idx); [~,im]=max(a(idx)); yc=a(idx(im));
        key=clean(accS{idx(1)}); nm=key;
        if isKey(m_nam,key), nm=m_nam(key); if numel(nm)>22, nm=nm(1:22); end, end
        text(xc,min(1.1,yc+0.06),nm,'Rotation',40,'FontSize',7, ...
             'Color',[.6 0 .1],'Interpreter','none');
    end
end
xlabel('residue (voxel), grouped by protein  \rightarrow  (same order in all panels)');

%% ============ Part 2: full pairwise overlap ============
TOPset=cell(nA,1); Pset=cell(nA,1);
for sN=1:nA
    a=abs(Zs(ATOMS(sN),:)); [~,o]=sort(a,'descend'); t=o(1:min(TOPN,N));
    TOPset{sN}=t; Pset{sN}=unique(pidxS(t));
end
fprintf('\nPairwise overlap on top-%d residues  (residue J / protein J):\n', TOPN);
fprintf('%8s',''); for j=1:nA, fprintf('   atom%-3d',ATOMS(j)); end; fprintf('\n');
for i=1:nA
    fprintf('atom%-4d',ATOMS(i));
    for j=1:nA
        if i==j, fprintf('     --     '); continue; end
        rj=numel(intersect(TOPset{i},TOPset{j}))/numel(union(TOPset{i},TOPset{j}));
        pj=numel(intersect(Pset{i},Pset{j}))/numel(union(Pset{i},Pset{j}));
        fprintf('  %.2f/%.2f',rj,pj);
    end
    fprintf('\n');
end
fprintf('(low everywhere = distinct cysteine sub-networks; high = redundant)\n');

%% ============ Part 2b: cysteine atoms vs DISULFIDE annotation (FIX 2) ======
% Does each cysteine atom preferentially capture disulfide-bonded vs non-disulfide cys?
dcol = find(~cellfun('isempty',regexpi(lnames,'disulf','once')),1);
if isempty(dcol)
    fprintf('\n[Part 2b] no disulfide label column found; skipping functional split.\n');
else
    isDisulf = labS(:,dcol)>0;                      % residue-level disulfide annotation
    nDisCys  = sum(isDisulf(isCys));
    cysBg    = mean(isDisulf(isCys)) + 1e-9;
    fprintf('\n=== CYSTEINE ATOMS vs DISULFIDE STATUS (label col %d) ===\n', dcol);
    fprintf('Background: %.1f%% of cysteines are disulfide-annotated (n=%d)\n', 100*cysBg, nDisCys);
    if nDisCys < 10
        fprintf('** WARNING: only %d disulfide-annotated cysteines -- split test UNDERPOWERED; do not over-read. **\n', nDisCys);
    end
    fprintf('%8s  %10s  %10s  %12s\n','atom','top-cys','%disulf','enr_vs_cysbg');
    for sN=1:nA
        a=abs(Zs(ATOMS(sN),:)); [~,o]=sort(a,'descend'); top=o(1:min(TOPN,N));
        topCys = top(isCys(top));
        if isempty(topCys), continue; end
        pD = mean(isDisulf(topCys));
        fprintf('%8d  %10d  %9.1f%%  %11.2fx\n', ATOMS(sN), numel(topCys), 100*pD, pD/cysBg);
    end
    fprintf('(>1x = disulfide-biased; <1x = disulfide-depleted; spread across atoms = functional partition)\n');
end

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

%% ============ Part 3: layering test (FIX 6: mean, not median) ==============
focal=ATOMS(1);
aF=abs(Zs(focal,:)); [~,oF]=sort(aF,'descend'); topF=oF(1:min(TOPN,N));
structCys=topF(isCys(topF));
cysOther =setdiff(find(isCys),structCys);
nonCys   =find(~isCys);
others=ATOMS(2:end);
fprintf('\nLayering test: do atom-%d structural cysteines (n=%d) load on the others?\n', ...
        focal, numel(structCys));
figure('Color','w','Position',[80 80 max(560*numel(others),560) 360]);
for sN=1:numel(others)
    at=others(sN); a=abs(Zs(at,:));
    mS=mean(a(structCys)); mC=mean(a(cysOther)); mN=mean(a(nonCys));   % FIX: mean
    fprintf('  atom %d:  mean|act| struct-cys=%.4f | other-cys=%.4f | non-cys=%.4f\n', ...
            at,mS,mC,mN);
    if mS>2*mN
        fprintf('     -> ABOVE non-cys baseline: identity+structure LAYERED on the same residue.\n');
    else
        fprintf('     -> not above non-cys baseline: separate pools, not layered.\n');
    end
    subplot(1,numel(others),sN); hold on;
    grp={{structCys,[.85 .1 .1]},{cysOther,[.1 .4 .85]},{nonCys,[.5 .5 .5]}};
    for g=1:3
        v=sort(a(grp{g}{1})); y=(1:numel(v))/max(numel(v),1);
        stairs(v,y,'Color',grp{g}{2},'LineWidth',1.3);
    end
    xlabel(sprintf('|act| on atom %d',at)); ylabel('cum. frac'); ylim([0 1]);
    title(sprintf('atom %d',at));
    legend({sprintf('atom-%d struct-cys',focal),'other cys','non-cys'}, ...
           'Location','southeast','FontSize',7);
end

% decoder cosines among cysteine atoms (REPORTED ONLY -- not used as evidence;
% sits at the ~sqrt(2/(pi*d)) random-vector floor and is consistent with absorption)
Dn = Wdec(:,ATOMS); Dn = Dn./sqrt(sum(Dn.^2,1)+eps); Ccos = abs(Dn'*Dn);
fprintf('\nSAE decoder cosines |cos(D_i,D_j)| among cysteine atoms (reference only):\n');
fprintf('%8s',''); for j=1:nA, fprintf('   atom%-3d',ATOMS(j)); end; fprintf('\n');
for i=1:nA
    fprintf('atom%-4d',ATOMS(i));
    for j=1:nA
        if i==j, fprintf('     --   '); else, fprintf('   %.2f  ',Ccos(i,j)); end
    end
    fprintf('\n');
end
fprintf('(NOTE: in d=%d, random unit vectors have |cos|~%.3f -- treat as reference, NOT evidence)\n', ...
        d, sqrt(2/(pi*d)));

%% ============ Part 4: ALL-ATOM CENSUS ======================================
ENR_THR=3;
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

ord_atoms=[];
for j=1:20
    idx=find(isIdentity & domAA==AA(j)); [~,sidx]=sort(domEnr(idx),'descend');
    ord_atoms=[ord_atoms; idx(sidx)];
end
df=find(~isIdentity); [~,sd]=sort(domEnr(df),'descend'); ord_atoms=[ord_atoms; df(sd)];

M=abs(Zs(ord_atoms,ord)); M=M./(max(M,[],2)+eps);
figure('Color','w','Position',[60 50 1300 820]);
imagesc(M.^0.5); colormap(flipud(gray)); cb=colorbar; cb.Label.String='|act| (row-norm, sqrt)';
ylab=cell(K,1);
for r=1:K, at=ord_atoms(r);
    if isIdentity(at), ylab{r}=sprintf('%s  a%d  %0.0fx',domAA(at),at,domEnr(at));
    else,             ylab{r}=sprintf('-   a%d',at); end
end
set(gca,'YTick',1:K,'YTickLabel',ylab,'FontSize',6,'TickLength',[0 0]);
xlabel('residue (voxel), grouped by protein  \rightarrow'); ylabel('atom (dominant amino acid)');
nIdent=sum(isIdentity);
title(sprintf('All %d atoms labeled by dominant amino acid  [%d identity / %d diffuse]',K,nIdent,K-nIdent));
yline(nIdent+0.5,'r-','LineWidth',1.2);

fprintf('\n=== ATOM CENSUS (K=%d, enr>=%g x) ===\n',K,ENR_THR);
fprintf('Identity atoms: %d  |  Diffuse: %d\n', nIdent, K-nIdent);
fprintf('\n AA  #atoms  meanEnr  bestEnr  meanPurity\n'); covered=0;
for j=1:20
    idx=find(isIdentity & domAA==AA(j)); if isempty(idx), continue; end; covered=covered+1;
    fprintf('  %s   %3d    %5.1fx   %5.1fx     %.2f\n', AA(j),numel(idx),mean(domEnr(idx)),max(domEnr(idx)),mean(domPur(idx)));
end
fprintf('\nAmino acids with a clean identity atom: %d of 20\n', covered);

%% ============ Part 5: STRUCTURE CENSUS vs Swiss-Prot (FIX 5: label SAE) ====
STR_THR=3; MINPOS=5; PUR_THR=0.25;
domCon=zeros(K,1); conEnr=zeros(K,1); conPur=zeros(K,1); conCnt=zeros(K,1);
for kk=1:K
    a=abs(Zs(kk,:)); [~,o]=sort(a,'descend'); top=o(1:min(TOPN,N));
    cnt=sum(labS(top,:)>0,1); fr=cnt/numel(top); en=fr./labBg;
    [conEnr(kk),ic]=max(en); domCon(kk)=ic; conPur(kk)=fr(ic); conCnt(kk)=cnt(ic);
end
isStruct = conEnr>=STR_THR & conCnt>=MINPOS & conPur>=PUR_THR;
fprintf('\n=== STRUCTURE CENSUS (SAE, K=%d, enr>=%gx, >=%d hits, purity>=%.2f) ===\n', K, STR_THR, MINPOS, PUR_THR);
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
        title(sprintf('SAE atom %d  ->  %s   (enr %.0fx, purity %.2f)', at, lnames{domCon(at)}, conEnr(at), conPur(at)));
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

%% ============ Part 6 + 7: FRACTIONATION SWEEP + NULL (FUSED, FIX 3+4) ======
% Verdict is driven by the chance-partition NULL (the authority), with a
% purity floor. Cosine is reported for reference only and never gates anything.
fprintf('\n=== FRACTIONATION SWEEP + CHANCE-PARTITION NULL (K=%d, layer=%d) ===\n', K, LAYER);
fprintf('%3s %7s %9s %10s %7s %8s %9s   %s\n', ...
        'AA','#atoms','obsResJ','nullResJ','z','meanPur','meanCos','verdict');

fracTable = nan(20,7);  % [#atoms obsJ nullMean z meanPur meanCos passNull]
for j = 1:20
    idx = find(isIdentity & domAA==AA(j));
    nAj = numel(idx);
    if nAj < 1, continue; end
    meanPur = mean(domPur(idx));

    if nAj == 1
        fprintf('%3s %7d %9s %10s %7s %8.2f %9s   single atom (no partition)\n', ...
                AA(j), nAj, '--','--','--', meanPur, '--');
        fracTable(j,:) = [nAj NaN NaN NaN meanPur NaN 0];
        continue;
    end

    % top-TOPN residue sets per atom + decoder cosines
    TS=cell(nAj,1); pool=[];
    for u=1:nAj, a=abs(Zs(idx(u),:)); [~,o]=sort(a,'descend'); TS{u}=o(1:min(TOPN,N)); pool=[pool TS{u}]; end
    Dn=Wdec(:,idx); Dn=Dn./sqrt(sum(Dn.^2,1)+eps); C=abs(Dn'*Dn); meanCos=mean(C(triu(true(nAj),1)));

    % observed mean pairwise residue-Jaccard
    obs=[]; for p=1:nAj, for q=p+1:nAj
        obs(end+1)=numel(intersect(TS{p},TS{q}))/numel(union(TS{p},TS{q})); %#ok
    end, end
    obsJ=mean(obs);

    % chance-partition null: repartition this AA's own residue pool at random
    upool=unique(pool); Tn=min(TOPN,N); nullJ=zeros(NULL_NPERM,1);
    for b=1:NULL_NPERM
        sets=cell(nAj,1);
        for u=1:nAj, r=randperm(numel(upool),min(Tn,numel(upool))); sets{u}=upool(r); end
        nj=[]; for p=1:nAj, for q=p+1:nAj
            nj(end+1)=numel(intersect(sets{p},sets{q}))/numel(union(sets{p},sets{q})); %#ok
        end, end
        nullJ(b)=mean(nj);
    end
    z=(obsJ-mean(nullJ))/(std(nullJ)+eps);
    passNull = obsJ < mean(nullJ)-2*std(nullJ);   % observed disjointness BELOW chance

    if passNull && meanPur>=PUR_MIN
        verd='FRACTIONATES (below-chance disjoint, pure)';
    elseif passNull
        verd=sprintf('below-chance disjoint but LOW PURITY (%.2f) -> noise', meanPur);
    else
        verd='not below chance -> redundant / abundance artifact';
    end
    fprintf('%3s %7d %9.3f %10.3f %7.1f %8.2f %9.3f   %s\n', ...
            AA(j), nAj, obsJ, mean(nullJ), z, meanPur, meanCos, verd);
    fracTable(j,:) = [nAj obsJ mean(nullJ) z meanPur meanCos passNull];
end

% Headline count: AAs with >=2 atoms, below-chance disjoint, AND pure
isFrac = fracTable(:,1)>=2 & fracTable(:,7)==1 & fracTable(:,5)>=PUR_MIN;
nFrac  = sum(isFrac);
fprintf('\n>>> %d of 20 amino acids fractionate: >=2 atoms, residue-disjointness BELOW chance (z<-2), purity>=%.2f\n', ...
        nFrac, PUR_MIN);
fprintf('>>> Headline count is null-authoritative; cosine is reference-only and excluded from the criterion.\n');
fracAAs = AA(isFrac);
fprintf('>>> Amino acids: %s\n', fracAAs);

%% ============ CSV export for the 3D structure figure (Python reads atomN_all.csv) ============
% Export the atoms you want to render in 3D. For the disulfide structure panel,
% use the disulfide-biased SAE atom (e.g. 32). Add TM / DNA atoms if you want those too.
EXPORT_ATOMS = [32 1 196];   % disulfide, transmembrane, DNA-binding (edit to your atom indices)

accClean = cellfun(@(s) clean(s), accResS, 'UniformOutput', false);
for ee = 1:numel(EXPORT_ATOMS)
    AT = EXPORT_ATOMS(ee);
    a  = abs(Zs(AT,:));
    Texp = table(string(accClean(:)), rposS(:), a(:), ...
                 'VariableNames', {'accession','resnum','activation'});
    fn = sprintf('atom%d_SAE.csv', AT);
    writetable(Texp, fn);
    fprintf('wrote %s\n', fn);
end
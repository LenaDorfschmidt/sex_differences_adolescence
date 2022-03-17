clear all

rng('default')

sample = 'regr_FD_by_sex';

if strcmp(sample,'main')
    outpath = '../../Results/gene_decoding/';
    inpath = '../../Results/maturational_index/';
elseif strcmp(sample,'CV_stratified')
    outpath = '../../Results/replication/CV_stratified/gene_decoding/';
    inpath = '../../Results/replication/CV_stratified/maturational_index/';
elseif strcmp(sample,'glob.fc.regr')
    outpath = '../../Results/replication/glob.fc.regr/gene_decoding/';
    inpath = '../../Results/replication/glob.fc.regr/maturational_index/';
elseif strcmp(sample,'GSR')
    outpath = '../../Results/replication/nspn.gsr/gene_decoding/';
    inpath = '../../Results/replication/nspn.gsr/maturational_index/';
elseif strcmp(sample,'motionmatched')
    outpath = '../../Results/replication/motionmatched/gene_decoding/';
    inpath = '../../Results/replication/motionmatched/maturational_index/';
elseif strcmp(sample,'regr_FD_by_sex')
    outpath = '../../Results/replication/regr_FD_by_sex/gene_decoding/';
    inpath = '../../Results/replication/regr_FD_by_sex/maturational_index/';
end

if ~isdir(outpath)
    mkdir outpath;
end

MI = readtable([inpath 'diff.mi.all.txt'],'TreatAsEmpty',{'NA'},'ReadVariableNames',0);

%mi = mean([table2array(MI(17:196,1)),table2array(MI(197:376,1))],2,'omitnan');
mi = table2array(MI(17:196,1));
 
load('~/Documents/Education/Cambridge/Datasets/AllenBrainAtlas/Petra/AHBAdata.mat')
parcelExpression = parcelExpression(:,2:end);

idx_mi = not(ismissing(mi));
parcelExpression = parcelExpression(idx_mi,:);
mi = mi(idx_mi);

idx_gene = find(not(isnan(sum(parcelExpression,2))));

keep_idx = find(idx_mi);
keep_idx = keep_idx(idx_gene);

%pls1_out = horzcat(PLS.XS(:,1), keep_idx);


X = parcelExpression(idx_gene,:);
X = zscore(X')';
Y = mi(idx_gene);
Y = zscore(Y);

%save_AHBA(tmp,:) = X;
%csvwrite('~/Documents/Education/Cambridge/fMRI_NSPN/Results/gene_decoding/AHBA_PLS_X_data.csv',save_AHBA)
%writetable(cell2table(nm),'~/Documents/Education/Cambridge/fMRI_NSPN/Results/gene_decoding/AHBA_gene_names.csv')

dim=2;
numPerm=1000;

[PLS.XL,PLS.YL,PLS.XS,PLS.YS,PLS.BETA,PLS.PCTVAR,PLS.MSE,PLS.stats]=plsregress(X,Y,dim,'CV',2);

temp=cumsum(100*PLS.PCTVAR(2,1:dim)); %percentage of variance explained in Y by each PLS
Rsquared = temp(dim); %Equal to sum, variance explained by both component

%assess significance of PLS result
for j=1:numPerm %10000
    j
    order=randperm(size(Y,1));
    Yp=Y(order,:);
    [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(X,Yp,dim,'CV',2);
    PCTVARr_all(:,j)=PCTVARr(2,:);
    temp=cumsum(100*PCTVARr(2,1:dim));
    Rsq(j) = temp(dim);
end

for i=1:dim
   v_real=sum(PLS.PCTVAR(2,1:i));
   v_rand=sum(PCTVARr_all(1:i,:),1);
   p_values(i)=sum(v_real<v_rand)/numel(v_rand);
   %fig_raf_hist(v_rand,min(v_rand):0.1:max(v_rand),['PLS ' num2str(i)], 'Probaility', v_real);
   figure
   hold on
   hist(v_rand,30);
    plot(v_real,20,'.r','MarkerSize',15);
    set(gca,'Fontsize',34)
   xlabel(['PLS ' num2str(i)]);
   ylabel('Probability');
end

figure
bar(cumsum(PLS.PCTVAR'),'EdgeColor','none','FaceAlpha',0.7)
xlabel('PLS Component', 'fontsize',14)
ylabel('Variance Explained', 'fontsize',14)
title('Variance Explained in X and Y','fontsize',18)
legend('Variance Explained in X','Variance Explained in Y')

[rho,pval] = corr(Y, PLS.XS);

nm=probeInformation.GeneSymbol;

%save stats
p_values
%Variance explained by PLS component, 
%first row, explained in X (genes matrix), second explained in Y (regional
%vals) this is the one you want
myStats=[PLS.PCTVAR]; 
path_output = [outpath, 'plsout'];
[s1 s2 s3]=mkdir(path_output);
csvwrite([path_output '/PLS_stats.csv'],myStats);
save([path_output '/Y.mat'],'Y');
save([path_output '/PCTVARr_all.mat'],'PCTVARr_all');
save([path_output '/PLS.mat'],'PLS');
save([path_output '/X.mat'],'X');
save([path_output '/p_values.mat'],'p_values');
save([path_output '/dim.mat'],'dim');
save([path_output '/nm.mat'],'nm');


tab = table(probeInformation.GeneSymbol, probeInformation.EntrezID);

maxl = 11;
for g = 1:size(probeInformation.EntrezID,1)
    nchar = length(num2str(probeInformation.EntrezID(g)));
    probeInformation.ENSG{g} = ['ENSG',repmat('0',1,(maxl-nchar)),num2str(probeInformation.EntrezID(g))]; 
end  
    
writetable(table(probeInformation.GeneSymbol, probeInformation.ENSG'),'~/Desktop/tmp.csv')  
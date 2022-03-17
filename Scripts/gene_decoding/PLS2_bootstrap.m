clear
path_data='../../Results/gene_decoding/plsout/';
%load('AHBAdata.mat')
%nm = nm(sum(isnan(parcelExpression),2)==0);

load([path_data 'X.mat']);
load([path_data 'PLS.mat']);
load([path_data 'Y.mat']);
load([path_data 'dim.mat']);
load([path_data 'nm.mat']);
load('~/Documents/Education/Cambridge/Datasets/AllenBrainAtlas/Petra/AHBAdata.mat')

%number of bootstrap iterations
bootnum=10000;
index=1:size(PLS.stats.W,1);

[R1,p1]=corr(PLS.XS,X);

% align PLS components with desired direction%
PLSw = []; roiindex = []; EntrezIDs = []; PLSids = {};
for ixs = 1:size(R1,1)
    Rlist = R1(ixs,:);
    if abs(max(Rlist))<abs(min(Rlist))
        PLS.stats.W(:,ixs) = -1*PLS.stats.W(:,ixs);
        PLS.XS(:,ixs) = -1*PLS.XS(:,ixs);
    end
    
    [PLSw(ixs,:),wSort(ixs,:)] = sort(PLS.stats.W(:,ixs),'descend');
    PLSids(ixs,:) = nm(wSort(ixs,:));
    geneindex(ixs,:) = index(wSort(ixs,:));
    EntrezIDs(ixs,:) = probeInformation.EntrezID(wSort(ixs,:));
end

%define variables for storing the (ordered) weights from all bootstrap runs
PLSweights1=[];

%start bootstrap
for i=1:bootnum
    i
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim,'CV',2); %perform PLS for resampled data
      
    PLSweights_ind=[];
    for iw=1:size(PLSw,1)
        temp=stats.W(:,iw);%extract PLS1 weights
        %CHANGEEE
        newW=temp(wSort(iw,:)); %order the newly obtained weights the same way as initial PLS 
   
       if corr(PLSw(iw,:)',newW)<0
          newW=-1*newW; 
       end
       %PLSweights_ind=[PLSweights_ind,newW]; 
       PLSweights(i,iw,:)=newW;
    end
end


%pval_rh=0.05;
for id=1:dim

    PLS_perm=squeeze(PLSweights(:,id,:));
    PLSsw=std(PLS_perm);
  
    %get bootstrap weights
    temp1_gen=PLSw(id,:)./PLSsw;
    PLS_real=repmat(PLSw(id,:),size(PLS_perm,1),1);
   
    temp1=temp1_gen;
    %order bootstrap weights (Z) and names of regions
    [Z1 ind1]=sort(temp1,'descend');
    PLSid=PLSids(id,ind1);
    geneindex_ind=geneindex(id,ind1);
    entrezID = EntrezIDs(id,ind1);

    tab_all = table(PLSid',categorical(entrezID)',geneindex_ind',Z1');
    tab_all.Properties.VariableNames = {'gene_name','entrez_ID','gene_index','z_uncorr'};
    writetable(tab_all,[path_data 'pls' num2str(id) '.csv'])
    
    %%%%%%%
    
%     weight=abs(Z1)/max(Z1);
%     fid1 = fopen([path_data 'PLS_pvals_geneWeights_' num2str(id) '.csv'],'w');
%     for i=1:length(geneindex)
%         if pval_corrected(i)<pval_rh
%             fprintf(fid1,'%s, %f\n', PLSid{i}, Z1(i));
%         end
%     end
%     fclose(fid1);
%     
%     
%     temp1=-temp1_gen;
%     %order bootstrap negative weights
%     [Z1 ind1]=sort(temp1,'descend');
%     PLSid=PLSids(id,ind1);
%     geneindex_ind=geneindex(id,ind1);
%     pvals=1-normcdf(Z1);
%     pval_corrected=mafdr(pvals,'BHFDR',true);
%     weight=Z1/max(Z1);
%     
%     %print out results
%     fid1 = fopen([path_data 'PLS_geneWeights_' num2str(id) '_neg.csv'],'w');
%     display(['Genes weighted ' path_data 'PLS_geneWeights_' num2str(id) '.csv'])
%     for i=1:length(geneindex)
%         fprintf(fid1,'%s, %d, %f, %f\n', PLSid{i}, geneindex_ind(i), Z1(i),pval_corrected(i) );
%     end
%     fclose(fid1);
%     
%         %print out results
%     fid1 = fopen([path_data 'PLS_pvals_geneWeights_' num2str(id) '_neg.csv'],'w');
%     display(['Genes weighted ' path_data 'PLS_geneWeights_' num2str(id) '.csv'])
%     for i=1:length(geneindex)
%         
%         if pval_corrected(i)<pval_rh
%             fprintf(fid1,'%s, %f\n', PLSid{i}, Z1(i));
%         end
%     end
%     fclose(fid1);
%     
%     
% 
%     temp1=abs(temp1_gen);
%     %order bootstrap negative weights
%     [Z1 ind1]=sort(temp1,'descend');
%     PLSid=PLSids(id,ind1);
%     geneindex_ind=geneindex(id,ind1);
%     pvals=(1-normcdf(Z1))*2;
%     pval_corrected=mafdr(pvals,'BHFDR',true);
%     
%      pvals_ordered=pval_corrected(ind1);
%      
%    
%     weight=Z1/max(Z1);
%     
%     %print out results
%     fid1 = fopen([path_data 'PLS_geneWeights_' num2str(id) '_abs.csv'],'w');
%     display(['Genes weighted ' path_data 'PLS_geneWeights_' num2str(id) '.csv'])
%     for i=1:length(geneindex)
%         fprintf(fid1,'%s, %d, %f, %f\n', PLSid{i}, geneindex_ind(i), Z1(i),PLSid(i));
%     end
%     fclose(fid1); 
%     
%         fid1 = fopen([path_data 'PLS_pvals_geneWeights_' num2str(id) '_abs.csv'],'w');
%     display(['Genes weighted ' path_data 'PLS_geneWeights_' num2str(id) '.csv'])
%     for i=1:length(geneindex)
%         
%         if pval_corrected(i)<pval_rh
%             fprintf(fid1,'%s, %f\n', PLSid{i}, Z1(i));
%         end
%     end
%     fclose(fid1); 
%     clear PLS_perm
end

%clear all
MI = readtable('~/Documents/Education/Cambridge/fMRI_NSPN/Results/gene_decoding/diff.mi.all.txt','TreatAsEmpty',{'NA'});
mi = mean([MI.x(17:196),MI.x(197:376)],2,'omitnan');

load('~/Documents/Education/Cambridge/Datasets/AllenBrainAtlas/Petra/AHBAdata.mat')
parcelExpression = parcelExpression(:,2:end);

idx_mi = not(ismissing(mi));
parcelExpression = parcelExpression(idx_mi,:);
mi = mi(idx_mi);

idx_gene = find(not(isnan(sum(parcelExpression,2))));

coordinates = readtable('/Users/Lena/Documents/Education/Cambridge/Datasets/HCP.Glasser.coordinates.csv','TreatAsEmpty','NA')
coordinates = coordinates(196:376,:);

coordinates = coordinates(idx_mi,:);
coordinates = coordinates(idx_gene,:);

for c=1:dim
    for a=1:3
        component = PLS.XS(:,c);
        coords = coordinates{:,(a+1)};
        [rho(c,a),pval(c,a)] = corr(component, coords)
    end
end

tmp = (1:180);
tmp = tmp(idx_mi);
tmp = tmp(idx_gene);

print_xs = repmat({'NA'}, 180,1);
print_xs(tmp) = num2cell(PLS.XS(:,1));
print_xs = [print_xs; print_xs];

print_xs_inv = repmat({'NA'}, 180,1);
print_xs_inv(tmp) = num2cell(PLS.XS(:,1)*(-1));
print_xs_inv = [print_xs_inv; print_xs_inv];


writetable(cell2table(print_xs),'~/Desktop/b4p.pls.xs.txt','WriteVariableNames',0,'WriteRowNames',0)
writetable(cell2table(print_xs_inv),'~/Desktop/b4p.pls.xs.inverted.txt','WriteVariableNames',0,'WriteRowNames',0)

print_ys = repmat({'NA'}, 180,1);
print_ys(tmp) = num2cell(PLS.YS(:,1));
print_ys = [print_ys; print_ys];

print_ys_inv = repmat({'NA'}, 180,1);
print_ys_inv(tmp) = num2cell(PLS.YS(:,1)*(-1));
print_ys_inv = [print_ys_inv; print_ys_inv];

writetable(cell2table(print_ys),'~/Desktop/b4p.pls.ys.txt','WriteVariableNames',0,'WriteRowNames',0)
writetable(cell2table(print_ys_inv),'~/Desktop/b4p.pls.ys.inverted.txt','WriteVariableNames',0,'WriteRowNames',0)



CORT = find(strcmp(probeInformation.GeneSymbol,'CORT'));
NPY = find(strcmp(probeInformation.GeneSymbol,'NPY'));
SST = find(strcmp(probeInformation.GeneSymbol,'SST'));

CORT_PLS = find(strcmp(probeInformation.GeneSymbol(geneindex(1,:)),'CORT'))
NPY_PLS = find(strcmp(probeInformation.GeneSymbol(geneindex(1,:)),'NPY'))
SST_PLS = find(strcmp(probeInformation.GeneSymbol(geneindex(1,:)),'SST'))

corr(mi,parcelExpression(:,CORT),'rows','complete')
corr(mi,parcelExpression(:,NPY),'rows','complete')
corr(mi,parcelExpression(:,SST),'rows','complete')


print_cort = repmat({'NA'}, 180,1);
print_cort(tmp) = num2cell(parcelExpression(idx_gene,CORT));
print_cort = [print_cort; print_cort];
writetable(cell2table(print_cort),'~/Desktop/AHBA_expression_CORT.txt','WriteVariableNames',0,'WriteRowNames',0)

print_npy = repmat({'NA'}, 180,1);
print_npy(tmp) = num2cell(parcelExpression(idx_gene,NPY));
print_npy = [print_npy; print_npy];
writetable(cell2table(print_npy),'~/Desktop/AHBA_expression_NPY.txt','WriteVariableNames',0,'WriteRowNames',0)

print_sst = repmat({'NA'}, 180,1);
print_sst(tmp) = num2cell(parcelExpression(idx_gene,SST));
print_sst = [print_sst; print_sst];
writetable(cell2table(print_sst),'~/Desktop/AHBA_expression_SST.txt','WriteVariableNames',0,'WriteRowNames',0)


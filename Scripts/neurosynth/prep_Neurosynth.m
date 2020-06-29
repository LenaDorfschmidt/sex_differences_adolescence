% %% Adjust HCP parcellation labels
%%% OUTDATED/NOT USED!
% V = MRIread('~/Desktop/prep_neurosynth/HCPMMP1_on_MNI152_ICBM2009a_nlin.nii'); 
% Y.vol(:,round(size(Y.vol,2)/2):end,:) = Y.vol(:,round(size(Y.vol,2)/2):end,:)+180;
% [r,c,v] = ind2sub(size(Y.vol),find(Y.vol == 180));
% mask = c<round(size(Y.vol,2)/2);
% Y.vol(Y.vol==180) = 0;
% Y.vol(r(mask), c(mask),v(mask)) = 180;
% 
% MRIwrite(Y, '~/Desktop/prep_neurosynth/HCP_MNI_152.nii.gz','float')

% %% ADD Freesurfer Subcortical ROIs to HCP
% freesurfer = MRIread('~/Desktop/prep_neurosynth/mni_aparc+aseg.nii');
% HCP = MRIread('~/Desktop/prep_neurosynth/HCP_mapped_to_correct_MNI.nii.gz');
% 
% V = HCP;
% 
% subc_ids = [10,11,12,13,17,18,26,28,49,50,51,52,53,54,58,60];
% new_ids = (361:376);
% for subc = 1:16
%     cur_label = subc_ids(subc);
%     V.vol(find(freesurfer.vol==cur_label)) = new_ids(subc);
% end
% 
% MRIwrite(V, '~/Desktop/prep_neurosynth/HCP_freesurfer_MNI_152.nii.gz','float')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%       MI Difference        %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HCP = MRIread('~/Desktop/prep_neurosynth/uchange_parcellation_to_MNI_quick.nii');
V = HCP;

MI_raw = table2array(readtable("~/Documents/Education/Cambridge/nspn_fmri/diff_cc_mat_idx.txt", 'ReadVariableNames', 0));
f_MI = table2array(readtable("~/Documents/Education/Cambridge/nspn_fmri/female/cc_mat_idx.txt",'Format','%*s%f'));
m_MI = table2array(readtable("~/Documents/Education/Cambridge/nspn_fmri/male/cc_mat_idx.txt",'Format','%*s%f'));
MI_sorted = arrange_frantisek_kirstie(MI_raw);
MI_sorted(isnan(MI_sorted)) = 0;
MI = vertcat(MI_sorted, MI_raw(1:16));

for label = 1:sum(unique(V.vol)~=0)
    V.vol(find(V.vol==label)) = MI(label);
end

MRIwrite(V, '~/Desktop/prep_neurosynth/neurosynth_MI_zeros.nii.gz','float')

nanV = V;
nanV.vol(find(nanV.vol==0)) = NaN;
MRIwrite(nanV, '~/Desktop/prep_neurosynth/neurosynth_MI_nan.nii.gz','float')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%      Disruptive ROIs       %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disr = HCP;

disruptive_binary = table2array(readtable('/Users/Lena/Documents/Education/Cambridge/nspn_fmri/perm_shuffle/MI_diff_disruptive_incl_subc.txt', 'ReadVariableNames', 0));
disruptive_binary(strcmp(disruptive_binary,'NA')) = {[NaN]}
sorted_disruptive = arrange_frantisek_kirstie(disruptive_binary);
full_disruptive = [sorted_disruptive; disruptive_binary(1:16,1)];
disruptive_delta_MI = MI;
disruptive_delta_MI(cellfun(@isnan, full_disruptive)) = 0;

for label = 1:sum(unique(disr.vol)~=0)
    disr.vol(find(disr.vol==label)) = disruptive_delta_MI(label);
end

MRIwrite(disr, '~/Desktop/prep_neurosynth/neurosynth_MI_disruptive.nii.gz','float')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%      Conservative ROIs       %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cons = HCP;

conservative_binary = table2array(readtable('/Users/Lena/Documents/Education/Cambridge/nspn_fmri/perm_shuffle/MI_diff_conservative_incl_subc.txt', 'ReadVariableNames', 0));
conservative_binary(strcmp(conservative_binary,'NA')) = {[NaN]}
sorted_conservative = arrange_frantisek_kirstie(conservative_binary);
full_conservative = [sorted_conservative; conservative_binary(1:16,1)];
conservative_delta_MI = MI;
conservative_delta_MI(cellfun(@isnan, full_conservative)) = 0;

for label = 1:sum(unique(cons.vol)~=0)
    cons.vol(find(cons.vol==label)) = conservative_delta_MI(label);
end

MRIwrite(disr, '~/Desktop/prep_neurosynth/neurosynth_MI_conservative.nii.gz','float')


% labels = unique(V.vol);
% cortical_ids = find(unique(V.vol)>1000);
% for roi = 1:360
%     cur_label = labels(cortical_ids(roi));
%     V.vol(find(V.vol==cur_label)) = MI_sorted(roi);
% end
%     
% subc_ids = [10,11,12,13,17,18,26,28,49,50,51,52,53,54,58,60];
% for subc = 1:16
%     cur_label = subc_ids(subc);
%     V.vol(find(V.vol==cur_label)) = MI(subc);
% end
% 
% zerosV = V; nanV = V;
% 
% uvol = unique(V.vol);
% uvol = uvol(uvol>=2);
% for roi = 1:length(uvol)
%     cur_label = uvol(roi);
%     zerosV.V(find(zerosV.V==cur_label)) = 0;
% end
% 
% MRIwrite(zerosV, '~/Desktop/MI_zeros.nii.gz','float')
% tmp = MRIread('~/Desktop/MI_zeros.nii.gz')
% 
% incl = ismember(labels,subc_ids);
% for del = 1:sum(unique(nanV.vol)<1000)
%     if incl(del) == 0
%         nanV.vol(find(nanV.vol == labels(del))) = NaN; 
%     end
% end
% 
% MRIwrite(nanV, '~/Desktop/MI_nan.nii.gz','float')
% tmp = MRIread('~/Desktop/MI_nan.nii.gz');
% 

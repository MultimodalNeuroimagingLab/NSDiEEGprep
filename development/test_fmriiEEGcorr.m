
localDataPath = setLocalDataPath(1); 






%% slow correlations

fmri_beta = nsd_fmri.beta_avgR;

el_nr = find(ismember(all_channels.name,'ROC2'));
t_id = tt>-.05 & tt<.85;
tt_temp = tt(t_id);
thisBB = squeeze(bb_1000(el_nr,t_id,:));

corr_out = zeros(length(tt_temp),size(fmri_beta,1));
for kk = 1:size(fmri_beta,1)
    if mod(kk,500)==0, disp([int2str(kk) ' of ' int2str(size(fmri_beta,1))]), end
    R = faster_corr_mtrx([thisBB' fmri_beta(kk,:)']);
    corr_out(:,kk) = R(1:end-1,end);
end

%%

[u,s,v] = svd(corr_out','econ');

figure,plot(cumsum(diag(s).^2)./sum(diag(s).^2))

%% render SVD 

% render and plot SVD

% Load fs average surfaces
gI_L = gifti(fullfile(localDataPath.fs,'fsaverage','inflated.L.surf.gii'));
gI_R = gifti(fullfile(localDataPath.fs,'fsaverage','inflated.R.surf.gii'));
sulcal_labels_L = read_curv(fullfile(localDataPath.fs,'fsaverage','surf','lh.sulc'));
sulcal_labels_R = read_curv(fullfile(localDataPath.fs,'fsaverage','surf','rh.sulc'));

for hh = 2

    if hh==1
        hemi = 'l';
        g = gI_L;
        sulcal_labels = sulcal_labels_L;
        views_plot = {[40,-30],[-45,-10],[-90,20]};
    elseif hh==2
        hemi = 'r';
        g = gI_R;
        sulcal_labels = sulcal_labels_R;
        views_plot = {[-40,-30],[45,-10],[90,20]};
    end

    % load the color map
    cmap = cmapsign4(200);

    % make a plot with electrode dots
    for vv = 1:length(views_plot)
        v_d = [views_plot{vv}(1),views_plot{vv}(2)];
        
        for pc_ind = 1:5
            figure
            vert_amp = u(:,pc_ind);
            vert_amp = vert_amp/max(abs(vert_amp)); % set max to 1
            vert_amp(abs(vert_amp)<.2) = NaN;
            vert_amp = 100+round(vert_amp*100); % multiply by 100 and add 100
    
            tH = ieeg_RenderGiftiLabels(g,vert_amp,cmap,[],sulcal_labels);
            ieeg_viewLight(v_d(1),v_d(2)) % change viewing angle   
            set(gcf,'PaperPositionMode','auto')
            print('-dpng','-r300',fullfile(localDataPath.input,'derivatives','zeeshan',['sub-' sub_label],...
                ['svd_sub-' sub_label '_' all_channels.name{el_nr} '_pc' int2str(pc_ind) '_' int2str(v_d(1)) '_' int2str(v_d(2)) '_' hemi]))
        end
    end
%     close all
end

%%

figure('Position',[0 0 200 600])
for kk = 1:5
    subplot(5,1,kk)
    plot(tt_temp,v(:,kk))
    yline(0)
end



function [vert_label, area_label, cmap] = loadAtlas(freesurferPath, hemi, atlasName)

% function to load Brain Atlas for fsaverage (currently)
%
% input:
%     freesurferPath: local path to FreeSurfer directory
%     hemi (l/r): left/right hemisphere 
%     atlasName (Wang/BensonV/BensonE/visf/HCP)
%     
% output:
%     vert_label: vertex labels
%     area_label: surface labels corresponding to freesurfer vertices
%     cmap: colormap for different areas
%     
%
% Example usage:
%   [vert_label, area_label, cmap] = loadAtlas(freesurferPath, 'r', 'BensonV')
% 
% ZQ 2023

if isempty(atlasName)
    print('Choose atlas (Wang/BensonV/BensonE/visf/HCP)');
else
    ky = atlasName;
end

switch ky

%-------------------------------
    case 'Wang'
        
    % Surface labels - Wang visual ROI
    % L Wang L, et.al. Probabilistic Maps of Visual Topography in Human Cortex. Cereb Cortex (2015) 

    surface_labels = MRIread(fullfile(freesurferPath, 'fsaverage','surf',...
                    [hemi 'h.wang15_mplbl.mgz']));
    vert_label = surface_labels.vol(:);

    area_label = { 'V1v', 'V1d', ...
                   'V2v', 'V2d', ...
                   'V3v', 'V3d', ...
                   'hV4', ...
                   'VO1', 'VO2', ...
                   'PHC1','PHC2',...    % parahippocampal cortex
                   'TO2', 'TO1', ... 
                   'LO2', 'LO1', ...
                   'V3B', 'V3A', ...
                   'IPS0','IPS1','IPS2','IPS3','IPS4','IPS5', ...
                   'SPL1','FEF'};       % superior parietal lobule

    cmap = [...
        0.9438    0.3910    0.2668
        0.1974    0.5129    0.7403
        0.5978    0.8408    0.6445
        0.3686    0.3098    0.6353
        0.9955    0.8227    0.4828
        0.8417    0.9409    0.6096
        0.6196    0.0039    0.2588
        0.8000    0.8000    0.4000
        0.3539    0.7295    0.6562
        0.9877    0.6154    0.3391
        0.8206    0.2239    0.3094
        0.9000    0.6000    1.0000
        0.3438    0.3910    0.3668
        0.6974    0.5129    0.8403
        0.1978    0.8408    0.7445
        0.3955    0.8227    0.5828
        0.9686    0.3098    0.7353
        0.2417    0.9409    0.7096
        0.2196    0.0039    0.7588
        0.2000    0.8000    0.5000
        0.9539    0.7295    0.7562
        0.2877    0.6154    0.4391
        0.2206    0.2239    0.4094
        0.9000    1.0000    0.7000
        0.2000    0.5000    0.2000];

%-------------------------------
    
    case 'BensonV'

    % Surface labels - Benson visual ROI
    % NC Benson & J Winawer. Bayesian analysis of retinotopic maps. eLife (2018) 

    surface_labels = MRIread(fullfile(freesurferPath, 'fsaverage', 'surf',...
        [hemi 'h.benson14_varea.mgz']));
    vert_label = surface_labels.vol(:);

    area_label = { 'V1', ...            % primary visual  
                   'V2', ...
                   'V3', ...
                   'hV4', ...           % ventral
                   'VO1', 'VO2', ...
                   'LO1', 'LO2', ...    % lateral
                   'TO1', 'TO2', ...    % temporal
                   'V3b', 'V3a'};       % dorsal

    cmap   =  [255   0   0; 
               255 128   0; 
               255 255   0; 
                 0 128 255;
                 0  76 153; 153  76   0;
                 0 102  51; 153 153   0;  
               255  51 153; 102   0 204;
                 0 255 255;   0 255   0]./255;

 %-------------------------------            
 
    case 'visf'

    % Surface labels - Probabilistic functional O-T visual (visf) atlas - KGS
    % M Rosenke, et al. A Probabilistic Functional Atlas of Human Occipito-Temporal Visual Cortex. Cerebral Cortex (2021)

    [~, vert_label, temp] = read_annotation(fullfile(freesurferPath, 'fsaverage','label',...
                            [hemi 'h.visfAtlas.annot']));

    % Labels are numbered in some way (nothing from 1-16).. fixing this
    for ii = 1 : size( temp.table, 1)   % 1 is 0/unknown/nolabel
        vert_label( find(vert_label == temp.table( ii, 5))) = ii - 1;
    end

    cmap = temp.table( 2:end, 1:3) ./ 255;
    area_label = temp.struct_names(2:end);
    
%-------------------------------

    case 'BensonE'
        
    % surface labels - Benson eccentricity
    % NC Benson & J Winawer. Bayesian analysis of retinotopic maps. eLife (2018) 

    surface_labels = MRIread(fullfile(freesurferPath, 'fsaverage', 'surf',...
        [hemi 'h.benson14_eccen.mgz']));
    vert_label = surface_labels.vol(:);

    % load the color map -logscale
    cmap = [];
    temp = logspace(0,3,90);
    temp1 = temp(end:-1:1);
    cm = jet(10^3);
    for ii = 1:90 
        cmap(ii, :) = cm(round(temp1(ii)), :);
    end

    area_label = [1 : max(vert_label)];
    
%-------------------------------

    case 'HCP'
        
    % Surface labels - HCP MMP1 Atlas
    % M Glasser, et al. A multi-modal parcellation of human cerebral cortex. Nature (2016)

    [~, vert_label, temp] = read_annotation(fullfile(freesurferPath, 'fsaverage','label',...
                              [hemi 'h.HCP-MMP1.annot']));

    % Labels are numbered in some way (nothing from 1-180).. fixing this
    for ii = 1 : size( temp.table, 1)
        vert_label( find(vert_label == temp.table( ii, 5))) = ii - 1;
    end

    cmap = temp.table( 2:end, 1:3) ./ 255;
    area_label = temp.struct_names(2:end);
            
        
    otherwise
        warning('Unexpected input. Expected: Wang/BensonV/BensonE/visf/HCP')
        
end

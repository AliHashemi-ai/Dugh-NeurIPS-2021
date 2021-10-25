function allplots_cortex4(sa, data_in, colorlimits, cm, unit, smooth, printfolder, varargin)

if length(data_in) == size(sa.cortex75K.vc, 1)
  data = data_in;
end

if isfield(sa.cortex75K, 'in_cort')
  if length(data_in) == length(sa.cortex75K.in_cort)
    data = nan*ones(size(sa.cortex75K.vc, 1), 1);
    data(sa.cortex75K.in_cort) = data_in;
  end
end

if length(data_in) == length(sa.cortex10K.in_from_cortex75K)
  data = data_in(sa.cortex10K.in_to_cortex75K_geod);
end

if isfield(sa.cortex10K, 'in_cort')
  if length(data_in) == length(sa.cortex10K.in_cort)
    data_ = nan*ones(length(sa.cortex10K.in_from_cortex75K), 1);
    data_(sa.cortex10K.in_cort) = data_in;
    data = data_(sa.cortex10K.in_to_cortex75K_geod);
  end
end

if length(data_in) == length(sa.cortex5K.in_from_cortex75K)
  data = data_in(sa.cortex5K.in_to_cortex75K_geod);
end

if isfield(sa.cortex5K, 'in_cort')
  if length(data_in) == length(sa.cortex5K.in_cort)
    data_ = nan*ones(length(sa.cortex5K.in_from_cortex75K), 1);
    data_(sa.cortex5K.in_cort) = data_in;
    data = data_(sa.cortex5K.in_to_cortex75K_geod);
  end
end

if length(data_in) == length(sa.cortex2K.in_from_cortex75K)
  data = data_in(sa.cortex2K.in_to_cortex75K_geod);
end

if isfield(sa.cortex2K, 'in_cort')
  if length(data_in) == length(sa.cortex2K.in_cort)
    data_ = nan*ones(length(sa.cortex2K.in_from_cortex75K), 1);
    data_(sa.cortex2K.in_cort) = data_in;
    data = data_(sa.cortex2K.in_to_cortex75K_geod);
  end
end

if length(data_in) == length(sa.cortex1K.in_from_cortex75K)
  data = data_in(sa.cortex1K.in_to_cortex75K_geod);
end

if isfield(sa.cortex1K, 'in_cort')
  if length(data_in) == length(sa.cortex1K.in_cort)
    data_ = nan*ones(length(sa.cortex1K.in_from_cortex75K), 1);
    data_(sa.cortex1K.in_cort) = data_in;
    data = data_(sa.cortex1K.in_to_cortex75K_geod);
  end
end


set(0,'DefaultFigureColor',[1 1 1])
printfolder = [printfolder '/'];
mkdir(printfolder)

res = '150';

if smooth 
    vc = sa.cortex75K.vc_smooth;
    sm = '_smooth';
else
    vc = sa.cortex75K.vc;
    sm = '';
end

surface_pars = struct('alpha_const', 1, 'colormap', cm, 'colorlimits', colorlimits, ...
  'showdirections', 0, 'colorbars', 0, 'dipnames', [], 'mymarkersize', 15, 'directions', [0 0 1 1 1 1], ...
  'printcbar', 1, 'userticks', []);

if length(varargin) > 0
   varargin1 = varargin{1};
else
    varargin1 = {};
end

if length(varargin) > 1
    input_pars = varargin{2};
    finames = fieldnames(input_pars);
    for ifi = 1:length(finames)
      surface_pars = setfield(surface_pars, finames{ifi}, getfield(input_pars, finames{ifi}));
    end
end
    
surface_pars.myviewdir = [-1 0 0];
figure; showsurface3(vc, sa.cortex75K.tri_left, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_left'], ['-r' num2str(res)], '-a2', '-transparent'); 

figure; showsurface3(vc, sa.cortex75K.tri_right, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_right_inner'], ['-r' num2str(res)], '-a2', '-transparent'); 


surface_pars.myviewdir = [1 0 0];

figure; showsurface3(vc, sa.cortex75K.tri_right, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_right'], ['-r' num2str(res)], '-a2', '-transparent'); 

figure; showsurface3(vc, sa.cortex75K.tri_left, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_left_inner'], ['-r' num2str(res)], '-a2', '-transparent'); 

surface_pars.myviewdir = [-1e-10 0 1];
surface_pars.directions = [1 0 1 1 0 0];

figure; showsurface3(vc, sa.cortex75K.tri, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_top'], ['-r' num2str(res)], '-a2', '-transparent'); 


surface_pars.myviewdir = [0 0 1];

figure; showsurface3(vc, sa.cortex75K.tri, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_top_upright'], ['-r' num2str(res)], '-a2'); 


surface_pars.myviewdir = [-1e-10 0 -1];
% 
figure; showsurface3(vc, sa.cortex75K.tri, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_bottom'], ['-r' num2str(res)], '-a2', '-transparent'); 


surface_pars.myviewdir = [0 1e-10 -1];

figure; showsurface3(vc, sa.cortex75K.tri, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_bottom_upright'], ['-r' num2str(res)], '-a2'); 

if isfield(surface_pars, 'printcbar') && surface_pars.printcbar
  figure; 
  hf = imagesc(randn(5)); colormap(cm)
  set(gca, 'clim', colorlimits, 'position', [0.1 0.1 0.6 0.8], 'visible', 'off')
  set(hf, 'visible', 'off')
  cb = colorbar; 
  set(cb, 'fontsize', 30)
  if ~isempty(surface_pars.userticks)
    set(cb, 'xtick', sort([colorlimits, surface_pars.userticks]))
  end
  ylabel(cb, unit)
  export_fig([printfolder 'cortex_cbar'], ['-r' num2str(res)], '-a2', '-transparent')  
end

% set(0,'DefaultFigureColor','remove')
% figure; 
% hf = imagesc(randn(5)); colormap(cm)
% set(gca, 'clim', colorlimits, 'position', [0.1 0.1 0.6 0.8], 'visible', 'off')
% set(hf, 'visible', 'off')
% cb = colorbar; 
% set(cb, 'fontsize', 30)
% ylabel(cb, unit)
% export_fig([printfolder 'slices_cbar'], ['-r' num2str(res)], '-a2')  

% set(0,'DefaultFigureColor',[1 1 1])
% 
% figure; showmri_transp3(sa.mri, struct('orientation', 'axial', ...
%     'trcut', 0.5, 'colormaps', {{cm}}, 'colorbars', 0, 'colorlimits', colorlimits), [sa.cortex75K.vc(sa.cortex5K.in_from_cortex75K, :) data(sa.cortex5K.in_from_cortex75K)], varargin1{:})
% export_fig([printfolder 'slices_axial'], '-r300', '-a2'); 
% 
% figure; showmri_transp3(sa.mri, struct('orientation', 'sagittal', ...
%     'trcut', 0.5, 'colormaps', {{cm}}, 'colorbars', 0, 'colorlimits', colorlimits), [sa.cortex75K.vc(sa.cortex5K.in_from_cortex75K, :) data(sa.cortex5K.in_from_cortex75K)], varargin1{:})
% export_fig([printfolder 'slices_sagittal'], '-r300', '-a2'); 
% 
% figure; showmri_transp3(sa.mri, struct('orientation', 'coronal', ...
%     'trcut', 0.5, 'colormaps', {{cm}}, 'colorbars', 0, 'colorlimits', colorlimits), [sa.cortex75K.vc(sa.cortex5K.in_from_cortex75K, :) data(sa.cortex5K.in_from_cortex75K)], varargin1{:})
% export_fig([printfolder 'slices_coronal'], '-r300', '-a2'); 

% close all












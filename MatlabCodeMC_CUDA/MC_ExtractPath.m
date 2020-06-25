function [path] = MC_ExtractPath(Sim, num_path, dpath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MC_ExtractPath.m
% 
% [path]=MC_ExtractPath(Sim,num_path,dpath)
% 
% This routine bins optical pathlength into num_path bins of dpath width.
% Optical pathlength is considered as L*refractive_index.
% The function returns also the mean optical pathlength of each bin.
% 
% Sim:      Simulation structure obtained by MC_ReadOut.m.
% num_path: number of pathlength
% dpath:    width of pathlength bin
%
% path:     dimension: num_path x (N_LAYERS+1).
%           column 1: number of photons in each bin.
%           columns(2->N_LAYERS): mean optical path of each bin.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 




%% READ SIM STRUCTURE %%
N_LAYERS = Sim.Sample.N_layers;
%N_DETECTORS=Sim.Detection.N_detectors;
Data = squeeze(Sim.Data(1, :, :));
%N_PHOTONS_MAX=Sim.Detection.N_photons_max;
Nref = Sim.Sample.prop(:, 1);

sampl_path  = (0:dpath:(num_path-1)*dpath)';
path = zeros(num_path, N_LAYERS+1);
%% BIN SIMULATION DATA %%
    
length = Data*Nref;  %cammini nei layers
[path(:, 1), bin] = histc(length, sampl_path);
% patch if happen that a bin value is zero: i put the photon at the end
if sum(bin==0)>0
    bin(bin==0) = numel(sampl_path);
end
for i=1:N_LAYERS
%     for ipath = 1:1:num_path
%         if (path(ipath, 1)>0)
%             path(ipath, i+1) = Nref(i).*mean(Data(bin==ipath, i));
%         end
%     end
         path(:,i+1)=accumarray(bin,length,[num_path 1],@mean);
         
     
  
end
%% patch...check
path(1,1) = 0;
path_2=path(:,2);
dummy=find(path_2==0);
path_2(dummy)=sampl_path(dummy)+dpath/2;
path(:,2)=path_2;


return;



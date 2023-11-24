function [time, counts, stdev, dmua, dmus,zmax,zmean] = MC_run(NPHOTONS,TMAX,RorT,radius,Rho,opt,mua,musp,n_chan,dt,SEED)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MC_run.m
% 
% [time counts stdev, dmua, dmus, Zmax, Zmean] = MC_run(NPHOTONS,TMAX,RorT,radius,Rho,opt,mua,musp,n_chan,dt)
% 
% This routine extracts run a single MC simulation.
% 
% NPHOTONS: number of photons to be received
% TMAX:     maximum propagation time
% RorT:     'R': reflectance, 'T': transmittance     
% radius:   radius of the cylinder
% Rho:      [rmin rmax] of the receiver
% opt:      each row for each layer:
%           opt = [1.0 0 MUS g thick;
%                  .. .. .. .. .. ]
% mua:      absorption vector (mua=[1 X Num_Layers])
% musp:     reduced scattering for rescaling [0 for disabling the feature]
% n_chan:   number of temporal windows (n_chan = 0 --> CW)
% dt:       width of the temporal window
% SEED:     integer number for seeding the random number generator. If not
%           provided the seed is random
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%   Run a single simulation
%
% NPHOTONS = 1e4;
% TMAX = 5000;      % ps
% Rho = [1.45 1.55];% cm
% mua = 0;
% g = 0.8;
% MUS = 10/(1-g);  % cm-1
% thick = 4;        % cm
% opt = [1.0 0 MUS g thick];

PLOT = 0;
MC_FuncGenInputFile('MCsingle.mci','MCsingle',NPHOTONS,TMAX,RorT,radius,Rho,opt);
if nargin < 11
    !~/Documents/MC_CUDA/CUDAMCML MCsingle.mci
else
    strcmd = ['~/Documents/MC_CUDA/CUDAMCML MCsingle.mci -S',num2str(SEED)];
    system(strcmd);
    %!~/Documents/MC_CUDA/CUDAMCML MCsingle.mci -S13131313
end

Sim = MC_ReadOut('MCsingle_1.mco');
[time, counts, stdev, dmua, dmus,zmax,zmean] = MC_ExtractSimulation(Sim,n_chan,dt,mua,musp,PLOT);
end


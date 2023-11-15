function [time, counts, stdev, dmua, dmus,zmax,zmean] = MC_run(NPHOTONS,TMAX,RorT,diam,Rho,opt,mua,musp,n_chan,dt)
%
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
MC_FuncGenInputFile('MCsingle.mci','MCsingle',NPHOTONS,TMAX,RorT,diam,Rho,opt);
!~/Documents/MC_CUDA/CUDAMCML MCsingle.mci
%!~/Documents/MC_CUDA/CUDAMCML MCsingle.mci -S13131313

Sim = MC_ReadOut('MCsingle_1.mco');
[time, counts, stdev, dmua, dmus,zmax,zmean] = MC_ExtractSimulation(Sim,n_chan,dt,mua,musp,PLOT);
end


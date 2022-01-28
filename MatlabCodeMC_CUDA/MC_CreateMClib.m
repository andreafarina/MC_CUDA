clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%% INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Files
FILE_OUT='VASCOVID.mclib';
FILE_MC_DIR='./';%'~/Documents/Code/MC_CUDA_Dev/Simulations/Thermal/';
FILE_MC_PREF='VASCOVID_';
FILE_MC_EXT='mco';

%% Geometric series
%q=sqrt(2);
q = 1.2;
MUS0=1;  %cm-1 reduced scattering coefficient
NUM_MUS=19;

%% Temporal windows
NUM_PATH = 8000; %10000;
dpath = 0.03; %cm/chan

%% Scaling factor
%  new_thick=lib_thick * k
% K=0.36/0.3;
%% Write Header
fid=fopen([FILE_MC_DIR FILE_OUT], 'wb');
fwrite(fid, NUM_MUS, 'int');
fwrite(fid, MUS0, 'double');
fwrite(fid, q, 'double');
fwrite(fid, NUM_PATH, 'int');
tic
% Theor=zeros(NUM_PATH, NUM_DOMAIN+1, NUM_MUS);

%% Read Simulation and Write on McLib File
config = '';
for imus=1:NUM_MUS
    imus
    FILE = [FILE_MC_DIR FILE_MC_PREF num2str(imus) '.' FILE_MC_EXT];
    Sim = MC_ReadOut(FILE);
    [path] = MC_ExtractPath(Sim, NUM_PATH, dpath);
    if(imus==1) 
        N_LAYERS = Sim.Sample.N_layers;
        r_min = Sim.Detection.det(1);
        r_max = Sim.Detection.det(2);
        if(r_min==0) 
            rho = 0;
        else
            rho = (r_min+r_max)/2;
        end
        if(Sim.Detection.Kind == 'R')
            next = Sim.Sample.n_up;
            config = 0;
        else
            next = Sim.Sample.n_down;
            config = 1;
        end
        fwrite(fid, N_LAYERS, 'int');
        fwrite(fid, config, 'int');
        fwrite(fid, next, 'double');
        fwrite(fid, rho, 'double');
        for idom=1:N_LAYERS
            fwrite(fid, Sim.Sample.prop(idom, 1), 'double'); %refindex
            fwrite(fid, Sim.Sample.prop(idom, 3), 'double'); %g factor
            fwrite(fid, Sim.Sample.prop(idom, 4), 'double'); %thickness
        end
    end
    for idom=1:N_LAYERS+1
        fwrite(fid, path(:, idom), 'float');
%        Theor(:, idom+1, imus) = path(:, idom);
    end
end
   
toc
fclose(fid);


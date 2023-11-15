function Sim=MC_ReadOut(FILE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MC_ReadOut.m
% 
% Sim=MC_ReadOut(FILE)
%
% This routine reads the rough data output FILE .mco from the MC_detectors.
% Data are scored in a structure Sim.
% 3D-array is stored in Sim.Data(i,j,k)
% 
% i -> N_DETECTORS
% j -> N_LAYERS
% k -> N_PHOTONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %% LOAD LIBRARY %%

fid=fopen(FILE,'r');
display('Simulation Data')

HEADER_DIM = fread(fid,1,'short');
VERSION = fread(fid,1,'uint8=>char');

N_UP = fread(fid,1,'float');
N_LAYERS=fread(fid,1,'short');

lay_prop=zeros(N_LAYERS,4);
display(sprintf('Layer\tn\tmus\tg\tthick\t'))
for i=1:N_LAYERS
    lay_prop(i,:)=fread(fid,4,'float');
    display(sprintf('%d\t%.2f\t%.2f\t%.2f\t%.2f',i,lay_prop(i,:))) 
end
N_DOWN=fread(fid,1,'float');

RT = fread(fid,1,'uint8=>char');
if (RT=='R')
    display('Reflectance')
else
    display('Transmittance')
end

N_DETECTORS=fread(fid,1,'short');
detectors=fread(fid,2*N_DETECTORS,'float');
for i=1:N_DETECTORS
    det_str=sprintf('%.2f\t%.2f',detectors(2*i-1),detectors(2*i));
    display(['r' num2str(i) '=' det_str])
end
%N_PHOTONS_MAX=fread(fid,1,'uint64');
%display(['Photons max per receiver=' num2str(N_PHOTONS_MAX)]);

N_PHOT_DETECTED = fread(fid, N_DETECTORS, 'uint64');
N_PHOT_LAUNCHED = fread(fid, N_DETECTORS, 'uint64');

for i=1:N_DETECTORS
    display(['Detector ' num2str(i) ': received photons = ' num2str(N_PHOT_DETECTED(i))...
        '; launched photons = ' num2str(N_PHOT_LAUNCHED(i))]);
end

% LETTURA CAMMINI %%%%%%

fseek(fid,HEADER_DIM,-1);

Data=zeros(N_DETECTORS, N_LAYERS, N_PHOT_DETECTED);
 for i=1:N_DETECTORS
     for j=1:N_LAYERS
        Data(i, j, :)=fread(fid, N_PHOT_DETECTED, 'float');
         
     end
 end
% LETTURA numeo interazioni dis cattering Kappa
Kappa = zeros(N_DETECTORS, N_LAYERS, N_PHOT_DETECTED);
 for i=1:N_DETECTORS
     for j=1:N_LAYERS
        Kappa(i, j, :)=fread(fid, N_PHOT_DETECTED, 'short');
         
     end
 end

% LETTURA Zmax
Zmax = zeros(N_DETECTORS,N_PHOT_DETECTED);
for i=1:N_DETECTORS
    Zmax(i, :)=fread(fid, N_PHOT_DETECTED, 'float');
end

% LETTURA sumZ somma degli z per il calcolo di <Z>
% LETTURA Zmax
sumZ = zeros(N_DETECTORS,N_PHOT_DETECTED);
for i=1:N_DETECTORS
    sumZ(i, :)=fread(fid, N_PHOT_DETECTED, 'float');
end
 %pippo=fread(fid,N_PHOTONS*N_LAYERS*N_DETECTORS,'float');

%N_PHOTONS_LAUNCHED=fread(fid,1,'uint64');
%display(['Launched photons=' num2str(N_PHOTONS_LAUNCHED)]);
%figure(7),plot(squeeze(Data))
%save('y3','Data');
fclose(fid);


%% CREATE STRUCTURE %%
%% Sim.File %%
Sim.File.Name=FILE;
Sim.File.Header_dim=HEADER_DIM;
Sim.File.Ver=VERSION;
%% Sim.Sample %%
Sim.Sample.N_layers=N_LAYERS;
Sim.Sample.prop=lay_prop;
Sim.Sample.n_up=N_UP;
Sim.Sample.n_down=N_DOWN;
%% Sim.Detection %%
Sim.Detection.N_detectors=N_DETECTORS;
Sim.Detection.Kind=RT;
Sim.Detection.det=detectors;
%Sim.Detection.N_photons_max=N_PHOTONS_MAX;
Sim.Detection.N_photons_det=N_PHOT_DETECTED;
Sim.Detection.N_photons_launched=N_PHOT_LAUNCHED;

%% Sim.Data %%
Sim.Data = Data;
Sim.Kappa = Kappa;
Sim.Zmax = Zmax;
Sim.sumZ = sumZ;

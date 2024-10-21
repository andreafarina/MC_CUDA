function [time, counts, stdev, dmua, dmus, Zmax, Zmean]=MC_ExtractSimulation(Sim,n_chan,dt,mua,musp,PLOT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MC_ExtractSimulation.m
% 
% [time counts stdev, dmua, dmus, Zmax, Zmean] = MC_ExtractSimulation(Sim,n_chan,dt,mua,musp,PLOT)
% 
% This routine extracts TR profiles, standard errors, derivatives,
% maximum and avarege penetration depths from Simulation Data.
% 
% Sim:      Simulation structure obtained by MC_ReadOut.m.
% n_chan:   number of temporal windows (if n_chan = 1 -> CW)
% dt:       width of the temporal window
% mua:      absorption vector (mua=[1 X Num_Layers])
% musp:     apply scaling rules [0 for disabling the feature, musprime!]
% PLOT:     flag for plotting output if =0 the plot is not shown
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 




%% READ SIM STRUCTURE %%
N_LAYERS = Sim.Sample.N_layers;
N_DETECTORS = Sim.Detection.N_detectors;
Data = Sim.Data;
Kappa = Sim.Kappa;
ZetaMax = Sim.Zmax;
ZetaMean = Sim.sumZ./reshape(sum(Kappa,2),N_DETECTORS,[]);
lay_prop = Sim.Sample.prop;
mus0 = lay_prop(:,2)';
g = lay_prop(:,3)';

if ((N_LAYERS)~=numel(mua))
    disp('Insert absorption for each layers as array');
    time = -1;counts = -1;stdev = -1;
    return
end

C = 0.0299792458; %cm/ps

if n_chan == 1
    % CW
    time = [0,dt];
else
    time = (0:dt:(n_chan-1)*dt);
end
refind = lay_prop(:,1);
v = C*ones(1,N_LAYERS)./refind';

Simul = zeros(N_DETECTORS,n_chan);

%% BIN SIMULATION DATA %%

for i=1:N_DETECTORS
    
    N_PHOTONS_DET = Sim.Detection.N_photons_det(i);
    length = reshape(Data(i,:,:),N_LAYERS,N_PHOTONS_DET);  %cammini nei layers
    k = reshape(Kappa(i,:,:),N_LAYERS,N_PHOTONS_DET);  %scattering events nei layers
    zmax = reshape(ZetaMax(i,:),1,N_PHOTONS_DET);
    zmean = reshape(ZetaMean(i,:),1,N_PHOTONS_DET);
    weight_abs = exp(-mua*length);
     if ((nargin>4)&&prod(musp~=0))
            mus = musp./(1-g);
            logws = log(mus./mus0) * k;
            logws2 = -(mus-mus0) * length;
            weight_sca = exp(logws + logws2);
        else
            weight_sca = 1;
     end
     
        
    weight = weight_abs.* weight_sca;
    timep = (1./v)*length;
    
    [n_bin,bins] = histc(timep,time);
    %     if ~all(bins)
    %         disp('Increase the number of temporal windows to include all data');
    %         return;
    %     end
    
    if ~all(bins)
        weight(bins==0)=[];
        length(:,bins==0)=[];
        k(:,bins==0)=[];
        zmax(:,bins==0) = [];
        zmean(:,bins==0) = [];
        bins(bins==0)=[];
    end
    
    if n_chan == 1 %CW
        time(1) = [];
    end
    
    Sumweight(i,:) = accumarray(bins',weight,[numel(time) 1],@sum);
    Sumweight2(i,:) = accumarray(bins',weight.^2,[numel(time) 1],@sum);
    Zmax(i,:) = accumarray(bins',weight.*zmax,[numel(time) 1],@sum)'./Sumweight(i,:);
    Zmean(i,:) = accumarray(bins',weight.*zmean,[numel(time) 1],@sum)'./Sumweight(i,:);
    
end

time = time+dt/2;

%% FLUENCE cm-2*ps-1 %%
N_PHOTONS_LAUNCHED = Sim.Detection.N_photons_launched;
R_DET = Sim.Detection.det;

Simul2 = zeros(N_DETECTORS,n_chan);

for i = 1:N_DETECTORS
    fact = pi*dt*(R_DET(2*i).^2-R_DET(2*i-1).^2)*N_PHOTONS_LAUNCHED(i);
    
    Simul(i,:) = Sumweight(i,:)./fact;
    Simul2(i,:) = Sumweight2(i,:)./(fact.^2);
    
    % this is used only for the plot legend
    if R_DET(2*i-1) == 0
        R_AV = 0;
        DR = R_DET(2*i);
    else
        R_AV = (R_DET(2*i)+R_DET(2*i-1))/2;
        DR = (R_DET(2*i)-R_DET(2*i-1))/2;
    end
end

counts = Simul;
    
%% Standard deviation
if nargout > 2
% metodo corretto
stdev = sqrt((Simul2.^2 - Simul.^2)./(n_bin-1));
% metodo approssimato con Bernoulli Poisson
%stdev1 = counts./sqrt(n_bin);
stdev(isnan(stdev)) = 0;
end

%% calculation of derivative (the scattering scaling is not considered)
if nargout > 3
    meanLength = zeros(N_LAYERS,n_chan);
    meanKmus = zeros(N_LAYERS,n_chan);
    % dmua
    l_weight = repmat(weight,N_LAYERS,1).* length;
    for i = 1:N_LAYERS
        meanLength(i,:) = accumarray(bins',l_weight(i,:),[numel(time) 1],@sum);
    end
    meanLength = meanLength./Sumweight;
    dmua = -meanLength.*counts;
    % dmus
    Kmus_weight = repmat(weight,N_LAYERS,1).* (k./mus0');
    for i = 1:N_LAYERS
        meanKmus(i,:) = accumarray(bins',Kmus_weight(i,:),[numel(time) 1],@sum);
    end
    meanKmus = meanKmus./Sumweight;
    dmus = meanKmus.*counts + dmua;
    dmua(isnan(dmua)) = 0;
    dmus(isnan(dmus)) = 0;
end

%% Plot
if (nargin > 5)&&(PLOT==1)
    stringa_leg=['r = ' num2str(R_AV,'%.2f')  ' +/- ' num2str(DR,'%.2f') ' cm'];
    m_legend(i,:)=stringa_leg;
    figure,semilogy(time,counts,'.'),grid
    y_label_str='Fluence [cm^{-2}ps^{-1}]';
    xlabel('time [ps]'), ylabel (y_label_str)
    legend(m_legend);
    figure,plot(time,[Zmax',Zmean',Zmax'/2]),grid
    xlabel('time [ps]'), ylabel ('<Zmax> [cm]'),
    legend('<Zmax>','<Zmean>','<Zmax>/2')
else
    return
end

return;



function [time counts stdev]=MC_ExtractSimulation(Sim,n_chan,dt,mua,N_TOT,PLOT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MC_ExtractSimulation.m
% 
% [time counts stdev]=MC_ExtractSimulation(Sim,n_chan,dt,mua,N_TOT,PLOT)
% 
% This routine extracts TR profiles from Simulation Data.
% 
% Sim:      Simulation structure obtained by MC_ReadOut.m.
% n_chan:   number of temporal windows
% dt:       width of the temporal window
% mua:      absorption vector (mua=[1 X Num_Layers]
% N_TOT:    total counts (Area) if =0 then it uses the effective photons
% PLOT:     flag for plotting output if =0 the plot is not shown
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 




%% READ SIM STRUCTURE %%
N_LAYERS=Sim.Sample.N_layers;
N_DETECTORS=Sim.Detection.N_detectors;
Data=Sim.Data;
lay_prop=Sim.Sample.prop;

if ((N_LAYERS)~=numel(mua))
    disp('Insert absorption for each layers as array');
    time=-1;counts=-1;stdev=-1;
    return
end

C=0.0299792458; %cm/ps

time=(0:dt:(n_chan-1)*dt);%+dt/2;
refind=lay_prop(:,1);
v=C*ones(1,N_LAYERS)./refind';

Simul=zeros(N_DETECTORS,n_chan);

%% BIN SIMULATION DATA %%

for i=1:N_DETECTORS
    
    N_PHOTONS_DET=Sim.Detection.N_photons_det(i);
    length=reshape(Data(i,:,:),N_LAYERS,N_PHOTONS_DET);  %cammini nei layers
    % length=length(:,1:N_PHOTONS_DET);
    weight=exp(-mua*length);
    timep=(1./v)*length;
    
    [n_bin,bins]=histc(timep,time);
    %     if ~all(bins)
    %         disp('Increase the number of temporal windows to include all data');
    %         return;
    %     end
    
    if ~all(bins)
        weight(bins==0)=[];
        bins(bins==0)=[];
    end
    
    %tic
    Simul(i,:)=accumarray(bins',weight,[numel(time) 1],@sum);
    
    %     for k=1:N_PHOTONS_DET
    %         chan=bins(k);
    %         if chan>0
    %             Simul(i,chan)=Simul(i,chan)+weight(k);
    %         end
    %     end
    %     toc
    
end

time=time+dt/2;

%% FLUENCE cm-2*ps-1 %%
N_PHOTONS_LAUNCHED=Sim.Detection.N_photons_launched;
R_DET=Sim.Detection.det;
for i=1:N_DETECTORS
    k=pi*dt*(R_DET(2*i).^2-R_DET(2*i-1).^2)*N_PHOTONS_LAUNCHED(i);
    %k=1;
    Simul(i,:)=Simul(i,:)./k;
    if R_DET(2*i-1)==0
        R_AV=0;
        DR=R_DET(2*i);
    else
        R_AV=(R_DET(2*i)+R_DET(2*i-1))/2;
        DR=(R_DET(2*i)-R_DET(2*i-1))/2;
    end
end

%% Normalization
if (N_TOT==0)
    counts=Simul;
else
    S=sum(Simul,2);
    counts=Simul./repmat(S,1,n_chan)*N_TOT;
    %semilogy(time,Counts,'.'),grid
    %return;
end
%% Standar deviation
stdev=counts./sqrt(n_bin);
stdev(isnan(stdev))=0;
%% Plot

if (PLOT==1)
    stringa_leg=['r = ' num2str(R_AV,'%.2f')  ' +/- ' num2str(DR,'%.2f') ' cm'];
    m_legend(i,:)=stringa_leg;
    figure,semilogy(time,counts,'.'),grid
    if (N_TOT==0)
        y_label_str='Fluence [cm^{-2}ps^{-1}]';
    else
        y_label_str='Photons';
    end
    xlabel('time [ps]'), ylabel (y_label_str)
    legend(m_legend);
else
    return
end

return;



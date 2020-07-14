function [time TOFs]=MC_PartialTOF(Sim,n_chan,dt,mua,kind)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MC_PartialTOF.m
% 
% [time TOFs]=MC_PartialTOF(Sim,n_chan,dt,mua)
% 
% This routine calculates the mean time of flight (TOF) spent in each layer
% for n_chan temporal windows large dt. It is possible also to add
% absorption (mua).
% 
% Sim:      Simulation structure obtained by MC_ReadOut.m.
% n_chan:   number of temporal windows
% dt:       width of the temporal window
% mua:      absorption vector (mua=[1 X Num_Layers]
% kind:     'time' or 'path' (default 'time')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 




%% READ SIM STRUCTURE %%
N_LAYERS=Sim.Sample.N_layers;
N_DETECTORS=Sim.Detection.N_detectors;
Data=Sim.Data;
lay_prop=Sim.Sample.prop;

if ((N_LAYERS)~=numel(mua))
    disp('Insert absorption for each layers as array');
    time=-1;TOFs=-1;
    return
end

if nargin==4 kind='time';end;

C=0.0299792458; %cm/ps

time=(0:dt:(n_chan-1)*dt);%+dt/2;
refind=lay_prop(:,1);
v=C*ones(1,N_LAYERS)./refind';

TOFs=zeros(N_LAYERS,n_chan);
Simul=zeros(1,n_chan);
%% BIN SIMULATION DATA %%

for i=1:1
    
    N_PHOTONS_DET=Sim.Detection.N_photons_det(i);
    length=reshape(Data(i,:,:),N_LAYERS,N_PHOTONS_DET);  %cammini nei layers
    weight=exp(-mua*length);
    %rep_weight=repmat(weight,[N_LAYERS 1]);
    timep=(1./v)*length;
    
    [~,bins]=histc(timep,time);
    if ~all(bins)
        weight(bins==0)=[];
        length(:,bins==0)=[];
        bins(bins==0)=[];
    end
    
    Simul(i,:)=accumarray(bins',weight,[numel(time) 1],@sum);
    l_weight=length.*repmat(weight,N_LAYERS,[]);
    for i=1:N_LAYERS
        TOFs(i,:)=accumarray(bins',l_weight(i,:),[numel(time) 1],@sum);
    end
%     tic;
%     for k=1:N_PHOTONS_DET
%         chan=bins(k);
%         if chan>0
%             TOFs(:,chan)=TOFs(:,chan)+length(:,k).*weight(:,k);
%             %Simul(i,chan)=Simul(i,chan)+weight(k);
%         end
%     end
%  toc;
   
end

TOFs=TOFs./repmat(Simul,[N_LAYERS 1]);
if kind=='time'
    TOFs=TOFs./repmat(v',[1 n_chan]);
    y_string='Time of flight (ps)';
else
    y_string='Path length (cm)';
end

time=time+dt/2;

%hold on
array_lay=1:N_LAYERS;
 str_legend = strtrim(cellstr(int2str(array_lay.')));
plot(time,TOFs),xlabel('Time (ps)'),ylabel(y_string),legend(str_legend),grid


%% FLUENCE cm-2*ps-1 %%
% N_PHOTONS_LAUNCHED=Sim.Detection.N_photons_launched;
% R_DET=Sim.Detection.det;
% for i=1:N_DETECTORS
%     k=pi*dt*(R_DET(2*i).^2-R_DET(2*i-1).^2)*N_PHOTONS_LAUNCHED(i);
%     %k=1;
%     Counts(i,:)=Counts(i,:)./k;
%      if R_DET(2*i-1)==0
%          R_AV=0;
%          DR=R_DET(2*i);
%      else
%          R_AV=(R_DET(2*i)+R_DET(2*i-1))/2;
%          DR=(R_DET(2*i)-R_DET(2*i-1))/2;
%      end
% end   
%     stringa_leg=['r = ' num2str(R_AV,'%.2f')  ' ± ' num2str(DR,'%.2f') ' cm'];
%     m_legend(i,:)=stringa_leg;
% end
% semilogy(time,Counts,'.'),grid
% 
%  xlabel('time [ps]'), ylabel ('Fluence [cm^{-2}ps^{-1}]')
%  legend(m_legend);

return;



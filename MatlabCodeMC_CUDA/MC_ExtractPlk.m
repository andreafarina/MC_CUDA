function [l,k,plk] = MC_ExtractPlk(Sim,NbinK,dk,NbinL,dl,norm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MC_ExtractPlk.m
% 
% [l,k,plk] = MC_ExtractPlk(sim,NbinK,dk,NbinL,dl,normalized)
% 
% This routine extracts the plkMC  from Simulation Data.
% 
% Sim:      Simulation structure obtained by MC_ReadOut.m.
% NbinK:    number of K bins
% dk:       width of a K bin
% NbinL:    number of L bins
% dl:       width of a L bin
% norm:     flag for normalizing the plkMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



L = squeeze(sum(Sim.Data));
K = squeeze(sum(Sim.Kappa));
%Lbin = linspace(0,max(L),NbinL + 1);
%Kbin = linspace(0,max(K),NbinK + 1);
Lbin = (0:NbinL) * dl;
Kbin = (0:NbinK) * dk;

[~,~,~,binEll,binK] = histcounts2(L,K,Lbin,Kbin);
if ~all(binEll)
   %weight(bins==0)=[];
   binEll(binEll==0)=[];
   binK(binK==0)=[];
end
plk = accumarray([binEll,binK],1,[numel(Lbin),numel(Kbin)] - 1);

%dl = Lbin(2) - Lbin(1);
%dk = Kbin(2) - Kbin(1);
l = Lbin + dl/2;
k = Kbin + dk/2;
l(end) = []; k(end) = [];
plk = plk';
plkCW = sum(plk,2);
if norm
    area = trapz(k,plk,1);
else
    area = 1;
end
plk = plk./area;

figure,subplot(1,2,1),imagesc(l,k,plk),xlabel('L(cm)'),ylabel('K')
subplot(1,2,2),plot(k,plk,k,plkCW,'r'),xlabel('k')

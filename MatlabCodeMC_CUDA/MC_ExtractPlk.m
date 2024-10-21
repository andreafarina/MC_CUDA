function [l,k,plk] = MC_ExtractPlk(Sim,dk,dl,norm)
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



L = squeeze(sum(Sim.Data,2));
K = squeeze(sum(Sim.Kappa,2));

Lbin = (0:dl:(max(L)+dl));
Kbin = (0:max(K)) + dk/2;

[plk,~,~,binEll,binK] = histcounts2(L,K,Lbin,Kbin);
if ~all(binEll)
   %weight(bins==0)=[];
   binEll(binEll==0)=[];
   binK(binK==0)=[];
end

l = Lbin + dl/2;
k = Kbin + dk/2;
l(end) = []; k(end) = [];
plk = plk';
plkCW = sum(plk,2) * dl;

if norm
    area = trapz(k,plk,1);
    plkCW = plkCW./sum(plkCW);
else
    area = 1;
end
plk = plk./area;

figure,subplot(1,2,1),imagesc(l,k,plk),xlabel('L(cm)'),ylabel('K')
subplot(1,2,2),plot(k,plk,k,plkCW,'k'),xlabel('k')

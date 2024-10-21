%% Creazione del file di simulazione %%
FILE = 'Bone05.mci';
fid = fopen(FILE,'w');
VER = 1.2;

%% Data 
FILE_OUT_PREFIX = 'Simulations/VAMSHI/Bone05';%'/media/DATA/andreafarina/MC/Pancreas';%'Simulations/Abdomen/Abdomen_3cm';
EXT_OUT = 'mco';
NRUN = 54;
NPHOTONS = 1e6;
TMAX = 2.0345*4096;%2000 %ps
RorT = 'R';
Radius = 2.5;%2.54/2;
%Rho=[0.45 0.55;0.95 1.05;1.45 1.55;1.95 2.05];
Rho = 0.5 + [-0.05 0.05];
NDET = size(Rho,1);



n_up = 1.;%refractive index up
n_down = 1.;%refractive index down

MUSp0 = 0.5;
g = 0.8; %anisotropy factor
%q=sqrt(2);
q = 1.1;
MUS0 = MUSp0/(1-g);

%opt0=[1.33 0 MUS0 g 1.0
%    1.54 0 MUS0 g 1.0];
opt0 = [1.40 0 MUS0 g 1.2];

NLAY = size(opt0,1);

%% creazione header
buffer='############################################## \n\n';
fprintf(fid,buffer);
buffer=sprintf('%.1f\t\t\t\t\t\t\t\t\t\t\t# file version \n',VER);
fprintf(fid,buffer);
buffer=[num2str(NRUN) '\t\t\t\t\t\t\t\t\t\t\t' '# number of runs \n\n'];
fprintf(fid,buffer);
%% specify data for run n
for i=1:NRUN
    buffer=['#### SPECIFY DATA FOR RUN ' num2str(i) '\n'];
    fprintf(fid,buffer);
    
    buffer='##InParm\t\t\t\t\t\t\t\t\t';
    comment='# Input parameters. cm is used.\n';
    fprintf(fid,[buffer comment]);
    
    buffer=[FILE_OUT_PREFIX '_' num2str(i) '.' EXT_OUT '\t' 'A'];
    comment='\t\t\t\t\t\t\t\t# output file name, NO ASCII.\n';
    fprintf(fid,[buffer comment]);
    
    buffer=num2str(NPHOTONS);
    comment='\t\t\t\t\t\t\t\t\t\t# Number of photons per detectors\n';
    fprintf(fid,[buffer comment]);
    
    buffer=num2str(TMAX);
    comment='\t\t\t\t\t\t\t\t\t\t# Time Max (ps)\n';
    fprintf(fid,[buffer comment]);
    
    buffer=RorT;
    comment='\t\t\t\t\t\t\t\t\t\t\t# R=Reflectance, T=Transmittance\n';
    fprintf(fid,[buffer comment]);
    
    buffer=sprintf('%.4f',Radius);
    comment='\t\t\t\t\t\t\t\t\t\t# Radius of the cylinder\n\n';
    fprintf(fid,[buffer comment]);
    
    buffer=num2str(NDET);
    comment='\t\t\t\t\t\t\t\t\t\t\t# Numberes of detectors\n';
    fprintf(fid,[buffer comment]);
    
    comment='#r1\tr2\t\t\t\t\t\t\t\t\t\t# One line for each detector (cm)\n';
    fprintf(fid,comment);
    
    for j=1:NDET
        buffer=sprintf('%f\t%f\n',Rho(j,1),Rho(j,2));
        fprintf(fid,buffer);
    end
    
    fprintf(fid,'\n\n');
    
    %% OPT PARAMETERS
    buffer=num2str(NLAY);
    comment='\t\t\t\t\t\t\t\t\t\t\t# Number of layers\n';
    fprintf(fid,[buffer comment]);
    comment='#n\t\tmua\t\tmus\t\t\tg\t\td\t\t# One line for each layer\n';
    fprintf(fid,comment);
   
    buffer=sprintf('%.2f',n_up);
    comment='\t\t\t\t\t\t\t\t\t\t# n for medium above and lateral \n';
    fprintf(fid,[buffer comment]);
    
    for k=1:NLAY
        buffer=sprintf('%.2f\t%.2f\t%f\t%.2f\t%.3f\t# layer %i\n',...
        opt0(k,1),opt0(k,2),opt0(k,3)*q^(i-1),opt0(k,4),opt0(k,5),k');
        %opt0(k,1),opt0(k,2),opt0(k,3),opt0(k,4),opt0(k,5)*i,k');
        fprintf(fid,buffer);
    end
    
    buffer=sprintf('%.2f',n_down);
    comment='\t\t\t\t\t\t\t\t\t\t# n for medium below\n';
    fprintf(fid,[buffer comment]);
    
    fprintf(fid,'\n\n');
        
    
    
    
    
    
    
    
end
    
    











fclose(fid)


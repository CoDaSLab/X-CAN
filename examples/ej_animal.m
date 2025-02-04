%% XCAN example of animal data
% J. Camacho, E. Acar, M. Rasmussen, R. Bro. Cross-product Penalized 
% Component Analysis (XCAN), Submitted to Chemometrics and Intelligent 
% Laboratory Systems, 2019.
%
% Needs the MEDA Toolbox (v1.6) and the XCAN (path should be properly set)
%
% If you use these data, please add a reference to the below paper:
%
% R. Bro, E.E. Papalexakis, E. Acar, N.D. Sidiropoulos. Coclustering-
% a useful tool for chemometrics. Journal of Chemometrics 26(6):256-263
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 04/Feb/2025
%
% Copyright (C) 2025  University of Granada, Granada
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

close all
clear
clc


%% Data viz

[Xa,txt] = xlsread('data/Animaldata.xlsx');
vars_l = txt(1,2:end);
obs_l =  txt(2:end,1);

% PCA scores and loadings
scoresPca(Xa,'PCs',1:2,'Preprocessing',2,'ObsLabel',obs_l);
loadingsPca(Xa,'PCs',1:2,'Preprocessing',2,'VarsLabel',vars_l);

[X,m] = preprocess2D(Xa,'Preprocessing',2);
Xa = X+ones(size(X,1),1)*m;

%% Cross-product (XP) matrices

% XP matrices with no thresholding

XXt = crossprod(X');
plotMap(XXt);
ylabel('XXt','FontSize',20)

XtX = crossprod(X);
plotMap(XtX);
ylabel('XtX','FontSize',20)

% XP matrices after thresholding

thr = 0.5; 
XXt = crossprod(X');
r = find((XXt)<thr);
XXt(r) = 0;
plotMap(XXt);
ylabel('XXt','FontSize',20)

thr = 0.5; 
XtX = crossprod(X);
r = find(abs(XtX)<thr);
XtX(r) = 0;
plotMap(XtX);
ylabel('XtX','FontSize',20)


%% XCAN: X-thresholding

L=0.01;
lambdas = [L;L];
pcs = [9];

Data = X;

for i=1:size(lambdas,2)
    
   [XP, XT] = xcan(Data,1:pcs(i),XtX,lambdas(1,i),XXt,lambdas(2,i));     

   varX = trace(Data'*Data);
   for j = 1:pcs(i)
       sPjind = find(abs(XT(:,j))==max(abs(XT(:,j))),1);
       sPj = sign(XT(sPjind,j));
       
       Xp = (XT(:,j)*pinv(XT(:,j)'*XT(:,j))*XT(:,j)')*Data*(XP(:,j)*pinv(XP(:,j)'*XP(:,j))*XP(:,j)');
       varP = trace(Xp'*Xp);
        
       f = figure, subplot(2,1,1),
       title(sprintf('XC: %i, variance: %.2f%%',j,100*(varP)/varX),'FontSize',20)
       hold on,  bar(sPj*XT(:,j)), a=get(f,'Children'); a.XTick = 1:length(obs_l); a.XTickLabel =  obs_l; a.XTickLabelRotation = 45; box on
       
       ylabel('Scores','FontSize',16), subplot(2,1,2), bar(sPj*XP(:,j)), ylabel('Loadings','FontSize',16),
       a=get(f,'Children');  a(1).XTick = 1:length(vars_l); a(1).XTickLabel =  vars_l; a(1).XTickLabelRotation = 45; box on
       axis([0.5,length(XP(:,j))+0.5,-1.1,1.1]),     
    end
        
   E = Data - XT*XP';        
   varE = trace(E'*E);
   ResidualVar = 100*(varE)/varX 
end




clc;clear;close all
warning off

addpath('./src')
addpath('./Evaluation')

%% Load Data
n=96        % Main-Features 1)colon=19 2) LSCC=60 3) KRCCC=4 4)GBM=5 5)BIC=38

%load Gene 
%load matlab1
%load matlab2
%load matlab3
%load label
%load Methy 
%load Mirna
[BRCAGE, BRCAMETH, BRCAMIR, SUR]=BRCA;
%% Normalize data 
%[I1,I2,I3, Label] = Data_Simulation(0.1,1,1);
  %  name = 'synthetic1';
%[t1] = Standard_Normalization(I1');
%[t2] = Standard_Normalization(I2');
%[t3] = Standard_Normalization(I3');
%clear I1  I2 I3
% Gene
%[X1,X2,X3, Label] = Data_Simulation(0.2,2,1);
 %   name = 'synthetic1';
   % X1 = Standard_Normalization(X1');
    %X2 = Standard_Normalization(X2');
%% Z-score and PCA 
%[I1,C1] = pca_whit(In1',m);
[I1] =  pca_whit(BRCAGE',n);
% Methy
%[I2,C2] = pca_whit(In2',m);
[I2] =  pca_whit(BRCAMETH',n);
% Mirna
%[I3,C3] = pca_whit(In3',m);
[I3] =  pca_whit(BRCAMIR',n);
%% construct the graphs for each data type
options = [];
options.NeighborMode = 'KNN';
options.k = 30; 
options.WeightMode = 'HeatKernel';
options.t = 100;

%% Construct W
W1  = full(constructW(I1,  options));
W2  = full(constructW(I2, options));
W3  = full(constructW(I3,options));

W_Total=cat(3,W1,W2,W3);

%% Clustering % set the number of clusters
%for K=5
%[idx_COMBINE,Vk, Lnew] = sc_ml(W_Total, K, 0.3);

%Tdisp('Clustering = ')
%disp(K)
%Wx=Vk*Vk';
%outs = SpectralClustering(Wx,K);

%ConcordanceMatrix = Concordance_Network_NMI({Wx,W1,W2,W3},K)

%displayimagesc(ConcordanceMatrix)
%title(['Concordance Matrix at K = ',num2str(K)])

%displayClusters(Wx,outs)
%title(['K = ',num2str(K)])

%[K1, K2, K12,K22] = Estimate_Number_of_Clusters_given_graph(Wx, [1 5:5:15]);
%fprintf('The best number of clusters according to eigengap is %d\n', K1);
%fprintf('The best number of clusters according to rotation cost is %d\n', K2);

%disp('*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\')
%end
K = 5; 
[idx_COMBINE,Vk, Lnew] = sc_ml(W_Total, K, 0.5);

displayClusters(Vk*Vk',idx_COMBINE);

%[nmi_value,ACC,f,p,r,Purity,AR,RI,MI,HI,MIhat]= Cluster_Evaluation(idx_COMBINE,Label) 
%ALI=1;
% 
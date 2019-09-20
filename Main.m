%Jezabel Guzman
%Abril, 2016

%This script executes the proposed method from my thesis project
N=20;    %number of total videos
Natt=15+1; %number of attacks
Nk=5;   %number of keyframes per video (max)
Metrics(1:Natt)=struct('TP', 0, 'TN', 0, 'FP', 0, 'FN', 0,'TPR', 0, 'PPV', 0, 'F1', 0, 'F05', 0);
FP=cell(N,Natt+1);
RetrVideo(1:N,1:Natt)=struct('ID', 0, 'Kf', 0, 'SimTh', 0, 'SimCC', 0, 'SimORB', 0, 'SimSSM', 0, 'SimML', 0);

for v=1:N
   %OFFLINE STAGE
   [V,A] = ReadVideoInDataSet(v,DB);%Read, process and segment the video 'v' in database DB
   FP{v,1} = FPExtraction(V,A);       % Extracttion of fingerprints
   %ONLINE STAGE
   for att=2:Natt
       [Vatt,Aatt] = AttackVideo(V,A);
       FP{v,att} = FPExtraction(V,A);
   end  
end

%Searching and Matching process
for v=1:N
    for att=1:Natt
       RetrVideo(v,att)=SearchAndMatch(FP{v,att}, FP{:,1}); 
       Metrics(Natt)=MatchingComparison(Metrics(Natt));
    end
end

Metrics=GetMetrics(Metrics);
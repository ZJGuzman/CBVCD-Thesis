function [ FP,LSHmp, LSHmp_Index ] = FPExtraction( V,A )
%FPEXTRACTION Extracts the set of fingerprints Th(V), CC(V), ORB(V), SSM(A)
%   V: is the set of visual keyframes
%   A: is the set of audio segment of the respective visual keyframes
    
    matcherORB = cv.DescriptorMatcher('FlannBased',...
        'Index',  {'LSH', 'TableNumber', 10, 'KeySize',10, 'MultiProbeLevel',2},...
        'Search', {'Sorted', true});
    nameMatFile='GuzmanDB570.mat';
    LUT_DBindexORB=0;

    [M,N,K]=size(V);
    Al=length(A);
    FP(1:K)=struct('Vk', [], 'Th', [], 'CC', [], 'ORB', [], 'SSM', []);
    
    for k=1:K
       FP.Vk=V(:,:,k);
       if k==1 %to print execution times
           tic;
           [FP.Th,GV]=GetTh(V(:,:,k));
           disp('Time to get Th: ');
           toc;
           tic;
           FP.CC=GetCC(V(:,:,k));
           disp('Time to get CC: ');
           toc;
           tic;
           FP.ORB=GetORB(GV);
           disp('Time to get ORB: ');
           toc;
           tic;
           FP.SSM=GetSSM(A(:,k));
           disp('Time to get SSM: ');
           toc;
       else
           [FP.Th,GV]=GetTh(V(:,:,k));
           FP.CC=GetCC(V(:,:,k));
           FP.ORB=GetORB(GV);
           FP.SSM=GetSSM(A(:,k));
       end
       %ORB index
    end
    
    function [Th,GV]=GetTh(Kc)
        %grayscale conversion
        GKs=rgb2gray(Kc);
        GV=normalizeG(GKs);
        
                figure; imshow(Im2);
        Im3=0.5*Im2(:,1:m_r,:)+0.5*(fliplr(Im2(:,m_r+1:m_r*2,:)));
        figure; imshow(Im3);

        %rgb2bin 
        %figure; imshow(gray2bin(rgb2gray(Im1),'psk',256));

        GVf=rgb2gray(Kc);
         Kf_T=imresize(kf_ref, [m_r m_r*2]);
        Kf=0.5*Kf_T(:,1:m_r)+0.5*(fliplr(Kf_T(:,m_r+1:m_r*2)));
        
    end
    function Gn=normalizeG(G)
        %normalization of G
        % block division (4 blocks) 
        [x,y]=size(G);
        b11=double(G(1:floor(x/2), 1:floor(y/2)));
        b12=double(G(1:floor(x/2), floor(y/2)+1:y));
        b21=double(G(floor(x/2)+1:x, 1:floor(y/2)));
        b22=double(G(floor(x/2)+1:x, floor(y/2)+1:y));

        B11=uint8(255*mat2gray(b11));
        B12=uint8(255*mat2gray(b12));
        B21=uint8(255*mat2gray(b21));
        B22=uint8(255*mat2gray(b22));
        Gn=[B11 B12; B21 B22];

    end
    
    tic;
        T=FoldFrame(Tbig,m);
        tm1=toc+tm1;
        tic;
        T_Th=imresize(T, [30 30]);
        tm2=toc+tm2;
        tic;
        %T_DCT=DCTbyBlocks(T,m,wi);
        T_CC=GetFpCC(TC);
        tm3=toc+tm3;
        tic;
        T_ORB=getORB(T); 
        tm4=toc+tm4;
        if i==1
            idxK=idxK+1;
            T1=Tbig;
            FP(idxK).TS_ini=i; %save number of initial and ending frames (video segment) of the keyframe
            FP(idxK).TS_end=i+window-1;
            FP(idxK).Kf=T;
            FP(idxK).Th=T_Th;
            FP(idxK).CC=T_CC;
            FP(idxK).ORB=T_ORB;
            Times.Downs=tm0;
            Times.Fold=tm1;
            Times.Th=tm2;
            Times.CC=tm3;
            Times.ORB=tm4;
        else
            mergeF=false;
            if corr2(Tbig,T1)>= Theta
                C1=corr2(T_Th, FP(idxK).Th);
                C2=1-distBin(T_CC, FP(idxK).CC);
                if C1>=0.9 && C2>=0.9
                    mergeF=true;
                elseif C1>= Theta && C2 >=Theta
                    mergeF=true;
                    if ~isempty(T_ORB)&&  ~isempty(FP(idxK).ORB)
                        indexPairs= matchFeatures(T_ORB.Features, FP(idxK).ORB.Features);
                        d_ORB=length(indexPairs)/min(length(T_ORB.Features),length(FP(idxK).ORB.Features));
                        if d_ORB<Theta
                            mergeF=false;
                        end
                        
                    end
                end
            end
            if ~mergeF
                idxK=idxK+1;
                FP(idxK).TS_ini=i; %save number of initial and ending frames (video segment) of the keyframe
                FP(idxK).ORB=T_ORB;
                FP(idxK).TS_end=i+window-1;
                FP(idxK).Kf=T;
                FP(idxK).Th=T_Th;
                FP(idxK).CC=T_CC;
                        
                T1=Tbig;
            else   %merge keyframes (re-arrange timestamp and get the better FP
                chg=false;
                if ~isempty(T_ORB) 
                    if isempty(FP(idxK).ORB)
                        chg=true;
                    else
                        if length(T_ORB.Features)>length(FP(idxK).ORB.Features)
                           chg=true; 
                        end
                    end
                end
                if chg %change all keep TS_ini
                    FP(idxK).ORB=T_ORB;
                    FP(idxK).TS_end=i+window-1;
                    FP(idxK).Kf=T;
                    FP(idxK).Th=T_Th;
                    FP(idxK).CC=T_CC;
                else
                    FP(idxK).TS_end=i+window-1; 
                end
            end
            
           
        end
    end
    ttm=toc;
    Times.TotalPerVideo=ttm;
    Times.TotalPerframe=ttm/TFrames;  
    Times.KfExtr=ttm-tm0-tm1-tm2-tm3-tm4;
    
    function [Kf, KfColor,ti]=GetKeyframe(Vd,VdColor,win)
        tic;
        for j=2:win
            Kf=w*Vd(:,:,j-1)+(1-w)*Vd(:,:,j);
            Vd(:,:,j)=Kf;
        end
        for j=2:win
            KfColor=w*VdColor(:,:,:,j-1)+(1-w)*VdColor(:,:,:,j);
            VdColor(:,:,:,j)=KfColor;
        end
        ti=toc;
    end %Function DownsizeVideo        
    

    function Kf=FoldFrame(kf_ref,m_r)
        Kf_T=imresize(kf_ref, [m_r m_r*2]);
        Kf=0.5*Kf_T(:,1:m_r)+0.5*(fliplr(Kf_T(:,m_r+1:m_r*2)));
    end %Function 

    function [FpCC]=GetFpCC(TC)
        %divide image in 16x16 blosck that represents the average of RGB values       
        FpCC=[];
        Frame=imresize (TC, [16 16], 'method','bicubic','Colormap', 'original'); %bicubic is the weighted average of 4x4 neighbors

        %% this section generates the Color Correlation Histogram, 
        %the method only needs the cardinality of the set of pixels in each case
        %Case #1: Rxy>Gxy>Bxy
        %Case #2: Rxy>Bxy>Gxy
        %Case #3: Gxy>Rxy>Bxy
        %Case #4: Gxy>Bxy>Rxy
        %Case #5: Bxy>Rxy>Gxy
        %Case #6: Bxy>Gxy>Rxy
        %The above color correlation divides the RGB color cube into six subspaces, thus being regarded as a special kind of color histogram 
        %(called color correlation histogram). It is well known that the color histogram [32] is one of the commonly used features for representing 
        %a color image or video. Com- pared with other shape-based or texture-based features, the color histogram is very fast to compute and is 
        %flexible in terms of storage. Up until now, many features derived from the color histogram have been successfully applied to video retrieval,
        %segmentation, and identification [33]–[35]. However, few methods use a single feature from the color histogram for video sequence matching. 
        %Here, we must note that since the color correlation histogram is a global feature without any spatial information about the image, 
        %the robustness against most common content-preserving operations is expected to be much better. Please refer to our analysis in Section II-C.
        %The discriminability of the color correlation histogram may be weak for a single image or video frame. 
        %However, it is shown that our feature is very promising for video sequence matching when the length of the video is increased 
        
        C(6) = zeros;
        %revisar la correlacion y quitar los casos donde RGB sean iguales
        for ic=1:16
            for jc=1:16
                if (Frame(ic,jc,1) >= Frame(ic,jc,2)) && (Frame(ic,jc,2) >= Frame(ic,jc,3)) && (Frame(ic,jc,1) ~= Frame(ic,jc,3))
                    C(1)=C(1)+1;
                elseif (Frame(ic,jc,1) >= Frame(ic,jc,3)) && (Frame(ic,jc,3) >= Frame(ic,jc,2)) && (Frame(ic,jc,1) ~= Frame(ic,jc,2))
                    C(2)=C(2)+1;
                elseif (Frame(ic,jc,2) >= Frame(ic,jc,1)) && (Frame(ic,jc,1) >= Frame(ic,jc,3)) && (Frame(ic,jc,2) ~= Frame(ic,jc,3))
                    C(3)=C(3)+1;
                elseif (Frame(ic,jc,2) >= Frame(ic,jc,3)) && (Frame(ic,jc,3) >= Frame(ic,jc,1)) && (Frame(ic,jc,2) ~= Frame(ic,jc,1))
                    C(4)=C(4)+1;
                elseif (Frame(ic,jc,3) >= Frame(ic,jc,1)) && (Frame(ic,jc,1) >= Frame(ic,jc,2)) && (Frame(ic,jc,3) ~= Frame(ic,jc,2))
                    C(5)=C(5)+1;
                elseif (Frame(ic,jc,3) >= Frame(ic,jc,2)) && (Frame(ic,jc,2) >= Frame(ic,jc,1)) && (Frame(ic,jc,3) >= Frame(ic,jc,1))
                    C(6)=C(6)+1;
                end
            end
        end   
       %%Obtains the percentage of the color correlation
        total=sum(C);
        for ic=1:5
            FHaux=dec2bin(100*C(ic)/total,7); %pecentage in binary form
            FpCC=[FpCC FHaux]; %#ok<*AGROW>
        end
        FpCC=logical(FpCC-'0');
    end %Function GetFpCC
    %%%%%%%%%%%%%%%%%%%%%%%
   
%     function s=DCTbyBlocks(Kf_r,m,w)
%         % DCT by blocks
%          bk=0;
%          %B=zeros(w*2,w*2); %is the number of traslaped blocks (16-1)*(16-1)=225
%          s(1:450)=false; %2 bits per block
%          for bi=1:w:m-w 
%             for bj=1:w:m-w
%                 bk=bk+1;
%                 B= Kf_r(bi:bi+2*w-1,bj:bj+2*w-1); %overlaped block
%                 DCT=dct2(B);    %2D-DCT
%                 AC1=DCT(1,2);    %1st AC coefficient
%                 AC2=DCT(2,1);      %2nd AC coefficient
%                 if AC1>=0         %hash
%                     s(2*bk-1)=true; %otherwise false, already initialized
%                 end
%                 if AC2>=0
%                     s(2*bk)=true; %otherwise false, already initialized
%                 end
%             end
%          end
%     end %Function DCTbyBlocks

    function ORB=getORB(Kf_r)
        
        [~,ORBdescriptor]=cv.ORB(Kf_r, 'NFeatures', 500);
        if ~isempty(ORBdescriptor)
            ORB= binaryFeatures(ORBdescriptor);
        else
            ORB=[];
        end
    end %Function getORB

end %Function KeyframeExtr


end


function [ V,A ] = ReadVideoInDataSet( v,DB,nfps,asr,seglength,M,N,sf, kM,ThSize, ThType)
%READVIDEOINDATASET read, process and segment a video, returns the keyframes per video 
%  1. Read a video number 'v' in dataset number DB


if nargin<3
    nfps=10; %downsampling frame rate
    asr=8240; %audio downsampling
    seglength=1; %segment length in seconds
    kM=1 ; %method for keyframe extraction, std=1 average
    M=320; %keyframe width
    N=320; %keyframe height
    sf=1.3;
    ThSize=[60,60]; %size of thumbnail 
    ThType=1;%std = 1 grayscale

end
    
switch DB
    case 1 %Hollywood 2
        Filepath='D:\Testing time\Video dataset Hollywood2\AVIClipsScenes\';
        NumVideo=['0000' num2str(v)];
        NumVideo=NumVideo(length(NumVideo)-4:length(NumVideo));
        VideoName= [Filepath  'sceneclipautoautotrain' NumVideo '.avi'];
    case 2 %OpenArchive
        
    case 3 %ReTRiEVED
        Filepath='D:\Testing time\Retrieved dataset\';
        VideoName= [Filepath  num2str(v) '.mp4'];
    case 4 %TRECVID1
        
    case 5 % TRECVID2
        
    otherwise %VCDB
end
        

    % 1. Read video 
    videoFReader = vision.VideoFileReader;
    videoFReader.Filename = VideoName;
    videoFReader.ImageColorSpace = 'RGB';
    videoFReader.VideoOutputDataType = 'uint8';
    videoFReader.AudioOutputPort = true;
   
    VideoInfo=info(videoFReader); 
    FRVal= VideoInfo.VideoFrameRate;
   
    bufferLength= FRVal*seglength; %2 seconds
    FNum=0;
    bufferNum=0;
   
    
%     
%     v = VideoReader(VideoName);
%     FNum=0;
%     fps=floor(v.FrameRate);
%     dwin=floor(fps/nfps); %to downsample to nfps 
%     Duration=v.Duration;
%     VideoFrames=zeros(m,n,floor(Duration)*nfps, 'uint8');
%     VideoFramesColor=zeros(m,n,3,floor(Duration)*nfps, 'uint8');
%     i=1;
%     
%     while hasFrame(v)
%         i=i+1;
%         f=readFrame(v);
%         if mod(i,dwin)==0
%             FNum=FNum+1;
%             
%             VideoFramesColor(:,:,:,FNum) = imresize(f, [m n]);
%             VideoFrames(:,:,FNum) = rgb2gray(VideoFramesColor(:,:,:,FNum));
%             
%             
%         end

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
    function [Imout]=ScaleResize(Im,M,N,sf)
           
        [i, j, c]=size(Im);
        Im1=imresize(Im, sf); %scale factor
        [i1, j1, c1]=size(Im1);
        i2=floor((i1-i)/2);
        j2=floor((j1-j)/2);
        Imout=imresize(Im1(i2:i1-i2,j2:j1-j2,:), [M N]);
    end
end


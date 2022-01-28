%% scripting

[vidfile,viddir]=uigetfile();
mp4file=fullfile(viddir,vidfile);

vid=VideoReader(mp4file);

randofind=randi([100 vid.NumFrames-100]);

testframe=read(vid,randofind,'native');

AInv = imcomplement(im2double(testframe));
BInv = imreducehaze(AInv,.8); %, 'ContrastEnhancement','none');
B = imcomplement(BInv);
figure; montage({testframe,B});

figure;
B = imlocalbrighten(im2double(testframe),'AlphaBlend',true);
montage({testframe,B})



figure;
subplot(2,2,1);
image(B);
B2=B;
B2(:,:,1)=imadjust(B2(:,:,1));
subplot(2,2,2); image(B2);
B2(:,:,2)=imadjust(B2(:,:,2));
subplot(2,2,3); image(B2);
B2(:,:,3)=imadjust(B2(:,:,3));
subplot(2,2,4); image(B2);

% smooth this out a bit?

B2=imgaussfilt(B,2);
figure; montage({B,B2});

PSF = fspecial('motion',2,2);
Idouble = im2double(B);
imshow(Idouble)
title('Blurred Image')
wnr1 = deconvwnr(Idouble,PSF,.1);
imshow(wnr1)
title('Restored Blurred Image')
%{

% this doesnt realyl work, its fine for greyscale though
shadow_lab = rgb2lab(testframe);
max_luminosity = 100;
L = shadow_lab(:,:,1)/max_luminosity;
shadow_imadjust = shadow_lab;
shadow_imadjust(:,:,1) = imadjust(L)*max_luminosity;
shadow_imadjust = lab2rgb(shadow_imadjust);

shadow_histeq = shadow_lab;
shadow_histeq(:,:,1) = histeq(L)*max_luminosity;
shadow_histeq = lab2rgb(shadow_histeq);

shadow_adapthisteq = shadow_lab;
shadow_adapthisteq(:,:,1) = adapthisteq(L)*max_luminosity;
shadow_adapthisteq = lab2rgb(shadow_adapthisteq);

%}


figure
montage({testframe,shadow_imadjust,shadow_histeq,shadow_adapthisteq, B},'Size',[1 5])
title("Original Image and Enhanced Images using imadjust, histeq, adapthisteq, and reducehaze")
%%
% TAKES FOREVER, NOT USEFUL


%[vidfile,viddir]=uigetfile();

ratdir=uigetdir();

vidfiles=getAllFiles(ratdir,'.mp4');


for i=1:length(vidfiles)
%mp4file=fullfile(viddir,vidfile);

mp4file=vidfiles{i};

vid=VideoReader(mp4file);

[a,b,c]=fileparts(fullfile(vid.Path,vid.Name));

outavifile=fullfile(a,[b '-bright']);

ovid = VideoWriter(outavifile);
%ovid = VideoWriter(outavifile,'MPEG-4');
open(ovid)
%c1 = onCleanup(@() close(ovid));


% build waitbar
h = waitbar(0,sprintf('Starting %s Please wait... (%02.0f%%)',b,0));
%c2 = onCleanup(@(x,y) close(h));

updateintv = max(floor(vid.NumFrames/100),1);


% go frame by frame
for ii = 1:vid.NumFrames
    % Read the frame
    frame = read(vid,ii);
    
    % Apply image filter
    frame2=imlocalbrighten(im2double(frame),.8,'AlphaBlend',true);
    writeVideo(ovid,frame2);
    
    if(mod(ii,updateintv)==0)
        if(ishandle(h))
            waitbar(ii/vid.NumFrames,h,sprintf('Processing %s Please wait... (%02.0f%%)',b,ii*100/vid.NumFrames));
        end
    end
end

close(ovid);
close(h);
end
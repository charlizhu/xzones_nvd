[FileName,PathName,FilterIndex] = uigetfile('.mat');
disp(['Opening....' PathName FileName])
load([PathName '\' FileName]);
disp(['Loaded: ' PathName FileName])

v = VideoWriter([PathName 'tmVid.mp4'], 'MPEG-4');
v.FrameRate = 1/(frame_period_ms*1e-3);
open(v);

i=1;

numFrames = size(F,2);
while (i <= numFrames && ~isempty(F(i).cdata))
    writeVideo(v,F(i));
    i=i+1;
    if(mod(i/numFrames*100,10)==0)
        fprintf('%d percent complete...\n',i/numFrames*100)
    end
end
close(v)
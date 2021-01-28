vidObj = VideoWriter('tm1500.mp4','MPEG-4');
vidObj.FrameRate = 40;
open(vidObj);
for n=1:10
    disp(['File #', num2str(n)]);
    matFilename = ['tm1500_', num2str(n,'%03d'),'.mat'];
    load(matFilename, 'F1');
    writeVideo(vidObj, F1); disp('part1'); F1 = [];
    load(matFilename, 'F2');
    writeVideo(vidObj, F2); disp('part2'); F2 = [];
    load(matFilename, 'F3');
    writeVideo(vidObj, F3); disp('part3'); F3 = [];
    load(matFilename, 'F4');
    writeVideo(vidObj, F4); disp('part4'); F4 = [];
    load(matFilename, 'F5');
    writeVideo(vidObj, F5); disp('part5'); F5 = [];
    load(matFilename, 'F6');
    writeVideo(vidObj, F6); disp('part6'); F6 = [];
    load(matFilename, 'F7');
    writeVideo(vidObj, F7); disp('part7'); F7 = [];
    load(matFilename, 'F8');
    writeVideo(vidObj, F8); disp('part8'); F8 = [];
    load(matFilename, 'F9');
    writeVideo(vidObj, F9); disp('part9'); F9 = [];
    load(matFilename, 'F10');
    writeVideo(vidObj, F10); disp('part10'); F10 = [];
end
close(vidObj);
disp('done');

vehStruct = struct('id', 0, 'lane', 0, 'type', '[]', 'width', 0, 'length', 0, 'rectangle', 0, 'entryTime', 0, 'tick', 0, 'posxy', '[]', 'velxy', '[]', 'accxy', '[]');
vHistAll = repmat(vehStruct, 4, 1000);

for m=1:10
    filePath = 'C:\Radar\gtrack\test\vectors\urv20\';
    fileName = ['tm1500_', num2str(m,'%03d'),'.mat'];
    load(fileName, 'cHist', 'vHist');
    
    for n=1:4
        for k = 1:vHist(n).num
            vh = vHist(n).veh{k};
            p = vh.id;
            vHistAll(n,p).id = vh.id;
            vHistAll(n,p).lane = vh.lane;
            vHistAll(n,p).type = vh.type;
            vHistAll(n,p).width = vh.width;
            vHistAll(n,p).length = vh.length;
            vHistAll(n,p).rectangle = vh.rectangle;
            vHistAll(n,p).entryTime = vh.entryTime;
            vHistAll(n,p).tick = vh.tick;
            vHistAll(n,p).posxy = vh.posxyHistory;
            vHistAll(n,p).velxy = vh.velxyHistory;
            vHistAll(n,p).accxy = vh.accxyHistory;
        end
    end    
    fout = fopen([filePath,'pointCloud1500_',num2str(m,'%03d'),'.dat'], 'w');
    for t=1:1500
        fprintf(fout, 'NUM=%d\n', cHist(t).num);
        if(cHist(t).num > 0)
            for n=1:cHist(t).num
                fprintf(fout, '%f,%f,%f,%f,%f,%f,%u\n', ...
                    cHist(t).u(1,n),cHist(t).u(2,n),cHist(t).u(3,n), ...
                    cHist(t).uvar(1,n),cHist(t).uvar(2,n),cHist(t).uvar(3,n),...
                    cHist(t).vid(n));
            end
        end
    end
    fclose(fout);
end
save([filePath,'vHistAll.mat'], 'vHistAll');



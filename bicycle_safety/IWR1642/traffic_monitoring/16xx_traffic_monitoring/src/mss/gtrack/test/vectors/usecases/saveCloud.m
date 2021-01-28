for m=1:10
    filePath = 'C:\Radar\gtrack\test\vectors\urv100\';
    fileName = [filePath,'tm1500_', num2str(m,'%03d'),'.mat'];
    load(fileName, 'cHist', 'vHist');
    vehStruct = struct('num', 0, 'veh', []);
    vhist = repmat(vehStruct,1, 4);
    
    for n=1:4
        for k = 1:vHist(n).num
            vhist(n).num = vHist(n).num;
            vhist(n).veh{k} = [];
            
            vh = vHist(n).veh{k};
            vhist(n).veh{k}.id = vh.id;
            vhist(n).veh{k}.lane = vh.lane;
            vhist(n).veh{k}.type = vh.type;
            vhist(n).veh{k}.width = vh.width;
            vhist(n).veh{k}.length = vh.length;
            vhist(n).veh{k}.rectangle = vh.rectangle;
            vhist(n).veh{k}.entryTime = vh.entryTime;
            vhist(n).veh{k}.tick = vh.tick;
            vhist(n).veh{k}.posxy = vh.posxyHistory;
            vhist(n).veh{k}.velxy = vh.posxyHistory;
            vhist(n).veh{k}.accxy = vh.posxyHistory;
        end
    end
    save([filePath,'vhist_',num2str(m,'%03d'),'.mat'], 'vhist');
    
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



function makePlots(car)
    %disp(car);
%     for i = 1:length(car.TID)
%         data = num2str([car.TID(i),car.x(i),car.y(i),car.vx(i),car.vy(i),length(car.TID)]);
%         disp(data)
%         disp(ischar(data))
%         savedArray = strcat(data,'\n');
%         outputFile = fopen('hookdata2.txt', 'a+');
%         fprintf(outputFile,savedArray);
%         fclose('all')
%     end

% outputFile = fopen('connect2.txt', 'a+');
% disp(car.TID);
% data = num2str(car.TID);
% savedArray = strcat(data,'\n');
% fprintf(outputFile,savedArray);
% fclose('all');

outputFile = fopen('rear_righthook.txt','a+');
%disp(car);
data = num2str(car);
savedArray = strcat(data,'\n');
fprintf(outputFile,savedArray);
fclose(outputFile);
end
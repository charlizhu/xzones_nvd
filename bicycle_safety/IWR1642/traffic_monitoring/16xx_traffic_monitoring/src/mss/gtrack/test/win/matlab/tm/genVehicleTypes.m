function vtype = genVehicleTypes(aNum, prob)
    tLanes = size(prob,1);
    vtype = cell(tLanes,1);
    for n = 1:tLanes
        R = rand(aNum(n),1);
        c = cumsum(prob(n,:));
        vtype{n} = arrayfun(@(x) find(x <= c, 1, 'first'), R);
    end
end


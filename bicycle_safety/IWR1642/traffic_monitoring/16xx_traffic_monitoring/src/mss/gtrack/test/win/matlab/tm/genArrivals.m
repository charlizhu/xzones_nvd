function ind = genArrivals(time, dt, lambda)
    N = time/dt; % number of simulation steps
    ind = zeros(N,length(lambda)); 
    event=zeros(N,length(lambda)); % array recording at each step if a "packet" arrived. 
    % initialize it to zeros
    R=rand(size(event)); % generate a random array (with elements in [0,1]) of the same size as "event"
    for n=1:length(lambda)
        ind(:,n) = R(:,n)<(lambda(n)*dt);
    end
end
Fs = 1/0.03

fileID = fopen('AWR_Frontal_hook_data.txt','r');
formatSpec = '%f %f %f %f %f %f';
sizeA = [6 Inf];
A = fscanf(fileID,formatSpec, sizeA);
x = A(1,:);
y = A(2,:);
t = [0:0.03:183*0.03-0.03];
L = length(y)
x_fft = fft(x);
y_fft = fft(y);
figure()
subplot(2, 2, 1);
plot(t, x)
subplot(2, 2, 2);
plot(t, y)
f = Fs*(0:(L -1))/L
subplot(2, 2, 3);
plot(f, x_fft)
subplot(2, 2, 4);
plot(f, y_fft)

figure();
title("lowpass")
cleanedx = lowpass(x, 9, Fs);
subplot(2, 2, 1);
plot(t, cleanedx)
title("cleaned x");
subplot(2, 2, 2);
plot(t, x)

cleanedy = lowpass(y, 9, Fs);
subplot(2, 2, 3);

plot(t, cleanedy)
title("cleaned y");
subplot(2, 2, 4);
plot(t, y)


%butterworth filter 
figure();
title("butterworth")
[b,a] = butter(5,[0.01 15]/(Fs/2))
butterx = filter(b,a,x);
buttery = filter(b,a,y);
cov_butterx = cov(butterx)
cov_buttery = cov(buttery)

subplot(2, 2, 1);
plot(t, butterx)
title("buttered x");
subplot(2, 2, 2);
plot(t, x)

subplot(2, 2, 3);
plot(t, buttery)
title("buttered y");
subplot(2, 2, 4);
plot(t, y)

%x-y plots
figure()
subplot(3, 1, 1);
plot (x,y)
title("original-xy")
subplot(3, 1, 2);
plot (cleanedx,cleanedy)
title("default_lowpass-xy")
subplot(3, 1, 3);
plot (butterx,buttery)
title("buttered-xy")

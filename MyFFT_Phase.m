function [Phase] = MyFFT_Phase(K)
%% MyFFT 
%This m is used to separate the torque
%Output the corresponding frequency, magnitude and phase
%Input Variable is K
[N M] = size(K);
%N is the length of the signal
%M is the sampling times
%One electrical period is the 1hz
Fs = N;
T = 1/Fs;%Sampling period
t = ((0:N-1)*T)';%Time vector
H = 1./t;
Y = fft(K);
%Magnitude
P2 = abs(Y/N);
P1 = P2(1:N/2+1,:);
P1(2:N/2,:) = 2*P2(2:N/2,:);
Magnitude = P1;
%Phase
Ph = angle(Y);
Phase = Ph(1:N/2+1,:)*180/pi; 
%P1 is the single side spectrum while P2 is the double side spectrum
f = Fs*(0:(N/2))/N;
W = 2*pi*f;
% figure
% subplot(2,1,1)
% bar(f,Magnitude)
% subplot(2,1,2)
% bar(f,Phase)
% Getting the deviation function
for j =1:1:M
for i = 1:1:24
SignalK(i,:) = -W(i)*Magnitude(i,j)*sin(W(i)*t+Phase(i,j)*pi/180)/(2*pi);
end
Dsignal(:,j) = sum(SignalK)';
end
end
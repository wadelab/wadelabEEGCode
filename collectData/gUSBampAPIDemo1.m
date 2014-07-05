
ai = analoginput('guadaq',1);
addchannel(ai,1:16);
set(ai,'SampleRate',256,'SamplesPerTrigger',512);
start(ai)
while strcmp(ai.running,'On')==1
end
data=getdata(ai,ai.SamplesAvailable);
plot(data(256:end,1));
xlabel('Samples');
ylabel('Signal [Volt]');
delete(ai);
clear ai
figure(2);
surf(data);
shading interp;

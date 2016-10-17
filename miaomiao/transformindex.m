function output = transformindex(input)
% output = transformindex(input)
% fixes the binary inputs for the ANT EEG amplifier because the pins are in a different order from the ViewPixx
% desired numbers must be <256
% DHB 18/8/14
% ARW broke it out and saved it separately 06/06/2016

truebits = 2.^(2:2:24);
dn = dec2bin(input,length(truebits));
output = 0;
for m = 1:length(truebits)
    output = output + truebits(m)*str2num(dn(end-m+1));
end


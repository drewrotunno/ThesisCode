clear all; close all;

load('all_hf_data.mat');

figure('position', [47         293        1800         450], 'renderer', 'painter');
hold on

plot(efreq/1e3, es11, 'color', [.5,.5,.5]);
plot(zfreq/1e3, zs11, 'b');
plot(ufreq/1e3, us11, 'r');

lgd=legend({'Evap','Z','U'});

ylim([-30,0]);

xlabel('Frequency (GHz)')
ylabel('Power Reflected (S11, dB)')

box on
grid on

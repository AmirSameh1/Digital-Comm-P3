clc; clear; close all;

%% 1. Initialization
num_bits = 1.2e6; % Multiple of 3 and 4 for 8PSK and 16-QAM
data_bits = randi([0 1], 1, num_bits);
EbNo_dB  = 0:1:15; 
EbNo_lin = 10.^(EbNo_dB/10);

% Pre-allocate results
BER_8PSK  = zeros(1, length(EbNo_dB));
BER_16QAM = zeros(1, length(EbNo_dB));

%% 2. 8PSK Section (Based on Figure 2)
k_8PSK = 3;
% Mapping based on the 8PSK constellation image:
% 000 -> 1, 001 -> exp(j*pi/4), 011 -> j, 010 -> exp(j*3pi/4), 
% 110 -> -1, 111 -> exp(j*5pi/4), 101 -> -j, 100 -> exp(j*7pi/4)
Map_8PSK = [exp(1j*0), exp(1j*pi/4), exp(1j*3*pi/4), exp(1j*pi/2), ...
            exp(1j*7*pi/4), exp(1j*3*pi/2), exp(1j*pi), exp(1j*5*pi/4)];

bit_groups8 = reshape(data_bits, k_8PSK, [])';
sym_id8 = bit_groups8(:,1)*4 + bit_groups8(:,2)*2 + bit_groups8(:,3);
Si_8 = Map_8PSK(sym_id8 + 1);

for i = 1:length(EbNo_dB)
    % AWGN Channel [cite: 21]
    sigma2 = 1/(2 * k_8PSK * EbNo_lin(i));
    noise = sqrt(sigma2) * (randn(1,length(Si_8)) + 1j*randn(1,length(Si_8)));
    Xi = Si_8 + noise;
    
    % Demapper: Nearest Neighbor [cite: 24]
    [~, dec_idx] = min(abs(Xi.' - Map_8PSK), [], 2);
    dec_idx = dec_idx.' - 1;
    
    % Bit Recovery
    rec_bits = zeros(1, num_bits);
    rec_bits(1:3:end) = floor(dec_idx/4);
    rec_bits(2:3:end) = floor(mod(dec_idx, 4)/2);
    rec_bits(3:3:end) = mod(dec_idx, 2);
    BER_8PSK(i) = sum(data_bits ~= rec_bits) / num_bits;
end

%% 3. 16-QAM Section (Based on Figure 2)
k_16QAM = 4;
% Coordinates from the 16-QAM image (b0 b1 b2 b3)
% Maps bits to grid coordinates +/-1, +/-3
Map_16QAM = [(-3-3j), (-3-1j), (-3+3j), (-3+1j), ... % 0000 to 0011
             (-1-3j), (-1-1j), (-1+3j), (-1+1j), ... % 0100 to 0111
             ( 3-3j), ( 3-1j), ( 3+3j), ( 3+1j), ... % 1000 to 1011
             ( 1-3j), ( 1-1j), ( 1+3j), ( 1+1j)];    % 1100 to 1111

bit_groups16 = reshape(data_bits, k_16QAM, [])';
sym_id16 = bit_groups16(:,1)*8 + bit_groups16(:,2)*4 + bit_groups16(:,3)*2 + bit_groups16(:,4);
Si_16 = Map_16QAM(sym_id16 + 1);

for i = 1:length(EbNo_dB)
    sigma2 = 1/(2 * k_16QAM * EbNo_lin(i));
    noise = sqrt(sigma2) * (randn(1,length(Si_16)) + 1j*randn(1,length(Si_16)));
    Xi = Si_16 + noise;
    
    [~, dec_idx] = min(abs(Xi.' - Map_16QAM), [], 2);
    dec_idx = dec_idx.' - 1;
    
    rec_bits = zeros(1, num_bits);
    rec_bits(1:4:end) = floor(dec_idx/8);
    rec_bits(2:4:end) = floor(mod(dec_idx, 8)/4);
    rec_bits(3:4:end) = floor(mod(dec_idx, 4)/2);
    rec_bits(4:4:end) = mod(dec_idx, 2);
    BER_16QAM(i) = sum(data_bits ~= rec_bits) / num_bits;
end

%% 4. Theoretical Calculations 
% Theoretical 8PSK BER
th_8PSK = (2/k_8PSK)*0.5*erfc(sqrt(k_8PSK*EbNo_lin)*sin(pi/8));
% Theoretical 16-QAM BER (assuming normalized power for comparison)
th_16QAM = (3/8)*erfc(sqrt(0.4*EbNo_lin)); 

%% 5. Plotting Results 
figure('Name','BER vs Eb/No - 8PSK and 16-QAM','Color','w','Position',[100 100 800 500]);
semilogy(EbNo_dB, BER_8PSK, 'g-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(EbNo_dB, th_8PSK, 'k--', 'LineWidth', 1.2);
semilogy(EbNo_dB, BER_16QAM, 'm-s', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(EbNo_dB, th_16QAM, 'r:', 'LineWidth', 1.2);

grid on; ylim([1e-5 1]); xlim([0 15]);
xlabel('E_b/N_0 (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('BER vs E_b/N_0: 8PSK and 16-QAM Simulation', 'FontSize', 13);
legend('8PSK Simulated', '8PSK Theoretical', '16-QAM Simulated', '16-QAM Theoretical', ...
       'Location', 'southwest');
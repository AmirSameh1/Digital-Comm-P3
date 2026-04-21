clc; clear; close all;

num_bits    = 1e6;
data_bits = randi([0 1], 1, num_bits);
EbNo_dB     = -2:1:18;      % Eb/No range in dB that will be ploted
EbNo_lin    = 10.^(EbNo_dB/10); % Eb/No range in linear unit so we could use it to generate noise
BER_BPSK  = zeros(1, length(EbNo_dB));

%%  ---------------- BPSK -----------------

Si = 2*data_bits-1;  % mapping the data into 1 ,-1 
k_BPSK =1;  %bits ber sympol in BPSK so we could find prober E/No for noise generation

for i=1:length(EbNo_dB)

    sigma2 = 1/(2 * k_BPSK * EbNo_lin(i));  % noise variance per dimension
    noise = sqrt(sigma2) * randn(1, num_bits);
    Xi = Si + noise;
    received_data = zeros(1, num_bits);

     for j=1:length(Xi) % demapper to recover bits from recieved sympol with threshold 0
        if(Xi(j)>0)
            received_data(j) =1;
        end
     end

    BER_BPSK(i) = sum(data_bits ~= received_data) / num_bits;
end



%% ----Plot results vs theoritical results

th_BER_BPSK  = 0.5*erfc(sqrt(EbNo_lin)); % theoritical for BPSK

figure('Name','BER vs Eb/No - All Modulations','Color','w','Position',[100 100 900 600]);

semilogy(EbNo_dB, BER_BPSK,  'b-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(EbNo_dB, th_BER_BPSK,  'r--',  'LineWidth', 1.5);

grid on; ylim([1e-5 1]);
xlabel('E_b/N_0 (dB)', 'FontSize', 13);
ylabel('BER',          'FontSize', 13);
title('BER vs E_b/N_0 — Simulated (solid+marker) vs Theoretical (dashed)', 'FontSize', 13);

legend('BPSK sim', ...
       'BPSK theory', ...
       'Location','southwest', 'FontSize', 10);



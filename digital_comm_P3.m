clc; clear; close all;

num_bits    = 1e6;
data_bits = randi([0 1], 1, num_bits);
EbNo_dB     = -2:1:18;      % Eb/No range in dB that will be ploted
EbNo_lin    = 10.^(EbNo_dB/10); % Eb/No range in linear unit so we could use it to generate noise
BER_BPSK  = zeros(1, length(EbNo_dB));
BER_QPSK  = zeros(2, length(EbNo_dB));

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

    BER_BPSK(i) = sum(data_bits ~= received_data) / num_bits; % compare each bit in each array
end

th_BER_BPSK  = 0.5*erfc(sqrt(EbNo_lin)); % theoritical for BPSK

%%  ---------------- QPSK -----------------

k_QPSK = 2;%bits ber sympol in QPSK
QPSK_maps = (1/sqrt(2)) *  [(-1-1j), (-1+1j), (1-1j), (1+1j);     % 00  01  10  11
             (-1-1j), (-1+1j), (1+1j), (1-1j)];    % 00  01  10  11

for x=1:2
    for i=1:length(EbNo_dB)

        bit_pairs = reshape(data_bits, 2, []);
        bit_pairs = bit_pairs';
        sym_id = bit_pairs(:,1)*2 + bit_pairs(:,2);
        Si = QPSK_maps(x, sym_id + 1);

        % AWGN
        sigma2 = 1/(2 * k_QPSK * EbNo_lin(i));
        noise = sqrt(sigma2) * (randn(1,length(Si)) + 1j*randn(1,length(Si)));
        Xi = Si + noise;

        % Demapper: nearest neighbor
        [~, dec_idx] = min(abs(Xi.' - QPSK_maps(x,:)), [], 2);
         dec_idx = dec_idx.' - 1;   % 0..3
    
        % Convert back to bits
         received_data = zeros(1, num_bits);
         received_data(1:2:end) = floor(dec_idx/2);
         received_data(2:2:end) = mod(dec_idx, 2);
         BER_QPSK(x,i) = sum(data_bits ~= received_data) / num_bits; % compare each bit in each array
    end
end
th_BER_QPSK1 = 0.5*erfc(sqrt(EbNo_lin)); % theoritical for BPSK (Gray)
th_BER_QPSK2 = erfc(sqrt(EbNo_lin)); % theoritical for QPSK (non- Gray)


%% ----Plot results vs theoritical results

%---BPSK
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


%---QPSK_map1
figure('Name','BER vs Eb/No - All Modulations','Color','w','Position',[100 100 900 600]);
semilogy(EbNo_dB, BER_QPSK(1,:),  'b-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(EbNo_dB, th_BER_QPSK1,  'r--',  'LineWidth', 1.5);
grid on; ylim([1e-5 1]);
xlabel('E_b/N_0 (dB)', 'FontSize', 13);
ylabel('BER',          'FontSize', 13);
title('BER vs E_b/N_0', 'FontSize', 13);
legend('QPSK (Gray) sim', ...
       'QPSK (Gray) theory', ...
       'Location','southwest', 'FontSize', 10);


       
%---QPSK_map2
figure('Name','BER vs Eb/No - All Modulations','Color','w','Position',[100 100 900 600]);
semilogy(EbNo_dB, BER_QPSK(2,:),  'b-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(EbNo_dB, th_BER_QPSK2,  'r--',  'LineWidth', 1.5);
grid on; ylim([1e-5 1]);
xlabel('E_b/N_0 (dB)', 'FontSize', 13);
ylabel('BER',          'FontSize', 13);
title('BER vs E_b/N_0', 'FontSize', 13);
legend('QPSK (non-Gray) sim', ...
       'QPSK (non-Gray) theory', ...
       'Location','southwest', 'FontSize', 10);

%---QPSK_map1&2
figure('Name','BER vs Eb/No - All Modulations','Color','w','Position',[100 100 900 600]);
semilogy(EbNo_dB, BER_QPSK(1,:),  'b-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(EbNo_dB, BER_QPSK(2,:),  'r-o',  'LineWidth', 1.5);
grid on; ylim([1e-5 1]);
xlabel('E_b/N_0 (dB)', 'FontSize', 13);
ylabel('BER',          'FontSize', 13);
title('BER vs E_b/N_0', 'FontSize', 13);
legend('QPSK (Gray) sim', ...
       'QPSK (non-Gray) sim', ...
       'Location','southwest', 'FontSize', 10);
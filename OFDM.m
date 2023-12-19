%This is the code that automatically adapts to a modulation technique based on the user input. The file OFDM2.m can be used to test for custom SNR and modulation techniques.
close all
clc

subcarriers = 64;
cyclic_prefix = 16;

%Change the value of snr_db from range 0-20. This is the typical possible range of SNR for regular modems.
snr_db = input('Enter SNR value within range 0-20:-  ');

fprintf('Enter 1 for AWGN Channel. Enter 2 for AWGN + Fading Channel. \n')
flag = input('Enter your choice:-  ');

if flag == 2
    fprintf('Enter 1 for Multipath Fading. Enter 2 for Rayleigh fading. Enter 3 for Rician fading. Enter 4 for Nakagami-m fading. \n')
    fading_type = input('Enter type of fading to be simulated:-  ');
    if fading_type == 1
        channel_taps = input('Enter the number of channel taps(number of reflections):-  ');
    elseif fading_type == 3
        K = input('Enter the value of K factor for a specific scenario:-  ');
    elseif fading_type == 4
        m = input('Enter Nakagami-m variable:- ');
    end
end

%Input any image file below. It should be a .jpg or .png type.
image = imread('SnapdragonStadium.jpg');

binary_form = dec2bin(image(:))';
binary_form = binary_form(:);

%BPSK = 1; QPSK = 2; 8PSK = 3; 16QAM = 4; 32QAM = 5; 64QAM = 6
if flag == 2
    if snr_db <= 0 || snr_db < 8
        bits = 1;
    elseif snr_db <= 8 || snr_db < 12
        bits = 2;
    elseif snr_db <= 12 || snr_db < 14
        bits = 3;
    elseif snr_db <= 14 || snr_db < 16
        bits = 4;
    elseif snr_db <= 16 || snr_db < 18
        bits = 5;
    elseif snr_db <= 18 || snr_db <= 20
        bits = 6;
    else
        fprintf('Impractical Value of SNR')
    end
elseif flag == 1
    if snr_db <= 0 || snr_db < 5
        bits = 1;
    elseif snr_db <= 5 || snr_db < 8
        bits = 2;
    elseif snr_db <= 8 || snr_db < 12
        bits = 3;
    elseif snr_db <= 12 || snr_db < 14
        bits = 4;
    elseif snr_db <= 14 || snr_db < 16
        bits = 5;
    elseif snr_db <= 16 || snr_db <= 20
        bits = 6;
    else
        fprintf('Impractical Value of SNR')
    end
else
    fprintf('Wrong Input')
end

missing_bits = mod((bits - mod(length(binary_form) , bits)), bits);
null_padding = repmat('0', missing_bits, 1);
padded_bits = [binary_form;null_padding];
matrix = reshape(padded_bits, bits, (length(padded_bits)/bits))';
matrix_decimal = bin2dec(matrix);

if bits == 1
    index = 2^(bits - 1);
    theta = 0 : (pi/index) : ((2*pi) - (pi/index));
    I = cos(theta);
    Q = sin(theta);
    OFDM_symbols = (I + 1i*Q);
elseif bits == 2 || bits == 3
    index = 2^(bits - 1);
    theta = 0 : (pi/index) : ((2*pi) - (pi/index));
    I = cos(theta + pi/4);
    Q = sin(theta + pi/4);
    OFDM_symbols = (I + 1i*Q);
elseif bits == 4 || bits == 6
    index = sqrt(2^bits);
    I = repmat(linspace(-1,1,index), index, 1);
    Q = repmat(linspace(-1,1,index)', 1, index);
    OFDM_symbols = I(:) + Q(:)*1i;
elseif bits == 5
    index = sqrt(2^bits + 4);
    I = repmat(linspace(-1,1,index),index,1);
    Q = repmat(linspace(-1,1,index)',1,index);
    OFDM_symbols = (I(:) + Q(:)*1i);
    OFDM_symbols = OFDM_symbols([2:5 7:30 32:35]);
end

modulated_symbols = OFDM_symbols(matrix_decimal + 1);

missing_subcarriers = mod((subcarriers - mod(length(modulated_symbols), subcarriers)), subcarriers);
padded_subcarriers = [modulated_symbols;zeros(missing_subcarriers,1)];
subcarrier_matrix = reshape(padded_subcarriers, subcarriers, (length(padded_subcarriers)/subcarriers));
ifft_subcarriers = ifft(subcarrier_matrix);

padded_cyclic_prefix = [ifft_subcarriers((end - cyclic_prefix +1 : end), :);ifft_subcarriers];
transmitted_signal = padded_cyclic_prefix(:);

signal_power = mean(abs(transmitted_signal.^2));
snr = 10^(snr_db/10);
noise_power = signal_power/snr;

AWGN_real = normrnd(0,sqrt(noise_power/2),size(real(transmitted_signal)));
AWGN_imaginary = normrnd(0,sqrt(noise_power/2),size(imag(transmitted_signal)))*1i;
noisy_signal = transmitted_signal + AWGN_real + AWGN_imaginary;

if flag == 2
    if fading_type == 1
        fading = exp(-(0:channel_taps));
        fading = fading/norm(fading);
        noisy_signal = conv(noisy_signal, fading ,'same');
        
    elseif fading_type == 2
        fading = randn(size(noisy_signal)) + 1i * randn(size(noisy_signal));
        noisy_signal = noisy_signal .* fading;

    elseif fading_type == 3
        LOS_component = sqrt(K / (K + 1)) * exp(1i * pi * randn(size(noisy_signal)));
        scattered_component = sqrt(1 / (K + 1)) * (randn(size(noisy_signal)) + 1i * randn(size(noisy_signal)));
        fading = LOS_component + scattered_component;
        noisy_signal = (LOS_component + scattered_component) .* noisy_signal;

    elseif fading_type == 4
        block_size = 1000;
        num_blocks = ceil(length(noisy_signal) / block_size);
    
        for block = 1:num_blocks
            start_idx = (block - 1) * block_size + 1;
            end_idx = min(block * block_size, length(noisy_signal));
            fading_block = sqrt(random('gamma', m/2, 2, 1, end_idx - start_idx + 1)); 
            noisy_signal(start_idx:end_idx) = noisy_signal(start_idx:end_idx) .* fading_block.';
        end
    end
end

without_prefix = reshape(noisy_signal, (subcarriers + cyclic_prefix), (length(noisy_signal)/(subcarriers + cyclic_prefix)));
noisy_signal = without_prefix((cyclic_prefix +1 : end), :);

noisy_signal = fft(noisy_signal);

if flag ==2
    null_subcarriers = 1:16;
    channel_estimate = mean(noisy_signal(:, null_subcarriers), 2);
    noisy_signal = noisy_signal ./ repmat(channel_estimate, 1, size(noisy_signal, 2));
end

noisy_signal = noisy_signal(:);
noisy_signal = noisy_signal(1 : end - missing_subcarriers);

symbols = [real(OFDM_symbols) imag(OFDM_symbols)];
if size(OFDM_symbols, 2) > 1
    symbols = [real(OFDM_symbols(:)) imag(OFDM_symbols(:))];
end
received_symbols = knnsearch(symbols, [real(noisy_signal) imag(noisy_signal)]) - 1;

if bits == 1
    modulation_order = 2^1;
elseif bits == 2 || bits == 3
    modulation_order = 2^2;
elseif bits == 4 || bits == 6
    modulation_order = sqrt(2^bits);
elseif bits == 5
    modulation_order = sqrt(2^bits + 4);
end

received_symbols_binary = dec2bin(received_symbols, log2(modulation_order));
received_image_binary = reshape(received_symbols_binary', numel(received_symbols_binary), 1);
received_image_binary = received_image_binary(1 : end - missing_bits);
binary_form_adjusted = binary_form(1 : end - missing_bits);
num_errors = sum(received_image_binary ~= binary_form_adjusted);
BER = num_errors / (length(binary_form_adjusted));

received_image = reshape(received_image_binary, 8, (numel(received_image_binary)/8));
received_image = uint8(bin2dec(received_image'));
received_image = reshape(received_image,size(image));

if bits == 1
    type = 'BPSK Modulation is used';
    constellation = 'BPSK Signal Constellation';
elseif bits == 2
    type = 'QPSK Modulation is used';
    constellation = 'QPSK Signal Constellation';
elseif bits == 3
    type = '8-PSK Modulation is used';
    constellation = '8-PSK Signal Constellation';
elseif bits == 4
    type = '16-QAM Modulation is used';
    constellation = '16-QAM Signal Constellation';
elseif bits == 5
    type = '32-QAM Modulation is used';
    constellation = '32-QAM Signal Constellation';
elseif bits == 6
    type = '64-QAM Modulation is used';
    constellation = '64-QAM Signal Constellation';
end

figure(1)
sgtitle(type);

subplot(2,1,1);
imshow(image);
title('Original Image');

subplot(2,1,2);
imshow(received_image);
title(sprintf('Received Image\n \\rmBER: %.2g',BER));

figure(2)
sgtitle(constellation);

subplot(2,1,1);
plot(modulated_symbols,'x','linewidth',2,'markersize',10);
xlim([-2 2]);
ylim([-2 2]);
xlabel('In phase')
ylabel('Qudrature')
title('Transmitted Signal Constellation');
grid on;

subplot(2,1,2);
plot(noisy_signal(1:2000:end),'x','markersize',1);
xlim([-4 4]);
ylim([-4 4]);
xlabel('In phase')
ylabel('Qudrature')
title(sprintf('Received Signal Constellation\n \\rmSNR: %d', snr_db));
grid on;
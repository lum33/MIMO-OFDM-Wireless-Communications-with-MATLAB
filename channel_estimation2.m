% channel_estimation2.m 
% Plotting the MSE vs SNR(dB)
% for LS/DFT Channel Estimation with linear/spline interpolation
Nfft=32; 
Ng=Nfft/8;                                                                  % length of cyclic prefix
Nofdm=Nfft+Ng;                                                              % total number of a OFDM symbol
Nsym=100;                                                                   % Number of OFDM symbol
Nps=4;                                                                      % ps = pilot spacing
Np=Nfft/Nps;                                                                % p = number of pilots per OFDM symbol|Np=8 
Nbps=4;                                                                     % bps = bits per modulated symbol
M=2^Nbps;                                                                   % Number of possible modulated symbol

% DELETED CODE 
% mod_object = modem.qammod('M',M, 'SymbolOrder','gray');
% demod_object = modem.qamdemod('M',M, 'SymbolOrder','gray');

Es=1;                                                                       % Signal energy 
A=sqrt(3/2/(M-1)*Es);                                                       % QAM normalization factor
% SNR = 30; 
sq2=sqrt(2); 

MSE_LSl = zeros(1,30);                                                          % Initial of MSE LS linear interpolation vectors
MSE_LSs = zeros(1,30);                                                          % Initial of MSE LS linear interpolation vectors
MSE_MMSE = zeros(1,30);                                                          % Initial of MSE LS linear interpolation vectors

nose = 0;                                                                   % Initial of number of symbol errors
% For each OFDM symbol(在Nofdm個OFDM symbol中有Nfft個subcarrier)
for nsym=1:Nsym
    Xp = 2*(randn(1,Np)>0)-1;                                               % Pilot sequence generation
    msgint=randi(1,Nfft-Np,M);                                              % bit generation|Nfft-Np=24,M=16
    % Adjusted part
    % Data = A*modulate(mod_object,msgint);
    Data = qammod(msgint,M)*A;
    ip = 0;                                                                 % Initial of index of pilot position
    pilot_loc = [];                                                         % Initial of pilot location
    
    % For each subcarrier, do the following:
    for k=1:Nfft
        if mod(k,Nps)==1                                                    % If the subcarrier is a pilot subcarrier
            X(k)=Xp(floor(k/Nps)+1);                                        % Assign the corresponding pilot value to the subcarrier.
            pilot_loc=[pilot_loc k];                                        % Store the index of the pilot subcarrier.
            ip = ip+1;
        else
            X(k) = Data(k-ip);                                              % Assign the corresponding modulated data symbol to the subcarrier.
        end
    end

    x = ifft(X,Nfft); xt = [x(Nfft-Ng+1:Nfft) x];                          	% IFFT and add CP
    h = [(randn+j*randn) (randn+j*randn)/2];                               	% A (2-tap) channel | 用於模擬multipath fading Model
    H = fft(h,Nfft); ch_length=length(h);                                  	% True channel and its length
    H_power_dB = 10*log10(abs(H.*conj(H)));                                	% True channel power in dB
    y_channel = conv(xt,h);                                                	% Channel path (convolution)

    for SNR=1:30
        yt = awgn(y_channel,SNR,'measured');                                    
        y = yt(Ng+1:Nofdm); Y = fft(y);                                         % Remove CP and FFT
    
        % 前半部為處理OFDM的mod和demod,後半部為channel estimation
        for m=1:3
            if m==1 
                H_est = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,'linear');
                method='LS-linear';                                             % LS estimation with linear interpolation
            elseif m==2
                H_est = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,'spline');
                method='LS-spline';                                         	% LS estimation with spline interpolation
            else
                H_est = MMSE_CE(Y,Xp,pilot_loc,Nfft,Nps,h,SNR);
                method='MMSE';                                                  % MMSE estimation
            end
        
            % Calculate the estimated channel power spectrum
            H_est_power_dB = 10*log10(abs(H_est.*conj(H_est)));
            h_est = ifft(H_est); h_DFT = h_est(1:ch_length);
            if m==1 
                MSE_LSl(SNR) = MSE_LSl(SNR) + (H-H_est)*(H-H_est)';
            elseif m==2 
                MSE_LSs(SNR) = MSE_LSs(SNR) + (H-H_est)*(H-H_est)';
            else
                MSE_MMSE(SNR) = MSE_MMSE(SNR) + (H-H_est)*(H-H_est)';
            end
        end
    end
    
    Y_eq = Y./H_est; 
    ip = 0;
    
    % Extract the data from received signal
    Data_extracted = [];
    for k=1:Nfft
        if mod(k,Nps)==1
            ip=ip+1; 
        else
            Data_extracted(k-ip)=Y_eq(k);
        end
    end
    % Adjusted part
    % msg_detected = demodulate(demod_object,Data_extracted/A);
    msg_detected = qamdemod(Data_extracted'/A,M);                           % Demodulate the extracted data
    nose = nose + sum(msg_detected~=msgint);                                % calculate the number of errors
    % MSEs = MSE/(Nfft*Nsym);                                                 % calculate the MSE
end

% Plotting the MSE vs SNR(dB)
x = linspace(1,30,30);
z1 = MSE_LSl(x);
plot(x,z1,'--ro');

hold on;
z2 = MSE_LSs(x);
plot(x,z2,'--go');
hold off;

hold on;
z3 = MSE_MMSE(x);
plot(x,z3,'--bo');
hold off;

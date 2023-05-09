% channel_estimation2.m 
% Aim :
%       1. LS/MMSE Channel Estimation
%       2. Plotting the MSE vs SNR(dB)
Nfft=32; 
Ng=Nfft/8;                                                                  % length of cyclic prefix
Nofdm=Nfft+Ng;                                                              % total number of a OFDM symbol
Nsym=100;                                                                   % Number of OFDM symbol
Nps=4;                                                                      % ps = pilot spacing
Np=Nfft/Nps;                                                                % p = number of pilots per OFDM symbol|Np=8 
Nbps=4;                                                                     % bps = bits per modulated symbol
M=2^Nbps;                                                                   % Number of possible modulated symbol

Es=1;                                                                       % Signal energy 
A=sqrt(3/2/(M-1)*Es);                                                       % QAM normalization factor
sq2=sqrt(2); 
% Initialization of MSE matrix
MSE_LSl = zeros(100,32);                                                      
MSE_LSs = zeros(100,32);                                                      
MSE_MMSE = zeros(100,32);
% Initialization of average MSE vectors
ave_MSE1=zeros(1,32);
ave_MSE2=zeros(1,32);
ave_MSE3=zeros(1,32);

nose = 0;                                                                   % Initial of number of symbol errors
for sample=1:100
    % For each OFDM symbol, do the following
    for nsym=1:Nsym
        Xp = 2*(randn(1,Np)>0)-1;                                               % Pilot sequence generation
        msgint=randi(1,Nfft-Np,M);                                              % bit generation|Nfft-Np=24,M=16
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

        for SNR=1:32                                                            % channel est. in the SNR range(-1, 30)
            yt = awgn(y_channel,SNR-2,'measured');                                    
            y = yt(Ng+1:Nofdm); Y = fft(y);                                     % Remove CP and FFT
            for m=1:3
                if m==1 
                    H_est = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,'linear');
                    method='LS-linear';                                         % LS estimation with linear interpolation
                elseif m==2
                    H_est = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,'spline');
                    method='LS-spline';                                         % LS estimation with spline interpolation
                else
                    H_est = MMSE_CE(Y,Xp,pilot_loc,Nfft,Nps,h,SNR);
                    method='MMSE';                                              % MMSE estimation
                end

                % Calculate the estimated channel power spectrum
                H_est_power_dB = 10*log10(abs(H_est.*conj(H_est)));
                h_est = ifft(H_est); h_DFT = h_est(1:ch_length);

                % Calculate MSE 
                if m==1 
                    MSE_LSl(sample,SNR) = MSE_LSl(sample,SNR) + (H-H_est)*(H-H_est)';
                elseif m==2 
                    MSE_LSs(sample,SNR) = MSE_LSs(sample,SNR) + (H-H_est)*(H-H_est)';
                else
                    MSE_MMSE(sample,SNR) = MSE_MMSE(sample,SNR) + (H-H_est)*(H-H_est)';
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
        msg_detected = qamdemod(Data_extracted'/A,M);                           % Demodulate the extracted data
        nose = nose + sum(msg_detected~=msgint);                                % calculate the number of errors
        % MSEs = MSE/(Nfft*Nsym);                                               
    end
end


% Calculate average MSE
for sample=1:100
    for SNR=1:32
        ave_MSE1(SNR) = mean(MSE_LSl(sample,SNR))/(Nfft*Nsym);
        ave_MSE2(SNR) = mean(MSE_LSs(sample,SNR))/(Nfft*Nsym);
        ave_MSE3(SNR) = mean(MSE_MMSE(sample,SNR))/(Nfft*Nsym);
    end
end

% Plotting the MSE vs SNR(dB)
x = linspace(-1,30,32);
y1 = ave_MSE1(x+2);
plot(x,y1,'--ro');

hold on;
y2 = ave_MSE2(x+2);
plot(x,y2,'--go');
hold off;

hold on;
y3 = ave_MSE3(x+2);
plot(x,y3,'--bo');
hold off;

xlabel('SNR(dB)')
ylabel('average MSE')
title('average MSE v.s. SNR(dB) graph')
legend('LS-linear','LS-spline','MMSE');

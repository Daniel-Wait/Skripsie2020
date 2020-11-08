%--------------------------------%
%             Setup              %
%--------------------------------%

fs = 44100;
Ts = 1/fs;

dsamp = 2;
Tb = 20e-3;
samps = Tb*fs;

n = 1;
delta_f = n/(2*Tb);

bw_ch = 100;
sep_f = 200;

f1 = 2200;
f0 = 1600;
f1 = [f1-sep_f f1];
f0 = [f0 f0+sep_f];

L = 2048
N = 4096
K = 3

bits = [0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1];

%Import Data from Table
T = table2array(teensydata);
T_len = size(T, 1)

%Display Teensy UPOLA Buffer Output
corr_out = T(T_len - K*N + 1 : T_len);
figure('Name', 'Teensy UPOLA Output')
plot(corr_out)

%Generate Template
tx = -1*cdma_chirp_fsk(bits, samps/2, f0*dsamp, f1*dsamp, fs, dsamp*bw_ch);
figure('Name', 'MATLAB TX Message')
plot(tx)

%Separate Downsampled Float Data and Display 
rx = T(T_len - 2*K*N + 1 : T_len - K*N);
rx = reshape(rx, 1, K*N);
figure('Name', 'Teensy RX Downsampled F32 Message')
plot(rx)

%Display Auto-Corr of Template
[rtt,lags] = xcorr(tx, tx);
figure('Name', 'MATLAB xcorr self')
plot(lags, rtt)

%Display Expected xCorr of Template & Downsampled Float Data
[rtr,lags] = xcorr(rx, tx);
figure('Name', 'MATLAB xcorr output')
plot(lags, rtr)

%--------------------------------%
%MATLAB script of UPOLA Algorithm%
%--------------------------------%

%Partitioned Filter FFTS
mf_coef= flip(tx);
mf_coef_n = zeros(K,N);
mf_fft_n = zeros(K,N);

for j  = 1:K
    start = (j-1)*L;
    stops = start+L;
    if (stops > size(mf_coef,2))
        stops = size(mf_coef,2);
    end
    mf_coef_n(j,:) = [mf_coef((start+1):stops) zeros(1, N-(stops-start))];
    mf_fft_n(j,:) = fft( mf_coef_n(j,:) );
end


%Zeropad & Read RX Data into arrays
xf = zeros(K*4,N);
xn = zeros(K*4,N);
yf = zeros(K*3,N);
yn = zeros(K*3,N);
ola = zeros(1,N+12*L);

rx = [rx zeros(1, L*3);]

for i = 1:K*3
  cnt = (i-1)*L+1;
  xn(i+K,:) = [rx(cnt:(cnt+L-1)) zeros(1, N-L)];
  xf(i+K,:) = fft( xn(i+K,:) );
end

%Freq. Domain Convolution
for i = K+1:K*4
    for j = 1:K
        yf(i-K,:) = yf(i-K,:) + xf(i-j,:).*mf_fft_n(j,:);
    end
    yn(i-K,:) = ifft( yf(i-K,:) );
end


%Overlap Management
ola(1:L) = yn(1, 1:L);
ola(L+1:N) = yn(1, L+1:N);

for i = 1:K*3
  ola( i*L+1 : N+(i-1)*L ) = ola( i*L+1 : N+(i-1)*L ) + yn(i, 1:N-L);
  ola( N+(i-1)*L+1 : (i+1)*L ) = yn(i, N-L+1:L);
  ola( (i+1)*L+1 : N+i*L) = yn(i, L+1:N);
end

figure('Name', 'MATLAB Block Convolver Output')
plot(ola)


%--------------------------------%
%      Modultaion Function       %
%--------------------------------%

function mn = cdma_chirp_fsk(bitseq, spb, f_0, f_1, F, df)
    T = 1/F;
    nb = size(bitseq, 2);
    num = 1:spb;
    mn = zeros(1, spb*nb);
    mn_len = 1:size(mn,2);
    chirp = (df/spb).*(0:(spb-1));
    cnt0 = 0;
    cnt1 = 0;
    
    for i = 1:nb
        if (bitseq(i) == 1)
            fi = f_1( mod(cnt1, size(f_1, 2))+1 );
            fi = fi + chirp - df/2;
            bit = sin(2*pi*T*(fi.*num));
            cnt1 = cnt1 + 1;
        end
        if (bitseq(i) == 0)
            fz = f_0( mod(cnt0,  size(f_0, 2))+1 );
            fz = fz - chirp + df/2;
            bit = sin(2*pi*T*(fz.*num));
            cnt0 = cnt0 + 1;
        end
        mn((i-1)*spb+1:i*spb) = bit;
    end
end
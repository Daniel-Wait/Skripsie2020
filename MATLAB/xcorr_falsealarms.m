fs = 44100;
Ts = 1/fs;

Tb = 20e-3;
samps = Tb*fs; 

delta_f = 1/(2*Tb);


f1 = delta_f*80;
f0 = delta_f*40;


fi = [f1-250, f1];
fz = [f0, f0+250];

dsamp = 2;
blocks = 3;
L = 2048;
N = L*2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fhss_codes12 = [
 [0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1,],
 [1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1,],
 [1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0,],
 [0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1,],
 [0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0,],
 [0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1,],
 [1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0,],
 [0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0,],
 [1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1,],
 [1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0,],
 [0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1,],
 [1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1,],
 [0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0,],
 ]

BITS = 12;
NUM = size(fhss_codes12, 1);

bit_begin = 1;
bit_end = bit_begin + BITS -1; 
trgcnt_begin = bit_end + 1;
trgcnt_end = trgcnt_begin + NUM -1;
xcorrmax_begin = trgcnt_end + 1;
xcorrmax_end = xcorrmax_begin + NUM -1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sat = @(x, delta) min(max(x/delta, -1), 1);

B = zeros(NUM, xcorrmax_end);

B(:, bit_begin:bit_end) = fhss_codes12;

% rv = cdma_chirp_fsk( B( 2, bit_begin : bit_end ), samps, fz, fi, fs, delta_f*4);
% rv = downsample(rv, dsamp);
% rv_polluted = sat(awgn(rv, SNRdb),1);
%nz = 558;
%rv_polluted = 0.7*rv + [zeros(1, nz) 0.3*rv(1, nz:5292-1)]
%rv_polluted = awgn(rv_polluted, SNRdb);
%rv_polluted = sat(rv_polluted,1);
%figure('Name', 'RX Message')
%plot(rv_polluted)

SNRdb = 5;
SNR = 10^(SNRdb/10);

thresh = ((BITS+5)/2);
thresh = 8.5;

noise_stddev = 1/SNR;
xcorr_thresh = thresh*(samps/(2*dsamp));

for j = 1:NUM
    j
    template = cdma_chirp_fsk(B( j, bit_begin : bit_end ), samps/dsamp, dsamp*fz, dsamp*fi, fs, dsamp*delta_f*4 );
    %template = chirp_fsk(B( j, bit_begin : bit_end ), samps/dsamp, dsamp*f0, dsamp*f1, fs, dsamp*delta_f*4);
    %template = fsk(B( j, bit_begin : bit_end ), samps/dsamp, dsamp*f0, dsamp*f1, fs);    
    
    for i = 1:NUM
        rx = cdma_chirp_fsk( B( i, bit_begin : bit_end ), samps, fz, fi, fs, delta_f*4);
        %rx = chirp_fsk(B( i, bit_begin : bit_end ), samps, f0, f1, fs, delta_f*4);
        %rx = fsk( B( i, bit_begin : bit_end ), samps, f0, f1, fs);
        rx = downsample(rx, dsamp);

        for n = 1:2000      
            nz = 558;
            %rx_polluted = 0.6*rx;
            rx_polluted = 0.7*rx + [zeros(1, nz) 0.4*rx(1, nz:5292-1)];
            rx_polluted = awgn(rx_polluted, SNRdb);
            
            rx_polluted = sat(rx_polluted,1);

            rtr = xcorr( rx_polluted, template );
            
            if ( max(rtr) >= xcorr_thresh )
                B(j, trgcnt_begin + i-1) = B(j, trgcnt_begin + i-1) + 1;
            end
            
            if ( max(rtr) > B(j, xcorrmax_begin+ i-1) )
                B(j, xcorrmax_begin+ i-1) = max(rtr);
            end           
          
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mn = fsk(bitseq, spb, f_0, f_1, F)
    T = 1/F;
    nb = size(bitseq, 2);
    num = 1:spb;
    mn = zeros(1, spb*nb);
    mn_len = 1:size(mn,2);
    
    for i = 1:nb
        if (bitseq(i) == 1) 
            bit = sin(2*f_1*pi*num*T);
        end
        if (bitseq(i) == 0)
            bit = sin(2*f_0*pi*num*T);
        end
        mn((i-1)*spb+1:i*spb) = bit;
    end
 
    %figure('Name', 'Modulated Signal')
    %plot(mn)
end


function mn = chirp_fsk(bitseq, spb, f_0, f_1, F, df)
    T = 1/F;
    nb = size(bitseq, 2);
    num = 1:spb;
    mn = zeros(1, spb*nb);
    mn_len = 1:size(mn,2);
    chirp = (df/spb).*(0:(spb-1));
    f_0 = f_0 - chirp;
    f_1 = f_1 + chirp;
   
    
    for i = 1:nb
        if (bitseq(i) == 1) 
            bit = sin(2*pi*T*(f_1.*num));
        end
        if (bitseq(i) == 0)
            bit = sin(2*pi*T*(f_0.*num));
        end
        mn((i-1)*spb+1:i*spb) = bit;
    end
end


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
            fz = fz - chirp - df/2;
            bit = sin(2*pi*T*(fz.*num));
            cnt0 = cnt0 + 1;
        end
        mn((i-1)*spb+1:i*spb) = bit;
    end
end
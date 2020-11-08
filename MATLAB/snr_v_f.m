f_in = linspace(1, 3.5, 10000);
snr_out = arrayfun(@snr_est, f_in);
figure('Name', 'SNR vs Freq.')
plot(f_in, snr_out)

function SNR = snr_est(f)
    %Init. Vars
    r = 40;         %distance m         
    k = 2;          %Spreading factor

    s = 0.5;        %shipping activity
    wkmh = 22;
    w = wkmh/3.6;   %wind speed m/s

    TVR = 110;      %hydrophone transmitter dB re uPa.m/V
    Vpk = 15;       %projector input Volts


    %Receiver Bandwidth
    BWrx = 600;       

    %Signal Level
    Vin = Vpk/sqrt(2);
    SL = TVR + 20*log10(Vin)      


    %Ambient Noise
    Nt_db  = 17 - 30*log10(f); %turbulence
    Ns_db  = 40 + 20*(s-0.5) + 26*log10(f)- 60*log10(f +0.03); %shipping
    Nw_db  = 50 + 7.5*(w^(1/2)) + 20*log10(f) - 40*log10(f +0.4); %wind
    Nth_db = -15 + 20*log10(f); %thermal

    N = 10^(Nt_db/10) + 10^(Ns_db/10) + 10^(Nw_db/10) + 10^(Nth_db/10);
    N_db = 10*log10(N)


    %Thorp
    af_dB = 0.11*(f^2)/(1+(f^2)) + 44*(f^2)/(4100+(f^2)) + 2.75*(10^(-4))*(f^2) + 0.003 % a_f = absorption coefficient

    %Loss dB
    PL = k*10*log10(r) + (r/1000)*(af_dB) ;

    %SNR estimation
    NL = ( N_db + 10*log10(BWrx) );
    SNR = SL - PL - NL;
end

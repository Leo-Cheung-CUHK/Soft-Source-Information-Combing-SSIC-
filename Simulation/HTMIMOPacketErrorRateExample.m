% Create a format configuration object for a 2-by-2 HT transmission
cfgHT = wlanHTConfig;
cfgHT.ChannelBandwidth = 'CBW20'; % 20 MHz channel bandwidth
cfgHT.NumTransmitAntennas = 1;    % 2 transmit antennas
cfgHT.NumSpaceTimeStreams = 1;    % 2 space-time streams
cfgHT.PSDULength = 20;          % PSDU length in bytes
cfgHT.MCS = 5;                   % 2 spatial streams, 64-QAM rate-5/6
cfgHT.ChannelCoding = 'LDPC';      % BCC channel coding

%% Simulation Parameters
snr = 8:1:18;

maxNumPEs = 100; % The maximum number of packet errors at an SNR point
maxNumPackets = 10000; % The maximum number of packets at an SNR point

%% Set the remaining variables for the simulation.
% Get the baseband sampling rate
fs = wlanSampleRate(cfgHT);

% Get the OFDM info
ofdmInfo = wlanHTOFDMInfo('HT-Data',cfgHT);

% Set the sampling rate of the channel
tgnChannel.SampleRate = fs;

% Indices for accessing each field within the time-domain packet
ind = wlanFieldIndices(cfgHT);

%% Initialization phase for our method
K = 7;
%%% derive a_{i,k}
[a_total,S_total] = DetMatRep(127,K);
L = 16;

%% Processing SNR Points
S = numel(snr);

CONV_SeedbitErrorRate = zeros(S,1);
SD_SeedbitErrorRate = zeros(S,1);

CONV_bitErrorRate = zeros(S,1);
CONV_packetErrorRate = zeros(S,1);

HRSX_bitErrorRate = zeros(S,1);
HRSX_packetErrorRate = zeros(S,1);

SRSX_bitErrorRate = zeros(S,1);
SRSX_packetErrorRate = zeros(S,1);

%parfor i = 1:S % Use 'parfor' to speed up the simulation
for i = 1:S % Use 'for' to debug the simulation
    % Set random substream index per iteration to ensure that each
    % iteration uses a repeatable set of random numbers
    stream = RandStream('combRecursive','Seed',0);
    stream.Substream = i;
    RandStream.setGlobalStream(stream);

    % Account for noise energy in nulls so the SNR is defined per
    % active subcarrier
    packetSNR = snr(i);

    % Loop to simulate multiple packets
    CONV_SeednumBitErrors = 0;
    SD_SeedbitError = 0;

    CONV_numBitErrors = 0;
    CONV_numPacketErrors = 0;
    HRSX_numBitErrors = 0;
    HRSX_numPacketErrors = 0;
    SRSX_numBitErrors = 0;
    SRSX_numPacketErrors = 0;

    n = 1; % Index of packet transmitted
    numPacketErrors = 0;
    while numPacketErrors<=maxNumPEs && n<=maxNumPackets
        scramInit = randi([1,127]);

        % Generate a packet waveform
        txPSDU = randi([0 1],cfgHT.PSDULength*8,1); % PSDULength in bytes

        %%% note that we still need extra 127-16 = 111 zeros bits
        if L>16
            txPSDU(1:(L-16)+1) = 0;
        end

        if (L>16)
            validate_txPSDU = txPSDU((L-16)+1:end);
        else
            validate_txPSDU = txPSDU;
        end

        tx = wlanWaveformGenerator(txPSDU,cfgHT,'ScramblerInitialization',scramInit);
        rx = awgn(tx,packetSNR,'measured');

        %%
        % Packet detect and determine coarse packet offset
        coarsePktOffset = wlanPacketDetect(tx,cfgHT.ChannelBandwidth);
        if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
            continue; % Go to next loop iteration
        end

        % Extract L-STF and perform coarse frequency offset correction
        lstf = tx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
        coarseFreqOff = wlanCoarseCFOEstimate(lstf,cfgHT.ChannelBandwidth);
        tx = helperFrequencyOffset(tx,fs,-coarseFreqOff);

        % Extract the non-HT fields and determine fine packet offset
        nonhtfields = tx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
        finePktOffset = wlanSymbolTimingEstimate(nonhtfields,...
            cfgHT.ChannelBandwidth);

        % Determine final packet offset
        pktOffset = coarsePktOffset+finePktOffset;

        % If packet detected outwith the range of expected delays from the
        % channel modeling; packet error
        if pktOffset>15
            continue; % Go to next loop iteration
        end

        % Extract L-LTF and perform fine frequency offset correction
        lltf = tx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
        fineFreqOff = wlanFineCFOEstimate(lltf,cfgHT.ChannelBandwidth);
        tx = helperFrequencyOffset(tx,fs,-fineFreqOff);

        % Extract HT-LTF samples from the waveform, demodulate and perform
        % channel estimation
        htltf = tx(pktOffset+(ind.HTLTF(1):ind.HTLTF(2)),:);
        htltfDemod = wlanHTLTFDemodulate(htltf,cfgHT);
        chanEst = wlanHTLTFChannelEstimate(htltfDemod,cfgHT);

        % Extract HT Data samples from the waveform
        htdata = tx(pktOffset+(ind.HTData(1):ind.HTData(2)),:);

        % Estimate the noise power in HT data field
        nVarHT = htNoiseEstimate(htdata,chanEst,cfgHT);

        % Recover the transmitted scrambled PSDU in HT Data
        txPSDU_scrambled = wlanHTCodedDataRecover(htdata,chanEst,nVarHT,cfgHT);

        % Conventional Descraming
        [~,Tx_scramInitBits,~] = wlanHTDescramer(double(txPSDU_scrambled < 0),cfgHT.PSDULength);

        %%
        % Pass the waveform through the TGn channel model
        % reset(tgnChannel); % Reset channel for different realization
        % rx = tgnChannel(tx);

        % Packet detect and determine coarse packet offset
        coarsePktOffset = wlanPacketDetect(rx,cfgHT.ChannelBandwidth);
        if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
            numPacketErrors = numPacketErrors+1;
            continue; % Go to next loop iteration
        end

        % Extract L-STF and perform coarse frequency offset correction
        lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
        coarseFreqOff = wlanCoarseCFOEstimate(lstf,cfgHT.ChannelBandwidth);
        rx = helperFrequencyOffset(rx,fs,-coarseFreqOff);

        % Extract the non-HT fields and determine fine packet offset
        nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
        finePktOffset = wlanSymbolTimingEstimate(nonhtfields,...
            cfgHT.ChannelBandwidth);

        % Determine final packet offset
        pktOffset = coarsePktOffset+finePktOffset;

        % If packet detected outwith the range of expected delays from the
        % channel modeling; packet error
        if pktOffset>15
            numPacketErrors = numPacketErrors+1;
            continue; % Go to next loop iteration
        end

        % Extract L-LTF and perform fine frequency offset correction
        lltf = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
        fineFreqOff = wlanFineCFOEstimate(lltf,cfgHT.ChannelBandwidth);
        rx = helperFrequencyOffset(rx,fs,-fineFreqOff);

        % Extract HT-LTF samples from the waveform, demodulate and perform
        % channel estimation
        htltf = rx(pktOffset+(ind.HTLTF(1):ind.HTLTF(2)),:);
        htltfDemod = wlanHTLTFDemodulate(htltf,cfgHT);
        chanEst = wlanHTLTFChannelEstimate(htltfDemod,cfgHT);

        % Extract HT Data samples from the waveform
        htdata = rx(pktOffset+(ind.HTData(1):ind.HTData(2)),:);

        % Estimate the noise power in HT data field
        nVarHT = htNoiseEstimate(htdata,chanEst,cfgHT);

        % Recover the transmitted scrambled PSDU in HT Data
        rxPSDU_scrambled = wlanHTCodedDataRecover(htdata,chanEst,nVarHT,cfgHT);

        %%% Conventional Receiver
        % Hard-decision making
        decBits = double(rxPSDU_scrambled < 0);
        % Conventional Descraming
        [rxPSDU,Rx_scramInitBits,~] = wlanHTDescramer(decBits,cfgHT.PSDULength);

        % Determine if any bits are in error, i.e. a packet error
        CONV_SeedbitError = biterr(Tx_scramInitBits,Rx_scramInitBits);
        CONV_bitError     = biterr(validate_txPSDU,rxPSDU(L-16+1:end));
        CONV_packetError  = any(CONV_bitError);

        CONV_SeednumBitErrors = CONV_SeednumBitErrors+CONV_SeedbitError;
        CONV_numBitErrors     = CONV_numBitErrors+CONV_bitError;
        CONV_numPacketErrors  = CONV_numPacketErrors+CONV_packetError;

        %% Our Method
        yllr = double(rxPSDU_scrambled); % resevese it for our method
        l_total = yllr(1:L); % extract the first 127 LLR
        % derive f(r_1,r_2,...,r_7)
        [f_total,r_total] = DeriveF(a_total,l_total,K,L);
        %% hard-descrambling soft-decoding (HDSC)
        [~,~,Rx_scramInitBits,txBits,errorBits,pass] = Hard_NonHTDescramer_llr(cfgHT.PSDULength,f_total,r_total,yllr,L,txPSDU);
        SD_SeedbitError = SD_SeedbitError + biterr(Tx_scramInitBits,Rx_scramInitBits);
        if pass
            if errorBits~=0
                HRSX_numBitErrors    = HRSX_numBitErrors + errorBits;
                HRSX_numPacketErrors = HRSX_numPacketErrors + 1;
            end
        end
        %% Soft-descrambling soft-decoding (SDSC)
        [~,~,errorBits,pass] = Soft_NonHTDescramer_llr(cfgHT.PSDULength,a_total,f_total,r_total,yllr,L,txPSDU);
        if pass
            if errorBits~=0
                SRSX_numBitErrors    = SRSX_numBitErrors + errorBits;
                SRSX_numPacketErrors = SRSX_numPacketErrors + 1;
            end
        end
        n = n+1;
    end

    %% CONV
    % Calculate Seed's bit error rate (BER) at SNR point
    CONV_SeedbitErrorRate(i) = CONV_SeednumBitErrors/((n-1)*7);
    disp(['SNR ' num2str(snr(i))...
        ' completed after '  num2str(n-1) ' packets,'...
        'Seed BER: ' num2str(CONV_SeedbitErrorRate(i))]);

    % Calculate Seed's bit error rate (BER) at SNR point
    SD_SeedbitErrorRate(i) = SD_SeedbitError/((n-1)*7);
    disp(['SNR ' num2str(snr(i))...
        ' completed after '  num2str(n-1) ' packets,'...
        'Seed BER: ' num2str(SD_SeedbitErrorRate(i))]);

    % Calculate bit error rate (BER) at SNR point
    CONV_bitErrorRate(i) = CONV_numBitErrors/[(n-1)*length(validate_txPSDU)];
    disp(['SNR ' num2str(snr(i))...
        ' completed after '  num2str(n-1) ' packets,'...
        ' BER: ' num2str(CONV_bitErrorRate(i))]);

    % Calculate packet error rate (PER) at SNR point
    CONV_packetErrorRate(i) = CONV_numPacketErrors/(n-1);
    disp(['SNR ' num2str(snr(i))...
        ' completed after '  num2str(n-1) ' packets,'...
        ' PER: ' num2str(CONV_packetErrorRate(i))]);

    %% HRSX
    % Calculate bit error rate (BER) at SNR point
    HRSX_bitErrorRate(i) = HRSX_numBitErrors/((n-1)*length(validate_txPSDU));
    disp(['SNR ' num2str(snr(i))...
        ' completed after '  num2str(n-1) ' packets,'...
        ' BER: ' num2str(HRSX_bitErrorRate(i))]);

    % Calculate packet error rate (PER) at SNR point
    HRSX_packetErrorRate(i) = HRSX_numPacketErrors/(n-1);
    disp(['SNR ' num2str(snr(i))...
        ' completed after '  num2str(n-1) ' packets,'...
        ' PER: ' num2str(HRSX_packetErrorRate(i))]);

    % Calculate bit error rate (BER) at SNR point
    SRSX_bitErrorRate(i) = SRSX_numBitErrors/((n-1)*length(validate_txPSDU));
    disp(['SNR ' num2str(snr(i))...
        ' completed after '  num2str(n-1) ' packets,'...
        ' BER: ' num2str(SRSX_bitErrorRate(i))]);

    % Calculate packet error rate (PER) at SNR point
    SRSX_packetErrorRate(i) = SRSX_numPacketErrors/(n-1);
    disp(['SNR ' num2str(snr(i))...
        ' completed after '  num2str(n-1) ' packets,'...
        ' PER: ' num2str(SRSX_packetErrorRate(i))]);
end
%% Plot Bit Error Rate vs SNR Results
figure;
semilogy(snr,CONV_SeedbitErrorRate,'--');
hold on;
semilogy(snr,SD_SeedbitErrorRate,'-o');
grid on;
xlabel('SNR [dB]');
ylabel('Seed BER');
legend('CONV','SD')

%% Plot Bit Error Rate vs SNR Results
figure;
semilogy(snr,CONV_bitErrorRate,'--');
hold on;
semilogy(snr,HRSX_bitErrorRate,'-o');
hold on;
semilogy(snr,SRSX_bitErrorRate,'-*');

grid on;
xlabel('SNR [dB]');
ylabel('BER');
legend('CONV','HRSX','SRSX')

%% Plot Packet Error Rate vs SNR Results
figure;
semilogy(snr,CONV_packetErrorRate,'--');
hold on;
semilogy(snr,HRSX_packetErrorRate,'-o');
hold on;
semilogy(snr,SRSX_packetErrorRate,'-*');
grid on;
xlabel('SNR [dB]');
ylabel('PER');
legend('CONV','HRSX','SRSX')

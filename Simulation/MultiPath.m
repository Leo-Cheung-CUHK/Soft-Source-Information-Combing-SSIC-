clear all

% Create a format configuration object for a 2-by-2 HT transmission
cfgHT = wlanHTConfig;
cfgHT.ChannelBandwidth = 'CBW20'; % 20 MHz channel bandwidth
cfgHT.NumTransmitAntennas = 1;    % 2 transmit antennas
cfgHT.NumSpaceTimeStreams = 1;    % 2 space-time streams
cfgHT.PSDULength = 100;          % PSDU length in bytes
cfgHT.MCS = 5;                   % 2 spatial streams, 64-QAM rate-5/6
cfgHT.ChannelCoding = 'LDPC';      % BCC channel coding

%% Simulation Parameters
% basic_range = 8:1:18;
% Path_Number = 4;
% SNR_vector = zeros(Path_Number,length(basic_range));

% for path_index = 1:Path_Number
%     SNR_vector(path_index,:) = basic_range+(path_index-1)*0.5;
% end

basic_range = 9.5:1:19.5;
Path_Number = 4;
SNR_vector = zeros(Path_Number,length(basic_range));

for path_index = 1:Path_Number
    SNR_vector(path_index,:) = basic_range;
end

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
L = 127;

%% Processing SNR Points
S = numel(basic_range);

CONV_SeedbitErrorRate     = zeros(Path_Number,S);
SD_SeedbitErrorRate       = zeros(Path_Number,S);
CONV_bitErrorRate         = zeros(Path_Number,S);
CONV_packetErrorRate      = zeros(Path_Number,S);
HRSX_bitErrorRate         = zeros(Path_Number,S);
HRSX_packetErrorRate      = zeros(Path_Number,S);
SRSX_bitErrorRate         = zeros(Path_Number,S);
SRSX_packetErrorRate      = zeros(Path_Number,S);

SSIC_CONV_bitErrorRate    = zeros(S,1);
SSIC_CONV_packetErrorRate = zeros(S,1);
SSIC_HRSX_bitErrorRate    = zeros(S,1);
SSIC_HRSX_packetErrorRate = zeros(S,1);
SSIC_SRSX_bitErrorRate    = zeros(S,1);
SSIC_SRSX_packetErrorRate = zeros(S,1);

parfor i = 1:S % Use 'for' to debug the simulation
    % Set random substream index per iteration to ensure that each
    % iteration uses a repeatable set of random numbers
    stream = RandStream('combRecursive','Seed',0);
    stream.Substream = i;
    RandStream.setGlobalStream(stream);

    % Loop to simulate multiple packets
    CONV_SeednumBitErrors = zeros(Path_Number,1);
    SD_SeedbitError       = zeros(Path_Number,1);

    CONV_numBitErrors     = zeros(Path_Number,1);
    CONV_numPacketErrors  = zeros(Path_Number,1);
    HRSX_numBitErrors     = zeros(Path_Number,1);
    HRSX_numPacketErrors  = zeros(Path_Number,1);
    SRSX_numBitErrors     = zeros(Path_Number,1);
    SRSX_numPacketErrors  = zeros(Path_Number,1);

    SSIC_CONV_numBitErrors    = 0;
    SSIC_CONV_numPacketErrors = 0;
    SSIC_HRSX_numBitErrors    = 0;
    SSIC_HRSX_numPacketErrors = 0;
    SSIC_SRSX_numBitErrors    = 0;
    SSIC_SRSX_numPacketErrors = 0;

    TX_PACKETS = zeros(Path_Number,1); % Index of packet transmitted
    numPacketErrors = zeros(Path_Number,1); 

    SSIC_CONV_TX_PACKETS = 0; % Index of packet transmitted
    SSIC_HRSX_TX_PACKETS = 0; % Index of packet transmitted
    SSIC_SRSX_TX_PACKETS = 0; % Index of packet transmitted

    while numPacketErrors(1)<=maxNumPEs && TX_PACKETS(1)<=maxNumPackets

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

        CONV_xLLR_total = 0;
        HRSX_xLLR_total = 0;
        SRSX_xLLR_total = 0;

        for path_index = 1:Path_Number

            scramInit = randi([1,127]);

            tx = wlanWaveformGenerator(txPSDU,cfgHT,'ScramblerInitialization',scramInit);
            rx = awgn(tx,SNR_vector(path_index,i),'measured');

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
                numPacketErrors(path_index) = numPacketErrors(path_index)+1;
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
                numPacketErrors(path_index) = numPacketErrors(path_index)+1;
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
            [rxPSDU,Rx_scramInitBits,I] = wlanHTDescramer(decBits,cfgHT.PSDULength);

            conv_xLLR = zeros(size(rxPSDU_scrambled));
            for d = 1:length(conv_xLLR)
                d_cycle = mod(d-1,127)+1;
                if I(d_cycle) == 0
                    conv_xLLR(d) = rxPSDU_scrambled(d);
                else
                    conv_xLLR(d) = -rxPSDU_scrambled(d);
                end
            end
            CONV_xLLR_total = CONV_xLLR_total + conv_xLLR(L+1:end);

            % Determine if any bits are in error, i.e. a packet error
            CONV_SeedbitError = biterr(Tx_scramInitBits,Rx_scramInitBits);
            CONV_bitError     = biterr(validate_txPSDU,rxPSDU(L-16+1:end));
            CONV_packetError  = any(CONV_bitError);

            CONV_SeednumBitErrors(path_index) = CONV_SeednumBitErrors(path_index)+CONV_SeedbitError;
            CONV_numBitErrors(path_index)     = CONV_numBitErrors(path_index)+CONV_bitError;
            CONV_numPacketErrors(path_index)  = CONV_numPacketErrors(path_index)+CONV_packetError;
            %% Our Method
            yllr = double(rxPSDU_scrambled); % resevese it for our method
            l_total = yllr(1:L); % extract the first 127 LLR
            % derive f(r_1,r_2,...,r_7)
            [f_total,r_total] = DeriveF(a_total,l_total,K,L);

            %% hard-descrambling soft-decoding (HDSC)
            [hard_llr,~,Rx_scramInitBits,txBits,errorBits,pass] = Hard_NonHTDescramer_llr(cfgHT.PSDULength,f_total,r_total,yllr,L,txPSDU);
            SD_SeedbitError(path_index) = SD_SeedbitError(path_index) + biterr(Tx_scramInitBits,Rx_scramInitBits);
            if pass
                if errorBits~=0
                    HRSX_numBitErrors(path_index)    = HRSX_numBitErrors(path_index) + errorBits;
                    HRSX_numPacketErrors(path_index) = HRSX_numPacketErrors(path_index) + 1;
                end
                HRSX_xLLR_total = HRSX_xLLR_total + hard_llr;
            end
            %% Soft-descrambling soft-decoding (SDSC)
            [soft_llr,~,errorBits,pass] = Soft_NonHTDescramer_llr(cfgHT.PSDULength,a_total,f_total,r_total,yllr,L,txPSDU);
            if pass
                if errorBits~=0
                    SRSX_numBitErrors(path_index)    = SRSX_numBitErrors(path_index) + errorBits;
                    SRSX_numPacketErrors(path_index) = SRSX_numPacketErrors(path_index) + 1;
                end
                SRSX_xLLR_total = SRSX_xLLR_total + soft_llr;
            end
            %%
            TX_PACKETS(path_index) = TX_PACKETS(path_index)+1;
        end

        %% SSIC 
        if (sum(CONV_xLLR_total) ~= 0)
            SSIC_CONV_TX_PACKETS = SSIC_CONV_TX_PACKETS+1;

            recBits = (CONV_xLLR_total < 0);
            errorBits = biterr(validate_txPSDU,double(recBits));
            if (errorBits ~=0)
                SSIC_CONV_numBitErrors    = SSIC_CONV_numBitErrors + errorBits;
                SSIC_CONV_numPacketErrors = SSIC_CONV_numPacketErrors + 1;
            end
        end

        if (sum(HRSX_xLLR_total) ~= 0)
            SSIC_HRSX_TX_PACKETS = SSIC_HRSX_TX_PACKETS+1;

            recBits = (HRSX_xLLR_total<0);
            errorBits = biterr(validate_txPSDU,double(recBits));
            if (errorBits~=0)
               SSIC_HRSX_numBitErrors    = SSIC_HRSX_numBitErrors + errorBits;
               SSIC_HRSX_numPacketErrors = SSIC_HRSX_numPacketErrors + 1;
            end
        end

        if (sum(SRSX_xLLR_total) ~= 0)
            SSIC_SRSX_TX_PACKETS = SSIC_SRSX_TX_PACKETS+1;

            recBits = (SRSX_xLLR_total<0);
            errorBits = biterr(validate_txPSDU,double(recBits));
            if (errorBits~=0)
               SSIC_SRSX_numBitErrors    = SSIC_SRSX_numBitErrors + errorBits;
               SSIC_SRSX_numPacketErrors = SSIC_SRSX_numPacketErrors + 1;
            end
        end
    end

    %% Seed
    % Calculate Seed's bit error rate (BER) at SNR point
    CONV_SeedbitErrorRate(:,i) = CONV_SeednumBitErrors./(TX_PACKETS*7);
    % Calculate Seed's bit error rate (BER) at SNR point
    SD_SeedbitErrorRate(:,i) = SD_SeedbitError./(TX_PACKETS*7);
    %% CONV
    % Calculate bit error rate (BER) at SNR point
    CONV_bitErrorRate(:,i) = CONV_numBitErrors./(TX_PACKETS*length(validate_txPSDU));
    % Calculate packet error rate (PER) at SNR point
    CONV_packetErrorRate(:,i) = CONV_numPacketErrors./TX_PACKETS;
    %% HRSX
    % Calculate bit error rate (BER) at SNR point
    HRSX_bitErrorRate(:,i) = HRSX_numBitErrors./(TX_PACKETS*length(validate_txPSDU));
    % Calculate packet error rate (PER) at SNR point
    HRSX_packetErrorRate(:,i) = HRSX_numPacketErrors./TX_PACKETS;
    %% SRSX
    % Calculate bit error rate (BER) at SNR point
    SRSX_bitErrorRate(:,i) = SRSX_numBitErrors./(TX_PACKETS*length(validate_txPSDU));
    % Calculate packet error rate (PER) at SNR point
    SRSX_packetErrorRate(:,i) = SRSX_numPacketErrors./TX_PACKETS;

    %% SSIC CONV
    % Calculate bit error rate (BER) at SNR point
    SSIC_CONV_bitErrorRate(i) = SSIC_CONV_numBitErrors/(SSIC_CONV_TX_PACKETS*length(validate_txPSDU));
    % Calculate packet error rate (PER) at SNR point
    SSIC_CONV_packetErrorRate(i) = SSIC_CONV_numPacketErrors/SSIC_CONV_TX_PACKETS;
    %% SSIC HRSX
    % Calculate bit error rate (BER) at SNR point
    SSIC_HRSX_bitErrorRate(i) = SSIC_HRSX_numBitErrors/(SSIC_HRSX_TX_PACKETS*length(validate_txPSDU));
    % Calculate packet error rate (PER) at SNR point
    SSIC_HRSX_packetErrorRate(i) = SSIC_HRSX_numPacketErrors/SSIC_HRSX_TX_PACKETS;
    %% SSIC SRSX
    % Calculate bit error rate (BER) at SNR point
    SSIC_SRSX_bitErrorRate(i) = SSIC_SRSX_numBitErrors/(SSIC_SRSX_TX_PACKETS*length(validate_txPSDU));
    % Calculate packet error rate (PER) at SNR point
    SSIC_SRSX_packetErrorRate(i) = SSIC_SRSX_numPacketErrors/SSIC_SRSX_TX_PACKETS;
end

%% Plot Bit Error Rate vs SNR Results
figure;
for path_index = 1:Path_Number
    semilogy(basic_range,CONV_SeedbitErrorRate(path_index,:),'--');
    hold on;
    semilogy(basic_range,SD_SeedbitErrorRate(path_index,:),'-o');
end
grid on;
xlabel('SNR [dB]');
ylabel('Seed BER');
legend('CONV 1','SD 1','CONV 2','SD 2','CONV 3','SD 3','CONV 4','SD 4')

%% Plot Bit Error Rate vs SNR Results
figure;
for path_index = 1:Path_Number
    semilogy(basic_range,CONV_bitErrorRate(path_index,:),'--');
    hold on;
    semilogy(basic_range,HRSX_bitErrorRate(path_index,:),'-o');
    hold on;
    semilogy(basic_range,SRSX_bitErrorRate(path_index,:),'-*');
end

hold on;
semilogy(basic_range,SSIC_CONV_bitErrorRate,'--');
hold on;
semilogy(basic_range,SSIC_HRSX_bitErrorRate,'-o');
hold on;
semilogy(basic_range,SSIC_SRSX_bitErrorRate,'-*');

grid on;
xlabel('SNR [dB]');
ylabel('BER');
legend('CONV 1','HRSX 1','SRSX 1','CONV 2','HRSX 2','SRSX 2',...
    'CONV 3','HRSX 3','SRSX 3','CONV 4','HRSX 4','SRSX 4','SSIC Naive','SSIC HRSX','SSIC SRSX')

%% Plot Packet Error Rate vs SNR Results
figure;
for path_index = 1:Path_Number
    semilogy(basic_range,CONV_packetErrorRate(path_index,:),'--');
    hold on;
    semilogy(basic_range,HRSX_packetErrorRate(path_index,:),'-o');
    hold on;
    semilogy(basic_range,SRSX_packetErrorRate(path_index,:),'-*');
end

hold on;
semilogy(basic_range,SSIC_CONV_packetErrorRate,'--');
hold on;
semilogy(basic_range,SSIC_HRSX_packetErrorRate,'-o');
hold on;
semilogy(basic_range,SSIC_SRSX_packetErrorRate,'-*');

grid on;
xlabel('SNR [dB]');
ylabel('PER');
legend('CONV 1','HRSX 1','SRSX 1','CONV 2','HRSX 2','SRSX 2',...
    'CONV 3','HRSX 3','SRSX 3','CONV 4','HRSX 4','SRSX 4','SSIC Naive','SSIC HRSX','SSIC SRSX')

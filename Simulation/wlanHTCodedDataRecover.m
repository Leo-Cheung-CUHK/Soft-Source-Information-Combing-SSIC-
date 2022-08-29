function [descramIn, eqDataSym, varargout] = wlanHTCodedDataRecover( ...
    rx, chanEst, noiseVarEst, cfgHT, varargin)
%wlanHTDataRecover Recover information bits from HT-Data field signal
%
%   BITS = wlanHTDataRecover(RX, CHANEST, NOISEVAREST, CFGHT) recovers the
%   information bits in the HT-Data field for a HT-Mixed format
%   transmission.
%
%   BITS is an int8 column vector of length 8*CFGHT.PSDULength containing
%   the recovered information bits.
%
%   RX is the received time-domain HT-Data field signal. It is a Ns-by-Nr
%   matrix of real or complex values, where Ns represents the number of
%   time-domain samples in the HT-Data field and Nr represents the number
%   of receive antennas. Ns can be greater than the HT-Data field length;
%   in this case additional samples at the end of RX are not used.
%
%   CHANEST is the estimated channel at data and pilot subcarriers based on
%   the HT-LTF. It is a real or complex array of size Nst-by-Nsts-by-Nr,
%   where Nst represents the total number of occupied subcarriers.
%
%   NOISEVAREST is the noise variance estimate. It is a real, nonnegative
%   scalar.
%
%   CFGHT is the format configuration object of type <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, which
%   specifies the parameters for the HT-Mixed format.
%
%   BITS = wlanHTDataRecover(..., NAME, VALUE) specifies additional
%   name-value pair arguments described below. When a name-value pair is
%   not specified, its default value is used.
%
%   'OFDMSymbolOffset'          OFDM symbol sampling offset. Specify the
%                               OFDMSymbolOffset as a fraction of the
%                               cyclic prefix (CP) length for every OFDM
%                               symbol, as a double precision, real scalar
%                               between 0 and 1, inclusive. The OFDM
%                               demodulation is performed based on Nfft
%                               samples following the offset position,
%                               where Nfft denotes the FFT length. The
%                               default value of this property is 0.75,
%                               which means the offset is three quarters of
%                               the CP length.
%
%   'EqualizationMethod'        Specify the equalization method as one of
%                               'MMSE' | 'ZF'. 'MMSE' indicates that the
%                               receiver uses a minimum mean square error
%                               equalizer. 'ZF' indicates that the receiver
%                               uses a zero-forcing equalizer. The default
%                               value of this property is 'MMSE'.
%
%   'PilotPhaseTracking'        Specify the pilot phase tracking performed
%                               as one of 'PreEQ' | 'None'. 'PreEQ' pilot
%                               phase tracking estimates and corrects a
%                               common phase offset across all subcarriers
%                               and receive antennas for each received OFDM
%                               symbol before equalization. 'None'
%                               indicates that pilot phase tracking does
%                               not occur. The default is 'PreEQ'.
%
%   'LDPCDecodingMethod'        Specify the LDPC decoding algorithm as one
%                               of these values:
%                               - 'bp'            : Belief propagation (BP)
%                               - 'layered-bp'    : Layered BP
%                               - 'norm-min-sum'  : Normalized min-sum
%                               - 'offset-min-sum': Offset min-sum
%                               The default is 'bp'.
%
%   'MinSumScalingFactor'       Specify the scaling factor for normalized
%                               min-sum LDPC decoding algorithm as a scalar
%                               in the interval (0,1]. This argument
%                               applies only when you set
%                               LDPCDecodingMethod to 'norm-min-sum'. The
%                               default is 0.75.
%
%   'MinSumOffset'              Specify the offset for offset min-sum LDPC
%                               decoding algorithm as a finite real scalar
%                               greater than or equal to 0. This argument
%                               applies only when you set
%                               LDPCDecodingMethod to 'offset-min-sum'. The
%                               default is 0.5.
%
%   'MaximumLDPCIterationCount' Specify the maximum number of iterations in
%                               LDPC decoding as a positive scalar integer.
%                               This applies when you set the channel
%                               coding property of format configuration
%                               object of type <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a> to 'LDPC'.
%                               The default is 12.
%
%   'EarlyTermination'          To enable early termination of LDPC
%                               decoding, set this property to true. Early
%                               termination applies if all parity-checks
%                               are satisfied before reaching the number of
%                               iterations specified in the
%                               'MaximumLDPCIterationCount' input. To let
%                               the decoding process iterate for the number
%                               of iterations specified in the
%                               'MaximumLDPCIterationCount' input, set this
%                               argument to false. This applies when you
%                               set the channel coding property of format
%                               configuration object of type <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a>
%                               to 'LDPC'.The default is false.
%
%   [..., EQDATASYM, CPE] = wlanHTDataRecover(...) also returns the
%   equalized subcarriers and common phase error.
%
%   EQDATASYM is a complex Nsd-by-Nsym-by-Nss array containing the
%   equalized symbols at data carrying subcarriers. Nsd represents the
%   number of data subcarriers, Nsym represents the number of OFDM symbols
%   in the HT-Data field, and Nss represents the number of spatial streams.
%
%   CPE is a column vector of length Nsym containing the common phase error
%   between each received and expected OFDM symbol.
%
%   Example:
%   %  Recover a HT-Data field signal through a SISO AWGN channel using
%   %  ZF equalization.
%
%     cfgHT = wlanHTConfig('PSDULength', 1024);     % HT format configuration
%     txBits = randi([0 1], 8*cfgHT.PSDULength, 1); % Payload bits
%     txHSig = wlanHTData(txBits, cfgHT);           % Generate HT-Data signal
%
%     % Pass through an AWGN channel with noise variance of 1
%     rxHTSig = awgn(txHSig, 1, 1);
%
%     % Recover payload bits
%     rxBits = wlanHTDataRecover(rxHTSig, ones(56,1), 1, cfgHT, 'EqualizationMethod', 'ZF');
%
%     [numerr, ber] = biterr(rxBits, txBits);       % Compare bits
%     disp(ber)
%
%   See also wlanHTConfig, wlanHTData.

%   Copyright 2015-2021 The MathWorks, Inc.

%#codegen

narginchk(4, 20);
nargoutchk(0, 3);

% Calculate CPE if requested
if nargout>2
    calculateCPE = true;
else
    calculateCPE = false;
end

% HT configuration input self-validation
validateattributes(cfgHT, {'wlanHTConfig'}, {'scalar'}, mfilename, 'HT-Mixed format configuration object');
s = validateConfig(cfgHT, 'MCS');

% Validate rxHTData
validateattributes(rx, {'double'}, {'2d','finite'}, 'rxHTData', 'HT-Data field signal');
% Validate chanEst
validateattributes(chanEst, {'double'}, {'3d','finite'}, 'chanEst', 'channel estimates');
% Validate noiseVarEst
validateattributes(noiseVarEst, {'double'}, {'real','scalar','nonnegative','finite'}, 'noiseVarEst', 'noise variance estimate');

numSTS     = cfgHT.NumSpaceTimeStreams;
numRx      = size(rx, 2);
mcsTable   = wlan.internal.getRateTable(cfgHT);
numSS      = mcsTable.Nss;
numCBPSSI  = mcsTable.NCBPS/numSS;
numDBPS    = mcsTable.NDBPS;
rate       = mcsTable.Rate;
numOFDMSym = s.NumDataSymbols;
STBC       = numSTS - numSS;
mSTBC      = 1 + (STBC~=0);

% If PSDU is empty there is no data to return
if cfgHT.PSDULength == 0
    bits     = zeros(0, 1, 'int8');
    eqDataSym  = zeros(mcsTable.NSD, 0, numSS);
    if calculateCPE==true
        varargout{1} = []; % CPE
    end
    return;
end

% Validate and parse optional inputs
recParams = wlan.internal.parseOptionalInputs(mfilename, varargin{:});

% Get OFDM configuration
ofdmInfo = wlan.internal.vhtOFDMInfo('HT-Data', cfgHT.ChannelBandwidth, cfgHT.GuardInterval);

% Cross validate input
coder.internal.errorIf(size(chanEst, 1) ~= ofdmInfo.NumTones, 'wlan:wlanHTDataRecover:InvalidHTChanEst1D', ofdmInfo.NumTones);
coder.internal.errorIf(size(chanEst, 2) ~= numSTS, 'wlan:wlanHTDataRecover:InvalidHTChanEst2D', numSTS);
coder.internal.errorIf(size(chanEst, 3) ~= numRx, 'wlan:wlanHTDataRecover:InvalidHTChanEst3D');

% Cross-validation between inputs
minInputLen = numOFDMSym*(ofdmInfo.FFTLength+ofdmInfo.CPLength);
coder.internal.errorIf(size(rx, 1) < minInputLen, 'wlan:wlanHTDataRecover:ShortHTDataInput', minInputLen);

% OFDM demodulation with de-normalization and removing phase rotation per subcarrier
demod = wlan.internal.legacyOFDMDemodulate(rx(1:minInputLen,:), ofdmInfo, recParams.OFDMSymbolOffset, numSTS);

if calculateCPE==true || strcmp(recParams.PilotPhaseTracking, 'PreEQ')
    % Get reference pilots, from IEEE Std 802.11-2012, Eqn 20-58/59
    % For HT-MF, offset by 3 to allow for L-SIG and HT-SIG pilot symbols
    z = 3;
    refPilots = wlan.internal.htPilots(numOFDMSym, z, cfgHT.ChannelBandwidth, numSTS);

    % Estimate CPE and phase correct symbols
    cpe = wlan.internal.commonPhaseErrorEstimate(demod(ofdmInfo.PilotIndices,:,:), chanEst(ofdmInfo.PilotIndices,:,:), refPilots);
    if strcmp(recParams.PilotPhaseTracking, 'PreEQ')
        demod(ofdmInfo.DataIndices,:,:) = wlan.internal.commonPhaseErrorCorrect(demod(ofdmInfo.DataIndices,:,:), cpe);
    end
    if calculateCPE==true
        varargout{1} = cpe.'; % Permute to Nsym-by-1
    end
end

% Equalization
if numSS < numSTS
    [eqDataSym, csiData] = wlan.internal.wlanSTBCCombine(demod(ofdmInfo.DataIndices,:,:), chanEst(ofdmInfo.DataIndices,:,:), numSS, recParams.EqualizationMethod, noiseVarEst);
else
    [eqDataSym, csiData] = wlan.internal.wlanEqualize(demod(ofdmInfo.DataIndices,:,:), chanEst(ofdmInfo.DataIndices,:,:), recParams.EqualizationMethod, noiseVarEst);
end

% Constellation demapping
qamDemodOut = wlanConstellationDemap(eqDataSym, noiseVarEst, mcsTable.NBPSCS);

% Apply bit-wise CSI and concatenate OFDM symbols in the first dimension
qamDemodOut = bsxfun(@times, ...
                     reshape(qamDemodOut, mcsTable.NBPSCS, [], numOFDMSym, numSS), ...
                     reshape(csiData, 1, [], 1, numSS)); % [Nbpscs Nsd Nsym Nss]
qamDemodOut = reshape(qamDemodOut, [], numSS);

% BCC Deinterleaving
if strcmp(cfgHT.ChannelCoding,'BCC')
    deintlvrOut = wlanBCCDeinterleave(qamDemodOut, 'VHT', numCBPSSI, cfgHT.ChannelBandwidth);
else
    % Deinterleaving is not required for LDPC
    deintlvrOut = qamDemodOut;
end

% Stream deparsing
streamDeparserOut = wlanStreamDeparse(deintlvrOut, mcsTable.NES, mcsTable.NCBPS, mcsTable.NBPSCS);

% Channel decoding
if strcmp(cfgHT.ChannelCoding,'BCC')
    % BCC channel decoding
    htDataBits = wlanBCCDecode(streamDeparserOut, rate);
    % BCC decoder deparser
    descramIn = reshape(htDataBits.', [], 1);
else
    % LDPC Channel decoding
    numPLD = cfgHT.PSDULength*8 + 16; % Number of payload bits
                                      % LDPC decoding parameters, IEEE Std 802.11-2012, Section
                                      % 20.3.11.17.5.
    cfg = wlan.internal.getLDPCparameters(numDBPS, rate, mSTBC, numPLD);
    descramIn = wlanLDPCDecode(streamDeparserOut(:), cfg, ...
                                             recParams.algChoice, recParams.alphaBeta, recParams.MaximumLDPCIterationCount, recParams.EarlyTermination);
end

end

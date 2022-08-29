function [bits, eqDataSym, varargout] = wlanNonHTDataRecover( ...
    rxNonHTData, chanEst, noiseVarEst, cfgNonHT, varargin)
%wlanNonHTDataRecover Recover information bits from non-HT Data field signal
%
%   BITS = wlanNonHTDataRecover(RXNONHTDATA, CHANEST, NOISEVAREST,
%   CFGNONHT) recovers the information bits in the non-HT Data field for a
%   non-HT OFDM format transmission.
%
%   BITS is an int8 column vector of length 8*CFGNONHT.PSDULength
%   containing the recovered information bits.
%
%   RXNONHTDATA is the received time-domain non-HT Data field signal. It is
%   a Ns-by-Nr matrix of real or complex values, where Ns represents the
%   number of time-domain samples in the non-HT Data field and Nr represents
%   the number of receive antennas. Ns can be greater than the non-HT Data
%   field length; in this case redundant samples at the end of RXNONHTDATA
%   are not used.
%
%   CHANEST is the estimated channel at data and pilot subcarriers based on
%   the L-LTF. It is a real or complex array of size Nst-by-1-by-Nr, where
%   Nst represents the total number of occupied subcarriers. The singleton
%   dimension corresponds to the single transmitted stream in the L-LTF
%   which includes the combined cyclic shifts if multiple transmit antennas
%   are used.
%
%   NOISEVAREST is the noise variance estimate. It is a real, nonnegative
%   scalar.
%
%   CFGNONHT is the format configuration object of type <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a> 
%   that specifies the non-HT format parameters. Only OFDM modulation is
%   supported.
%
%   BITS = wlanNonHTDataRecover(..., NAME, VALUE) specifies additional
%   name-value pair arguments described below. When a name-value pair is
%   not specified, its default value is used.
%
%   'OFDMSymbolOffset'      OFDM symbol sampling offset. Specify the
%                           OFDMSymbolOffset as a fraction of the cyclic
%                           prefix (CP) length for every OFDM symbol, as a
%                           double precision, real scalar between 0 and 1,
%                           inclusive. The OFDM demodulation is performed
%                           based on Nfft samples following the offset
%                           position, where Nfft denotes the FFT length.
%                           The default value of this property is 0.75,
%                           which means the offset is three quarters of the
%                           CP length.
%
%   'EqualizationMethod'    Specify the equalization method as one of
%                           'MMSE' | 'ZF'. 'MMSE' indicates that the
%                           receiver uses a minimum mean square error
%                           equalizer. 'ZF' indicates that the receiver
%                           uses a zero-forcing equalizer. The default
%                           value of this property is 'MMSE'.
%
%   'PilotPhaseTracking'    Specify the pilot phase tracking performed as
%                           one of 'PreEQ' | 'None'. 'PreEQ' pilot phase
%                           tracking estimates and corrects a common phase
%                           offset across all subcarriers and receive
%                           antennas for each received OFDM symbol before
%                           equalization. 'None' indicates that pilot phase
%                           tracking does not occur. The default is 'PreEQ'.
%
%   [..., EQDATASYM, CPE, SCRAMINIT] = wlanNonHTDataRecover(...) also
%   returns the equalized subcarriers, common phase error, and recovered
%   scrambler initial state.
%
%   EQDATASYM is a complex 48-by-Nsym matrix containing the equalized
%   symbols at data carrying subcarriers. There are 48 data carrying
%   subcarriers in the non-HT Data field. Nsym represents the number of
%   OFDM symbols in the non-HT Data field.
%
%   CPE is a column vector of length Nsym containing the common phase error
%   between each received and expected OFDM symbol.
%
%   SCRAMINIT is an int8 scalar containing the recovered initial scrambler
%   state. The function maps the initial state bits X1 to X7, as specified
%   in IEEE 802.11-2016, Section 17.3.5.5 to SCRAMINIT, treating the
%   rightmost bit as most significant.
%
%   Example: 
%   %  Recover a non-HT Data field signal through a SISO AWGN channel
%   %  using ZF equalization.
%
%     cfgNonHT = wlanNonHTConfig('PSDULength', 1024);  % non-HT OFDM 
%     txBits = randi([0 1], 8*cfgNonHT.PSDULength, 1); % PSDU bits
%     tx = wlanNonHTData(txBits, cfgNonHT);       % non-HT Data field signal
%
%     % Add AWGN, with noise variance of 1
%     rx = awgn(tx, 1, 1);
%
%     % Recover PSDU bits 
%     rxBits = wlanNonHTDataRecover(rx, ones(52,1), 1, cfgNonHT, 'EqualizationMethod', 'ZF');
%
%     [numerr, ber] = biterr(rxBits, txBits); % Compare bits
%     disp(ber)
%
%   See also wlanNonHTConfig, wlanLLTFChannelEstimate, wlanNonHTData. 

%   Copyright 2015-2020 The MathWorks, Inc.

%#codegen

narginchk(4, 10);
nargoutchk(0, 4);

% Calculate CPE if requested
if nargout>2
    calculateCPE = true;
else
    calculateCPE = false;
end

% Non-HT configuration input self-validation
validateattributes(cfgNonHT, {'wlanNonHTConfig'}, {'scalar'}, mfilename, 'format configuration object');
% Only applicable for OFDM and DUP-OFDM modulations
coder.internal.errorIf(~strcmp(cfgNonHT.Modulation, 'OFDM'), 'wlan:shared:InvalidModulation');
s = validateConfig(cfgNonHT);

% Validate rxNonHTData
validateattributes(rxNonHTData, {'double'}, {'2d','finite'}, 'rxHTData', 'Non-HT OFDM Data field signal'); 
% Validate chanEst
validateattributes(chanEst, {'double'}, {'3d','finite'}, 'chanEst', 'channel estimates'); 
% Validate noiseVarEst
validateattributes(noiseVarEst, {'double'}, {'real','scalar','nonnegative','finite'}, 'noiseVarEst', 'noise variance estimate'); 

% Validate and parse optional inputs
recParams = wlan.internal.parseOptionalInputs(mfilename, varargin{:});

numRx = size(rxNonHTData, 2);

numOFDMSym = s.NumDataSymbols;

% Get OFDM configuration
[cfgOFDM,dataInd,pilotInd] = wlan.internal.wlanGetOFDMConfig(cfgNonHT.ChannelBandwidth, 'Long', 'Legacy');

% Cross validate inputs
numST = numel([dataInd; pilotInd]); % Total number of occupied subcarriers
coder.internal.errorIf(size(chanEst, 1) ~= numST, 'wlan:wlanNonHTDataRecover:InvalidNHTChanEst1D', numST);
coder.internal.errorIf(size(chanEst, 2) ~= 1, 'wlan:wlanNonHTDataRecover:InvalidNHTChanEst2D');
coder.internal.errorIf(size(chanEst, 3) ~= numRx, 'wlan:wlanNonHTDataRecover:InvalidNHTChanEst3D');

% Extract data and pilot subcarriers from channel estimate
chanEstData = chanEst(dataInd,:,:);
chanEstPilots = chanEst(pilotInd,:,:);

% Cross-validation between inputs
minInputLen = numOFDMSym*(cfgOFDM.FFTLength+cfgOFDM.CyclicPrefixLength);
coder.internal.errorIf(size(rxNonHTData, 1) < minInputLen, 'wlan:wlanNonHTDataRecover:ShortNHTDataInput', minInputLen);

% Processing 
% OFDM Demodulation 
[ofdmDemodData, ofdmDemodPilots] = wlan.internal.wlanOFDMDemodulate(rxNonHTData(1:minInputLen, :), cfgOFDM, recParams.OFDMSymbolOffset);

% Pilot phase tracking
if calculateCPE==true || strcmp(recParams.PilotPhaseTracking, 'PreEQ')
    % Get reference pilots, from IEEE Std 802.11-2012, Eqn 18-22
    z = 1; % Offset by 1 to account for L-SIG pilot symbol
    refPilots = wlan.internal.nonHTPilots(numOFDMSym, z, cfgNonHT.ChannelBandwidth);
    
    % Estimate CPE and phase correct symbols
    cpe = wlan.internal.commonPhaseErrorEstimate(ofdmDemodPilots, chanEstPilots, refPilots);
    if strcmp(recParams.PilotPhaseTracking, 'PreEQ')
        ofdmDemodData = wlan.internal.commonPhaseErrorCorrect(ofdmDemodData, cpe);
    end
    if calculateCPE==true
        varargout{1} = cpe.'; % Permute to Nsym-by-1
    end
end

% Merge num20 channel estimates and demodulated symbols together for the
% repeated subcarriers for data carrying subcarriers
NsdSeg = 48; % Number of subcarriers in 20 MHz segment
num20MHz = size(ofdmDemodData,1)/NsdSeg; % Number of 20 MHz subchannels
ofdmDataOutOne20MHz = coder.nullcopy(complex(zeros(NsdSeg, numOFDMSym, numRx*num20MHz))); % Preallocate for codegen
chanEstDataOne20MHz = coder.nullcopy(complex(zeros(NsdSeg, 1, numRx*num20MHz))); % Preallocate for codegen
[ofdmDataOutOne20MHz(:), chanEstDataOne20MHz(:)] = wlan.internal.mergeSubchannels(ofdmDemodData, chanEstData, num20MHz);

% Equalization
[eqDataSym, csiData] = wlan.internal.wlanEqualize(ofdmDataOutOne20MHz, chanEstDataOne20MHz, recParams.EqualizationMethod, noiseVarEst);

% Demap and decode
[bits,~] = wlanNonHTDataBitRecover(eqDataSym,noiseVarEst,csiData,cfgNonHT);

%varargout{2} = scraminit; % disable this 

end
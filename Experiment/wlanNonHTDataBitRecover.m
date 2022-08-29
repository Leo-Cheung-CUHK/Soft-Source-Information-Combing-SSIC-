function [psdu, scramInit] = wlanNonHTDataBitRecover(sym, noiseVarEst, varargin)
%wlanNonHTDataBitRecover Recover data bits from non-HT Data field
%   PSDU = wlanNonHTDataBitRecover(RX,NOISEVAREST,CFG) recovers the data
%   bits given the equalized Data field from a non-HT transmission, the
%   noise variance estimate, and the non-HT configuration object.
%
%   PSDU is an int8 column vector of length 8*CFG.PSDULength containing the
%   recovered information bits.
%
%   RX contains the demodulated and equalized Data field OFDM symbols,
%   specified as a 48-by-Nsym complex-valued matrix, where 48 is the number
%   of data subcarriers in the Data field and Nsym is the number of OFDM
%   symbols.
%
%   NOISEVAREST is the noise variance estimate, specified as a nonnegative
%   scalar.
%
%   CFG is the format configuration object of type <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a> which
%   specifies the parameters for the non-HT format. Only OFDM modulation
%   type is supported.
%
%   DATABITS = wlanNonHTDataBitRecover(...,CSI,CFG) uses the channel state
%   information to enhance the demapping of OFDM subcarriers. CSI is a
%   48-by-1 column vector of real values.
%
%   [...,SCRAMINIT] = wlanNonHTDataBitRecover(...) additionally returns the
%   recovered initial scrambler state as an int8 scalar. The function maps
%   the initial state bits X1 to X7 as specified in IEEE Std 802.11-2016,
%   Section 17.3.5.5 to SCRAMINIT, treating the leftmost bit as most
%   significant.
%
%   Example:
%   % Recover PSDU for a non-HT packet in AWGN.
%
%      cfg = wlanNonHTConfig('MCS',4);
%      
%      % Generate transmit waveform
%      txPSDU = randi([0 1],cfg.PSDULength*8,1,'int8');
%      txWaveform = wlanWaveformGenerator(txPSDU,cfg);
%      
%      % Add noise
%      snr = 30;
%      rxWaveform = awgn(txWaveform,snr);
%      
%      % Extract the data field
%      ind = wlanFieldIndices(cfg);
%      rxData = rxWaveform(ind.NonHTData(1): ind.NonHTData(2),:);
%      
%      % OFDM demodulate
%      demod = wlanNonHTOFDMDemodulate(rxData,'NonHT-Data',cfg);
%      
%      % Extract data subcarriers
%      info = wlanNonHTOFDMInfo('NonHT-Data',cfg);
%      rxSym = demod(info.DataIndices,:,:);
%      
%      % Recover data bits
%      csi = ones(size(rxSym,1),1); % Assume CSI estimate of all ones
%      nVar = 10^(-snr/10); % Noise variance
%      rxPSDU = wlanNonHTDataBitRecover(rxSym,nVar,csi,cfg);
%      
%      % Compare against original information bits
%      disp(isequal(txPSDU,rxPSDU));
%
%   See also wlanNonHTConfig, wlanNonHTOFDMDemodulate, wlanNonHTData. 

%   Copyright 2015-2020 The MathWorks, Inc.

%#codegen

narginchk(3,4);

[nsd,nsym,nss] = size(sym);
validateattributes(sym, {'double'}, {'finite','2d','nrows',48}, mfilename, 'SYM');

if isa(varargin{1},'wlanNonHTConfig')
    % wlanNonHTDataBitRecover(RX,NOISEVAREST,CFG)
    % If no CSI input is present then assume 1 for processing
    csi = ones(nsd,nss);
    cfg = varargin{1};
elseif nargin>3 && isa(varargin{2}, 'wlanNonHTConfig')
    % wlanNonHTDataBitRecover(RX,NOISEVAREST,CSI,CFG)
    csi = varargin{1};
    cfg = varargin{2};
    validateattributes(csi, {'double'}, {'real','finite','size',[48,1]}, mfilename, 'CSI');
else
    coder.internal.error('wlan:wlanNonHTDataBitRecover:InvalidSyntax');
end

% Non-HT configuration input self-validation
validateattributes(cfg, {'wlanNonHTConfig'}, {'scalar'}, mfilename, 'format configuration object');
% Only applicable for OFDM and DUP-OFDM modulations
coder.internal.errorIf(~strcmp(cfg.Modulation, 'OFDM'), 'wlan:shared:InvalidModulation');
s = validateConfig(cfg);
coder.internal.errorIf(nsym<s.NumDataSymbols, 'wlan:shared:IncorrectNumOFDMSym', s.NumDataSymbols, nsym);

% Validate noiseVarEst
validateattributes(noiseVarEst, {'double'}, {'real','scalar','nonnegative','finite'}, 'noiseVarEst', 'noise variance estimate'); 

mcsTable = wlan.internal.getRateTable(cfg);

% Constellation demapping
qamDemodOut = wlanConstellationDemap(sym, noiseVarEst, mcsTable.NBPSCS);

% Apply bit-wise CSI and concatenate OFDM symbols in the first dimension
qamDemodOut = bsxfun(@times, ...
              reshape(qamDemodOut, mcsTable.NBPSCS, [], nsym), ...
              reshape(csi, 1, [])); % [Nbpscs Nsd Nsym]
qamDemodOut = reshape(qamDemodOut, [], 1);

% Deinterleave
deintlvrOut = wlanBCCDeinterleave(qamDemodOut, 'Non-HT', mcsTable.NCBPS);

% Channel decoding
decBits = wlanBCCDecode(deintlvrOut, mcsTable.Rate);

% Derive initial state of the scrambler
scramSeqInit = decBits(1:7);
scramInitBits = wlan.internal.scramblerInitialState(scramSeqInit);

% Remove pad and tail bits, and descramble
if all(scramInitBits==0)
    % Scrambler initialization invalid (0), therefore do not descramble
    descramDataOut = decBits(1:(16+8*cfg.PSDULength));
else
    descramDataOut = wlanScramble(decBits(1:(16+8*cfg.PSDULength)), scramInitBits);
end

% Remove the 16 service bits
psdu = descramDataOut(17:end);

% Convert scrambler initialization bits to number
scramInit = bi2de(scramInitBits.', 'left-msb');
end
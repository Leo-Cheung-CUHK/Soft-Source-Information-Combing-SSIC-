function [y,numIterations,parityCheck] = wlanLDPCDecode(x,cfg,algChoice,alphaBeta,maxNumIter,earlyTermination)
%wlanLDPCDecode Low-Density-Parity-Check (LDPC) decoder
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = wlanLDPCDecode(X,CFG,ALGCHOICE,ALPHABETA,MAXNUMITER,EARLYTERMINATION)
%   decodes the input X for the specified rate, LDPC algorithm, and
%   options. Output Y is a hard decision decoded output of the information
%   bits.
%
%   Each column of X specifies the log-likelihood ratios of a codeword. The
%   number of rows in X should be equal to one of the four valid block
%   lengths: 648, 672, 1296, and 1944.
%
%   CFG should be a structure including the fields:
%   VecPayloadBits  - Number of payload bits within a codeword
%   Rate            - Coding rate
%   NumCW           - Number of LDPC codewords
%   LengthLDPC      - LDPC codeword length
%   VecShortenBits  - Vector of shortening bits in each codeword
%   VecPunctureBits - Vector of puncture bits in each codeword
%   VecRepeatBits   - Number of coded bits to be repeated
%
%   ALGCHOICE specifies the LDPC decoding algorithm as one of these values:
%     0 - Belief Propagation
%     1 - Layered Belief Propagation
%     2 - Layered Belief Propagation with Normalized Min-Sum approximation
%     3 - Layered Belief Propagation with Offset Min-Sum approximation
%
%   ALPHABETA specifies the scaling factor for Normalized Min-Sum
%   approximation or the offset factor for Offset Min-Sum approximation.
%   Its value is irrelevant for the other two LDPC algorithms but still
%   needed.
%
%   MAXNUMITER specifies the number of decoding iterations required to
%   decode the input X.
%
%   EARLYTERMINATION specifies if the conditions for an early termination
%   should be checked. If true, after each iteration, wlanLDPCDecode will
%   determine independently for each column of X if all parity checks are
%   satisfied, and will stop for column(s) with all parity checks
%   satisfied; otherwise, the function will execute for the maximum number
%   of iterations MAXNUMITER.
%
%   [Y, NUMITERATIONS] = wlanLDPCDecode(...) decodes the LDPC encoded data.
%   The function returns the actual number of LDPC decoding iterations, one
%   per codeword. The NUMITERATIONS is NUMCW-by-1 vector, where NUMCW is
%   the number of codewords.
%
%   [Y, NUMITERATIONS, PARITYCHECK] = wlanLDPCDecode(...) decodes the LDPC
%   encoded data and returns the parity check per codeword. The PARITYCHECK
%   is NUMINP-by-NUMCW matrix, where NUMINP is the number of information
%   bits within a codeword and NUMCW is the number of codewords.
%
%   See also wlanLDPCEncode, getLDPCparameters. 

%   Copyright 2016-2019 The MathWorks, Inc.

%#codegen

    numCW           = cfg.NumCW;
    lengthLDPC      = cfg.LengthLDPC;
    vecShortenBits  = cfg.VecShortenBits;
    vecPunctureBits = cfg.VecPunctureBits;
    vecRepeatBits   = cfg.VecRepeatBits;
    vecPayloadBits  = cfg.VecPayloadBits;
    rate            = cfg.Rate;

    % Initialize output
    y = coder.nullcopy(zeros(sum(vecPayloadBits),1,'int8'));
    offset = 0;
    depuncturedCW = coder.nullcopy(zeros(lengthLDPC,numCW));

    for idxCW = 1:numCW
        % Retrieve information bits
        inpBits = x(offset + (1:vecPayloadBits(idxCW)),1);
        % Size of the parity bits after puncturing
        pBlkSize = round(lengthLDPC*(1-rate)) - vecPunctureBits(idxCW);
        % Get parity bits
        parityBits = double(x(offset + vecPayloadBits(idxCW) + (1:pBlkSize),1));
        % Convert into LLRs
        shortenBits =  realmax * ones(vecShortenBits(idxCW),1);
        % Extra bits to compensate for puncturing
        extraBits = zeros(vecPunctureBits(idxCW),1);
        % Depunctured codeword
        depuncturedCW(:,idxCW) = [inpBits;shortenBits;parityBits;extraBits];
        offset = offset + vecPayloadBits(idxCW) + pBlkSize + vecRepeatBits(idxCW);
    end

    [out,numIterations,parityCheck] =  wlan.internal.ldpcDecodeCore(...
        depuncturedCW,rate,algChoice,alphaBeta,maxNumIter,earlyTermination,'soft decision');

    idx = 0;
    for idxCW = 1:numCW
        y(idx + (1:vecPayloadBits(idxCW))) = out(1:vecPayloadBits(idxCW),idxCW);
        idx = vecPayloadBits(idxCW) + idx;
    end

end

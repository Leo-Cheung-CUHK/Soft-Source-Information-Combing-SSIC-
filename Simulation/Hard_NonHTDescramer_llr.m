function [psdu_llr, scramInit,scramInitBits,txBits,errorBits,pass] = Hard_NonHTDescramer_llr(PSDULength,f_total,r_total,yllr,L,mpdu)

pass = 0;
txBits = 0;
errorBits = 0;

%% decode r_1,r_2,...,r_7 (Derive initial state of the scrambler)
[~,r_index]=max(f_total);
scramInitBits = r_total(r_index,:).';

yllr = yllr(1:(16+8*PSDULength));
scramblerInitBits = scramInitBits;

% Remove pad and tail bits, and descramble
if all(scramInitBits==0)
    % Scrambler initialization invalid (0), therefore do not descramble
    descramDataOut = yllr;
else
    buffSize = min(127,size(yllr,1));
    I = zeros(buffSize,1,'int8');
    
    % Scrambling sequence generated using generator polynomial
    for d = 1:buffSize
        I(d) = xor(scramblerInitBits(1),scramblerInitBits(4)); % x7 xor x4
        scramblerInitBits(1:end-1) = scramblerInitBits(2:end); % Left-shift
        scramblerInitBits(7) = I(d);                           % Update x1
    end

    scramblerSequence = repmat(I,ceil(size(yllr,1)/buffSize),1);
    
    scramblerSequenceSign = double(scramblerSequence*(-2)+1);
    descramLLR = yllr.*scramblerSequenceSign(1:size(yllr,1));
    %descramDataOut = double((descramLLR<0));
end


% Remove the 16 service bits + 111 zeros bits
if (L >= 16)
    psdu_llr = descramLLR(L+1:end);
    mpdu     = mpdu((L-16)+1:end);
else
    psdu_llr = descramLLR(17:end);
end

if length(mpdu) == length(psdu_llr)
    txBits    = length(psdu_llr);
    errorBits = biterr(mpdu,double((psdu_llr<0)));
    pass = 1;
end

% Convert scrambler initialization bits to number
scramInit = bi2de(scramInitBits.', 'left-msb');

end
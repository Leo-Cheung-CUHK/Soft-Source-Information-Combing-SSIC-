function [psdu_llr, scramInit,scramInitBits] = Hard_NonHTDescramer_llr(PSDULength,f_total,r_total,yllr)

%% decode r_1,r_2,...,r_7 (Derive initial state of the scrambler)
[~,r_index]=max(f_total);
scramInitBits = r_total(r_index,:).';

yllr = yllr(1:(16+8*PSDULength));
scramblerInitBits = scramInitBits;

% Remove pad and tail bits, and descramble
if all(scramInitBits==0)
    % Scrambler initialization invalid (0), therefore do not descramble
    descramLLR = yllr;
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

% Remove the 16 service bits 
psdu_llr = descramLLR(17:end);

% Convert scrambler initialization bits to number
scramInit = bi2de(scramInitBits.', 'left-msb');

end
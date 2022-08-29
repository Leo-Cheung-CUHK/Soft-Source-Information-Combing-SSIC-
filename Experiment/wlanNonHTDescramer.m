function [psdu, scramInit,I] = wlanNonHTDescramer(decBits,PSDULength)
I = zeros(127,1);
% Derive initial state of the scrambler
scramSeqInit = decBits(1:7);
scramInitBits = wlan.internal.scramblerInitialState(scramSeqInit);

% Remove pad and tail bits, and descramble
if all(scramInitBits==0)
    % Scrambler initialization invalid (0), therefore do not descramble
    descramDataOut = decBits(1:(16+8*PSDULength));
else
    [descramDataOut,I] = wlanScramble(decBits(1:(16+8*PSDULength)), scramInitBits);
end

% Remove the 16 service bits
psdu = descramDataOut(17:end);

% Convert scrambler initialization bits to number
scramInit = bi2de(scramInitBits.', 'left-msb');
end
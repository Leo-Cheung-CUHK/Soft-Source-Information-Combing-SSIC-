function validateScramblerInitialization(scramInit,inDSSSMode,isDMG,isNonHT,isWUR,isOFDM,cfgFormat)
    % Validate ScramblerInitialization
    if (inDSSSMode || isWUR)
        return % Not applicable if DSSS or WUR
    end
    if isDMG && any(scramInit~=93)
        if strcmp(phyType(cfgFormat),'Control')
            coder.internal.errorIf(any((scramInit<1) | (scramInit>15)),'wlan:wlanWaveformGenerator:InvalidScramblerInitialization','Control',1,15);
        elseif wlan.internal.isDMGExtendedMCS(cfgFormat.MCS)
            % At least one of the initialization bits must be non-zero,
            % therefore determine if the pseudorandom part can be 0
            % given the extended MCS and PSDU length.
            if all(wlan.internal.dmgExtendedMCSScramblerBits(cfgFormat)==0)
                minScramblerInit = 1; % Pseudorandom bits cannot be all zero
            else
                minScramblerInit = 0; % Pseudorandom bits can be all zero
            end
            coder.internal.errorIf(any((scramInit<minScramblerInit) | (scramInit>31)),'wlan:wlanWaveformGenerator:InvalidScramblerInitialization','SC extended MCS',minScramblerInit,31);
        else
            coder.internal.errorIf(any((scramInit<1) | (scramInit>127)),'wlan:wlanWaveformGenerator:InvalidScramblerInitialization','SC/OFDM',1,127);
        end
    else
        if isNonHT && isOFDM && ...
                any(strcmp(cfgFormat.ChannelBandwidth,{'CBW20','CBW40','CBW80','CBW160'})) ...
                && cfgFormat.SignalChannelBandwidth
            % Non-HT may include bandwidth signaling

            % Validate type
            validateattributes(scramInit,{'double','int8'}, ...
                               {'real','integer','2d','nonempty'},mfilename,'''ScramblerInitialization'' value');

            % Validate range
            range = scramblerRange(cfgFormat);
            minVal = range(1);
            maxVal = range(2);
            % Check for correct range
            if any((scramInit<minVal) | (scramInit>maxVal),'all')
                coder.internal.error('wlan:wlanWaveformGenerator:InvalidScramInitBWSignaling',minVal,maxVal);
            end
        else
            % Validate scrambler initialization
            validateattributes(scramInit,{'double','int8'},{'real','integer','2d','nonempty'},mfilename,'''ScramblerInitialization'' value');
            if any((scramInit<1) | (scramInit>127),'all')
                coder.internal.error('wlan:wlanWaveformGenerator:InvalidScramInit',1,127);
            end
        end
    end
end
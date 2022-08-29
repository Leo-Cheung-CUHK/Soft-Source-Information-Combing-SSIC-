function [macConfig, payload, status] = wlanMPDUDecode(mpdu, phyConfig, varargin)
%wlanMPDUDecode MPDU decoder
%   [MACCONFIG, PAYLOAD, STATUS] = wlanMPDUDecode(MPDU, PHYCONFIG)
%   validates the FCS and decodes the given MPDU. The decoding status is
%   returned along with a MAC frame configuration object and a payload.
%
%   MACCONFIG is the MAC frame configuration returned as an object of type
%   <a href="matlab:help('wlanMACFrameConfig')">wlanMACFrameConfig</a>.
%
%   PAYLOAD represents one or more MSDUs returned as a cell array
%   containing one or more character arrays, one for each MSDU. Each row in
%   the character array is the hexadecimal representation of an octet. For
%   all the MAC frames that do not contain data, PAYLOAD is returned as an
%   empty cell array.
%
%   STATUS represents the result of MPDU decoding, specified as an
%   enumeration value of type <a href="matlab:help('wlanMACDecodeStatus')">wlanMACDecodeStatus</a>. Any value of status
%   other than 'Success' (0) indicates that the MPDU decoding has failed.
%   If the decoding fails, the output MACCONFIG does not display any
%   properties as it may not be valid and the PAYLOAD is returned as
%   an empty cell array.
%   
%   MPDU represents the MAC protocol data unit, specified as one of the
%   following types:
%     - A binary vector representing MPDU bits
%     - A character vector representing octets in hexadecimal format.
%     - A string scalar representing octets in hexadecimal format.
%     - A numeric vector, where each element is in the range of [0 - 255]
%       inclusive, representing octets in decimal format.
%     - An n-by-2 character array, where each row represents an octet in
%       hexadecimal format.
%
%   PHYCONFIG is a format configuration object of type <a href="matlab:help('wlanHESUConfig')">wlanHESUConfig</a>,
%   <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>, or <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a>.
%
%   [MACCONFIG, PAYLOAD, STATUS] = wlanMPDUDecode(..., 'DataFormat', FRMT)
%   specifies the format, FRMT, of the input data MPDU as one of 'bits' or
%   'octets'. The default value is 'bits'. If it is specified as 'octets',
%   the MPDU can be a numeric vector representing octets in decimal format
%   or alternatively, it can be a character array or string scalar
%   representing octets in hexadecimal format. Otherwise, MPDU is a binary
%   vector.
%
%   Examples:
%
%   Example 1:
%   % Decode a data frame given in the form of octets
%
%       % Create a data frame
%       txFrameCfg = wlanMACFrameConfig('FrameType', 'QoS Data');
%       mpduOctets = wlanMACFrame(randi([0 255], 1, 40), txFrameCfg);
%
%       % PHY configuration
%       phyConfig = wlanNonHTConfig;
%
%       % Decode the data frame given in the form of octets
%       [rxFrameCfg, payload, status] = wlanMPDUDecode(mpduOctets, phyConfig, ...
%                                                   'DataFormat', 'octets');
%
%   Example 2:
%   % Decode a data frame given in the form of bits
%
%       % Create a data frame
%       txFrameCfg = wlanMACFrameConfig('FrameType', 'QoS Data');
%       mpduOctets = wlanMACFrame(randi([0 255], 1, 40), txFrameCfg);
%
%       % Convert the hexadecimal octets to bits
%       mpduBits = reshape(de2bi(hex2dec(mpduOctets), 8)', [], 1);
%
%       % PHY configuration
%       phyConfig = wlanNonHTConfig;
%
%       % Decode the data frame given in the form of bits
%       [rxFrameCfg, payload, status] = wlanMPDUDecode(mpduBits, phyConfig);
%
%   Example 3:
%   % Decode MPDUs extracted from an A-MPDU.
%
%       % Configure a VHT format A-MPDU
%       txFrameCfg = wlanMACFrameConfig('FrameType', 'QoS Data', ...
%                                       'FrameFormat', 'VHT');
%
%       % PHY configuration
%       phyCfg = wlanVHTConfig;
%
%       % Create random payload. 8 cell array elements correspond to
%       % 8 MSDUs.
%       payload = repmat({randi([0 255], 1, 40)}, 1, 8);
%
%       % Create an A-MPDU containing 8 MPDUs
%       ampdu = wlanMACFrame(payload, txFrameCfg, phyCfg);
%
%       % Decode the A-MPDU
%       [mpduList, delimiterFails, status] = wlanAMPDUDeaggregate(ampdu, ...
%                                                   phyCfg, ...
%                                                   'DataFormat', 'octets');
%
%       % Decode all the MPDUs in the mpduList
%       for i = 1:numel(mpduList)
%           if ~delimiterFails(i)
%               [macCfg, payload, status] = wlanMPDUDecode(mpduList{i}, ...
%                                                   phyCfg, ...
%                                                   'DataFormat', 'octets');
%           end
%       end
%
%   See also wlanAMPDUDeaggregate, wlanMACFrame, wlanMACFrameConfig,
%   wlanMACManagementConfig.

%   Copyright 2018-2020 The MathWorks, Inc.

%#codegen

narginchk(2, 4);
% For codegen: Assigning the return data types to the properties whose
% default values are of different data type.
mgmtConfig = wlanMACManagementConfig('Timestamp', uint64(0), 'AdditionalRates', {'1 Mbps'});
macConfig = wlanMACFrameConfig('ManagementConfig', mgmtConfig);
payload = cell(1, 0);
nvPair = varargin;
% Refer Table 9-19 in IEEE Std 802.11-2016 for frame length limits.
maxVHTorHEMPDULength = 11454;
maxHTorNonHTMMPDULength = 2304;

% Validate input arguments
[status, mpduColVector] = validateInputs(mpdu, phyConfig, nvPair);
if status ~= wlanMACDecodeStatus.Success
    return;
end

% Frame Format
switch class(phyConfig)
    case 'wlanNonHTConfig'
        macConfig.FrameFormat = 'Non-HT';
    case 'wlanHTConfig'
        macConfig.FrameFormat = 'HT-Mixed';
        macConfig.MPDUAggregation = phyConfig.AggregatedMPDU;
    case 'wlanVHTConfig'
        macConfig.FrameFormat = 'VHT';
        % VHT format frame is aggregated by default
        macConfig.MPDUAggregation = true;
    otherwise % wlanHESUConfig
        macConfig.FrameFormat = phyConfig.packetFormat;
        % HE format frame is aggregated by default
        macConfig.MPDUAggregation = true;
end

% Validate FCS
[status, mpdu] = checkFCS(mpduColVector);
if status ~= wlanMACDecodeStatus.Success
    macConfig.DecodeFailed = true;
    return;
end

% Position of frame bits within the MPDU
pos = 1;

% Frame control (16-bits)
frameControl = mpdu(pos : pos+15);
[macConfig, status] = decodeFrameControl(macConfig, frameControl, status);
if status ~= wlanMACDecodeStatus.Success
    macConfig.DecodeFailed = true;
    return;
end
pos = pos + 16;

% MPDU Length
mpduLength = numel(mpdu)/8;

% Validate MPDU length limit for VHT and HE formats
if any(strcmp(macConfig.FrameFormat, {'VHT', 'HE-SU', 'HE-EXT-SU'}))
    % Maximum length limit for an MPDU is 11454 octets for VHT or HE
    % formats.
    if mpduLength > maxVHTorHEMPDULength
        status = wlanMACDecodeStatus.MaxMPDULengthExceeded;
        return;
    end
end

% Validate MMPDU length limit for Non-HT and HT formats
if strcmp(macConfig.getType, 'Management') && any(strcmp(macConfig.FrameFormat, {'Non-HT', 'HT-Mixed'}))
    % Maximum length limit for a management frame is 2304 octets for Non-HT
    % or HT formats.
    if mpduLength > maxHTorNonHTMMPDULength
        status = wlanMACDecodeStatus.MaxMMPDULengthExceeded;
        return;
    end
end

% Decode rest of the frame based on frame-type
switch macConfig.getType
    case 'Control'
        [macConfig, status] = decodeControlFrame(macConfig, mpdu(pos:end), status);
        
    case 'Management'
        [macConfig, status] = decodeManagementFrame(macConfig, mpdu(pos:end), status);
        
    otherwise % Data
        [macConfig, payload, status] = decodeDataFrame(macConfig, mpdu(pos:end), status);
end

if status ~= wlanMACDecodeStatus.Success
    macConfig.DecodeFailed = true;
    return;
end
end

% Decodes control frame
function [macConfig, status] = decodeControlFrame(macConfig, mpduBits, status)
    pos = 1;
    minOctets = 14;
    numDataBits = numel(mpduBits);

    % The smallest control frame contains 14 octets. FCS (4 octets) and
    % Frame Control (2 octets) fields are already parsed. There should be
    % at least 8 octets more.
    if (numDataBits < 8*8)
        coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseFrame', minOctets, 'control');
        status = wlanMACDecodeStatus.NotEnoughData;
        return;
    end

    % Duration (16 bits) (0 to 14 bits are valid for duration, 15th bit is
    % reserved)
    macConfig.Duration = bi2deOptimized(mpduBits(pos : pos+14)');
    pos = pos + 16;

    % Address1 (48 bits)
    address = dec2hex(bi2deOptimized(reshape(mpduBits(pos : pos+47), 8, [])'), 2);
    macConfig.Address1 = reshape(address', 1, []);
    pos = pos + 48;

    % Address2 (48 bits)
    if ~any(strcmp(macConfig.FrameType, {'CTS', 'ACK'}))
        minOctets = minOctets + 6;
        if numDataBits >= (pos + 47)
            address = dec2hex(bi2deOptimized(reshape(mpduBits(pos : pos+47), 8, [])'), 2);
            macConfig.Address2 = reshape(address', 1, []);
            pos = pos + 48;
        else
            coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseField', 'Address2', minOctets);
            status = wlanMACDecodeStatus.NotEnoughData;
            return;
        end
    end

    if strcmp(macConfig.FrameType, 'Block Ack')
        % BA Control field (16 bits)
        minOctets = minOctets + 2;
        if numDataBits >= (pos + 15)
            % BA Ack Policy
            pos = pos + 1;

            % Check BA variant. Only compressed BA is supported.
            baVariant = bi2deOptimized(mpduBits(pos : pos + 3)');
            % Compressed BA variant code is 2
            if (baVariant ~= 2)
                coder.internal.warning('wlan:wlanMPDUDecode:UnsupportedBAVariant', baVariant);
                status = wlanMACDecodeStatus.UnsupportedBAVariant;
                return;
            end
            pos = pos + 4;

            % Skip reserved fields
            pos = pos + 7;

            % TID Info (Only TIDs [0 - 7] are supported)
            macConfig.TID = bi2deOptimized(mpduBits(pos : pos + 2)');
            pos = pos + 4;
        else
            coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseField', 'BA Control', minOctets);
            status = wlanMACDecodeStatus.NotEnoughData;
            return;
        end

        % BA Starting Sequence Control field (16 bits)
        minOctets = minOctets + 2;
        if numDataBits >= (pos + 15)
            % Fragment Number: If bit-3 is 0 and bits-(1,2) is 0, bitmap
            % size is 8 octets. If bit-3 is 0 and bits-(1,2) is 2, bitmap
            % size is 32 octets. Other combinations are reserved.
            if (double(mpduBits(pos + 3)) == 0) && (bi2deOptimized(mpduBits(pos+1 : pos+2)') == 0)
                bitmapSize = 8;
            elseif (double(mpduBits(pos + 3)) == 0) && (bi2deOptimized(mpduBits(pos+1 : pos+2)') == 2)
                bitmapSize = 32;
            else
                status = wlanMACDecodeStatus.UnknownBitmapSize;
                return;
            end
            pos = pos + 4;

            % Sequence Number
            macConfig.SequenceNumber = bi2deOptimized(mpduBits(pos : pos + 11)');
            pos = pos + 12;
        else
            coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseField', 'BA Starting Sequence Control', minOctets);
            status = wlanMACDecodeStatus.NotEnoughData;
            return;
        end

        % BA Bitmap field
        minOctets = minOctets + bitmapSize;
        if numDataBits >= (pos + bitmapSize*8-1)
            % BA Bitmap
            bitmapOctets = bi2deOptimized(reshape(mpduBits(pos : pos+bitmapSize*8-1), 8, [])');
            bitmapOctets = bitmapOctets(end : -1 : 1);
            macConfig.BlockAckBitmap = reshape(dec2hex(bitmapOctets, 2)', 1, []);
        else
            coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseField', 'BA Bitmap', minOctets);
            status = wlanMACDecodeStatus.NotEnoughData;
            return;
        end
    end
end

% Decodes data frame
function [macConfig, payload, status] = decodeDataFrame(macConfig, mpduBits, status)
    pos = 1;
    minOctets = 28;
    payload = cell(1, 0);
    numDataBits = numel(mpduBits);

    % A minimum data frame contains at least 28 octets. FCS (4 octets) and
    % Frame Control (2 octets) fields are already parsed. There should be
    % at least 22 octets more.
    if (numDataBits < 22*8)
        coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseFrame', minOctets, 'Data');
        status = wlanMACDecodeStatus.NotEnoughData;
        return;
    end

    % Duration (16 bits) (0 to 14 bits are valid for duration, 15th bit is
    % reserved)
    macConfig.Duration = bi2deOptimized(mpduBits(pos : pos+14)');
    pos = pos + 16;

    % Address1 (48 bits)
    address = dec2hex(bi2deOptimized(reshape(mpduBits(pos : pos+47), 8, [])'), 2);
    macConfig.Address1 = reshape(address', 1, []);
    pos = pos + 48;

    % Address2 (48 bits)
    address = dec2hex(bi2deOptimized(reshape(mpduBits(pos : pos+47), 8, [])'), 2);
    macConfig.Address2 = reshape(address', 1, []);
    pos = pos + 48;

    % Address3 (48 bits)
    address = dec2hex(bi2deOptimized(reshape(mpduBits(pos : pos+47), 8, [])'), 2);
    macConfig.Address3 = reshape(address', 1, []);
    pos = pos + 48;

    % Sequence Control (16 bits)
    sequenceControl = mpduBits(pos : pos+15)';
    pos = pos + 16;
    macConfig.SequenceNumber = bi2deOptimized(sequenceControl(5 : 16));

    % Skip Address4 field (48 bits)
    if macConfig.ToDS && macConfig.FromDS
        minOctets = minOctets + 6;

        % Skip Address4 field
        if numDataBits >= (pos + 47)
            pos = pos + 48;
        else
            coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseField', 'Address4', minOctets);
            status = wlanMACDecodeStatus.NotEnoughData;
            return;
        end
    end

    qosControl = zeros(1, 16);

    % Decode the fields applicable to QoS frames
    if any(strcmp(macConfig.FrameType, {'QoS Data', 'QoS Null'}))
        minOctets = minOctets + 2;

        % QoS Control (16 bits)
        if (numDataBits >= (pos + 15))
            qosControl = mpduBits(pos : pos+15)';
            pos = pos + 16;
            % TID
            tid = bi2deOptimized(qosControl(1 : 3));
            macConfig.TID = tid;

            % AckPolicy
            ackPolicy = bi2deOptimized(qosControl(6 : 7));
            switch ackPolicy
                case 0
                    macConfig.AckPolicy = 'Normal Ack/Implicit Block Ack Request';
                case 1
                    macConfig.AckPolicy = 'No Ack';
                case 2
                    macConfig.AckPolicy = 'No explicit acknowledgment/PSMP Ack/HTP Ack';
                otherwise
                    macConfig.AckPolicy = 'Block Ack';
            end

            % AMSDUPresent flag
            macConfig.MSDUAggregation = double(qosControl(8));
        else
            coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseField', 'QoS Control', minOctets);
            status = wlanMACDecodeStatus.NotEnoughData;
            return;
        end

        % HT Control (32 bits)
        if ~strcmp(macConfig.FrameFormat, 'Non-HT') && macConfig.HTControlPresent
            minOctets = minOctets + 4;

            if numDataBits >= (pos + 31)
                htControl = dec2hex(bi2deOptimized(mpduBits(pos : pos+31)'), 8);
                macConfig.HTControl = htControl;
                pos = pos + 32;
            else
                coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseField', 'HT Control', minOctets);
                status = wlanMACDecodeStatus.NotEnoughData;
                return;
            end
        end
    end

    % Skip Mesh Control field
    if strcmp(macConfig.FrameType, 'QoS Data')
        mcsControlPresent = double(qosControl(9));
        minOctets = minOctets + 6; % Mesh Flags + Mesh TTL + Mesh sequence number

        if mcsControlPresent
            if (numDataBits >= (pos + 47))
                % Extract address extension mode (2 bits) from Mesh flags
                % field (1 octet)
                addressExtMode = bi2deOptimized(mpduBits(pos : pos+1)');

                % Mesh Control = Mesh Flags (1 octet) + Mesh TTL (1 octet)
                % + Mesh Sequence Number (4 octets) + Mesh Address
                % Extension (0 or 6 or 12 octets)

                % Skip (Mesh Flags + Mesh TTL + Mesh sequence number)
                pos = pos + 8 + 8 + 4*8;

                % Skip Address Extension Mode:
                % 0 - indicates no extra addresses
                % 1 - indicates 1 extra address
                % 2 - indicates 2 extra addresses
                % 3 - reserved
                if (addressExtMode ~= 3)
                    minOctets = minOctets + addressExtMode*6;
                    pos = pos + addressExtMode*6*8;
                else
                    status = wlanMACDecodeStatus.UnknownAddressExtMode;
                    return;
                end

                % Not enough data to skip mesh control and extract the payload
                if (numDataBits < (pos - 1))
                    coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseField', 'Mesh Control', minOctets);
                    status = wlanMACDecodeStatus.NotEnoughData;
                    return;
                end
            else
                % Not enough data to skip mesh control and extract the payload
                coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseField', 'Mesh Control', minOctets);
                status = wlanMACDecodeStatus.NotEnoughData;
                return;
            end
        end
    end

    % A-MSDU subframe header length
    amsduSubframeHdrLen = 14;
    % Maximum length limits of MSDU and A-MSDU to be validated. Refer Table
    % 9-19 in IEEE Std 802.11-2016 for frame length limits.
    maxNonHTAMSDULength = 4065; % in octets
    maxHTAMSDULength = 7935; % in octets
    maxMSDULength = 2304; % in octets

    % Frame Body (variable)
    if strcmp(macConfig.FrameType, 'QoS Data') && macConfig.MSDUAggregation
        % A-MSDU length
        amsduLength = numel(mpduBits(pos:end))/8;

        % Validate A-MSDU length
        if strcmp(macConfig.FrameFormat, 'Non-HT') && (amsduLength > maxNonHTAMSDULength)
            % Max A-MSDU length for Non-HT format is 4065 octets
            coder.internal.warning('wlan:wlanMPDUDecode:MaxAMSDULengthExceeded', 'Non-HT', maxNonHTAMSDULength);
            status = wlanMACDecodeStatus.MaxAMSDULengthExceeded;
            return;
        elseif strcmp(macConfig.FrameFormat, 'HT-Mixed') && (amsduLength > maxHTAMSDULength)
            % Max A-MSDU length for HT format is 7935 octets
            coder.internal.warning('wlan:wlanMPDUDecode:MaxAMSDULengthExceeded', 'HT', maxHTAMSDULength);
            status = wlanMACDecodeStatus.MaxAMSDULengthExceeded;
            return;
        end

        % Decode the A-MSDU: Iterate over the MSDUs in the A-MSDU
        while numDataBits >= (pos + amsduSubframeHdrLen*8)
            % A-MSDU Destination Address (6 octets)
            address = dec2hex(bi2deOptimized(reshape(mpduBits(pos : pos+47), 8, [])'), 2);
            macConfig.AMSDUDestinationAddress = reshape(address', 1, []);
            pos = pos + 48;

            % A-MSDU Source Address (6 octets)
            address = dec2hex(bi2deOptimized(reshape(mpduBits(pos : pos+47), 8, [])'), 2);
            macConfig.AMSDUSourceAddress = reshape(address', 1, []);
            pos = pos + 48;

            % MSDU length (2 octets)
            msduLength = bi2deOptimized([mpduBits(pos+8 : pos+15)' mpduBits(pos : pos+7)']);
            pos = pos + 16;

            if msduLength
                % Validate MSDU length
                if msduLength > maxMSDULength
                    status = wlanMACDecodeStatus.MaxMSDULengthExceeded;
                    payload = cell(1, 0);
                    return;
                end

                % Extract the MSDUs
                if (numDataBits >= (pos + msduLength*8-1))
                    payload{end + 1} = dec2hex(bi2deOptimized(reshape(mpduBits(pos : pos + msduLength*8-1), 8, [])'), 2);
                    pos = pos + msduLength*8;
                else
                    status = wlanMACDecodeStatus.MalformedAMSDULength;
                    payload = cell(1, 0);
                    return;
                end
            end

            % Skip subframe padding
            pad = abs(mod((amsduSubframeHdrLen + msduLength), -4));
            if pad
                pos = pos + pad*8;
            end
        end

    elseif any(strcmp(macConfig.FrameType, {'Data', 'QoS Data'}))
        % MSDU length
        msduLength = numel(mpduBits(pos:end))/8;

        % Validate MSDU length
        if msduLength > maxMSDULength
            status = wlanMACDecodeStatus.MaxMSDULengthExceeded;
            return;
        end

        % Extract the MSDU
        if (numDataBits >= pos)
            payload{end + 1} = dec2hex(bi2deOptimized(reshape(mpduBits(pos:end), 8, [])'), 2);
        end
    end
end

% Decodes management frame
function [macConfig, status] = decodeManagementFrame(macConfig, mpduBits, status)
    pos = 1;
    minOctets = 28;
    numDataBits = numel(mpduBits);

    % The smallest management frame contains 28 octets. FCS (4 octets) and
    % Frame Control (2 octets) fields are already parsed. There should be
    % at least 22 octets more.
    if (numDataBits < 22*8)
        coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseFrame', minOctets, 'management');
        status = wlanMACDecodeStatus.NotEnoughData;
        return;
    end

    % Duration (16 bits) (0 to 14 bits are valid for duration, 15th bit is
    % reserved)
    macConfig.Duration = bi2deOptimized(mpduBits(pos : pos+14)');
    pos = pos + 16;

    % Address1 (48 bits)
    address = dec2hex(bi2deOptimized(reshape(mpduBits(pos : pos+47), 8, [])'), 2);
    macConfig.Address1 = reshape(address', 1, []);
    pos = pos + 48;

    % Address2 (48 bits)
    address = dec2hex(bi2deOptimized(reshape(mpduBits(pos : pos+47), 8, [])'), 2);
    macConfig.Address2 = reshape(address', 1, []);
    pos = pos + 48;

    % Address3 (48 bits)
    address = dec2hex(bi2deOptimized(reshape(mpduBits(pos : pos+47), 8, [])'), 2);
    macConfig.Address3 = reshape(address', 1, []);
    pos = pos + 48;

    % Sequence Control (16 bits)
    sequenceControl = mpduBits(pos : pos+15)';
    macConfig.SequenceNumber = bi2deOptimized(sequenceControl(5 : 16));
    pos = pos + 16;

    % Frame Body (variable) - Only Beacon frame is supported. So the
    % frame-body corresponds to the beacon frame.

    % Timestamp (64-bits)
    minOctets = minOctets + 8;
    if numDataBits >= (pos + 63)
        timestampLSB = uint64(bi2deOptimized(mpduBits(pos : pos + 31)'));
        timestampMSB = uint64(bi2deOptimized(mpduBits(pos + 32 : pos + 63)'));
        macConfig.ManagementConfig.Timestamp = bitor(timestampLSB, bitshift(timestampMSB, 32));
        pos = pos + 64;
    else
        coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseField', 'Timestamp', minOctets);
        status = wlanMACDecodeStatus.NotEnoughData;
        return;
    end

    % Beacon Interval (16-bits)
    minOctets = minOctets + 2;
    if numDataBits >= (pos + 15)
        macConfig.ManagementConfig.BeaconInterval = bi2deOptimized(mpduBits(pos : pos+15)');
        pos = pos + 16;
    else
        coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseField', 'Beacon Interval', minOctets);
        status = wlanMACDecodeStatus.NotEnoughData;
        return;
    end

    % Capability Information field (16-bits)
    minOctets = minOctets + 2;
    if numDataBits >= (pos + 15)
        macConfig.ManagementConfig.ESSCapability = double(mpduBits(pos));
        pos = pos + 1;
        macConfig.ManagementConfig.IBSSCapability = double(mpduBits(pos));
        pos = pos + 1;
        pos = pos + 1; % CF-Pollable
        pos = pos + 1; % CF-Poll Request
        macConfig.ManagementConfig.Privacy = double(mpduBits(pos));
        pos = pos + 1;
        macConfig.ManagementConfig.ShortPreamble = double(mpduBits(pos));
        pos = pos + 1;
        pos = pos + 1; % Reserved
        pos = pos + 1; % Reserved
        macConfig.ManagementConfig.SpectrumManagement = double(mpduBits(pos));
        pos = pos + 1;
        macConfig.ManagementConfig.QoSSupport = double(mpduBits(pos));
        pos = pos + 1;
        macConfig.ManagementConfig.ShortSlotTimeUsed = double(mpduBits(pos));
        pos = pos + 1;
        macConfig.ManagementConfig.APSDSupport = double(mpduBits(pos));
        pos = pos + 1;
        macConfig.ManagementConfig.RadioMeasurement = double(mpduBits(pos));
        pos = pos + 1;
        pos = pos + 1; % Reserved
        macConfig.ManagementConfig.DelayedBlockAckSupport = double(mpduBits(pos));
        pos = pos + 1;
        macConfig.ManagementConfig.ImmediateBlockAckSupport = double(mpduBits(pos));
        pos = pos + 1;
    else
        coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseField', 'Capability Information', minOctets);
        status = wlanMACDecodeStatus.NotEnoughData;
        return;
    end

    % Flags for tracking mandatory IEs
    ssidIEPresent = false;
    supportedRatesIEPresent = false;

    while numDataBits >= (pos + 15)
        % Element ID (1 octet)
        elementID = bi2deOptimized(mpduBits(pos : pos+7)');
        pos = pos + 8;

        % Element Length (1 octet)
        ieLength = bi2deOptimized(mpduBits(pos : pos+7)');
        pos = pos + 8;

        % Information field (variable)
        if (numDataBits >= pos+ieLength*8-1)
            elementIDExtension = 0;

            % Element ID Extension
            if elementID == 255
                elementIDExtension = bi2deOptimized(mpduBits(pos : pos+7)');
                % IE length includes the extension ID. Refer section
                % 9.4.2.1 of IEEE Std 802.11-2016
                ieLength = ieLength - 1;
                pos = pos + 8;
            end

            % IE information
            if ieLength > 0
                informationOctets = bi2deOptimized(reshape(mpduBits(pos : pos+ieLength*8-1), 8, [])');
                pos = pos + ieLength*8;
            else
                % If IE length is 0, return a 0-by-1 empty double array
                % that will be converted to a 1-by-0 char array in the
                % following steps.
                informationOctets = zeros(0, 1);
            end

            if elementID == 0
                % SSID must not exceed 32 octets
                if ieLength > 32
                    status = wlanMACDecodeStatus.MalformedSSID;
                    return;
                end

                % SSID IE (variable)
                macConfig.ManagementConfig.SSID = char(informationOctets)';
                ssidIEPresent = true;

            elseif elementID == 1
                % Number of supported rates must not exceed 8 rates
                if (ieLength < 1) || (ieLength > 8)
                    coder.internal.warning('wlan:wlanMPDUDecode:MalformedSupportedRatesIE', 'Beacon', ieLength);
                    status = wlanMACDecodeStatus.MalformedSupportedRatesIE;
                    return;
                end

                % Supported Rates IE (variable)
                rates = informationOctets;
                basicRates = cell(1, 0);
                additionalRates = cell(1, 0);
                for i = 1:numel(rates)
                    % Basic Rates: If the encoded rate has most significant
                    % bit set to 1, it is considered as a basic rate.
                    if (bitand(rates(i), 128) == 128)
                        datarate = getDataRate(bitand(rates(i), 127));
                        if ~isempty(datarate)
                            basicRates{end + 1} = datarate;
                        end
                    else % Additional Rates
                        datarate = getDataRate(bitand(rates(i), 127));
                        if ~isempty(datarate)
                            additionalRates{end + 1} = datarate;
                        end
                    end
                end
                macConfig.ManagementConfig.BasicRates = basicRates;
                macConfig.ManagementConfig.AdditionalRates = additionalRates;
                supportedRatesIEPresent = true;

            else
                % IE information
                information = reshape(dec2hex(informationOctets, 2)', 1, []);

                % Fill the IE information in the management configuration
                if elementID == 255
                    % IEs with element ID extension
                    macConfig.ManagementConfig = macConfig.ManagementConfig.addIE([elementID elementIDExtension], information);
                else
                    % IEs without element ID extension
                    macConfig.ManagementConfig = macConfig.ManagementConfig.addIE(elementID, information);
                end
            end
        else
            % Decoded IE length is found invalid. The specified length is more
            % than the remaining data.
            coder.internal.warning('wlan:wlanMPDUDecode:MalformedIELength', 'Beacon', elementID);
            status = wlanMACDecodeStatus.MalformedIELength;
            return;
        end
    end

    % SSID and Supported Rates IEs are mandatory for a beacon frame. If any
    % of these IEs are missing, consider it as a malformed packet.
    if ~ssidIEPresent || ~supportedRatesIEPresent
        status = wlanMACDecodeStatus.MissingMandatoryIEs;
        return;
    end
end

% Decodes frame control field
function [macConfig, status] = decodeFrameControl(macConfig, frameControl, status)
    pos = 1;

    % Protocol Version (2-bits)
    protocolVersion = bi2deOptimized(frameControl(pos : pos+1)');
    if (protocolVersion ~= 0)
        status = wlanMACDecodeStatus.InvalidProtocolVersion;
        return;
    end
    pos = pos + 2;

    % Type (2-bits)
    [type, status] = getType(frameControl(pos : pos+1), status);
    if status ~= wlanMACDecodeStatus.Success
        return;
    end
    pos = pos + 2;

    % Subtype (4-bits)
    [subtype, status]= getSubtype(frameControl(pos : pos+3), type, status);
    if status ~= wlanMACDecodeStatus.Success
        return;
    end
    macConfig.FrameType = subtype;
    pos = pos + 4;

    % ToDS (1-bit)
    macConfig.ToDS = double(frameControl(pos));
    pos = pos + 1;

    % FromDS (1-bit)
    macConfig.FromDS = double(frameControl(pos));
    pos = pos + 1;

    % MoreFragments (1-bit)
    pos = pos + 1;

    % Retransmission (1-bit)
    macConfig.Retransmission = double(frameControl(pos));
    pos = pos + 1;

    % Power Management (1-bit)
    macConfig.PowerManagement = double(frameControl(pos));
    pos = pos + 1;

    % More Data (1-bit)
    macConfig.MoreData = double(frameControl(pos));
    pos = pos + 1;

    % Protected Frame (1-bit)
    pos = pos + 1;

    % Order (1-bit)
    order = double(frameControl(pos));

    if any(strcmp(macConfig.FrameType, {'QoS Data', 'QoS Null'})) && ~strcmp(macConfig.FrameType, 'Non-HT')
        % +HTC (1-bit)
        macConfig.HTControlPresent = order;
    end
end

% Checks FCS
function [status, mpdu] = checkFCS(mpduWithFCS)
    persistent detect

    % CRC Detector object
    if isempty(detect)
        % Refer section 9.2.48 in IEEE Std 802.11-2016 for FCS calculation.
        detect = comm.CRCDetector([32 26 23 22 16 12 11 10 8 7 5 4 2 1 0], 'InitialConditions', 1, 'DirectMethod', true, 'FinalXOR', 1);
    end

    % Validate the FCS
    [mpdu, err] = detect(double(mpduWithFCS));
    mpdu = reshape(mpdu, [], 1);

    % Update status
    if err
        status = wlanMACDecodeStatus.FCSFailed;
    else
        status = wlanMACDecodeStatus.Success;
    end
end

% Return frame type code
function [type, status] = getType(typeCode, status)
    code = bi2deOptimized(typeCode');
    switch code
        case 0
            type = 'Management';
        case 1
            type = 'Control';
        case 2
            type = 'Data';
        otherwise
            % For codegen
            type = '';
            % Display a warning specifying the code of the unsupported frame
            % type and also mention the list of supported frame type codes.
            coder.internal.warning('wlan:wlanMPDUDecode:UnsupportedFrameType', code);
            status = wlanMACDecodeStatus.UnsupportedFrameType;
    end
end

% Return frame subtype code
function [subtype, status] = getSubtype(subtypeCode, type, status)
    code = bi2deOptimized(subtypeCode');

    if strcmp(type, 'Management')
        if (code == 8)
            subtype = 'Beacon';
        else
            % For codegen
            subtype = '';
            % Display a warning specifying the code of the unsupported
            % management subtype and also mention the list of supported
            % subtype codes.
            coder.internal.warning('wlan:wlanMPDUDecode:UnsupportedFrameSubtype', code, 'management', '8');
            status = wlanMACDecodeStatus.UnsupportedFrameSubtype;
        end
    elseif strcmp(type, 'Data')
        switch code
            case 0
                subtype = 'Data';
            case 4
                subtype = 'Null';
            case 8
                subtype = 'QoS Data';
            case 12
                subtype = 'QoS Null';
            otherwise
                % For codegen
                subtype = '';
                % Display a warning specifying the code of the unsupported
                % data subtype and also mention the list of supported
                % subtype codes.
                coder.internal.warning('wlan:wlanMPDUDecode:UnsupportedFrameSubtype', code, ...
                    'data', '0, 4, 8, and 12');
                status = wlanMACDecodeStatus.UnsupportedFrameSubtype;
        end
    else % Control
        switch code
            case 11
                subtype = 'RTS';
            case 12
                subtype = 'CTS';
            case 13
                subtype = 'ACK';
            case 9
                subtype = 'Block Ack';
            otherwise
                % For codegen
                subtype = '';
                % Display a warning specifying the code of the unsupported
                % control subtype and also mention the list of supported
                % subtype codes
                coder.internal.warning('wlan:wlanMPDUDecode:UnsupportedFrameSubtype', ...
                    code, 'control', '9, 11, 12, and 13');
                status = wlanMACDecodeStatus.UnsupportedFrameSubtype;
        end
    end
end

% Returns the data rate for the given code
function rate = getDataRate(code)
    % Refer Table-18.4 in Std IEEE 802.11-2016
    switch(code)
        case 2
            rate = '1 Mbps';
        case 4
            rate = '2 Mbps';
        case 11
            rate = '5.5 Mbps';
        case 12
            rate = '6 Mbps';
        case 18
            rate = '9 Mbps';
        case 22
            rate = '11 Mbps';
        case 24
            rate = '12 Mbps';
        case 36
            rate = '18 Mbps';
        case 48
            rate = '24 Mbps';
        case 72
            rate = '36 Mbps';
        case 96
            rate = '48 Mbps';
        case 108
            rate = '54 Mbps';
        otherwise
            % For codegen
            rate = '';
            coder.internal.warning('wlan:wlanMPDUDecode:UnknownRateReceived', code);
    end
end

% Validates the input arguments
function [status, mpduBits] = validateInputs(mpdu, phyConfig, nvPair)
    % Initialize
    status = wlanMACDecodeStatus.Success;
    mpduLength = numel(mpdu);
    isAMPDU = 0;

    % Validate the inputs
    dataFormat = wlan.internal.validateMACDecodeInputs(mpdu, phyConfig, nvPair, isAMPDU);

    if strcmpi(dataFormat, 'bits')
        % Validate MPDU given in the form of bits
        validateattributes(mpdu, {'logical', 'numeric'}, {'binary', 'vector'}, '', 'MPDU');
        coder.internal.errorIf((rem(mpduLength, 8) ~= 0), 'wlan:shared:InvalidDataSize');
        mpduBits = double(reshape(mpdu, [], 1));

    else % dataFormat == 'octets'
        % MPDU format must be in either hexadecimal or decimal octets
        validateattributes(mpdu, {'char', 'numeric', 'string'}, {}, mfilename, 'MPDU')

        if isnumeric(mpdu)
            validateattributes(mpdu, {'numeric'}, {'vector', 'integer', 'nonnegative', '<=', 255}, mfilename, 'MPDU');
            mpduBits = double(reshape(de2bi(mpdu, 8)', [], 1));

        else % char or string
            if ischar(mpdu)
                if isvector(mpdu)
                    % Convert row vector to column of octets.
                    columnOctets = reshape(mpdu, 2, [])';
                else
                    validateattributes(mpdu, {'char'}, {'2d', 'ncols', 2}, mfilename, 'MPDU', 1);
                    columnOctets = mpdu;
                end
            else % string
                validateattributes(mpdu, {'string'}, {'scalar'}, mfilename, 'MPDU')

                % Convert octets to char type
                columnOctets = reshape(char(mpdu), 2, [])';
            end

            % Validate hex-digits
            wlan.internal.validateHexOctets(columnOctets, 'MPDU');

            % Converting hexadecimal format octets to integer format
            decOctets = hex2dec(columnOctets);

            mpduBits = reshape(de2bi(decOctets, 8)', [], 1);
        end
    end

    % Validate minimum length required to decode an MPDU. FCS (4 octets)
    % and Frame Control (2-octets) fields are the basic fields needed to
    % decode an MPDU.
    if (numel(mpduBits)/8 < 6)
        status = wlanMACDecodeStatus.NotEnoughData;
        coder.internal.warning('wlan:wlanMPDUDecode:NotEnoughDataToParseMPDU');
    end
end

% Optimized bi2de taken from comms.internal.utilities
function dec = bi2deOptimized(bin)
    dec = comm.internal.utilities.bi2deRightMSB(double(bin), 2);
end
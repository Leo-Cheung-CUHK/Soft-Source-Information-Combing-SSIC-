    function numDataSym = getNumDataSymbols(pCfgNonHT)
    % Get number of OFDM data symbols
        mcsTable = wlan.internal.getRateTable(pCfgNonHT); 
        Ntail = 6; Nservice = 16;
        numDataSym = ceil((8*pCfgNonHT.PSDULength + Nservice + Ntail)/mcsTable.NDBPS);
    end    
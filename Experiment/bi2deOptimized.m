% Optimized bi2de taken from comms.internal.utilities
function dec = bi2deOptimized(bin)
    dec = comm.internal.utilities.bi2deRightMSB(double(bin), 2);
end
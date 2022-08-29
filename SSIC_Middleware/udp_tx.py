import socket
import time
import scipy.io

import sys

# sys.setdefaultencoding('cp1252')

client_IP   = '192.168.10.1'
client_port = 9999

server_IP   = '192.168.10.2'
server_port = 9999

s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
s.bind((client_IP, client_port))

# import matlab data
mat = scipy.io.loadmat('SN_string.mat')
SN_strings = mat['SN_strings']

tx_index = 0

while True:
    time.sleep(0.001)

    print( SN_strings[tx_index])
    print(len( SN_strings[tx_index]))

    tx_data = SN_strings[tx_index]
        
    s.sendto(tx_data.encode('ISO-8859-1'), (server_IP, server_port))

    tx_index = tx_index + 1
    if tx_index == len(SN_strings):
        tx_index = 0

s.close()
import pytun
import struct
import socket
import select
import queue
import os
from pytun import TunTapDevice

import random

def get_ip(addr):
    return '.'.join(map(str, addr))

def tcp_head( raw_data):
    (src_port, dest_port, sequence, acknowledgment, offset_reserved_flags) =  struct.unpack('! H H L L H', raw_data[:14])
    offset = (offset_reserved_flags >> 12) * 4
    flag_urg = (offset_reserved_flags & 32) >> 5
    flag_ack = (offset_reserved_flags & 16) >> 4
    flag_psh = (offset_reserved_flags & 8) >> 3
    flag_rst = (offset_reserved_flags & 4) >> 2
    flag_syn = (offset_reserved_flags & 2) >> 1
    flag_fin = offset_reserved_flags & 1
    data = raw_data[offset:]
    return src_port, dest_port, sequence, acknowledgment, flag_urg, flag_ack, flag_psh, flag_rst, flag_syn, flag_fin, data

def udp_head(raw_data):
    (src_port, dest_port, length, checksum) =  struct.unpack('! H H H H', raw_data[:8])
    data = raw_data[8:]
    return src_port, dest_port, length, checksum, data

def get_ip(addr):
    return '.'.join(map(str, addr))

def ipv4_head(raw_data):
    version_header_length = raw_data[0]
    version = version_header_length >> 4
    header_length = (version_header_length & 15) * 4
    ttl, proto,src, target = struct.unpack('! 8x B B 2x 4s 4s', raw_data[0:20])
    src = get_ip(src)
    target = get_ip(target)
    data = raw_data[header_length:]
    return version, header_length, ttl, proto, src, target, data

# Return properly formatted MAC address: (ie AA:BB:CC:DD:EE:FF)
def get_mac_addr(bytes_addr):
    bytes_str = map('{:02x}'.format, bytes_addr)
    result = ':'.join(bytes_str).upper()
    return result

def ethernet_head(raw_data):
    dest, src, prototype = struct.unpack('! 6s 6s H', raw_data[:14])
    dest_mac = get_mac_addr(dest)
    src_mac = get_mac_addr(src)
    proto = socket.htons(prototype)
    data = raw_data[14:]
    mac_header = raw_data[0:14]
    return dest_mac, src_mac, proto, data, mac_header

### Global Setting 
local_IP = '192.168.10.1'
peer_IP = '192.168.10.2'

trans_SN = 0
recv_SN_v = []

tap_MAC = b'\x00\x11\x22\x33\x44\x55'

# TWO NIC 
local_NIC_A = b'\x20\x7C\x8F\x9B\x56\xE7'
local_NIC_A_string = '20:7c:8f:9b:56:e7'
peer_NIC_A = b'\x20\x7C\x8F\x9B\x85\x78'
peer_NIC_A_string = '20:7c:8f:9b:85:78'

interface_A = 'wlp5s0'

# One TAP
tap = TunTapDevice('mytap', pytun.IFF_TAP | pytun.IFF_NO_PI)
# print (tap.name)

tap.addr = local_IP
tap.netmask = '255.255.255.0'

tap.hwaddr = tap_MAC

tap.mtu = 1500
tap.persist(True)
tap.up()
# print("tap created.")

os.system('sudo ip route add '+peer_IP+'/32 dev mytap')
os.system('sudo ip route add '+local_IP+'/32 dev mytap')

# Manually setup the ARP table to avoid auto-ARP process
os.system('sudo arp -i mytap -s '+peer_IP+' '+peer_NIC_A_string)
os.system('sudo arp -i mytap -s '+local_IP+' '+local_NIC_A_string)

ETH_P_ALL = 3
ETH_FRAME_LEN = 1514  # Max. octets in frame sans FCS

# RAW sockets
raw_socket_A = socket.socket(socket.AF_PACKET, socket.SOCK_RAW, socket.htons(ETH_P_ALL))
raw_socket_A.bind((interface_A, 0))
print("Raw socket is created.")

TX_SN = 0

try:
    while True:
        packet = tap.read(tap.mtu)
        eth_packet = ethernet_head(packet)

        if eth_packet[2] == 8:
            ipv4_packet = ipv4_head(eth_packet[3])
            if ipv4_packet[5] == peer_IP:
                # Use SSIC packet format 
                trans_SN = trans_SN + 1

                MAC_A_header = peer_NIC_A + local_NIC_A + eth_packet[4][12:14] 
                # packet_A = MAC_A_header  + struct.pack('H',TX_SN) + eth_packet[3]
                packet_A = MAC_A_header  + eth_packet[3]

                eth_packet = ethernet_head(packet_A)
                print('\nOutgoing Ethernet Frame to NIC A:')
                print('\tDestination: {}, Source: {}, Protocol: {}'.format(eth_packet[0], eth_packet[1],eth_packet[2]))

                raw_socket_A.sendall(packet_A)

        # Update sequnce number
        if trans_SN == 1000:
            trans_SN = 0

except KeyboardInterrupt:
    tap.down()
    tap.close()
    raw_socket_A.close()
    os.system('sudo ip link delete '+tap.name)
    print("Exit Successfully")

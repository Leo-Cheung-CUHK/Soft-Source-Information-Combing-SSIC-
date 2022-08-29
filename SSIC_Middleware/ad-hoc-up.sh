#!/bin/bash

# Author: Lihao ZHANG
# sudo apt update
# sudo apt install net-tools

if [ $# -ne 4 ]
  then
    echo "Please input NIC_name ch_number ip_addr as input parameter!"
    exit
fi

nic_name=$1
ch_number=$2
ip_addr=$3
cell=$4

echo $nic_name
echo $ch_number
echo $ip_addr
echo $cell

# stop the network manager -- To revise - lihao
sudo service network-manager stop
sudo rfkill unblock all

sudo ifconfig $nic_name down
sudo iwconfig $nic_name channel $ch_number
sudo iwconfig $nic_name mode ad-hoc

sudo ifconfig $nic_name up
sudo ip link set $nic_name down

if [[ $ch_number -gt 14 ]]
then
  sudo iwconfig $nic_name essid 'sdr-ad-hoc-5G'
else
  sudo iwconfig $nic_name essid 'sdr-ad-hoc-2.4G'
fi

sudo ip link set $nic_name up
sudo iwconfig $nic_name ap $cell

# modulation type
# if [[ $ch_number -gt 14 ]]
# then
#   sudo iwconfig $nic_name modulation 11a 
# else
#   sudo iwconfig $nic_name modulation 11g 
# fi

# To fix the 802.11 parameters-- To revise  - lihao
sudo iwconfig $nic_name rate 6M        # data rate
sudo iwconfig $nic_name retry 1        # retransmission number
sudo iwconfig $nic_name rts off
sudo iwconfig $nic_name frag off 
sudo iwconfig $nic_name key off
# sudo iwconfig $nic_name power off

sudo ifconfig $nic_name $ip_addr netmask 255.255.255.0
iwconfig $nic_name

#!/bin/bash

BIN_LIST[0]="MTA_HOME"
BIN_LIST[1]="TCOFFEE_BIN"
BIN_LIST[2]="CLUSTALW_BIN"
BIN_LIST[3]="CLUSTALO_BIN"
BIN_LIST[4]="NORMD_BIN"
BIN_LIST[5]="STRIKE_BIN"

#Clean previously declared paths
for i in `seq 0 5`;
	do
                sed -i "/${BIN_LIST[$i]}/ d" /home/$USER/.bashrc
        done

#Add bins to bashrc
echo "Info: Enter the path of MTA home. Ex: /home/username/mta"

echo -n "Path of ${BIN_LIST[0]}: "
read tmp_bin
echo "export ${BIN_LIST[0]}=$tmp_bin" >> /home/$USER/.bashrc

echo "Info: Enter the absolute path of the following binaries."

for i in `seq 1 5`;
	do
                echo -n "Path of ${BIN_LIST[$i]}: "
		read tmp_bin
		echo "export ${BIN_LIST[$i]}=$tmp_bin" >> /home/$USER/.bashrc
        done

echo "Manually edit the .bashrc if a problem occurs."

source /home/$USER/.bashrc



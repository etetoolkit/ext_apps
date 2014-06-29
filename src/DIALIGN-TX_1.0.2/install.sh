#! /bin/bash
more LICENSE.TXT
 echo "press any key" 
 read
 cd source
 make clean
if [ ! -e ../bin ]; then
	mkdir ../bin
fi
 make
 mv dialign-tx ../bin

if [ -e ../bin/dialign-tx ]; then
echo "installation succesful"
fi

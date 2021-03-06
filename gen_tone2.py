import random
import wave
import struct
import math

import sys

if len(sys.argv)!=4:
	print "%s sample_rate seconds freq" % sys.argv[0]
	sys.exit(1)

SAMPLE_FREQ=int(sys.argv[1])
SAMPLE_LEN=int(sys.argv[2])*SAMPLE_FREQ
OUT_FREQA=int(sys.argv[3])
noise_output = wave.open('%dk_%dk.wav' % (OUT_FREQA/1000,SAMPLE_FREQ/1000), 'w')
noise_output.setparams((2, 2, SAMPLE_FREQ, 0, 'NONE', 'not compressed'))

values = []

a=0
for i in range(0, SAMPLE_LEN):
	value = int(math.sin(2*float(i)*OUT_FREQA*math.pi/(SAMPLE_FREQ))*32767/2)
	scale = float(math.sin(math.pi*float(i)/SAMPLE_LEN))
	value = value*scale
        packed_value = struct.pack('h', value)
        values.append(packed_value)
        values.append(packed_value)
print SAMPLE_LEN,a

value_str = ''.join(values)
noise_output.writeframes(value_str)

noise_output.close()

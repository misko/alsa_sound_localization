import random
import wave
import struct
import math
import random
import sys

if len(sys.argv)!=5:
	print "%s sample_rate seconds freqA freqB" % sys.argv[0]
	sys.exit(1)

SAMPLE_FREQ=int(sys.argv[1])
SAMPLE_LEN=int(sys.argv[2])*SAMPLE_FREQ
OUT_FREQA=int(sys.argv[3])
OUT_FREQB=int(sys.argv[4])
noise_output = wave.open('phase_%dk.%dk_%dk.wav' % (OUT_FREQA/1000,OUT_FREQB/1000,SAMPLE_FREQ/1000), 'w')
noise_output.setparams((2, 2, SAMPLE_FREQ, 0, 'NONE', 'not compressed'))

values = []

a=0
i=0
d=1
while i<SAMPLE_LEN:
	l=min(SAMPLE_LEN-i,int(SAMPLE_LEN/100+int(random.random()*SAMPLE_LEN/100))) # int(random.random()*(max(SAMPLE_LEN-i,20000))+1000)
	#for i in range(0, SAMPLE_LEN):
	print l
	r=random.random()
	freq_from=OUT_FREQA*r+(1-r)*OUT_FREQB
	freq_to=OUT_FREQA*r+(1-r)*OUT_FREQB
	for j in range(0, l):
		t = float(j)/l
		value = int(math.sin(2*float(j)*((1-t)*freq_from+t*freq_to)*math.pi/(SAMPLE_FREQ))*32767/2)
		scale=1.0
		if t<0.1:
			scale = t/0.1
		elif t>0.9:
			scale = 1-(t-0.9)/0.1
		value = value*scale*scale
		packed_value = struct.pack('h', value)
		values.append(packed_value)
		values.append(packed_value)
	i+=l
print SAMPLE_LEN,a

value_str = ''.join(values)
noise_output.writeframes(value_str)

noise_output.close()

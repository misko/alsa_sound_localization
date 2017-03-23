import random
import wave
import struct
import math
import random
import sys

if len(sys.argv)!=4:
	print "%s sample_rate seconds freqA,freqB,freqC,freqD" % sys.argv[0]
	sys.exit(1)

SAMPLE_FREQ=int(sys.argv[1])
SAMPLE_LEN=int(sys.argv[2])*SAMPLE_FREQ
freqs=map(lambda x  : int(x) ,sys.argv[3].split(','))
noise_output = wave.open('phase.wav', 'w')
noise_output.setparams((2, 2, SAMPLE_FREQ, 0, 'NONE', 'not compressed'))

values = []

a=0
i=0
d=1
x=0
while i<SAMPLE_LEN:
	#l=min(SAMPLE_LEN-i,int(SAMPLE_LEN/100+int(random.random()*SAMPLE_LEN/100))) # int(random.random()*(max(SAMPLE_LEN-i,20000))+1000)
	l=SAMPLE_LEN
	#for i in range(0, SAMPLE_LEN):
	print l
	for j in range(0, l):
		t = float(j)/l
		value = 0
		for f in freqs:
			value += int(math.sin(2*float(j)*f*math.pi/(SAMPLE_FREQ))*32767/(len(freqs)))
		if value>x:
			x=value
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
print "MAX",x
print SAMPLE_LEN,a

value_str = ''.join(values)
noise_output.writeframes(value_str)

noise_output.close()

import random
import wave
import struct
import math

SAMPLE_LEN=11000
SAMPLE_FREQ=44100
OUT_FREQA=500
noise_output = wave.open('noise2.wav', 'w')
noise_output.setparams((2, 2, SAMPLE_FREQ, 0, 'NONE', 'not compressed'))

values = []

a=0
for i in range(0, SAMPLE_LEN):
        value = random.randint(-32767, 32767)
	value = int(math.sin(2*float(i)*OUT_FREQA*math.pi/(SAMPLE_FREQ))*32767/5)
	scale = float(math.sin(math.pi*float(i)/SAMPLE_LEN))
	value = value*scale
        packed_value = struct.pack('h', value)
        values.append(packed_value)
        values.append(packed_value)
print SAMPLE_LEN,a

value_str = ''.join(values)
noise_output.writeframes(value_str)

noise_output.close()

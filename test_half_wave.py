import sys
import time
import numpy as np
from datetime import datetime

sample_window=4096
RATE = 44100

def generate_tone(f,n_frames):
  return (np.sin(2*np.arange(n_frames)*np.pi*f/RATE))

audio = generate_tone(3000,sample_window)
audio -= np.average(audio)
audio /= np.var(audio)
fftData=abs(np.fft.rfft(audio))**2
freqs=np.fft.rfftfreq(len(audio),d=1.0/RATE)

print "Full wave"
for i in xrange(len(freqs)):
	if fftData[i]/fftData.sum()>0.1:
		print freqs[i],fftData[i]/fftData.sum()

audio = np.hstack((np.zeros(sample_window/2) ,generate_tone(3000,sample_window/2)))
audio -= np.average(audio)
audio /= np.var(audio)
fftData=abs(np.fft.rfft(audio))**2
freqs=np.fft.rfftfreq(len(audio),d=1.0/RATE)

print "Half wave"
for i in xrange(len(freqs)):
	if fftData[i]/fftData.sum()>0.1:
		print freqs[i],fftData[i]/fftData.sum()

#!/usr/bin/python -tt

import sys
import pyaudio
#from pylab import *
import time
import Queue
import threading
import numpy as np
#import matplotlib.animation as animation
#import matplotlib.pyplot as plt
import scipy.io.wavfile
import scipy.fftpack
from datetime import datetime
#import time
actual_steps=[]
actual_steps_bins=[]
RATE = 44100

sample_window=2048

Chunk = 1024
Format = pyaudio.paInt16

#Create a synchronized queue that can only hold one list
recordQueue = Queue.Queue(Chunk)

#Callback functions for playing back and recording audio in
#non-blocking mode. See reference 4
def playcallback(in_data, frame_count, time_info, status):
    out_data = transmitChirp.readframes(frame_count)
    return (out_data , pyaudio.paContinue)

def recordcallback(in_data, frame_count, time_info, status):
    recordQueue.put(np.fromstring(in_data, np.int8))
    return (in_data, pyaudio.paContinue)   

def showtimespectrum(audio):
    assert(len(audio)==sample_window)
    fftData=abs(np.fft.rfft(audio))**2
    freqs=np.fft.rfftfreq(len(audio),d=1.0/RATE)
    #print freqs[np.argmax(fftData)],np.argmax(fftData),actual_steps_bins,"MAX FREQ",freqs[-1]
    z=np.zeros(len(actual_steps))
    i=0
    for b in actual_steps_bins:
      z[i]=fftData[b]
      i+=1
    if z.sum()/fftData.sum()>0.1:
	    s=[]
	    for k in z:
		s.append("%0.2f" % (k/z.sum()))
            print datetime.now().strftime("%Y-%m-%d %H:%M %S.%f") ,",".join(s)


#def play_ladder(p,l=1,steps=[260,500,1000,2000,3000,17000]): #,2000,4000]):
def play_ladder(p,l=3): #,2000,4000]):
    steps=actual_steps
    stream = p.open(#format = p.get_format_from_width(1), 
                    channels = 1, 
                    rate = RATE, 
		    format=pyaudio.paUInt8,
                    output = True)
    while True:
      for f in steps:
        n_frames = int(RATE*l) # sample_window #int(RATE*l)
        stream.write((np.sin(2*np.arange(n_frames)*np.pi*f/RATE)*127+128).astype(np.uint8))
    stream.stop_stream()
    stream.close()


def main():
    if len(sys.argv)!=2:
	print "%s input_index[>=0for record]" % sys.argv[0]
	sys.exit(1)
    input_index=int(sys.argv[1])
    steps=[1000,3000,4000]
    freqs=np.fft.rfftfreq(sample_window,d=1.0/RATE)
    for f in steps:
      z=np.floor(f*sample_window*(1.0/RATE))
      actual_steps.append(freqs[z])
      actual_steps_bins.append(z)
    steps=actual_steps
	
    p = pyaudio.PyAudio()

    if input_index>=0:
	    recordstream = p.open(#format=p.get_format_from_width(SampWidth),
				  #channels=Channels,
				   format=pyaudio.paUInt8,
				  channels=1,
				  #rate=Rate,
				  rate=RATE,
				  input=True,
					input_device_index = input_index,
				  output=False,
				  stream_callback=recordcallback)
	    recordstream.start_stream()    
	    
	    while recordstream.is_active():                
		#Get the recorded data from the synchronized priority queue
		recordData = []
		while len(recordData)!=sample_window:
		    recordData.extend(recordQueue.get())
		    recordQueue.task_done()            
		showtimespectrum(recordData)
	    recordQueue.join()
	    recordstream.stop_stream()
	    recordstream.close()
    else:
	play_ladder(p)
    #Create a thread to chirp continuously
    #t1 = threading.Thread(target=play_ladder, args=(p,))
    #t1.setDaemon(True)
    #t1.start()
    #Create a "thread" to record continuosly
        
    #t1.join()
        
    
    p.terminate()

if __name__ == '__main__':
    main()

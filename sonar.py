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
response_steps=[]
response_steps_bins=[]
RATE = 44100

#sample_window=2048
sample_window=4096

Chunk = 1024
Format = pyaudio.paInt16

#Create a synchronized queue that can only hold one list
recordQueue = Queue.Queue(Chunk)
playQueue = Queue.Queue(Chunk)

#Callback functions for playing back and recording audio in
#non-blocking mode. See reference 4
def playcallback(in_data, frame_count, time_info, status):
    out_data = transmitChirp.readframes(frame_count)
    return (out_data , pyaudio.paContinue)

def recordcallback(in_data, frame_count, time_info, status):
    recordQueue.put(np.fromstring(in_data, np.int8))
    return (in_data, pyaudio.paContinue)   

def get_bin_from_FFT(fftData):
    z=np.zeros(len(actual_steps))
    i=0
    print len(fftData)
    for b in actual_steps_bins:
      assert(len(fftData) in (sample_window/2+1,sample_window/4+1))
      if len(fftData)==sample_window/4+1:
        b/=2
      z[i]=fftData[b]
      i+=1
    return z

def showtimespectrum(audio):
    assert(len(audio)==sample_window)
    fftData=abs(np.fft.rfft(audio))**2
    freqs=np.fft.rfftfreq(len(audio),d=1.0/RATE)
    #print freqs[np.argmax(fftData)],np.argmax(fftData),actual_steps_bins,"MAX FREQ",freqs[-1]
    z=get_bin_from_FFT(fftData)
    if z.sum()/fftData.sum()>0.1:
	    #l=[]
            #for i in xrange(len(actual_steps_bins)):
	    #	l.append((z[i]/z.sum(),i))
	    #l.sort(reverse=True)

	    #tones=[l[0][1],l[1][1]]
	    #if tones[0]>tones[1] and tones[0]!=len(actual_steps_bins)-1:
	    #	t=tones[0]
	    #	tones[0]=tones[1]
	    #	tones[1]=t
	    #print "testing,",tones 
	    
	    #do smaller FFTs with the main freq
	    if False:
		    for x in xrange(1,len(audio)-sample_window/2,sample_window/16):
			    fftData=abs(np.fft.rfft(audio[x:x+sample_window/2]))**2
			    freqs=np.fft.rfftfreq(sample_window/2,d=1.0/RATE)
			    zz=get_bin_from_FFT(fftData)
			    #print freqs[np.argmax(fftData)],np.argmax(fftData),actual_steps_bins,"MAX FREQ",freqs[-1]
			    print "WTF",zz/zz.sum()

	    #tones=[ [l[0][1],(l[0][1]+1)%len(actual_steps_bins) ] , [ (l[0][1]-1)%len(actual_steps_bins) , l[0][1]]] 
	    #for ft, tt in tones:
	    #	f_ft  = actual_steps[ft]
	    #	f_tt  = actual_steps[tt]
	    #	print f_ft,f_tt
	    #	a = np.hstack( ( generate_tone(f_ft,sample_window/8) , generate_tone(f_tt,sample_window/8 ) ) )
    	    #	zz=np.correlate(a, audio)*1.0
	    #	print zz.max(),zz.argmax()
	    s=[]
	    for k in z:
		s.append("%0.2f" % (k/z.sum()))
            print datetime.now().strftime("%Y-%m-%d %H:%M %S.%f") ,",".join(s)


#def play_ladder(p,l=1,steps=[260,500,1000,2000,3000,17000]): #,2000,4000]):
def play(p):#,l=3): #,2000,4000]):
    stream = p.open(#format = p.get_format_from_width(1), 
                    channels = 1, 
                    rate = RATE, 
		    format=pyaudio.paUInt8,
                    output = True)
    while True:
	d=playQueue.get()
        stream.write(d) 
    stream.stop_stream()
    stream.close()


def generate_tone(f,n_frames):
  return (np.sin(2*np.arange(n_frames)*np.pi*f/RATE)*127+128).astype(np.uint8)

def main():
    #a= np.hstack((generate_tone(10000,100),generate_tone(15000,100)))
    #v= np.hstack((np.zeros(100),generate_tone(10000,200),generate_tone(15000,100),generate_tone(1000,200)))
    #a= np.ones(300)
    #b= np.hstack((np.ones(200),np.zeros(100)))
    #z=np.correlate(v, a)*1.0
    #z/=z.max()
    #i=0
    #for x in z:
    #	if x>0.99:
    #		print i,x
    #	i+=1
    #sys.exit(1)
    if len(sys.argv)!=2:
	print "%s input_index[>=0for record]" % sys.argv[0]
	sys.exit(1)
    input_index=int(sys.argv[1])
    freqs=np.fft.rfftfreq(sample_window,d=1.0/RATE)
    #z=np.floor(1000*sample_window*(1.0/RATE))
    steps=[300,1500,4500]
    for f in steps:
      z=np.floor(f*sample_window*(1.0/RATE))
      if z%2!=0:
	z+=1
      actual_steps.append(freqs[z])
      actual_steps_bins.append(z)
      response_steps.append(freqs[z+10])
      response_steps_bins.append(z+10)
    steps=actual_steps
	
    p = pyaudio.PyAudio()
    #Create a thread to chirp continuously
    t1 = threading.Thread(target=play, args=(p,))
    t1.setDaemon(True)
    t1.start()
    #Create a "thread" to record continuosly

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
        steps=actual_steps
        while True:
          for f in steps:
            n_frames = sample_window #int(RATE*l) # sample_window #int(RATE*l)
	    playQueue.put(generate_tone(f,n_frames))
            #playQueue.put((np.sin(2*np.arange(n_frames)*np.pi*f/RATE)*127+128).astype(np.uint8))
        
    t1.join()
        
    
    p.terminate()

if __name__ == '__main__':
    main()

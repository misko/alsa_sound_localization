#!/usr/bin/python -tt

import sys
import pyaudio
import time
import Queue
import threading
import numpy as np
from datetime import datetime

sample_window=4096

blip_freq=900
blip_freq_bin=0
blip_length=sample_window
blop_freq=1500
blop_freq_bin=0
blop_length=sample_window
RATE = 44100

Format = pyaudio.paInt16

#Create a synchronized queue that can only hold one list
recordLock = threading.Lock()
recordQueue = Queue.Queue()
playQueue = Queue.Queue()
n_threads = 1
resultQueue = Queue.Queue()

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
    for b in actual_steps_bins:
      assert(len(fftData) in (sample_window/2+1,sample_window/4+1))
      if len(fftData)==sample_window/4+1:
        b/=2
      z[i]=fftData[b]
      i+=1
    return z

def analyze(audio,blip_blop):
    #blip_freq_bin,blop_freq_bin=blip_blop_freq_bins
    #assert(len(audio)==sample_window)
    fftData=abs(np.fft.rfft(audio))**2
    freqs=np.fft.rfftfreq(len(audio),d=1.0/RATE)
    #print freqs[np.argmax(fftData)],np.argmax(fftData),actual_steps_bins,"MAX FREQ",freqs[-1]
    #z=get_bin_from_FFT(fftData)
    z=fftData[blip_freq_bin]
    if blip_blop:
	z=fftData[blop_freq_bin]
    #print fftData[blip_freq_bin]/fftData.sum(), fftData[blop_freq_bin]/fftData.sum()
    if z/fftData.sum()>0.1:
            print datetime.now().strftime("%Y-%m-%d %H:%M %S.%f") , blip_blop
	    if not blip_blop: # we are blop and we heard blip
		playQueue.put(generate_tone(blop_freq,blop_length))
		


def play(p,blip_blop):#,l=3): #,2000,4000]):
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

def listen(x,blip_blop):
	while True:
		with recordLock:
			recordData = []
			while len(recordData)!=sample_window:
			    recordData.extend(recordQueue.get())
			    recordQueue.task_done() 
		#print "process",x
		analyze(recordData,blip_blop)

def generate_tone(f,n_frames):
  return (np.sin(2*np.arange(n_frames)*np.pi*f/RATE)*127+128).astype(np.uint8)

def round_freq(f):
    freqs=np.fft.rfftfreq(sample_window,d=1.0/RATE)
    f_bin=np.floor(f*sample_window*(1.0/RATE))
    if f_bin%2==1:
	f_bin+=1
    return freqs[f_bin],f_bin 

#def main():
if __name__ == '__main__':
    if len(sys.argv)!=3:
	print "%s input_index[>=0for record] blip[0]/blop[1]" % sys.argv[0]
	sys.exit(1)
    input_index=int(sys.argv[1])
    blip_blop=int(sys.argv[2])==0

    blip_freq,blip_freq_bin = round_freq(blip_freq)
    blop_freq,blop_freq_bin = round_freq(blop_freq)

    p = pyaudio.PyAudio()
    #Create a thread to chirp continuously
    emit_thread = threading.Thread(target=play, args=(p,blip_blop))
    emit_thread.setDaemon(True)
    emit_thread.start()
    #Create a "thread" to record continuosly

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

    #start listening threads
    listen_threads = []
    for x in range(n_threads):
	t = threading.Thread(target=listen,args=(x,blip_blop))
	t.setDaemon(True)
	t.start()
	listen_threads.append(t)

    if blip_blop: # play our blip!
	playQueue.put(generate_tone(blip_freq,blip_length))

    print "Online"
    while True:
	line=sys.stdin.readline().strip()
	if line=="X":
		break

    sys.exit(1)


    if input_index>=0:
	    print "START LISEN"
	    while True:
		line=sys.stdin.readline().strip()
		if line=="X":
			recordQueue.join()
			recordstream.stop_stream()
			recordstream.close()
    			p.terminate()
			sys.exit(1)
			break
		else:
			print "COMMAND %s not understood\n",line
    else:
        steps=actual_steps
        while True:
          for f in steps:
            n_frames = sample_window #int(RATE*l) # sample_window #int(RATE*l)
	    playQueue.put(generate_tone(f,n_frames))
            #playQueue.put((np.sin(2*np.arange(n_frames)*np.pi*f/RATE)*127+128).astype(np.uint8))
        
    t1.join()
        
    

    #main()

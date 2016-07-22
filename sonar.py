#!/usr/bin/python -tt

import pyaudio
import wave
from pylab import *
import time
import Queue
import threading
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import scipy.io.wavfile
import scipy.fftpack

actual_steps=[]
actual_steps_bins=[]
RATE = 44100

ChirpFile = "sound files/plain_chirp.wav"
transmitChirp = wave.open(ChirpFile, 'rb')
ChirpPeriod = 2

sample_window=2048

Chunk = 1024
Format = pyaudio.paInt16
SampWidth = transmitChirp.getsampwidth()
Channels = transmitChirp.getnchannels()
Rate = transmitChirp.getframerate()
#Time in seconds to display at one time on the waveform display
WindowTime = 1

#Create a synchronized queue that can only hold one list
recordQueue = Queue.Queue(Chunk)

#Callback functions for playing back and recording audio in
#non-blocking mode. See reference 4
def playcallback(in_data, frame_count, time_info, status):
    out_data = transmitChirp.readframes(frame_count)
    return (out_data , pyaudio.paContinue)

def recordcallback(in_data, frame_count, time_info, status):
    recordQueue.put(fromstring(in_data, 'Int16'))
    #print len(in_data)
    return (in_data, pyaudio.paContinue)   

#Require: Requires an exisitng pyAudio object which has set
#       up the port audio system
#Ensures: Plays the ChirpFile wave file over the speakers
#       if the speakers are available

#Note: look into threading this. I think the print delay is the only thing
#causing this to play properly
def chirp(p):
    while True:
        #Returns a tuple (nchannels, sampwidth, framerate, nframes, comptype, compname)
        audioParams = transmitChirp.getparams()
        playstream = p.open(format=p.get_format_from_width(audioParams[1]),
                        channels=audioParams[0],
                        rate=audioParams[2],
                        input=False,
                        output=True,
                        stream_callback=playcallback)
        audioTime = round(float(audioParams[-3]) / audioParams[2], 6)
        #print 'AudioTime:', audioTime
        #Keep just playing the chirp and showing the data
        playstream.start_stream()
        #Temporary delay. Use some function of chirp length or
        #travel time for audio to determine this after threading function
        time.sleep(ChirpPeriod)
        #pause(ChirpPeriod)
        playstream.stop_stream()
        playstream.close()
        #transmitChirp.setpos(0)
        transmitChirp.rewind()

#Requires: A start and stop number and the number of points n
#Ensures: Returns a list of floats of size n such that the first
#       item is the start number, the last item is the stop number
#Code from [5]
def myfrange(start, stop, n):
    L = [0.0] * n
    nm1 = n - 1
    nm1inv = 1.0 / nm1
    for i in range(n):
        L[i] = nm1inv * (start*(nm1 - i) + stop*i)
    return L

#Requires: A numpy array of non-zero length of audio data that has been recorded
#Ensures: Returns a numpy array result of the correlation with the ChirpFile

def findechoes(receivedAudio):
    #transmittedWave = wave.open(ChirpFile, 'rb')
    (sampRate, dataArray) = scipy.io.wavfile.read(ChirpFile)
    #print 'opened file at', sampRate
    numChunks = np.ceil(len(receivedAudio) / float(len(dataArray)))
    corrData = []
    for i in range(int(numChunks)):
        startIndex = len(dataArray) * i
        endIndex = startIndex + len(dataArray)
        corrData.extend(np.correlate(dataArray, 
                                     receivedAudio[startIndex:endIndex]))
    #print np.correlate(dataArray, receivedAudio[0:len(dataArray)+1], 'same')
    #print 'Corrdata', corrData
    #if len(corrData)>0:
    #	print "MAX",np.array(corrData).max()
    return corrData

#Using [1] as reference
def showtimespectrum(audio):
    #subplot(211)
    #plot(audio)
    title('Spectrogram of recorded audio')
    #subplot(212)    
    audioTime = float(len(audio)) / RATE
        
    fftData=abs(np.fft.rfft(audio))**2
    #print "LEN of FFT", len(fftData),"len input",len(audio)
    freqs=np.fft.rfftfreq(len(audio),d=1.0/RATE)
    print freqs[np.argmax(fftData)],np.argmax(fftData),actual_steps_bins,"MAX FREQ",freqs[-1]
    z=np.zeros(len(actual_steps))
    i=0
    for b in actual_steps_bins:
      z[i]=fftData[b]
      i+=1
    print z/z.sum()
    #Pxx, freqs, bins, im = spec]gram(audio, Fs = Rate, 
    #                                scale_by_freq=True, sides='default')
    #print len(audio),freqs.size,bins.size,Pxx.size
    #print freqs
    
    #xlim([0, audioTime])
    #corrResult = findechoes(audio)
    #subplot(212)
    #plt.plot(corrResult)
    #print "DONE"
    #pause(0.5)
    #plt.clf()
    #show()

def plotspectrum(i, line, ax, rStream):
    if rStream.is_active():
        recordData = []        
        for x in range(recordQueue.qsize()):
            recordData.extend(recordQueue.get())
            recordQueue.task_done()
            
            FFT = abs(scipy.fft(recordData))
            timestep = 1 / float(Rate)
            freqs = scipy.fftpack.fftfreq(len(recordData), timestep)
            line.set_xdata(freqs)
            line.set_ydata(20*scipy.log10(FFT))
            #ax.autoscale_view()
            ax.set_ylim(0, 200)
            ax.set_xlim(-Rate/2, Rate/2)
            #print 'DUUUDE. I am FFTING', freqs
    return line,

def showspectrum(rStream):
    #Pxx, freqs, bins, im = specgram([], Fs = Rate, scale_by_freq = True, sides = 'default')
    if rStream.is_active():
        recordData = []        
        for x in range(recordQueue.qsize()):
            recordData.extend(recordQueue.get())
            recordQueue.task_done()
            
            FFT = abs(scipy.fft(recordData))
            timestep = 1 / float(Rate)
            freqs = scipy.fftpack.fftfreq(len(recordData), timestep)
    return 
    #old code
    fig, ax = plt.subplots()
    line, = ax.plot([],[])
    aniF = animation.FuncAnimation(fig, plotspectrum, interval = 25, 
                                   blit = True, fargs = (line, ax, rStream))
    ax.clear()
    show()

def plotwave(i, line, ax, rStream):

    if rStream.is_active():
        recordData = []
        for x in range(recordQueue.qsize()):
            recordData.extend(recordQueue.get())
            recordQueue.task_done()

        old_x_data = line.get_xdata()
        old_y_data = line.get_ydata()
        old_sample_len = len(old_x_data)
        new_sample_len = len(recordData)
        #Length in seconds of old and new data                        
        old_len = old_sample_len / float(Rate)
        new_len = new_sample_len / float(Rate)
        
        #x-axis in seconds
        xmin, xmax = ax.get_xlim()
        
        #If old data is atleast the length of waterfall time then just
        #move the window forward. Otherwise keep adding data until atleast
        #window time
        if old_len >= WindowTime:
            ax.set_xlim(xmin+new_len, xmax+new_len)
            line.set_xdata(myfrange(xmin+new_len,xmax+new_len,old_sample_len))
            line.set_ydata(append(old_y_data[new_sample_len:], recordData))
        else:
            ax.set_xlim(xmin, xmax+new_len)
            #print 'Build window'
            line.set_xdata(myfrange(xmin, xmax+new_len, 
                                    old_sample_len + new_sample_len))
            line.set_ydata(append(old_y_data, recordData))
                  
    return line, 

def showwave(rStream):    
    #Pxx, freqs, t, fig = specgram([])
    #fig, ax = plt.specgram([1])
    fig, ax = plt.subplots()
    line, = ax.plot([], [])
    #Y Limit based on audio bit depth of 16 for now
    #
    ax.set_ylim(-1000, 1000)
    print "Brah"
    ani = animation.FuncAnimation(fig, plotwave, interval = 15, blit = True,
                                  fargs=(line, ax, rStream))
    ax.clear()
    plt.show()


#def play_ladder(p,l=1,steps=[260,500,1000,2000,3000,17000]): #,2000,4000]):
def play_ladder(p,l=3): #,2000,4000]):
    #actual_steps=[]
    #freqs=np.fft.rfftfreq(sample_window,d=1.0/RATE)
    #for f in steps:
    #  actual_steps.append(freqs[np.floor(f*sample_window*(1.0/RATE))])
    steps=actual_steps
    #r = 16000 #number of frames per second/frameset.      
    stream = p.open(#format = p.get_format_from_width(1), 
                    channels = 1, 
                    rate = RATE, 
		    format=pyaudio.paUInt8,
                    output = True)
    while True:
      for f in steps:
        n_frames = sample_window #int(RATE*l)
	#np.sin(2 * np.pi * frequency * samples)
        stream.write((np.sin(2*np.arange(n_frames)*np.pi*f/RATE)*127+128).astype(uint8))
        #d = []    
        #for x in xrange(n_frames): #TODO use numpy array
        # d.append(chr(int( math.sin(x*math.pi*f/r)*127+128)))
        #stream.write("".join(d))
    stream.stop_stream()
    stream.close()


def main():
    steps=[1000,3000,4000]
    freqs=np.fft.rfftfreq(sample_window,d=1.0/RATE)
    for f in steps:
      z=np.floor(f*sample_window*(1.0/RATE))
      actual_steps.append(freqs[z])
      actual_steps_bins.append(z)
    steps=actual_steps
	
    p = pyaudio.PyAudio()
    #Instantiate pyAudio
    #Open stream to record audio data from the mike
    #recordstream = p.open(format=p.get_format_from_width(SampWidth),
    recordstream = p.open(format=p.get_format_from_width(SampWidth),
                          #channels=Channels,
                          channels=1,
                          #rate=Rate,
                          rate=RATE,
                          input=True,
                          output=False,
                          stream_callback=recordcallback)
    #Create a thread to chirp continuously
    t1 = threading.Thread(target=play_ladder, args=(p,))
    t1.setDaemon(True)
    t1.start()
    #Create a "thread" to record continuosly
    recordstream.start_stream()    
    
    #Comment out one of the below two lines to show either the audio
    #waveform or the FFT spectrum
    #Keep both commented out if using the following for loop to show
    #timespectrum
    #showwave(recordstream)
    #showspectrum(reco1)
    
    #Old code to display fft spectrum versus time
    #Buggy and slow as hell but worth keeping
    while recordstream.is_active():                
        #Get the recorded data from the synchronized priority queue
        recordData = []
        #for i in range(recordQueue.qsize()):            
        #    recordData.extend(recordQueue.get())
        #    recordQueue.task_done()            
	#while len(recordData)!=8*4096:
	while len(recordData)!=sample_window:
	    #print(len(recordData))
            recordData.extend(recordQueue.get())
            recordQueue.task_done()            
	    #print(len(recordData),"@")
        #Display it
        showtimespectrum(recordData)
        #showtimewave(recordData)
        
    #Wait for the threads to end and cleanup
    t1.join()
        
    recordQueue.join()
    recordstream.stop_stream()
    recordstream.close()
    
    p.terminate()

if __name__ == '__main__':
    main()

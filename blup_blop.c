#include <complex.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <alsa/asoundlib.h>
#include <pthread.h>
#include <math.h>

#ifndef _STRUCT_TIMEVAL
#  define _STRUCT_TIMEVAL
#endif
#include <sys/time.h>

#include <stdlib.h>
#include <string.h>
#include <netdb.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <strings.h>

#define h_addr h_addr_list[0] /* for backward compatibility */

#define M_PI (3.14159265358979323846264338327950288)

long read_sample_window=256;
long fft_length=256;
double sensitivity=0.5;
struct timeval udp_time;
char * host = NULL;
char * device = NULL;
int portno = 8880;
#define NDELTAS 32
int ndeltas=0;
double deltas[NDELTAS];

snd_output_t *output = NULL;
int sampling_rate = 44100;

pthread_mutex_t lock_time = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lock_fftw = PTHREAD_MUTEX_INITIALIZER;

fftw_complex * signal_ext;
fftw_complex * out;

fftw_plan pa;

double target_freqA=0;
double target_freqB=0;
int target_freqA_bin=0;
int target_freqB_bin=0;
double * freq_bins=NULL;
double freq_step=0;

void init_fft() {
	signal_ext = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * fft_length);
	out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * fft_length);
	pa = fftw_plan_dft_1d(fft_length , signal_ext, out, FFTW_FORWARD, FFTW_ESTIMATE);
}
void close_fft() {
	fftw_destroy_plan(pa);
	fftw_free(signal_ext);
	fftw_free(out);
	fftw_cleanup();
}
void fft() {
	//only need to lock plan calls
	//the fftw_execute is thread-safe?
	fftw_execute(pa);
	return;
}

int16_t * add_tone(float freq, int n_frames, int16_t* buffer) {
	for (int i=0; i<n_frames; i++) {
		if (freq==0) {
			buffer[i]=0;
		} else {
			buffer[i]=sin(2*i*M_PI*freq/sampling_rate)*(pow(2,15));
		}
	}
	return buffer; 
}

int16_t * generate_tone(float freq, int n_frames) {
	int16_t * buffer = (int16_t *)malloc(sizeof(int16_t)*n_frames);
	if (buffer==NULL) {
		fprintf(stderr,"Failed to allocate space in generate_tone\n");
	}
	add_tone(freq,n_frames,buffer);
	return buffer; 
}

void ShortToReal(signed short* shrt,double* real,int siz) {
	int i;
	for(i = 0; i < siz; ++i) {
		real[i] = shrt[i]; // 32768.0;
	}
}

snd_pcm_t *capture_handle;
snd_pcm_hw_params_t *hw_params;
snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;

void init_record_audio(snd_pcm_t ** capture_handle, snd_pcm_hw_params_t ** hw_params , char * s) {
	fprintf(stderr,"Trying to open %s\n",s);
	int err;


	if ((err = snd_pcm_open (capture_handle, s, SND_PCM_STREAM_CAPTURE, 0)) < 0) {
		fprintf (stderr, "cannot open audio device(%s)\n", 
				snd_strerror (err));
		exit (1);
	}

	if ((err = snd_pcm_hw_params_malloc (hw_params)) < 0) {
		fprintf (stderr, "cannot allocate hardware parameter structure (%s)\n",
				snd_strerror (err));
		exit (1);
	}

	if ((err = snd_pcm_hw_params_any (*capture_handle, *hw_params)) < 0) {
		fprintf (stderr, "cannot initialize hardware parameter structure (%s)\n",
				snd_strerror (err));
		exit (1);
	}

	if ((err = snd_pcm_hw_params_set_access (*capture_handle, *hw_params, SND_PCM_ACCESS_RW_NONINTERLEAVED)) < 0) {
		fprintf (stderr, "cannot set access type (%s)\n",
				snd_strerror (err));
		exit (1);
	}

	if ((err = snd_pcm_hw_params_set_format (*capture_handle, *hw_params, format)) < 0) {
		fprintf (stderr, "cannot set sample format (%s)\n",
				snd_strerror (err));
		exit (1);
	}

	unsigned int rate=sampling_rate;	
	if ((err = snd_pcm_hw_params_set_rate_near (*capture_handle, *hw_params, &rate, 0)) < 0) {
		fprintf (stderr, "cannot set sample rate (%s)\n",
				snd_strerror (err));
		exit (1);
	}
	fprintf(stderr,"sampling_rate IS %u should be %d\n",rate,sampling_rate);
	assert(rate==sampling_rate);

	if ((err = snd_pcm_hw_params_set_channels (*capture_handle, *hw_params,1)) < 0) {
		fprintf (stderr, "cannot set channel count (%s)\n",
				snd_strerror (err));
		exit (1);
	}

	if ((err = snd_pcm_hw_params (*capture_handle, *hw_params)) < 0) {
		fprintf (stderr, "cannot set parameters (%s)\n",
				snd_strerror (err));
		exit (1);
	}

	snd_pcm_hw_params_free (*hw_params);

	if ((err = snd_pcm_prepare (*capture_handle)) < 0) {
		fprintf (stderr, "cannot prepare audio interface for use (%s)\n",
				snd_strerror (err));
		exit (1);
	}

	//try to lower latency
	snd_pcm_uframes_t buffer_size = 1024;
	snd_pcm_uframes_t period_size = 64;

	snd_pcm_hw_params_set_buffer_size_near (*capture_handle, *hw_params, &buffer_size); //how big the full buffe ris
	snd_pcm_hw_params_set_period_size_near (*capture_handle, *hw_params, &period_size, NULL); //this is how much is moved to CPU at a time
	fprintf(stderr,"REC %d %d\n",buffer_size,period_size);
}


void *thread_single_capture(void * x) {
	signed short * buffer = (signed short*)x;
	signed short * buffers[] = {buffer};
	int err;
	assert(capture_handle!=NULL);
	struct timeval start_capt;
	if(gettimeofday( &start_capt, 0 )) {
		fprintf(stderr,"FAILED TO get time...\n");
		exit(1);
	}

	if ((err = snd_pcm_readn (capture_handle, buffers, read_sample_window)) != read_sample_window) {
		fprintf (stderr, "read from audio interface failed (%s)\n",
				snd_strerror (err));
		exit (1);
	}
	struct timeval end_capt;
	if(gettimeofday( &end_capt, 0 )) {
		fprintf(stderr,"FAILED TO get time...\n");
		exit(1);
	}
	long delta_micro = 1000000 * (end_capt.tv_sec-start_capt.tv_sec) + (end_capt.tv_usec-start_capt.tv_usec); // - 1000000*((float)SAMPLES)/sampling_rate;
	fprintf(stderr,"%ld microseconds on overhead capture %ld\n",delta_micro, (1000000*read_sample_window)/sampling_rate );
	return NULL; 
}


int print_min() {
	int nsamples=ndeltas;
	if (nsamples>NDELTAS) {
		nsamples=NDELTAS;
	}
	double mn = deltas[0];
	for (int i=0; i<nsamples; i++) {
		if (mn>deltas[i]) {
			mn=deltas[i];
		}
	}	
	fprintf(stderr,"MIN %e\n",mn);
}

int print_avg() {
	//now copy over to complex buffer
	//compute mean and stddev
	int nsamples=ndeltas;
	if (nsamples>NDELTAS) {
		nsamples=NDELTAS;
	}
	double mean=0,var=0, stddev=0;
	for (int i=0; i<nsamples; i++) {
		mean+=deltas[i];
	}
	mean/=nsamples;
	for (int i=0; i<nsamples; i++) {
		var+=pow(deltas[i]-mean,2);
	}
	var/=nsamples;
	stddev=sqrt(var);
	fprintf(stderr,"MEAN %e STDDEV %e\n",mean,stddev);
	return 0;
}


void *thread_capture(void *threadarg) {
	//fprintf(stderr,"%0.3lf %d\n",target_freq,target_freq_bin);
	//lock the audio system
	signed short *buffer =  malloc(read_sample_window * snd_pcm_format_width(format) / 8 );
	fftw_complex * result = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * fft_length);
	fftw_complex * power_spectrum = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (fft_length/2+1));
	
	int *th_id = (int*)threadarg;

	while (1) {
		//get mutex for sound listening 
		pthread_mutex_lock(&lock);
		//pthread_create(threads+i,NULL,thread_single_capture,addys+2*i);
		thread_single_capture(buffer);
		pthread_mutex_unlock(&lock);

		//now copy over to complex buffer
		//compute mean and stddev
		double mean=0,var=0, stddev=0;
		for (int i=0; i<read_sample_window; i++) {
			mean+=buffer[i];
		}
		mean/=read_sample_window;
		for (int i=0; i<read_sample_window; i++) {
			buffer[i]-=mean;
		}
		for (int i=0; i<read_sample_window; i++) {
			var+=pow(buffer[i],2);
		}
		var/=read_sample_window;
		//fprintf(stderr,"MEAN %e STDDEV %e\n",mean,stddev);
		stddev=sqrt(var);

		for (int j=0; j<read_sample_window/fft_length; j++) {
			pthread_mutex_lock(&lock_fftw);
			for (int i=0; i<read_sample_window; i++) {
				signal_ext[i]=buffer[j*fft_length+i]/stddev;
			}
			fft();
			memcpy(result,out,sizeof(fftw_complex)*read_sample_window);
			pthread_mutex_unlock(&lock_fftw);


			power_spectrum[0] = abs(result[0]*result[0]); //DC component
			power_spectrum[read_sample_window/2] = abs(result[read_sample_window/2]*result[read_sample_window/2]);  /* Nyquist freq. */
			for (int k = 1; k <read_sample_window/2; ++k)  /* (k < N/2 rounded up) */
				power_spectrum[k] = abs(result[k]*result[k] + result[read_sample_window-k]*result[read_sample_window-k]);

			double normalize = 0;
			for (int k = 1; k<read_sample_window/2; ++k)
				normalize+=power_spectrum[k];
			//for (int k = target_freq_bin-5; k<target_freq_bin+5; ++k)
			//	normalize+=abs(power_spectrum[k]);
			//for (int i=0; i<=SAMPLES/2; i++) {
			//	fprintf(stderr,"%0.2fHz %0.3f,", freq_bins[i],creal(power_spectrum[i]/normalize));
			//}
			//fprintf(stderr,"%0.3lf %d %0.3lf\n",creal(power_spectrum[target_freq_bin]/normalize),target_freq_bin,normalize);
			//fprintf(stderr,"\n");
			//for (int i=0; i<SAMPLES; i++) {
			//	fprintf(stderr,"%0.3f %0.3f,", creal(result[i]), cimag(result[i]));
			//}

			/*int mxi=1;
			  float mx=creal(power_spectrum[mxi]);
			  for (int k=1; k<SAMPLES/2; k++) {
			  float v = creal(power_spectrum[k]);
			  if (v>mx) {
			  mxi=k;
			  mx=v;
			  }
			  }
			  fprintf(stderr, "MAX BIN %0.3f %0.3f %d\n",mx/normalize,freq_bins[mxi],mxi); */


			/*

			   Slide a quarter window across and try to find out the maximum values, then where
			   the maximum values start and taper off in front

			   if we cant find the start accurately because its too close to the front maybe we can save 
			   all blocks and have block IDs or maybe it doesnt matter significantly

			   OR

			   just have a bigger window and scan in quarter size sliding window,
			   that way will have a better chance of catching the signal start in one window that can be 			reasonably scanned

			   ALSO

			   measure time of start capture call and return time, to get an estimate on card delay
			   which could be variable!	

			   card_delay = (end - start) - samples/sampling_rate

			 */

			if (creal(power_spectrum[target_freqA_bin]/normalize)>sensitivity && udp_time.tv_usec>0) {
				int ret = pthread_mutex_trylock(&lock_time);
				if (ret==0) {
					struct timeval time;
					if(gettimeofday( &time, 0 )) {
						fprintf(stderr,"FAILED TO get time...\n");
						exit(1);
					}

					long delta_micro = 1000000 * (time.tv_sec-udp_time.tv_sec) + (time.tv_usec-udp_time.tv_usec) - 1000000*((float)read_sample_window)/sampling_rate;
					deltas[ndeltas++%NDELTAS]=delta_micro;
					//print_avg();
					print_min();
					//fprintf(stderr, "Delta %ld\n",delta_micro);
					udp_time.tv_usec=0;
					pthread_mutex_unlock(&lock_time);
				}
			}
		}	
	}
	return NULL;
}

void * udp_server(void * x) {
	struct sockaddr_in serveraddr; /* server's addr */
	struct sockaddr_in clientaddr; /* client addr */
	int buffer_length=512;
	char buf[buffer_length]; /* message buf */

	int sockfd = socket(AF_INET, SOCK_DGRAM, 0);
	if (sockfd < 0) {
		fprintf(stderr,"Failed to make socket\n");
		exit(1);
	} 
	int optval = 1;
	setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, 
			(const void *)&optval , sizeof(int));

	bzero((char *) &serveraddr, sizeof(serveraddr));
	serveraddr.sin_family = AF_INET;
	serveraddr.sin_addr.s_addr = htonl(INADDR_ANY);
	serveraddr.sin_port = htons((unsigned short)portno);

	if (bind(sockfd, (struct sockaddr *) &serveraddr, 
				sizeof(serveraddr)) < 0)  {
		fprintf(stderr,"Error on binding\n");
		exit(1);
	}

	int clientlen = sizeof(clientaddr);
	while (1) {
		bzero(buf, buffer_length);
		int n = recvfrom(sockfd, buf, buffer_length, 0,
				(struct sockaddr *) &clientaddr, &clientlen);
		if (n < 0) {
			fprintf(stderr,"Failed recvfrom\n");
			exit(1);
		}
	
		if(gettimeofday( &udp_time, 0 )) {
			fprintf(stderr,"Failed to get time\n");
			exit(1);
		}

		/*struct hostent * hostp = gethostbyaddr((const char *)&clientaddr.sin_addr.s_addr, 
				sizeof(clientaddr.sin_addr.s_addr), AF_INET);
		if (hostp == NULL) {
			fprintf(stderr,"Failed gethost\n");
			exit(1);
		}
		char * hostaddrp = inet_ntoa(clientaddr.sin_addr);
		if (hostaddrp == NULL) {
			fprintf(stderr,"Failed ntoa\n");
			exit(1);
		}
		printf("server received datagram from %s (%s)\n", 
				hostp->h_name, hostaddrp);
		printf("server received %d/%d bytes: %s\n", strlen(buf), n, buf);*/

		n = sendto(sockfd, buf, strlen(buf), 0, 
				(struct sockaddr *) &clientaddr, clientlen);
		if (n < 0)  {
			fprintf(stderr,"Error on sendto\n");
			exit(1);
		}
	}
}

int run_server() {
	init_fft();
	init_record_audio(&capture_handle,&hw_params,device); 

	//make a thread for UDP server
	pthread_t udp_thread; 
	pthread_create(&udp_thread,NULL,udp_server,NULL);


	//make the listen threads
	int nthreads = 4;
	int th_ids[nthreads];
	pthread_t threads[nthreads];
	for (int i=0; i<nthreads; i++) {
		th_ids[i]=i;
		pthread_create(threads+i,NULL,thread_capture,th_ids+i);
	}

	for (int i=0; i<nthreads; i++) {
		pthread_join(threads[i], NULL );
	}
	close_fft();
	snd_pcm_close (capture_handle); 
}

int run_client() {


	//open the socket
	int sockfd = socket(AF_INET, SOCK_DGRAM, 0);
	if (sockfd < 0) {
		fprintf(stderr,"failed to make sokcet\n");
		exit(0);
	}

	struct hostent *server = gethostbyname(host);
	if (server == NULL) {
		fprintf(stderr,"ERROR, no such host as %s\n", host);
		exit(0);
	}

	struct sockaddr_in serveraddr;
	bzero((char *) &serveraddr, sizeof(serveraddr));
	serveraddr.sin_family = AF_INET;
	bcopy((char *)server->h_addr, 
			(char *)&serveraddr.sin_addr.s_addr, server->h_length);
	serveraddr.sin_port = htons(portno);

	int buffer_length=512;
	char buf[buffer_length];
	bzero(buf, buffer_length);

	int serverlen = sizeof(serveraddr);

	//init the audio
	int err;
	snd_pcm_t *handle;
	if ((err = snd_pcm_open(&handle, device, SND_PCM_STREAM_PLAYBACK, 0)) < 0) {
		printf("Playback open error: %s\n", snd_strerror(err));
		exit(EXIT_FAILURE);
	}
	if ((err = snd_pcm_set_params(handle,
					//SND_PCM_FORMAT_U8,
					SND_PCM_FORMAT_S16,
					SND_PCM_ACCESS_RW_INTERLEAVED,
					1,
					sampling_rate,
					0,
					50000000)) < 0) {   /* 0.5sec */
					//0)) < 0) {   /* 0.5sec */
		printf("Playback open error: %s\n", snd_strerror(err));
		exit(EXIT_FAILURE);
	}
	snd_pcm_hw_params_t *hw_params;

	snd_pcm_hw_params_malloc (&hw_params);
	snd_pcm_hw_params_any (handle, hw_params);
	snd_pcm_hw_params_set_access (handle, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED);
	snd_pcm_hw_params_set_format (handle, hw_params, SND_PCM_FORMAT_S16);
	unsigned int rrate = sampling_rate;
	snd_pcm_hw_params_set_rate_near (handle, hw_params, &rrate, NULL);
	snd_pcm_hw_params_set_channels (handle, hw_params, 1);
	/* These values are pretty small, might be useful in
	   situations where latency is a dirty word. */
	long signal_period = sampling_rate/4; //0.5 signal period
	snd_pcm_uframes_t buffer_size = signal_period*1024;
	snd_pcm_uframes_t period_size = signal_period/2;
	
	fprintf(stderr,"PLAY buffer-size %ld, period-size %ld, signal-period %ld\n",buffer_size,period_size, signal_period);
	snd_pcm_hw_params_set_buffer_size_near (handle, hw_params, &buffer_size);
	snd_pcm_hw_params_set_period_size_near (handle, hw_params, &period_size, NULL);
	fprintf(stderr,"PLAY buffer-size %ld, period-size %ld, signal-period %ld\n",buffer_size,period_size, signal_period);
	signal_period = period_size; 
	fprintf(stderr,"PLAY buffer-size %ld, period-size %ld, signal-period %ld\n",buffer_size,period_size, signal_period);
	snd_pcm_sw_params_t *sw_params;
	snd_pcm_sw_params_malloc(&sw_params);
	snd_pcm_sw_params_current(handle, sw_params);
	//if (snd_pcm_sw_params_set_start_threshold(handle, sw_params, buffer_size - period_size)<0) {
	if (snd_pcm_sw_params_set_start_threshold(handle, sw_params, 0)<0) {
		fprintf(stderr,"Failed to set start threshold\n");
		exit(1);
	}
	//snd_pcm_sw_params_set_avail_min(handle, sw_params, period_size);


	int16_t * signal_buffer = (int16_t *)malloc(sizeof(int16_t)*signal_period);
	add_tone(target_freqA, fft_length, signal_buffer);
	add_tone(target_freqA, fft_length, signal_buffer+fft_length);
	add_tone(0, signal_period-fft_length*2, signal_buffer+fft_length*2);

	while (1) {
		//send packet
		int n = sendto(sockfd, buf, strlen(buf), 0, &serveraddr, serverlen);
		if (n < 0)  {
			fprintf(stderr,"Failed to send packets?\n");
			exit(1);
		}
		//play sound
		snd_pcm_sframes_t frames = snd_pcm_writei(handle, signal_buffer, signal_period);
		snd_pcm_start(handle);
		if (frames < 0)
			frames = snd_pcm_recover(handle, frames, 0);
		if (frames < 0) {
			printf("snd_pcm_writei failed: %s\n", snd_strerror(err));
		} else if (frames < signal_period) {
			printf("Short write (expected %li, wrote %li)\n", signal_period, frames);
		}
		snd_pcm_sframes_t delayp;
		snd_pcm_sframes_t availp;
		int ret= snd_pcm_avail_delay( handle, &availp, &delayp);
		//fprintf(stderr, "%d %ld %ld %ld\n",ret,delayp, availp, availp+delayp);
		//sleep for 
		//fprintf(stderr,"EMIT 1\n");
		struct timespec sleep_time,slept_time;
		sleep_time.tv_sec=signal_period/sampling_rate;
		sleep_time.tv_nsec=1e9*(((double)signal_period)/sampling_rate-sleep_time.tv_sec)/4;
		//fprintf(stderr,"Sleep for %ld %ld\n",sleep_time.tv_sec,sleep_time.tv_nsec);
		nanosleep(&sleep_time,&slept_time);
		//fprintf(stderr,"EMIT 2\n");
		//	long delta_micro = 1000000 * (time.tv_sec-udp_time.tv_sec) + (time.tv_usec-udp_time.tv_usec) ;
		//sleep(2);
	}

	free(signal_buffer);	

	snd_pcm_close(handle);

}

int main (int argc, char *argv[]) {
	if (argc!=9) {
		fprintf(stderr,"%s [0server/1client] DEV freq1 freq2 fft_length host port sensitivity\n",argv[0]);
		exit(1);
	}
	for (int i=0; i<NDELTAS; i++) {
		deltas[i]=0;
	}
	ndeltas=0;
	udp_time.tv_usec=0;
	int server=atoi(argv[1]);
	device = argv[2];
	float target_freqA_orig = atof(argv[3]);
	float target_freqB_orig = atof(argv[4]);
	fft_length=atoi(argv[5]);
	host=argv[6];
	portno=atoi(argv[7]);
	sensitivity=atof(argv[8]);

	freq_step = (sampling_rate/2)/(fft_length/2); // go up to the nyquist frequency, linearly spaced over SAMPLES/2+1
	freq_bins = (double*)malloc(sizeof(double)*(fft_length/2+1));
	for (int i=0; i<=fft_length/2; i++) {
		freq_bins[i]=i*freq_step;
	}
	target_freqA_bin = (int)(target_freqA_orig/freq_step);
	target_freqB_bin = (int)(target_freqB_orig/freq_step);
	target_freqA = target_freqA_bin*freq_step;
	target_freqB = target_freqB_bin*freq_step;
	fprintf(stderr,"Requested freqA for %0.3f, changed to %0.3f instead\n",target_freqA_orig,target_freqA);
	fprintf(stderr,"Requested freqB for %0.3f, changed to %0.3f instead\n",target_freqB_orig,target_freqB);

	if (server==0) {
		run_server();
	} else {
		run_client();
	}

	return 0; 
}

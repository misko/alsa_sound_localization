#include <complex.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <alsa/asoundlib.h>
#include <pthread.h>
#include <math.h>

#define MIC_DEVICE "plughw:%d,0"

#define RATE 44100
#define SAMPLES 64 //512 //4096 // ((int)(RATE*SECONDS)) 

pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lock_fftw = PTHREAD_MUTEX_INITIALIZER;

fftw_complex * signal_ext;
fftw_complex * out;

fftw_plan pa;

double listen_freq=0;
int listen_bin=0;
double * freq_bins=NULL;
double freq_step=0;

void init_fft() {
	signal_ext = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * SAMPLES);
	out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * SAMPLES);
	pa = fftw_plan_dft_1d(SAMPLES , signal_ext, out, FFTW_FORWARD, FFTW_ESTIMATE);
}
void close_fft() {
	fftw_destroy_plan(pa);
	fftw_free(signal_ext);
	fftw_free(out);
	fftw_cleanup();
}
void fft() {
	fftw_execute(pa);
	return;
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

void init_audio(snd_pcm_t ** capture_handle, snd_pcm_hw_params_t ** hw_params , char * s) {
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

	unsigned int rate=RATE;	
	if ((err = snd_pcm_hw_params_set_rate_near (*capture_handle, *hw_params, &rate, 0)) < 0) {
		fprintf (stderr, "cannot set sample rate (%s)\n",
				snd_strerror (err));
		exit (1);
	}
	fprintf(stderr,"RATE IS %u should be %d\n",rate,RATE);
	assert(rate==RATE);

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
}


void *thread_single_capture(void * x) {
	signed short * buffer = (signed short*)x;
	signed short * buffers[] = {buffer};
	int err;
	assert(capture_handle!=NULL);
	if ((err = snd_pcm_readn (capture_handle, buffers, SAMPLES)) != SAMPLES) {
		fprintf (stderr, "read from audio interface failed (%s)\n",
				snd_strerror (err));
		exit (1);
	}
	return NULL;
}


void *thread_capture(void *threadarg) {
	//lock the audio system
	signed short *buffer =  malloc(SAMPLES * snd_pcm_format_width(format) / 8 );
	fftw_complex * result = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * SAMPLES);
	fftw_complex * power_spectrum = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (SAMPLES/2+1));
	
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
		for (int i=0; i<SAMPLES; i++) {
			mean+=buffer[i];
		}
		mean/=SAMPLES;
		for (int i=0; i<SAMPLES; i++) {
			buffer[i]-=mean;
		}
		for (int i=0; i<SAMPLES; i++) {
			var+=pow(buffer[i],2);
		}
		var/=SAMPLES;
		//fprintf(stderr,"MEAN %e STDDEV %e\n",mean,stddev);
		stddev=sqrt(var);
		pthread_mutex_lock(&lock_fftw);
		for (int i=0; i<SAMPLES; i++) {
			signal_ext[i]=buffer[i]/stddev;
		}
		fft();
		memcpy(result,out,sizeof(fftw_complex)*SAMPLES);
		power_spectrum[0] = result[0]*result[0]; //DC component
		power_spectrum[SAMPLES/2] = result[SAMPLES/2]*result[SAMPLES/2];  /* Nyquist freq. */
		for (int k = 1; k <SAMPLES/2; ++k)  /* (k < N/2 rounded up) */
			power_spectrum[k] = result[k]*result[k] + result[SAMPLES-k]*result[SAMPLES-k];
		//for (int i=0; i<=SAMPLES/2; i++) {
		//	fprintf(stderr,"%0.2fHz %0.3f,", freq_bins[i],creal(power_spectrum[i]));
		//}
		//for (int i=0; i<SAMPLES; i++) {
		//	fprintf(stderr,"%0.3f %0.3f,", creal(result[i]), cimag(result[i]));
		//}
		fprintf(stderr,"%0.3f\n",power_spectrum[listen_bin]);
		pthread_mutex_unlock(&lock_fftw);
	}
	return NULL;
}

int main (int argc, char *argv[]) {
	if (argc!=3) {
		fprintf(stderr,"%s DEV listen_freq\n",argv[0]);
		exit(1);
	}
	char * dev = argv[1];
	float listen_freq_orig = atof(argv[2]);

	freq_step = RATE*1.0/SAMPLES;
	freq_bins = (double*)malloc(sizeof(double)*(SAMPLES/2+1));
	for (int i=0; i<=SAMPLES/2; i++) {
		freq_bins[i]=i*freq_step;
	}
	listen_bin = (int)(listen_freq_orig/freq_step);
	listen_freq = listen_bin*freq_step;
	fprintf(stderr,"Requested to listen for %0.3f, lisetning for %0.3f instead\n",listen_freq_orig,listen_freq);

	init_fft();
	init_audio(&capture_handle,&hw_params,dev); 

	int th_ids[] = {0,1};

	pthread_t threads[2];
	for (int i=0; i<2; i++) {
		pthread_create(threads+i,NULL,thread_capture,th_ids+i);
	}

	for (int i=0; i<2; i++) {
		pthread_join(threads[i], NULL );
	}
	close_fft();
	snd_pcm_close (capture_handle); 


	return 0; 
}

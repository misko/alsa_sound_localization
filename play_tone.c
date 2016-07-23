/*
 *  This extra small demo sends a random samples to your speakers.
 */

#include <alsa/asoundlib.h>
#include <math.h>
#define M_PI (3.14159265358979323846264338327950288)

static char *device = "default";                        /* playback device */

snd_output_t *output = NULL;

int sampling_rate = 48000;

int16_t * generate_tone(float freq, int n_frames) {
	int16_t * buffer = (int16_t *)malloc(sizeof(int16_t)*n_frames);
	for (int i=0; i<n_frames; i++) {
		//(np.sin(2*np.arange(n_frames)*np.pi*f/RATE)*32000).astype(np.int16)
		buffer[i]=sin(2*i*M_PI*freq/sampling_rate)*(pow(2,15));
	}
	return buffer; 
}

int main(int argc, char ** argv) {
	if (argc!=2) {
		fprintf(stderr,"%s frequency\n",argv[0]);
		exit(1);
	}
	float freq = atof(argv[1]);

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
					1,
					500000)) < 0) {   /* 0.5sec */
		printf("Playback open error: %s\n", snd_strerror(err));
		exit(EXIT_FAILURE);
	}


	int n_samples = sampling_rate*3;
	int16_t * tone = generate_tone(freq,n_samples);

	snd_pcm_sframes_t frames = snd_pcm_writei(handle, tone, n_samples);
	if (frames < 0)
		frames = snd_pcm_recover(handle, frames, 0);
	if (frames < 0) {
		printf("snd_pcm_writei failed: %s\n", snd_strerror(err));
	} else if (frames < n_samples) {
		printf("Short write (expected %li, wrote %li)\n", n_samples, frames);
	}

	free(tone);	

	snd_pcm_close(handle);
	return 0;
}

all: play_tone get_tone


play_tone: play_tone.c
	gcc -std=c99 -lfftw3 -lm -lasound -lpthread play_tone.c -o play_tone

get_tone: get_tone.c
	gcc -std=c99 -lfftw3 -lm -lasound -lpthread get_tone.c -o get_tone

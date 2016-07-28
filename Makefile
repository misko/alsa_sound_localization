all: play_tone get_tone blup_blop blippity


play_tone: play_tone.c
	gcc -std=c99 -lfftw3 -lm -lasound -lpthread play_tone.c -o play_tone

get_tone: get_tone.c
	gcc -std=c99 -lfftw3 -lm -lasound -lpthread get_tone.c -o get_tone

blup_blop: blup_blop.c
	gcc -std=c99 -lfftw3 -lm -lasound -lpthread blup_blop.c -o blup_blop

blippity: blippity.c
	gcc -std=c99 -lfftw3 -lm -lasound -lpthread blippity.c -o blippity

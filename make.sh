python gen_tone4.py 48000 3 500 1000
python gen_tone4.py 48000 3 19000 21000
python gen_tone4.py 48000 3 18000 21000
python gen_tone5.py 48000 3 19000 21000
ffmpeg -i 19k.21k_48k.wav -codec:a libmp3lame -q:a 0 cat.mp3
ffmpeg -i 16k.16k_48k.wav -codec:a libmp3lame -q:a 0 cat2.mp3
ffmpeg -i 17k.17k_48k.wav -codec:a libmp3lame -q:a 0 cat3.mp3
ffmpeg -i 18k.21k_48k.wav -codec:a libmp3lame -q:a 0 hi.mp3
ffmpeg -i 0k.1k_48k.wav  human.mp3
#python gen_tone2.py 192000 6 20000
#python gen_tone2.py 192000 6 30000
#python gen_tone2.py 192000 6 40000
#python gen_tone2.py 192000 6 50000
#python gen_tone2.py 192000 6 60000

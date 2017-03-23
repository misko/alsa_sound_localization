#!/bin/bash
if [ $# -ne 1 ]; then
	echo $0 wav_filename
	exit
fi
wav_filename=$1
filename="${wav_filename%.*}"

ffmpeg -i ${wav_filename} -codec:a libmp3lame -q:a 0 ${filename}.mp3

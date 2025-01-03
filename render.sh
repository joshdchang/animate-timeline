cd text-frames && ffmpeg -r 30 -f image2 -s 2400x1600 -start_number 1 -i frame-%d.png -vframes 8784 -vcodec libx264 -crf 25 -pix_fmt yuv420p ../movie.mp4

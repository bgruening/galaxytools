ffmpeg -ss 00:00:02 -to 00:00:04 -i SampleVideo_360x240_1mb.mp4 -c copy SampleVideo_360x240_cut.mp4

ffmpeg -i SampleVideo_360x240_cut.mp4 SampleVideo_360x240_cut.webm
ffmpeg -i SampleVideo_360x240_cut.mp4 SampleVideo_360x240_cut.flac
ffmpeg -i SampleVideo_360x240_cut.mp4 SampleVideo_360x240_cut.mp3
ffmpeg -i SampleVideo_360x240_cut.flac flac_to_ogg.ogg
ffmpeg -i SampleVideo_360x240_cut.mp4 SampleVideo_360x240_cut.mkv

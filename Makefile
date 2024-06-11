timer:
	g++ gfft_parallel.cpp gfft.cpp gfft_inplace.cpp timer.cpp utils.cpp parallel_natali.cpp FFT.cpp radix_parallel.cpp -o timer

run_parallel:
	g++ gfft_parallel.cpp gfft.cpp gfft_inplace.cpp run_parallel.cpp utils.cpp -o run_parallel

test:
	g++ gfft_parallel.cpp FFT.cpp radix_parallel.cpp gfft.cpp gfft_inplace.cpp test.cpp utils.cpp -o test

weather:
	g++ gfft_parallel.cpp gfft_inplace.cpp utils.cpp weather.cpp -o weather

weather_2radix:
	g++ fft.cpp gfft_parallel.cpp gfft_inplace.cpp utils.cpp weather_2_radix.cpp -o weather_2radix
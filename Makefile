timer:
	g++ general/gfft_parallel.cpp general/gfft.cpp general/gfft_inplace.cpp timer.cpp utils.cpp radix/FFT.cpp radix/radix_parallel.cpp -o timer

run_parallel:
	g++ general/gfft_parallel.cpp general/gfft.cpp general/gfft_inplace.cpp run_parallel.cpp utils.cpp -o run_parallel

test:
	g++ general/gfft_parallel.cpp radix/FFT.cpp radix/radix_parallel.cpp general/gfft.cpp general/gfft_inplace.cpp test.cpp utils.cpp -o test

weather:
	g++ general/gfft_parallel.cpp general/gfft_inplace.cpp utils.cpp weather.cpp -o weather

weather_2radix:
	g++ radix/FFT.cpp  utils.cpp weather_2_radix.cpp -o weather_2radix

all:
	make timer
	make run_parallel
	make test
	make weather
	make weather_2radix

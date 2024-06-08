timer:
	g++ gfft_parallel.cpp gfft.cpp gfft_inplace.cpp timer.cpp utils.cpp -o timer

run_parallel:
	g++ gfft_parallel.cpp gfft.cpp gfft_inplace.cpp run_parallel.cpp utils.cpp -o run_parallel
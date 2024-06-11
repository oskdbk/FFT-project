timer:
	g++ gfft_parallel.cpp gfft.cpp gfft_inplace.cpp timer.cpp utils.cpp parallel_natali.cpp -o timer

run_parallel:
	g++ gfft_parallel.cpp gfft.cpp gfft_inplace.cpp run_parallel.cpp utils.cpp -o run_parallel

test:
	g++ gfft_parallel.cpp gfft.cpp gfft_inplace.cpp test.cpp utils.cpp -o test
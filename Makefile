all:
	g++ FFT.cpp -o fft

general:
	g++ General_FFT_algo.cpp -o general

parallel:
	g++ gfft_parallel.cpp -o parallel
drawing:
	g++ -o application application.cpp ../general/gfft_parallel.cpp ../general/gfft.cpp ../general/gfft_inplace.cpp ../utils.cpp `pkg-config --cflags --libs opencv4` -lsfml-graphics -lsfml-window -lsfml-system

parallel:
	g++ -o application_parallel application_parallel.cpp ../general/gfft_parallel.cpp ../general/gfft.cpp ../general/gfft_inplace.cpp ../utils.cpp `pkg-config --cflags --libs opencv4` -lsfml-graphics -lsfml-window -lsfml-system

drawing:
	g++ -o application application.cpp ../gfft_parallel.cpp ../gfft.cpp ../gfft_inplace.cpp ../utils.cpp `pkg-config --cflags --libs opencv4` -lsfml-graphics -lsfml-window -lsfml-system

parallel:
	g++ -o application_parallel application_parallel.cpp ../gfft_parallel.cpp ../gfft.cpp ../gfft_inplace.cpp ../utils.cpp `pkg-config --cflags --libs opencv4` -lsfml-graphics -lsfml-window -lsfml-system

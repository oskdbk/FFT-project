# Application
## How to run
You can compile the code for the non-parallel version by ```make drawing``` or just ```make``` and then run by ```./application``` from inside the ```drawing``` folder.
```
cd drawing
make drawing
./application
```

You can compile the code for the parallel version by ```make parallel``` and then run by ```./application_parallel``` from inside the ```drawing``` folder.
```
cd drawing
make parallel
./application_parallel
```

Dependencies: code relies on OpenCV and SFML. You can get it from terminal on Linux (Ubuntu/Debian) by:
```
sudo apt install libopencv-dev
sudo apt install libsfml-dev
```

<div align="center">
  
## Image of PI
Number of contours: 3, Number of points: 581
#### Initial image
<img src="drawing/pictures/pi_try.png" alt="PI" style="height:200px;">

#### Result image

<img src="drawing/results/pi_result.png" alt="PI result" style="height:200px;">

## Image of Psyduck
Number of contours: 20, Number of points: 4825
#### Initial image
<img src="drawing/pictures/example.png" alt="Psyduck" style="height:200px;">

#### Result image
With connecting lines

<img src="drawing/results/psyduck_result_np.png" alt="Psyduck result" style="height:200px;">

Without connecting lines

<img src="drawing/results/psyduck.png" alt="Psyduck result" style="height:200px;">

## Image of Gojo
Number of contours: 123, Number of points: 2844
#### Initial image
<img src="drawing/pictures/gojo.png" alt="Gojo" style="height:200px;">

#### Result image
With connecting lines

<img src="drawing/results/gojo_result_np.png" alt="Gojo result" style="height:200px;">


Without connecting lines

<img src="drawing/results/gojo_result.png" alt="Gojo result" style="height:200px;">
</div>

#include "general/gfft.h"
#include "general/gfft_parallel.h"
#include "general/gfft_inplace.h"
#include "utils.h"
#include "radix/FFT.h"
#include "radix/radix_parallel.h"
#include <vector>

using namespace std;


#ifndef NO_MAIN
int main(int argc, char* argv[]){
    // This just generates a vector of length 15 of (1, 0)
    // vector<complex_num> P(15, complex_num(1, 0));
    vector<complex_num> P;
    int iters = 5;
    int max_threads = 8;
    if (argc == 1){
        P = Read_CSV("weather_data/Weather_data.csv");
    } else {
        int n = stoi(argv[1]);
        P =  std::vector<complex_num>(40320, complex_num(1.0, 0.0));
    }
    if (argc > 2){
        max_threads = stoi(argv[2]);
    }
    if (argc > 3){
        iters = stoi(argv[3]);
    }
    int n = P.size();
    
    //cout << "Running with n = " << n << ", max_threads = " << max_threads << " and iters = " << iters << endl;
    

    vector<int> lengths = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900};
    int step = 1000;
    int min_length = 1000;
    int max_length = 40000;
    for (int i = min_length; i <= max_length; i+=step){
        lengths.push_back(i);
    }
    iters = 10;
    // print the lengths vector in the format of a python list
    cout << "lengths = []";
    for (int i = 0; i < lengths.size(); i++){
        cout << lengths[i];
        if (i < lengths.size() - 1)
            cout << ", ";
    }
    cout << "];" << endl;

    return 0;



    vector<int> threads = {2, 4, 8, 12, 16, 32};
    vector<vector<double>> times;
    vector<string> names;
    vector<vector<double>> times_radix;
    vector<string> names_radix;
    for (int i = 0; i < lengths.size(); i++){
        // flush cout
        cout  << "\r" << "Running with n = " << lengths[i] << " (" << i +1 << "/" << lengths.size() << ")";
        cout.flush();
        vector<complex_num> P(lengths[i], complex_num(1.0, 0.0));
        vector<double> row;
        row.push_back(measure_time(GeneralFFT, iters, P, false));
        if (i == 0)
            names.push_back("Sequential");
        row.push_back(measure_time(GeneralFFT_inplace, iters, P, false, false));
        if (i == 0)
            names.push_back("Inplace");
        for (int j = 0; j < threads.size(); j++){
            double a = measure_time(GeneralFFT_Parallel, iters, P, threads[j], false);
            row.push_back(a);
            if (i == 0)
                names.push_back("Parallel_" + to_string(threads[j]));
        }
        times.push_back(row);

        row.clear();
        row.push_back(measure_time(Radix2FFT, iters, P));
        if (i == 0)
            names_radix.push_back("Radix2");
        
        for (int j = 0; j < threads.size(); j++){
            double a = measure_time(Radix2FFT_parallel, iters, P, threads[j], false);
            row.push_back(a);
            if (i == 0)
                names_radix.push_back("Radix2_Parallel_" + to_string(threads[j]));
        }
        times_radix.push_back(row);
    }
    // reverse rows and columns in times
    vector<vector<double>> reversed_times(times[0].size(), vector<double>(times.size()));
    for (int i = 0; i < times.size(); i++) {
        for (int j = 0; j < times[i].size(); j++) {
            reversed_times[j][i] = times[i][j];
        }
    }
    times = reversed_times;

    // reverse rows and columns in times_radix
    vector<vector<double>> reversed_times_radix(times_radix[0].size(), vector<double>(times_radix.size()));
    for (int i = 0; i < times_radix.size(); i++) {
        for (int j = 0; j < times_radix[i].size(); j++) {
            reversed_times_radix[j][i] = times_radix[i][j];
        }
    }
    times_radix = reversed_times_radix;


    // save to csv
    Write_CSV_Columns(times, names, "timing/general_times.csv");
    Write_CSV_Columns(times_radix, names_radix, "timing/radix_times.csv");
    return 0;


}
#endif
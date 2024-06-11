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
    int iters = 25;
    int max_threads = 8;
    if (argc == 1){
        P = Read_CSV("Weather_data.csv");
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
    

    vector<int> lengths = {100, 1000, 10000, 25000, 50000, 100000};
    vector<int> threads = {1, 2, 4, 8, 12, 16, 32};
    vector<vector<double>> times;
    vector<string> names;
    vector<vector<double>> times_radix;
    vector<string> names_radix;
    for (int i = 0; i < lengths.size(); i++){
        cout << "Running with n = " << lengths[i] << endl;
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


    // print the table
    cout << "Length\t";
    for (int i = 0; i < names.size(); i++){
        cout << names[i] << "\t";
    }
    cout << endl;
    for (int i = 0; i < lengths.size(); i++){
        cout << lengths[i] << "\t";
        for (int j = 0; j < times[i].size(); j++){
            cout << times[i][j] << "\t";
        }
        cout << endl;
    }


    // print the table for radix
    cout << "Length\t";
    for (int i = 0; i < names_radix.size(); i++){
        cout << names_radix[i] << "\t";
    }
    cout << endl;
    for (int i = 0; i < lengths.size(); i++){
        cout << lengths[i] << "\t";
        for (int j = 0; j < times_radix[i].size(); j++){
            cout << times_radix[i][j] << "\t";
        }
        cout << endl;
    }
    return 0;


}
#endif
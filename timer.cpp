#include "gfft.h"
#include "gfft_parallel.h"
#include "gfft_inplace.h"
#include "utils.h"
#include "parallel_natali.h"
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
    for (int i = 0; i < lengths.size(); i++){
        cout << "Running with n = " << lengths[i] << endl;
        vector<complex_num> P(lengths[i], complex_num(1.0, 0.0));
        vector<double> row;
        row.push_back(measure_time(GeneralFFT, iters, P, false));
        if (i == 0)
            names.push_back("Sequential");
        row.push_back(measure_time(GeneralFFT_inplace, iters, P, false));
        if (i == 0)
            names.push_back("Inplace");
        for (int j = 0; j < threads.size(); j++){
            double a = measure_time(GeneralFFT_Parallel, iters, P, threads[j]);
            row.push_back(a);
            if (i == 0)
                names.push_back("Parallel_" + to_string(threads[j]));
        }
        row.push_back(measure_time(FFT, iters, P));
        if (i == 0)
            names.push_back("Natali");
        times.push_back(row);
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
    return 0;


    // measure time of both versions

    // sequential version
    clock_t start, end;
    double avg_time = 0;
    for (int i = 0; i < iters; i++)
    {
        start = clock();
        vector<complex_num> P_star = GeneralFFT(P, false);
        end = clock();
        avg_time += double(end - start);
    }
    avg_time /= iters;

    cout << "Sequential version took " <<  avg_time << " ticks" << endl;

    // inplace version
    avg_time = 0;
    for (int i = 0; i < iters; i++)
    {
        start = clock();
        vector<complex_num> P_star_inplace = GeneralFFT_inplace(P, false);
        end = clock();
        avg_time += double(end - start);
    }
    avg_time /= iters;

    cout << "Inplace version took " << avg_time << " ticks" << endl;
    
    // try with different number of threads
    for(int i = 1; i <= max_threads; i++){
        avg_time = 0;
        for(int j = 0; j < iters; j++){
            start = clock();
            vector<complex_num> P_star_parallel = GeneralFFT_Parallel(P, i);
            end = clock();
            avg_time += double(end - start);
        }
        avg_time /= iters;
        cout << "Parallel version with " << i << " threads took " << avg_time << " ticks" << endl;
    }

    // Natali version
    avg_time = 0;
    for (int i = 0; i < iters; i++)
    {
        start = clock();
        vector<complex_num> P_star_natali = FFT(P);
        end = clock();
        avg_time += double(end - start);
    }
    avg_time /= iters;

    cout << "Natali version took " << avg_time << " ticks" << endl;
    
    vector<complex_num> P_star = GeneralFFT(P, false);
    vector<complex_num> P_star_parallel = GeneralFFT_Parallel(P, max_threads);
    vector<complex_num> P_star_inplace = GeneralFFT_inplace(P, false);
    vector<complex_num> P_star_natali = FFT(P);
    is_same_vector(P_star, P_star_parallel, "Sequential", "Parallel", true);
    is_same_vector(P_star, P_star_inplace, "Sequential", "Inplace", true);
    is_same_vector(P_star, P_star_natali, "Sequential", "Natali", true);
    return 0;
}
#endif
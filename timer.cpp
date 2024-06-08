#include "gfft.h"
#include "gfft_parallel.h"
#include "gfft_inplace.h"
#include "utils.h"
#include <vector>

using namespace std;


#ifndef NO_MAIN
int main(int argc, char* argv[]){
    // This just generates a vector of length 15 of (1, 0)
    // vector<complex_num> P(15, complex_num(1, 0));
    vector<complex_num> P;
    int iters = 5;
    int max_threads = 24;
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
    
    cout << "Running with n = " << n << ", max_threads = " << max_threads << " and iters = " << iters << endl;
      

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


    vector<complex_num> P_star = GeneralFFT(P, false);
    vector<complex_num> P_star_parallel = GeneralFFT_Parallel(P, max_threads);
    vector<complex_num> P_star_inplace = GeneralFFT_inplace(P, false);
    size_t t = P_star.size();
    for(int i = 0; i < t; i++){
        if(norm(P_star[i] - P_star_parallel[i]) > 1e-6){
            cout << "Parallel differs at index " << i << endl;
            cout << P_star[i] << " != " << P_star_parallel[i] << endl;
            break;
        } 
    }
    cout << "Parallel version is the same as sequential" << endl;

    return 0;
}
#endif
#include <vector>
#include "gfft_parallel.h"
#include "utils.h"

using namespace std;

int main(int argc, char* argv[]){
    vector<complex_num> P;
    int num_threads = 24;
    if (argc == 1){
        cout << "No input file provided. Using default file Weather_data.csv" << endl;
        P = Read_CSV("Weather_data.csv");
    } else {
        char* filename = argv[1];
        P = Read_CSV(filename);
    }
    if (argc > 2){
        num_threads = stoi(argv[2]);
    }

    int n = P.size();
    cout << "Running with n = " << n << " and num_threads = " << num_threads << endl;

    // measure time of parallel version
    clock_t start, end;
    start = clock();
    vector<complex_num> P_star_parallel = GeneralFFT_Parallel(P, num_threads);
    end = clock();
    double time = double(end - start);
    cout << "Parallel version took " << time << " ticks (" << time / CLOCKS_PER_SEC << "s)" << endl;
    // PRT1(P_star_parallel, "Result:");

}
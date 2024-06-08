#include "gfft.h"
#include "gfft_parallel.h"
#include "gfft_inplace.h"
#include <vector>

using namespace std;


#ifndef NO_MAIN
int main(){
    // This just generates a vector of length 15 of (1, 0)
    // vector<complex_num> P(15, complex_num(1, 0));

    vector<complex_num> P = Read_CSV("Weather_data.csv");
    int n = P.size();
    // copy P to Q
    vector<complex_num> Q(P.size());
    for(int i = 0; i < P.size(); i++){
        Q[i] = P[i];
    }
    // cout<<P.size()<<endl;
    

    // measure time of both versions

    // sequential version
    clock_t start, end;
    double avg_time = 0;
    for (int i = 0; i < 30; i++)
    {
        start = clock();
        vector<complex_num> P_star = GeneralFFT(P, false);
        end = clock();
        avg_time += double(end - start);
    }
    avg_time /= 30;

    cout << "Sequential version took " <<  avg_time << " ticks" << endl;

    // inplace version
    avg_time = 0;
    for (int i = 0; i < 30; i++)
    {
        start = clock();
        vector<complex_num> P_star_inplace = GeneralFFT_inplace(Q, false);
        end = clock();
        avg_time += double(end - start);
    }
    avg_time /= 30;

    cout << "Inplace version took " << avg_time << " ticks" << endl;
    
    // try with different number of threads
    for(int i = 1; i <= 30; i++){
        if (n % i != 0){
            continue;
        }
        avg_time = 0;
        for(int j = 0; j < 30; j++){
            start = clock();
            vector<complex_num> P_star_parallel = GeneralFFT_Parallel(Q, i);
            end = clock();
            avg_time += double(end - start);
        }
        avg_time /= 30;
        cout << "Parallel version with " << i << " threads took " << avg_time << " ticks" << endl;
    }


    vector<complex_num> P_star = GeneralFFT(P, false);
    vector<complex_num> P_star_parallel = GeneralFFT_Parallel(Q, 20);
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
#include "utils.h"
#include "general/gfft_parallel.h"
#include <vector>

using namespace std;

int main(){
    vector<complex_num> P = Read_CSV("Weather_data.csv");
    vector<complex_num> Res = GeneralFFT_Parallel(P, 8, false);
    int keep = (double) Res.size() * 0.005;
    vector<int> keeps{1,2,4,8, 16};
    vector<string> headers;
    vector<vector<complex_num>> primes;

    for (int keep: keeps)
    {
        vector<complex_num> Res_copy = Res;
        keep_largest_n(Res_copy, keep);
        primes.push_back(GeneralFFT_Parallel(Res_copy, 8, true));
        headers.push_back(to_string(keep));
    }
    Write_CSV_Columns(primes, headers, "weather_data/weather_cleaned_multiple.csv");

    
}
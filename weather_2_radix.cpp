#include "utils.h"
#include "gfft_parallel.h"
#include <vector>
#include "FFT.h"

using namespace std;

int main(){
    vector<complex_num> P = Read_CSV("Weather_data.csv");
    int n = P.size();
    while(P.size()< 2048){
        P.push_back(complex_num(0.0,0.0));
    }
    vector<complex_num> P_star = Radix2FFT(P);
    vector<int> keeps{1,2,4,8,16};
    vector<string> headers;
    vector<vector<complex_num>> primes;

    for (int keep: keeps)
    {
        vector<complex_num> Res_copy = P_star;
        keep_largest_n(Res_copy, keep);
        vector<complex_num> prime = InverseRadix2FFT(Res_copy);
        while(prime.size() > n){
            prime.pop_back();
        }
        primes.push_back(prime);
        headers.push_back(to_string(keep));
    }
    Write_CSV_Columns(primes, headers, "weather_cleaned_radix.csv");

    
}
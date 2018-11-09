#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
#include <numeric>
#include <cmath>

#define func auto

func nChoosek(uint64_t n, uint64_t k) -> uint64_t
{
    if (k > n)
        return 0ULL;

    if (k * 2ULL > n)
        k = n - k;

    if (k == 0ULL)
        return 1;

    uint64_t result = n;
    uint64_t i = 1ULL;
    while( i != k ) {
        result *= (n - i);
        ++i;
        result /= i;
    }
    return result;
}


func nCk(uint64_t n, uint64_t k) -> double { // computes binomial coefficient
    if (k > n)
        return 0.0;

    if (k * 2ULL > n)
        k = n - k;

    if (k == 0.0)
        return 1.0;

    return exp(lgamma(n+1) - (lgamma(k+1) + lgamma(n-k+1)));
}


func PMF(uint64_t N, uint64_t K, uint64_t n, uint64_t k) -> double { // Hypergeometric distribution PMF
    return nCk(K, k) * nCk(N - K, n - k) / nCk(N, n);
}


func sample_hyper(int N, int K, int n,            // sampling from Hypergeometric distribution
                  std::mt19937_64&  rng,
                  std::vector<int>& nums) -> int {

    int rc{ 0 };
    for (int k = 0; k != n; ++k) {
        std::uniform_int_distribution<int> uni{ k, N-1 };
        auto s = uni(rng);
        std::swap(nums[k], nums[s]);

        rc += (nums[k] < K);
    }
    return rc;
}


func main() -> int {
    auto rng = std::mt19937_64{1234567891ULL};

    int N = 500; // compare with Wiki
    int K = 50;
    int n = 100;

    auto nums = std::vector<int>( N );
    std::iota(nums.begin(), nums.end(), 0); // array to shuffle, filled with 0, 1, 2, 3 ... sequence

    auto h = std::vector<int>( K, 0 ); // histogram

    int NT = 100000; // number of trials
    for(int k = 0; k != NT; ++k) {
        int q = sample_hyper(N, K, n, rng, nums);
        h[q] += 1; // count it
    }

    for (int k = 0; k != 20; ++k) { // and print it
        double e = double(h[k]) / double(NT);
        std::cout << k << "   " << e << "     " << PMF(N, K, n, k) << '\n';
    }

    std::cin >> n;
    return 0;
}

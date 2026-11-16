// Compute an integer upper bound on the reciprocal sum for Section 4 of
// Ochem and Zelinsky (arXiv:2607.19746), using their exclusions below 3,000,000.
// Build: g++ -std=c++17 -O3 -pthread QBoundPrefix.cpp -o QBoundPrefix
// Usage: QBoundPrefix [upper [threads [lower [residue-limit]]]]
// The second output is the upper bound multiplied by 10^18.
#include <algorithm>
#include <atomic>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

using u64 = std::uint64_t;
using u128 = __uint128_t;

u64 power_mod(u64 base, int exponent, int modulus) {
    u64 result = 1;
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result = result * base % modulus;
        }
        base = base * base % modulus;
        exponent /= 2;
    }
    return result;
}

int main(int argc, char** argv) {
    constexpr u64 seed_limit = 3000000;
    constexpr int block_size = 1 << 19;
    constexpr u64 scale = 1000000000000000000ULL;
    // This keeps the prime sieve index below 10^7 and all address arithmetic
    // and modular products within 64 bits.
    constexpr u64 upper_limit = 100000000000000ULL;

    u64 lower = seed_limit;
    u64 upper = 10000000000;
    u64 workers = 1;
    u64 residue_limit = 100000;
    u64* arguments[] = {&upper, &workers, &lower, &residue_limit};
    bool valid = argc <= 5;
    for (int i = 1; valid && i < argc; ++i) {
        const std::string_view argument(argv[i]);
        const auto end = argument.data() + argument.size();
        const auto parsed = std::from_chars(argument.data(), end, *arguments[i - 1]);
        valid = parsed.ec == std::errc{} && parsed.ptr == end;
    }
    if (!valid || lower < seed_limit || upper <= lower ||
        upper > upper_limit || workers < 1 || workers > 128 ||
        residue_limit > seed_limit) {
        std::cerr << "Usage: QBoundPrefix [upper [threads [lower [residue-limit]]]]\n";
        return 1;
    }
    int prime_limit = static_cast<int>(residue_limit);
    while (u64(prime_limit) * prime_limit < upper) {
        ++prime_limit;
    }

    std::vector<bool> composite(prime_limit + 1);
    std::vector<int> odd_primes;
    std::vector<std::pair<int, int>> forbidden_residues;
    for (int p = 3; p <= prime_limit; p += 2) {
        if (composite[p]) {
            continue;
        }
        if (u64(p) * p <= upper) {
            odd_primes.push_back(p);
        }
        if (p <= prime_limit / p) {
            for (int n = p * p; n <= prime_limit; n += 2 * p) {
                composite[n] = true;
            }
        }
        if (u64(p) <= residue_limit && p % 3 == 1) {
            int root = 1;
            for (int base = 2; root == 1; ++base) {
                root = static_cast<int>(power_mod(base, (p - 1) / 3, p));
            }
            forbidden_residues.emplace_back(p, root);
            forbidden_residues.emplace_back(p, p - 1 - root);
        }
    }

    constexpr u64 block_width = 2 * block_size;
    const u64 first_odd = (lower + 1) | 1;
    std::atomic<u64> next_block{first_odd};
    struct Total {
        u64 survivors = 0;
        u128 reciprocal_upper = 0;
    };
    std::vector<Total> totals(workers);
    auto sieve_blocks = [&](Total& total) {
        std::vector<std::uint8_t> marked(block_size);
        for (u64 lo = next_block.fetch_add(block_width, std::memory_order_relaxed);
             lo <= upper;
             lo = next_block.fetch_add(block_width, std::memory_order_relaxed)) {
            const auto size = static_cast<std::size_t>(
                std::min<u64>(block_size, (upper - lo) / 2 + 1));
            const u64 hi = lo + 2 * (size - 1);
            std::fill_n(marked.begin(), size, 0);

            // Odd multiples within the block advance by 2*p in value, p in index.
            for (int p : odd_primes) {
                const u64 square = u64(p) * p;
                if (square > hi) {
                    break;
                }
                const u64 first = std::max(square, ((lo + p - 1) / p) * p);
                const u64 first_odd = first % 2 == 1 ? first : first + p;
                for (u64 i = (first_odd - lo) / 2; i < size; i += p) {
                    marked[i] = 1;
                }
            }

            for (const auto& [p, root] : forbidden_residues) {
                const u64 first = lo + (root + p - lo % p) % p;
                const u64 first_odd = first % 2 == 1 ? first : first + p;
                for (u64 i = (first_odd - lo) / 2; i < size; i += p) {
                    marked[i] = 1;
                }
            }

            for (std::size_t i = 0; i < size; ++i) {
                if (marked[i]) {
                    continue;
                }
                ++total.survivors;
                const u64 denominator = lo + 2 * i - 1;
                total.reciprocal_upper += (scale + denominator - 1) / denominator;
            }
        }
    };
    std::vector<std::thread> threads;
    for (u64 i = 1; i < workers; ++i) {
        threads.emplace_back(sieve_blocks, std::ref(totals[i]));
    }
    sieve_blocks(totals.front());
    for (auto& thread : threads) {
        thread.join();
    }
    Total result;
    for (const auto& total : totals) {
        result.survivors += total.survivors;
        result.reciprocal_upper += total.reciprocal_upper;
    }
    if (result.reciprocal_upper > std::numeric_limits<u64>::max()) {
        std::cerr << "Scaled reciprocal sum exceeds the output range.\n";
        return 1;
    }
    std::cout << result.survivors << '\n'
              << static_cast<u64>(result.reciprocal_upper) << '\n';
}

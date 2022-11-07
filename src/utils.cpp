#ifndef UTILS_H
#define UTILS_H
#include <fmt/color.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <gmp.h>
#include <gmpxx.h>
#include <openssl/evp.h>
#include <stdio.h>

#include <chrono>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
#define BLAKE2s256_MD_SIZE 32
// Uses OpenSSL BLAKE2s hashing.

string get_gmp_version() {
    return fmt::format("GMP: {}.{}.{}", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
}

string get_openssl_version() {
    return fmt::format("OpenSSL: {}.{}.{}", OPENSSL_VERSION_MAJOR, OPENSSL_VERSION_MINOR, OPENSSL_VERSION_PATCH);
}

void bn2gmp(mpz_class &result, BIGNUM *bn_input) {
    string hexstr = BN_bn2hex(bn_input);
    mpz_set_str(result.get_mpz_t(), hexstr.c_str(), 16);
}

void generate_RSA_modulus(mpz_class &result, mpz_class &P, mpz_class &Q, uint bits) {
    BIGNUM *bn_r = BN_new();
    BIGNUM *bn_p = BN_new();
    BIGNUM *bn_q = BN_new();
    BN_CTX *ctx = BN_CTX_new();

    do {
        BN_generate_prime_ex(bn_p, bits / 2, 1, NULL, NULL, NULL);
        BN_generate_prime_ex(bn_q, bits - (bits / 2), 1, NULL, NULL, NULL);
        BN_mul(bn_r, bn_p, bn_q, ctx);
    } while (BN_num_bits(bn_r) != bits);

    bn2gmp(result, bn_r);
    bn2gmp(P, bn_p);
    bn2gmp(Q, bn_q);
}

string DigestToHexStr(const unsigned char *digest, const unsigned int *digest_len) {
    std::stringstream ss;
    for (unsigned int i = 0; i < *digest_len; i++) {
        ss << hex << (int)digest[i];
    }
    string hexstr = ss.str();
    return hexstr;
}

bool IsPowOf2(size_t m) {
    size_t flag = m & (m - 1);
    if (m > 0 && flag == 0) {
        return true;
    } else {
        return false;
    }
}

// Compute F([X_1, X_2, X_3, X_4]) = X_1 * X_2 * X_3 * X_4
void MulMPZVec(mpz_class &prod, const vector<mpz_class> &X) {
    size_t N = X.size();
    if (N == 0 || IsPowOf2(N) == false) {
        cout << "Exception: MulMPZVec: Not a power of two. Use MulMPZVecAlt instead." << endl;
        exit(0);
    }

    if (N == 1) {
        mpz_set(prod.get_mpz_t(), X[0].get_mpz_t());
    } else {
        vector<mpz_class> I;
        mpz_class temp;

        I.reserve(N / 2);

        for (size_t i = 1; i < X.size(); i += 2) {
            mpz_mul(temp.get_mpz_t(), X[i - 1].get_mpz_t(), X[i].get_mpz_t());
            I.push_back(temp);
        }
        MulMPZVec(prod, I);
    }
}

// Compute F(a, [X_1, X_2, X_3, X_4]) = a *X_1 + a * X_2 + a * X_3 + a* X_4
void ScalarMulMPZVec(vector<mpz_class> &Y, const mpz_class a, const vector<mpz_class> &X) {
    size_t N = X.size();
    Y.resize(N);

    for (size_t i = 0; i < N; i++) {
        mpz_class temp;
        mpz_mul(temp.get_mpz_t(), X[i].get_mpz_t(), a.get_mpz_t());
        mpz_set(Y[i].get_mpz_t(), temp.get_mpz_t());
    }
}

// Given a_1, a_2, a_3, a_4. Let A = a_1 * a_2 * a_3 * a_4.
// Computes Y_1 = A / a_1, Y_2 = A / a_2, Y_3 = A / a_3, Y_4 = A / a_4
void ComputeYs(vector<mpz_class> &Y, const mpz_class &product, const vector<mpz_class> &I) {
    size_t N = I.size();
    if (N == 0 || IsPowOf2(N) == false) {
        cout << "Exception: ComputeYs: Not a power of two" << endl;
        exit(0);
    }

    Y.reserve(N);
    if (N == 1) {
        Y.push_back(product);
    } else {
        size_t mid = N / 2;
        vector<mpz_class> left;
        vector<mpz_class> right;

        vector<mpz_class> I_left(I.begin(), I.begin() + mid);
        vector<mpz_class> I_right(I.begin() + mid, I.end());

        mpz_class product_left;
        mpz_class product_right;

        MulMPZVec(product_left, I_left);
        MulMPZVec(product_right, I_right);

        mpz_mul(product_left.get_mpz_t(), product_left.get_mpz_t(), product.get_mpz_t());
        mpz_mul(product_right.get_mpz_t(), product_right.get_mpz_t(), product.get_mpz_t());

        ComputeYs(left, product_right, I_left);
        ComputeYs(right, product_left, I_right);

        Y.insert(Y.end(), left.begin(), left.end());
        Y.insert(Y.end(), right.begin(), right.end());
    }
}

// Computes the bezout coefficients of multiple elements.
void mpz_gcdext_vec(mpz_class &g, vector<mpz_class> &bezouts, const vector<mpz_class> &X) {
    size_t N = X.size();
    if (N == 0 || N == 1 || IsPowOf2(N) == false) {
        cout << "Exception: mpz_gcdext_vec: Not a power of two" << endl;
        exit(0);
    }

    bezouts.reserve(N);
    if (N == 2) {
        bezouts.resize(N);
        mpz_gcdext(g.get_mpz_t(), bezouts[0].get_mpz_t(), bezouts[1].get_mpz_t(), X[0].get_mpz_t(), X[1].get_mpz_t());
    } else {
        size_t mid = N / 2;

        mpz_class gcd_left;
        mpz_class gcd_right;

        vector<mpz_class> bezouts_left;
        vector<mpz_class> bezouts_right;

        vector<mpz_class> X_left(X.begin(), X.begin() + mid);
        vector<mpz_class> X_right(X.begin() + mid, X.end());

        mpz_gcdext_vec(gcd_left, bezouts_left, X_left);
        mpz_gcdext_vec(gcd_right, bezouts_right, X_right);

        mpz_class a, b;
        mpz_gcdext(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(), gcd_left.get_mpz_t(), gcd_right.get_mpz_t());

        ScalarMulMPZVec(bezouts_left, a, bezouts_left);
        ScalarMulMPZVec(bezouts_right, b, bezouts_right);

        bezouts.insert(bezouts.end(), bezouts_left.begin(), bezouts_left.end());
        bezouts.insert(bezouts.end(), bezouts_right.begin(), bezouts_right.end());
    }
}

void MultiExpMPZ(mpz_class &result, const vector<mpz_class> &X, const vector<mpz_class> &exponents, const mpz_class &MODULUS) {
    size_t N = X.size();
    if (N == 0 || N != exponents.size()) {
        cout << "Exception: MultiExpMPZ: Vector is sized Zero." << endl;
        exit(1);
    }

    mpz_class product, temp;
    mpz_set_ui(product.get_mpz_t(), 1);

    for (size_t i = 0; i < N; i++) {
        mpz_powm(temp.get_mpz_t(), X[i].get_mpz_t(), exponents[i].get_mpz_t(), MODULUS.get_mpz_t());
        mpz_mul(product.get_mpz_t(), product.get_mpz_t(), temp.get_mpz_t());
        mpz_mod(product.get_mpz_t(), product.get_mpz_t(), MODULUS.get_mpz_t());
    }

    mpz_set(result.get_mpz_t(), product.get_mpz_t());
}

// Uses OpenSSL BLAKE2s hashing.
string HashStrVec(const vector<string> &input) {
    unsigned char *digest = new unsigned char[BLAKE2s256_MD_SIZE];
    unsigned int *digest_len = new unsigned int;

    size_t N = input.size();
    if (N == 0) {
        cout << "Exception: Hash: Vector is sized Zero." << endl;
        exit(1);
    }

    size_t data_len = 0;
    string a = std::accumulate(input.begin(), input.end(), std::string(""));
    data_len = a.length() + 1;

    char *data_char = (char *)malloc((data_len) * sizeof(char));
    strcpy(data_char, a.c_str());

    unsigned char *data = reinterpret_cast<unsigned char *>(data_char);

    EVP_Digest(data, data_len, digest, digest_len, EVP_blake2s256(), NULL);
    string hexstr = DigestToHexStr(digest, digest_len);
    return hexstr;
}

string HashMPZVec(const vector<mpz_class> &input) {
    vector<string> input_str(input.size());
    for (int i = 0; i < input.size(); i++) {
        input_str[i] = mpz_get_str(NULL, 10, input[i].get_mpz_t());
    }
    return HashStrVec(input_str);
}

// Performs 15 iterations of Miller-Rabin's test. GMP suggests that reps should be between 15-50.
void Hash_to_Prime(mpz_class &result, const mpz_class &input, int reps = 15) {
    unsigned char *digest = new unsigned char[BLAKE2s256_MD_SIZE];
    unsigned int *digest_len = new unsigned int;

    size_t data_len = mpz_sizeinbase(input.get_mpz_t(), 10);
    char *data_char = (char *)malloc((data_len) * sizeof(char));
    mpz_get_str(data_char, 10, input.get_mpz_t());
    unsigned char *data_t = reinterpret_cast<unsigned char *>(data_char);

    EVP_Digest(data_t, data_len, digest, digest_len, EVP_blake2s256(), NULL);
    string hexstr = DigestToHexStr(digest, digest_len);

    mpz_class p;
    mpz_set_str(p.get_mpz_t(), hexstr.c_str(), 16);

    int status = mpz_probab_prime_p(p.get_mpz_t(), reps);

    unsigned char *data;
    data = (unsigned char *)realloc(data, (*digest_len) * sizeof(unsigned char));
    if (data == NULL) {
        cout << "HASH to primes: NOT ENOUGH MEMORY" << endl;
        exit(1);
    }

    while (status == 0) {
        memcpy(data, digest, *digest_len);
        EVP_Digest(data, *digest_len, digest, digest_len, EVP_blake2s256(), NULL);
        hexstr = DigestToHexStr(digest, digest_len);     // Convert to Hex
        mpz_set_str(p.get_mpz_t(), hexstr.c_str(), 16);  // From Hex
        status = mpz_probab_prime_p(p.get_mpz_t(), reps);
    }
    mpz_set(result.get_mpz_t(), p.get_mpz_t());
    free(digest);
    free(data);
}

// Performs 15 iterations of Miller-Rabin's test. GMP suggests that reps should be between 15-50.
void Hash_to_Prime(mpz_class &result, const vector<mpz_class> &input, int reps = 15) {
    string hexstr = HashMPZVec(input);

    mpz_class p;
    mpz_set_str(p.get_mpz_t(), hexstr.c_str(), 16);

    int status = mpz_probab_prime_p(p.get_mpz_t(), reps);

    while (status == 0) {
        hexstr = HashStrVec(vector<string>{hexstr});     // Convert to Hex
        mpz_set_str(p.get_mpz_t(), hexstr.c_str(), 16);  // From Hex
        status = mpz_probab_prime_p(p.get_mpz_t(), reps);
    }
    mpz_set(result.get_mpz_t(), p.get_mpz_t());
}

// Will hash the inputs using the Blake256 hashing. Use the output as seed for PRNG.
// Will generate a pseudorandom number between 0 to 2^bits - 1
void Hash_to_Group(mpz_class &result, const vector<mpz_class> &input, const mpz_class &MODULUS) {
    mpz_class seed;
    string hexstr = "1" + HashMPZVec(input);
    mpz_set_str(seed.get_mpz_t(), hexstr.c_str(), 16);
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed(rstate, seed.get_mpz_t());
    mpz_urandomm(result.get_mpz_t(), rstate, MODULUS.get_mpz_t());
}

// Compute F([X_1, X_2, X_3, X_4]) MOD MODULUS = (X_1 * X_2 * X_3 * X_4) % MODULUS
// To do make it tree based approach so that A, B in A * B is of same size.
void MulMPZModVec(mpz_class &prod, const vector<mpz_class> &X, const mpz_class &MODULUS) {
    size_t N = X.size();
    if (N == 0) {
        cout << "Exception: MulMPZVec: Vector is sized Zero." << endl;
        exit(0);
    }

    mpz_set(prod.get_mpz_t(), X[0].get_mpz_t());
    mpz_mod(prod.get_mpz_t(), prod.get_mpz_t(), MODULUS.get_mpz_t());
    for (size_t i = 1; i < N; i++) {
        mpz_mul(prod.get_mpz_t(), prod.get_mpz_t(), X[i].get_mpz_t());
        mpz_mod(prod.get_mpz_t(), prod.get_mpz_t(), MODULUS.get_mpz_t());
    }
}

void FmtString(int index, string operation, std::chrono::nanoseconds t, int runs) {
    fmt::print("Benchmark\t{:>3}\t{:<20}\t{:>50.3f} us/op\t{:>3}\n", index, operation, chrono::duration<double, micro>(t).count(), runs);
}

void check_status(int index, string operation, bool status) {
    if (status == false) {
        auto str = fmt::format("Err:\t{:>2}\t{:<20}\t{:=^20}\n", index, operation, " FAIL ");
        fmt::print(fg(fmt::color::red), str);
    }
}

void PrintSize(int index, string variable, mpz_class t) {
    fmt::print("Size\t{:>3}\t{:<20}\t{:>20} bits\n", index, variable, mpz_sizeinbase(t.get_mpz_t(), 2));
}

void hash_to_prime(vector<mpz_class> &exponents, const uint &count, int reps = 15) {
    exponents.resize(count);

    for (size_t i = 0; i < exponents.size(); i++) {
        mpz_set_ui(exponents[i].get_mpz_t(), i);
    }

    auto t2 = chrono::high_resolution_clock::now();
    auto t3 = chrono::high_resolution_clock::now();
    auto t4 = t2 - t2;

    t2 = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < exponents.size(); i++) {
        Hash_to_Prime(exponents[i], exponents[i], reps);
    }
    t3 = chrono::high_resolution_clock::now();
    t4 = t3 - t2;
    t4 = t4 / (exponents.size());
    FmtString(100, fmt::format("Hash2Prime", reps), t4, 100);
}

void exponentiations(const mpz_class &base, const mpz_class &MODULUS, const vector<mpz_class> &exponents) {
    int N = exponents.size();
    N = 100 < N ? 100 : N;
    auto t2 = chrono::high_resolution_clock::now();
    auto t3 = chrono::high_resolution_clock::now();
    auto t4 = t2 - t2;

    mpz_class result;
    t2 = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < N; i++) {
        mpz_powm(result.get_mpz_t(), base.get_mpz_t(), exponents[i].get_mpz_t(), MODULUS.get_mpz_t());
    }
    t3 = chrono::high_resolution_clock::now();
    t4 = t3 - t2;
    t4 = t4 / (N);
    FmtString(N, "Exponentiations", t4, N);
}
#endif

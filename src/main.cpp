#include <gmp.h>
#include <gmpxx.h>
#include <openssl/evp.h>
#include <stdio.h>

#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#include "rsa.hpp"
using namespace std;

void RSA_test() {
    cout << get_gmp_version() << endl;
    cout << get_openssl_version() << endl;

    auto t2 = chrono::high_resolution_clock::now();
    auto t3 = chrono::high_resolution_clock::now();
    auto t4 = t2 - t2;

    uint security = 2048;
    uint runs = 1;

    // vector<uint> ell = {10, 12, 14, 16, 18};
    vector<uint> ell = {8};

    vector<mpz_class> exponents;
    // Generate all primes
    hash_to_prime(exponents, (1 << (ell[ell.size() - 1])), 15);


    bool status;
    Acc acc(security);
    acc.Setup();

    // Benchmark exponentiations
    exponentiations(acc.BASE, acc.MODULUS, exponents);

    cout << std::fixed << std::boolalpha;
    acc.Print();

    for (const auto &l : ell) {
        cout << " l = " << l << endl;
        mpz_class digest_X;
        mpz_class digest_XI;
        vector<mpz_class> memProofs_I;
        vector<pair<mpz_class, mpz_class> > nonmemProofs_I;

        mpz_class I_prod;
        mpz_class aggMemProof;
        pair<mpz_class, mpz_class> aggNonMemProof;
        mpz_class aggMemProofTemp, aggMemProofPoE;
        NiNonMemPi aggNonMemProofPoE;

        t2 = chrono::high_resolution_clock::now();
        size_t N = 1 << l;
        vector<mpz_class> data(exponents.begin(), exponents.begin() + N);
        cout << "Done copying X union I set." << endl;
        vector<mpz_class> X(data.begin(), data.begin() + (N / 2));
        cout << "Done copying X set." << endl;
        vector<mpz_class> I(data.begin() + (N / 2), data.begin() + N);
        cout << "Done copying I set." << endl;

        // =================================================================
        mpz_class product;
        t2 = chrono::high_resolution_clock::now();
        for (uint r = 0; r < (runs + 1); r++) {
            MulMPZVec(product, I);
        }
        t3 = chrono::high_resolution_clock::now();
        t4 = t3 - t2;
        t4 = t4 / (runs + 1);
        FmtString(l - 1, "MulMPZVec", t4, runs + 1);
        // =================================================================
        mpz_class result;
        t2 = chrono::high_resolution_clock::now();
        for (uint r = 0; r < runs; r++) {
            mpz_powm(result.get_mpz_t(), acc.BASE.get_mpz_t(), product.get_mpz_t(), acc.MODULUS.get_mpz_t());
        }
        t3 = chrono::high_resolution_clock::now();
        t4 = t3 - t2;
        t4 = t4 / runs;
        FmtString(l - 1, "LargeModExp", t4, runs);
        // =================================================================

        acc.CommitFake(digest_X, acc.BASE, X);
        cout << "Done commit." << endl;

        acc.MemProveFake(memProofs_I, digest_X, I);
        cout << "Done MemProveFake." << endl;
        mpz_powm(digest_XI.get_mpz_t(), memProofs_I[0].get_mpz_t(), I[0].get_mpz_t(), acc.MODULUS.get_mpz_t());

        acc.NonMemProveFake(nonmemProofs_I, X, I);
        cout << "Done NonMemProveFake." << endl;
        t3 = chrono::high_resolution_clock::now();
        t4 = t3 - t2;
        FmtString(l - 1, "DataGen", t4, 1);

        // =================================================================
        t2 = chrono::high_resolution_clock::now();
        for (uint r = 0; r < runs; r++) {
            acc.AggMemProve(aggMemProof, I_prod, I, memProofs_I);
        }
        t3 = chrono::high_resolution_clock::now();
        t4 = t3 - t2;
        t4 = t4 / runs;
        FmtString(l - 1, "AggMemProve", t4, runs);
        // =================================================================

        t2 = chrono::high_resolution_clock::now();
        for (uint r = 0; r < runs; r++) {
            status = acc.AggMemVerify(digest_XI, I, aggMemProof);
        }
        t3 = chrono::high_resolution_clock::now();
        t4 = t3 - t2;
        t4 = t4 / runs;
        FmtString(l - 1, "AggMemVerify", t4, runs);
        check_status(l - 1, "AggMemVerify", status);

        // =================================================================
        // mpz_set_si(I[I.size() - 1].get_mpz_t(), 20);
        t2 = chrono::high_resolution_clock::now();
        for (uint r = 0; r < runs; r++) {
            acc.AggMemProvePoE(aggMemProofTemp, aggMemProofPoE, I_prod, digest_XI, I, memProofs_I);
        }
        t3 = chrono::high_resolution_clock::now();
        t4 = t3 - t2;
        t4 = t4 / runs;
        FmtString(l - 1, "AggMemProvePoE", t4, runs);

        // =================================================================
        status = acc.AggMemVerify(digest_XI, I, aggMemProofTemp);
        check_status(l - 1, "AggMemVerifyPoENaive", status);
        t2 = chrono::high_resolution_clock::now();
        for (uint r = 0; r < runs; r++) {
            status = acc.AggMemVerifyPoE(digest_XI, I, aggMemProofTemp, aggMemProofPoE);
        }
        t3 = chrono::high_resolution_clock::now();
        t4 = t3 - t2;
        t4 = t4 / runs;
        FmtString(l - 1, "AggMemVerifyPoE", t4, runs);
        check_status(l - 1, "AggMemVerifyPoE", status);

        // =================================================================
        t2 = chrono::high_resolution_clock::now();
        for (uint r = 0; r < runs; r++) {
            acc.AggNonMemProve(aggNonMemProof, I_prod, digest_X, I, nonmemProofs_I);
        }
        t3 = chrono::high_resolution_clock::now();
        t4 = t3 - t2;
        t4 = t4 / runs;
        FmtString(l - 1, "AggNonMemProve", t4, runs);

        // =================================================================
        t2 = chrono::high_resolution_clock::now();
        for (uint r = 0; r < runs; r++) {
            status = acc.AggNonMemVerify(digest_X, I, aggNonMemProof);
        }
        t3 = chrono::high_resolution_clock::now();
        t4 = t3 - t2;
        t4 = t4 / runs;
        FmtString(l - 1, "AggNonMemVerify", t4, runs);
        check_status(l - 1, "AggNonMemVerify", status);

        // =================================================================
        t2 = chrono::high_resolution_clock::now();
        for (uint r = 0; r < runs; r++) {
            acc.AggNonMemProvePoE(aggNonMemProofPoE, I_prod, digest_X, I, nonmemProofs_I);
        }
        t3 = chrono::high_resolution_clock::now();
        t4 = t3 - t2;
        t4 = t4 / runs;
        FmtString(l - 1, "AggNonMemProvePoE", t4, runs);

        // =================================================================
        status = acc.AggNonMemVerify(digest_X, I, aggNonMemProofPoE.nonmem_pi);
        check_status(l - 1, "AggNonMemVerifyPoENaive", status);
        t2 = chrono::high_resolution_clock::now();
        for (uint r = 0; r < runs; r++) {
            status = acc.AggNonMemVerifyPoE(aggNonMemProofPoE, digest_X, I);  // {I.begin(), I.end()}
        }
        t3 = chrono::high_resolution_clock::now();
        t4 = t3 - t2;
        t4 = t4 / runs;
        FmtString(l - 1, "AggNonMemVerifyPoE", t4, runs);
        check_status(l - 1, "AggNonMemVerifyPoE", status);

        // =================================================================
        t2 = chrono::high_resolution_clock::now();
        for (uint r = 0; r < runs; r++) {
            acc.Commit(digest_X, acc.BASE, X);
        }
        t3 = chrono::high_resolution_clock::now();
        t4 = t3 - t2;
        t4 = t4 / runs;
        FmtString(l - 1, "Commit", t4, runs);

        // =================================================================
        status = true;
        t2 = chrono::high_resolution_clock::now();
        for (uint r = 0; r < runs; r++) {
            for (int k = 0; k < memProofs_I.size(); k++) {
                status = status & acc.MemVerifySingle(digest_XI, I[k], memProofs_I[k]);
            }
        }
        t3 = chrono::high_resolution_clock::now();
        t4 = t3 - t2;
        t4 = t4 / runs;
        FmtString(l - 1, "MemVerify", t4, runs);
        check_status(l - 1, "MemVerify", status);

        // =================================================================
        status = true;
        t2 = chrono::high_resolution_clock::now();
        for (uint r = 0; r < runs; r++) {
            for (int k = 0; k < nonmemProofs_I.size(); k++) {
                status = status & acc.NonMemVerifySingle(digest_X, I[k], nonmemProofs_I[k]);
            }
        }
        t3 = chrono::high_resolution_clock::now();
        t4 = t3 - t2;
        t4 = t4 / runs;
        FmtString(l - 1, "NonMemVerify", t4, runs);
        check_status(l - 1, "NonMemVerify", status);

        // =================================================================
        PrintSize(l - 1, "digest_XI", digest_XI);
        PrintSize(l - 1, "aggMemProof", aggMemProof);
        PrintSize(l - 1, "aggMemProofPoE", aggMemProofPoE);
        PrintSize(l - 1, "digest_X", digest_X);
        PrintSize(l - 1, "aggNonMemProof.first", aggNonMemProof.first);
        PrintSize(l - 1, "aggNonMemProof.second", aggNonMemProof.second);
        PrintSize(l - 1, "aggNonMemProofPoE.V", aggNonMemProofPoE.V);
        PrintSize(l - 1, "aggNonMemProofPoE.nonmem_pi.first", aggNonMemProofPoE.nonmem_pi.first);
        PrintSize(l - 1, "aggNonMemProofPoE.nonmem_pi.second", aggNonMemProofPoE.nonmem_pi.second);
        PrintSize(l - 1, "aggNonMemProofPoE.poke2_pi.z", aggNonMemProofPoE.poke2_pi.z);
        PrintSize(l - 1, "aggNonMemProofPoE.poke2_pi.Q", aggNonMemProofPoE.poke2_pi.Q);
        PrintSize(l - 1, "aggNonMemProofPoE.poke2_pi.r", aggNonMemProofPoE.poke2_pi.r);
        PrintSize(l - 1, "aggNonMemProofPoE.poe_pi", aggNonMemProofPoE.poe_pi);
        PrintSize(l - 1, "memProofs_I[0]", memProofs_I[0]);
        PrintSize(l - 1, "nonmemProofs_I[0].first", nonmemProofs_I[0].first);
        PrintSize(l - 1, "nonmemProofs_I[0].second", nonmemProofs_I[0].second);
        PrintSize(l - 1, "exponents[0]", exponents[0]);
    }
}

int main() {
    RSA_test();
    return 0;
}

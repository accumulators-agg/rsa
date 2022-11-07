#ifndef RSA_HPP
#define RSA_HPP

#include <gmp.h>
#include <gmpxx.h>
#include <openssl/evp.h>
#include <stdio.h>

#include <iostream>
#include <string>
#include <vector>

#include "utils.cpp"

using namespace std;

class NiPoKE2Pi {
   public:
    mpz_class z;
    mpz_class Q;
    mpz_class r;
};

class NiNonMemPi {
   public:
    mpz_class V;
    pair<mpz_class, mpz_class> nonmem_pi;
    NiPoKE2Pi poke2_pi;
    mpz_class poe_pi;
};

class Acc {
   private:
    mpz_class P;
    mpz_class Q;
    mpz_class PHI;

   public:
    uint bits;
    mpz_class MODULUS;
    mpz_class BASE;

    unsigned char digest[BLAKE2s256_MD_SIZE];
    unsigned int *digest_len;
    Acc(uint bits = 2048);
    void Setup();
    void Commit(mpz_class &digest, const mpz_class &base, const vector<mpz_class> &element);
    void MemProve(vector<mpz_class> &proofVec, const mpz_class &digest, const vector<mpz_class> &I);
    bool MemVerifySingle(const mpz_class &digest, const mpz_class &I, const mpz_class &proof);
    void NonMemProve(vector<pair<mpz_class, mpz_class> > &proofVec, const vector<mpz_class> &X, const vector<mpz_class> &I);
    bool NonMemVerifySingle(const mpz_class &digest, const mpz_class &y, const pair<mpz_class, mpz_class> &pi);
    void AggMemProve(mpz_class &proof, mpz_class &I_prod, const vector<mpz_class> &I, const vector<mpz_class> &proofs);
    void ShamirTrick(mpz_class &result, const mpz_class &wX, const mpz_class &wY, const mpz_class &X, const mpz_class &Y);
    void AggMemWit(mpz_class &result, mpz_class &xy, const mpz_class &wX, const mpz_class &wY, const mpz_class &X, const mpz_class &Y);
    bool AggMemVerify(const mpz_class &digest, const vector<mpz_class> &I, const mpz_class &proof);
    void AggNonMemProve(pair<mpz_class, mpz_class> &result, mpz_class &I_prod, const mpz_class &digest, const vector<mpz_class> &I, const vector<pair<mpz_class, mpz_class> > &proofs);
    bool AggNonMemVerify(const mpz_class &digest, const vector<mpz_class> &I, const pair<mpz_class, mpz_class> &proof);
    void BreakNonMemWit(vector<pair<mpz_class, mpz_class> > &proofVec, const mpz_class &digest, const pair<mpz_class, mpz_class> &proof, const vector<mpz_class> &I);
    void AggNonMemWit(pair<mpz_class, mpz_class> &result, const mpz_class &digest, const mpz_class &xy, const mpz_class &I, const pair<mpz_class, mpz_class> &pi_I, const mpz_class &J, const pair<mpz_class, mpz_class> &pi_J);
    void CommitFake(mpz_class &digest, const mpz_class &base, const vector<mpz_class> &elements);
    void MemProveFake(vector<mpz_class> &proofVec, const mpz_class &digest, const vector<mpz_class> &I);
    void NonMemProveFake(vector<pair<mpz_class, mpz_class> > &proofVec, const vector<mpz_class> &X, const vector<mpz_class> &I);
    void NiPoEProve(mpz_class &Q, const mpz_class &x, const mpz_class &u, const mpz_class &w);
    bool NiPoEVerify(const mpz_class &Q, const mpz_class &x, const mpz_class &u, const mpz_class &w);
    void NiPoKE2Prove(NiPoKE2Pi &pi, const mpz_class &x, const mpz_class &u, const mpz_class &w);
    bool NiPoKE2Verify(const NiPoKE2Pi &pi, const mpz_class &u, const mpz_class &w);
    void AggMemProvePoE(mpz_class &proof, mpz_class &poe_pi, mpz_class &I_prod, const mpz_class &digest, const vector<mpz_class> &I, const vector<mpz_class> &proofs);
    bool AggMemVerifyPoE(const mpz_class &digest, const vector<mpz_class> &I, const mpz_class &proof, const mpz_class &poe_pi);
    void AggNonMemProvePoE(NiNonMemPi &proof, mpz_class &I_prod, const mpz_class &digest, const vector<mpz_class> &I, const vector<pair<mpz_class, mpz_class> > &proofs);
    bool AggNonMemVerifyPoE(const NiNonMemPi pi, const mpz_class &digest, const vector<mpz_class> &I);
    void Print();
};

#include "rsa-agg.cpp"
#include "rsa-basic.cpp"
#include "rsa-fake.cpp"
#include "rsa-poe.cpp"
#endif

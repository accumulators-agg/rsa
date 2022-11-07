#ifndef RSA_AGG_CPP
#define RSA_AGG_CPP
void Acc::AggMemProve(mpz_class &proof, mpz_class &I_prod, const vector<mpz_class> &I, const vector<mpz_class> &proofs) {
    size_t N = I.size();
    if (N == 0 || IsPowOf2(N) == false) {
        cout << "Exception: ProveMem: Not a power of two" << endl;
        exit(0);
    }

    if (N != proofs.size()) {
        cout << "Exception: ProveMem: |I| != |Proofs|" << endl;
        exit(0);
    }

    if (N == 1) {
        mpz_set(proof.get_mpz_t(), proofs[0].get_mpz_t());
        mpz_set(I_prod.get_mpz_t(), I[0].get_mpz_t());
    } else {
        vector<mpz_class> agg_I;
        vector<mpz_class> agg_proofs;
        agg_I.reserve(N / 2);
        agg_proofs.reserve(N / 2);

        for (size_t i = 1; i < I.size(); i += 2) {
            mpz_class small_proof;
            mpz_class xy;
            AggMemWit(small_proof, xy, proofs[i - 1], proofs[i], I[i - 1], I[i]);
            agg_proofs.push_back(small_proof);
            agg_I.push_back(xy);
        }

        AggMemProve(proof, I_prod, agg_I, agg_proofs);
    }
}

void Acc::AggMemWit(mpz_class &result, mpz_class &xy, const mpz_class &wX, const mpz_class &wY, const mpz_class &X, const mpz_class &Y) {
    mpz_mul(xy.get_mpz_t(), X.get_mpz_t(), Y.get_mpz_t());  // Since X and Y are assumed to primes, it is fine.
    ShamirTrick(result, wX, wY, X, Y);
}

void Acc::ShamirTrick(mpz_class &result, const mpz_class &wX, const mpz_class &wY, const mpz_class &X, const mpz_class &Y) {
    mpz_class g, a, b;
    mpz_gcdext(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(), X.get_mpz_t(), Y.get_mpz_t());

    mpz_class r, p, q;
    mpz_powm(p.get_mpz_t(), wX.get_mpz_t(), b.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_powm(q.get_mpz_t(), wY.get_mpz_t(), a.get_mpz_t(), MODULUS.get_mpz_t());

    mpz_mul(r.get_mpz_t(), p.get_mpz_t(), q.get_mpz_t());
    mpz_mod(r.get_mpz_t(), r.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_set(result.get_mpz_t(), r.get_mpz_t());
}

bool Acc::AggMemVerify(const mpz_class &digest, const vector<mpz_class> &I, const mpz_class &proof) {
    size_t N = I.size();
    if (N == 0 || IsPowOf2(N) == false) {
        cout << "Exception: VerifyMemAgg: Not a power of two" << endl;
        exit(0);
    }

    mpz_class prod, result;
    MulMPZVec(prod, I);
    mpz_powm(result.get_mpz_t(), proof.get_mpz_t(), prod.get_mpz_t(), MODULUS.get_mpz_t());
    int status = mpz_cmp(result.get_mpz_t(), digest.get_mpz_t());
    return status == 0;
}

// From BBF19
void Acc::AggNonMemProve(pair<mpz_class, mpz_class> &result, mpz_class &I_prod, const mpz_class &digest,
                         const vector<mpz_class> &I, const vector<pair<mpz_class, mpz_class> > &proofs) {
    size_t N = I.size();
    if (N == 0 || IsPowOf2(N) == false) {
        cout << "Exception: ProveNonMemAgg: Not a power of two" << endl;
        exit(0);
    }

    if (proofs.size() != N) {
        cout << "Exception: ProveNonMemAgg: |I| != |Proofs|" << endl;
        exit(0);
    }

    if (N == 1) {
        mpz_set(result.first.get_mpz_t(), proofs[0].first.get_mpz_t());
        mpz_set(result.second.get_mpz_t(), proofs[0].second.get_mpz_t());
        mpz_set(I_prod.get_mpz_t(), I[0].get_mpz_t());
    } else {
        vector<mpz_class> agg_I;
        vector<pair<mpz_class, mpz_class> > agg_proofs;
        agg_I.reserve(N / 2);
        agg_proofs.reserve(N / 2);

        for (size_t i = 1; i < I.size(); i += 2) {
            pair<mpz_class, mpz_class> small_proof;
            mpz_class xy;
            mpz_mul(xy.get_mpz_t(), I[i - 1].get_mpz_t(), I[i].get_mpz_t());
            AggNonMemWit(small_proof, digest, xy, I[i - 1], proofs[i - 1], I[i], proofs[i]);

            agg_proofs.push_back(small_proof);
            agg_I.push_back(xy);
        }
        AggNonMemProve(result, I_prod, digest, agg_I, agg_proofs);
    }
}

bool Acc::AggNonMemVerify(const mpz_class &digest, const vector<mpz_class> &I, const pair<mpz_class, mpz_class> &proof) {
    mpz_class product, t1, t2, t;

    MulMPZVec(product, I);
    mpz_powm(t1.get_mpz_t(), digest.get_mpz_t(), proof.first.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_powm(t2.get_mpz_t(), proof.second.get_mpz_t(), product.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_mul(t.get_mpz_t(), t1.get_mpz_t(), t2.get_mpz_t());
    mpz_mod(t.get_mpz_t(), t.get_mpz_t(), MODULUS.get_mpz_t());

    return mpz_cmp(BASE.get_mpz_t(), t.get_mpz_t()) == 0;
}

void Acc::AggNonMemWit(pair<mpz_class, mpz_class> &result,
                       const mpz_class &digest, const mpz_class &xy,
                       const mpz_class &I, const pair<mpz_class, mpz_class> &pi_I,
                       const mpz_class &J, const pair<mpz_class, mpz_class> &pi_J) {
    mpz_class g, alpha, beta;
    mpz_gcdext(g.get_mpz_t(), alpha.get_mpz_t(), beta.get_mpz_t(), I.get_mpz_t(), J.get_mpz_t());
    int gcd = mpz_cmp_ui(g.get_mpz_t(), 1);
    if (gcd != 0) {
        cout << "Exception: AggNonMemWit: GCD was not 1!" << endl;
        exit(0);
    }
    mpz_class B, temp1, temp2;
    mpz_powm(temp1.get_mpz_t(), pi_I.second.get_mpz_t(), beta.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_powm(temp2.get_mpz_t(), pi_J.second.get_mpz_t(), alpha.get_mpz_t(), MODULUS.get_mpz_t());

    mpz_mul(B.get_mpz_t(), temp1.get_mpz_t(), temp2.get_mpz_t());
    mpz_mod(B.get_mpz_t(), B.get_mpz_t(), MODULUS.get_mpz_t());

    mpz_mul(temp1.get_mpz_t(), pi_I.first.get_mpz_t(), J.get_mpz_t());
    mpz_mul(temp1.get_mpz_t(), temp1.get_mpz_t(), beta.get_mpz_t());

    mpz_mul(temp2.get_mpz_t(), pi_J.first.get_mpz_t(), I.get_mpz_t());
    mpz_mul(temp2.get_mpz_t(), temp2.get_mpz_t(), alpha.get_mpz_t());

    mpz_class a, a_prime;
    mpz_add(a_prime.get_mpz_t(), temp1.get_mpz_t(), temp2.get_mpz_t());

    mpz_class temp;
    mpz_cdiv_qr(a_prime.get_mpz_t(), a.get_mpz_t(), a_prime.get_mpz_t(), xy.get_mpz_t());

    mpz_powm(temp.get_mpz_t(), digest.get_mpz_t(), a_prime.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_mul(temp.get_mpz_t(), temp.get_mpz_t(), B.get_mpz_t());
    mpz_mod(temp.get_mpz_t(), temp.get_mpz_t(), MODULUS.get_mpz_t());
    result = pair<mpz_class, mpz_class>{a, temp};
}

#endif

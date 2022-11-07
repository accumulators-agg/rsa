#ifndef RSA_FAKE_CPP
#define RSA_FAKE_CPP

void Acc::CommitFake(mpz_class &digest, const mpz_class &base, const vector<mpz_class> &elements) {
    if (elements.size() == 0) {
        mpz_set_ui(digest.get_mpz_t(), 1);
        return;
    }

    mpz_class exponent, result;
    MulMPZModVec(exponent, elements, PHI);

    mpz_powm(result.get_mpz_t(), base.get_mpz_t(), exponent.get_mpz_t(), MODULUS.get_mpz_t());

    mpz_set(digest.get_mpz_t(), result.get_mpz_t());
}

void Acc::MemProveFake(vector<mpz_class> &proofVec, const mpz_class &digest, const vector<mpz_class> &I) {
    size_t N = I.size();
    if (N == 0 || IsPowOf2(N) == false) {
        cout << "Exception: ProveMem: Not a power of two" << endl;
        exit(0);
    }

    proofVec.reserve(N);
    if (N == 1) {
        proofVec.push_back(digest);
    } else {
        size_t mid = N / 2;
        vector<mpz_class> left_proofs;
        vector<mpz_class> right_proofs;

        vector<mpz_class> I_left(I.begin(), I.begin() + mid);
        vector<mpz_class> I_right(I.begin() + mid, I.end());

        mpz_class digest_left;
        mpz_class digest_right;

        CommitFake(digest_left, digest, I_left);
        CommitFake(digest_right, digest, I_right);

        MemProveFake(left_proofs, digest_right, I_left);   // Note that LEFT side of the proofs require digest of the RIGHT.
        MemProveFake(right_proofs, digest_left, I_right);  // Note that RIGHT side of the proofs require digest of the LEFT.

        proofVec.insert(proofVec.end(), left_proofs.begin(), left_proofs.end());
        proofVec.insert(proofVec.end(), right_proofs.begin(), right_proofs.end());
    }
}

void Acc::NonMemProveFake(vector<pair<mpz_class, mpz_class> > &proofVec, const vector<mpz_class> &X, const vector<mpz_class> &I) {
    size_t N = I.size();
    if (N == 0 || IsPowOf2(N) == false) {
        cout << "Exception: NonProveMem: Not a power of two" << endl;
        exit(0);
    }
    proofVec.reserve(N);

    mpz_class S;
    MulMPZModVec(S, X, PHI);

    mpz_class a, b, g, B;
    for (size_t i = 0; i < N; i++) {
        mpz_gcdext(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(), S.get_mpz_t(), I[i].get_mpz_t());
        int gcd = mpz_cmp_ui(g.get_mpz_t(), 1);
        if (gcd != 0) {
            cout << "Exception: NonProveMem: GCD was not 1!" << endl;
            exit(0);
        }
        mpz_mod(b.get_mpz_t(), b.get_mpz_t(), PHI.get_mpz_t());
        mpz_powm(B.get_mpz_t(), BASE.get_mpz_t(), b.get_mpz_t(), MODULUS.get_mpz_t());
        pair<mpz_class, mpz_class> non_mem_proof{a, B};
        proofVec.push_back(non_mem_proof);
    }
}

#endif

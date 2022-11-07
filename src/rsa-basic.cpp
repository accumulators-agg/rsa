#ifndef RSA_BASIC_CPP
#define RSA_BASIC_CPP

// https://crypto.stackexchange.com/a/8692
// 2048 **approx** gives 128 bits of security
Acc::Acc(uint bits) {
    this->bits = bits;
    digest_len = new unsigned int;
}

void Acc::Setup() {
    generate_RSA_modulus(MODULUS, P, Q, bits);
    mpz_nextprime(BASE.get_mpz_t(), Q.get_mpz_t());
    mpz_mod(BASE.get_mpz_t(), BASE.get_mpz_t(), MODULUS.get_mpz_t());

    mpz_class P_SUB_ONE, Q_SUB_ONE;
    mpz_sub_ui(P_SUB_ONE.get_mpz_t(), P.get_mpz_t(), 1);
    mpz_sub_ui(Q_SUB_ONE.get_mpz_t(), Q.get_mpz_t(), 1);
    mpz_lcm(PHI.get_mpz_t(), P_SUB_ONE.get_mpz_t(), Q_SUB_ONE.get_mpz_t());
}


// Assumes that values are primes
void Acc::Commit(mpz_class &digest, const mpz_class &base, const vector<mpz_class> &elements) {
    if (elements.size() == 0) {
        mpz_set_ui(digest.get_mpz_t(), 1);
        return;
    }

    mpz_class prod, result;
    MulMPZVec(prod, elements);
    mpz_powm(result.get_mpz_t(), base.get_mpz_t(), prod.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_set(digest.get_mpz_t(), result.get_mpz_t());
}

void Acc::MemProve(vector<mpz_class> &proofVec, const mpz_class &digest, const vector<mpz_class> &I) {
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

        Commit(digest_left, digest, I_left);
        Commit(digest_right, digest, I_right);

        MemProve(left_proofs, digest_right, I_left);   // Note that LEFT side of the proofs require digest of the RIGHT.
        MemProve(right_proofs, digest_left, I_right);  // Note that RIGHT side of the proofs require digest of the LEFT.

        proofVec.insert(proofVec.end(), left_proofs.begin(), left_proofs.end());
        proofVec.insert(proofVec.end(), right_proofs.begin(), right_proofs.end());
    }
}

bool Acc::MemVerifySingle(const mpz_class &digest, const mpz_class &I, const mpz_class &proof) {
    mpz_class result;
    mpz_powm(result.get_mpz_t(), proof.get_mpz_t(), I.get_mpz_t(), MODULUS.get_mpz_t());

    int status = mpz_cmp(result.get_mpz_t(), digest.get_mpz_t());
    return status == 0;
}

void Acc::NonMemProve(vector<pair<mpz_class, mpz_class> > &proofVec, const vector<mpz_class> &X, const vector<mpz_class> &I) {
    size_t N = I.size();
    if (N == 0 || IsPowOf2(N) == false) {
        cout << "Exception: NonProveMem: Not a power of two" << endl;
        exit(0);
    }
    proofVec.reserve(N);

    mpz_class S;
    MulMPZVec(S, X);

    mpz_class a, b, g, B;
    for (size_t i = 0; i < N; i++) {
        mpz_gcdext(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(), S.get_mpz_t(), I[i].get_mpz_t());
        int gcd = mpz_cmp_ui(g.get_mpz_t(), 1);
        if (gcd != 0) {
            cout << "Exception: NonProveMem: GCD was not 1!" << endl;
            exit(0);
        }
        mpz_powm(B.get_mpz_t(), BASE.get_mpz_t(), b.get_mpz_t(), MODULUS.get_mpz_t());
        pair<mpz_class, mpz_class> non_mem_proof{a, B};
        proofVec.push_back(non_mem_proof);
    }
}

void Acc::BreakNonMemWit(vector<pair<mpz_class, mpz_class> > &proofVec, const mpz_class &digest,
                         const pair<mpz_class, mpz_class> &proof, const vector<mpz_class> &I) {
    size_t N = I.size();
    if (N == 0 || IsPowOf2(N) == false) {
        cout << "Exception: BreakNonMemWit: Not a power of two" << endl;
        exit(0);
    }
    proofVec.reserve(N);
    if (N == 1) {
        pair<mpz_class, mpz_class> temp;
        mpz_set(temp.first.get_mpz_t(), proof.first.get_mpz_t());
        mpz_set(temp.second.get_mpz_t(), proof.second.get_mpz_t());
        proofVec.push_back(temp);
    } else {
        size_t mid = N / 2;
        vector<mpz_class> I_left(I.begin(), I.begin() + mid);
        vector<mpz_class> I_right(I.begin() + mid, I.end());

        mpz_class product_left;
        mpz_class product_right;
        MulMPZVec(product_left, I_left);
        MulMPZVec(product_right, I_right);

        mpz_class a_left;
        mpz_class a_right;
        mpz_class a, B;
        mpz_set(a.get_mpz_t(), proof.first.get_mpz_t());
        mpz_set(B.get_mpz_t(), proof.second.get_mpz_t());

        mpz_mod(a_left.get_mpz_t(), a.get_mpz_t(), product_left.get_mpz_t());
        mpz_mod(a_right.get_mpz_t(), a.get_mpz_t(), product_right.get_mpz_t());

        mpz_class B_left;
        mpz_class B_right;

        mpz_class e_left, e_right;
        mpz_tdiv_q(e_left.get_mpz_t(), a.get_mpz_t(), product_left.get_mpz_t());
        mpz_tdiv_q(e_right.get_mpz_t(), a.get_mpz_t(), product_right.get_mpz_t());

        vector<mpz_class> base = {B, digest};
        vector<mpz_class> exponents1 = {product_right, e_left};
        vector<mpz_class> exponents2 = {product_left, e_right};
        MultiExpMPZ(B_left, base, exponents1, MODULUS);
        MultiExpMPZ(B_right, base, exponents2, MODULUS);

        vector<pair<mpz_class, mpz_class> > proofVec_left;
        vector<pair<mpz_class, mpz_class> > proofVec_right;

        pair<mpz_class, mpz_class> proof_left{a_left, B_left};
        pair<mpz_class, mpz_class> proof_right{a_right, B_right};

        BreakNonMemWit(proofVec_left, digest, proof_left, I_left);
        BreakNonMemWit(proofVec_right, digest, proof_right, I_right);

        proofVec.insert(proofVec.end(), proofVec_left.begin(), proofVec_left.end());
        proofVec.insert(proofVec.end(), proofVec_right.begin(), proofVec_right.end());
    }
}

bool Acc::NonMemVerifySingle(const mpz_class &digest, const mpz_class &y, const pair<mpz_class, mpz_class> &pi) {
    mpz_class p1, p2, result;
    mpz_powm(p1.get_mpz_t(), digest.get_mpz_t(), pi.first.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_powm(p2.get_mpz_t(), pi.second.get_mpz_t(), y.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_mul(result.get_mpz_t(), p1.get_mpz_t(), p2.get_mpz_t());
    mpz_mod(result.get_mpz_t(), result.get_mpz_t(), MODULUS.get_mpz_t());
    int status = mpz_cmp(BASE.get_mpz_t(), result.get_mpz_t());
    return status == 0;
}

void Acc::Print() {
    PrintSize(100, "BASE", BASE);
    PrintSize(100, "P", P);
    PrintSize(100, "Q", Q);
    PrintSize(100, "PHI", PHI);
    PrintSize(100, "MODULUS", MODULUS);
}
#endif

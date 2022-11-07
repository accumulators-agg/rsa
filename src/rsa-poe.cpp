#ifndef RSA_POE_CPP
#define RSA_POE_CPP

void Acc::NiPoEProve(mpz_class &Q, const mpz_class &x, const mpz_class &u, const mpz_class &w) {
    mpz_class l, q;
    Hash_to_Prime(l, vector<mpz_class>{x, u, w});

    mpz_tdiv_q(q.get_mpz_t(), x.get_mpz_t(), l.get_mpz_t());
    mpz_powm(Q.get_mpz_t(), u.get_mpz_t(), q.get_mpz_t(), MODULUS.get_mpz_t());
}

bool Acc::NiPoEVerify(const mpz_class &Q, const mpz_class &x, const mpz_class &u, const mpz_class &w) {
    mpz_class l, r;
    Hash_to_Prime(l, vector<mpz_class>{x, u, w});
    mpz_tdiv_r(r.get_mpz_t(), x.get_mpz_t(), l.get_mpz_t());
    mpz_class t1, t2, t;
    mpz_powm(t1.get_mpz_t(), Q.get_mpz_t(), l.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_powm(t2.get_mpz_t(), u.get_mpz_t(), r.get_mpz_t(), MODULUS.get_mpz_t());

    mpz_mul(t.get_mpz_t(), t1.get_mpz_t(), t2.get_mpz_t());
    mpz_mod(t.get_mpz_t(), t.get_mpz_t(), MODULUS.get_mpz_t());

    int status = mpz_cmp(t.get_mpz_t(), w.get_mpz_t());
    return status == 0;
}

void Acc::NiPoKE2Prove(NiPoKE2Pi &pi, const mpz_class &x, const mpz_class &u, const mpz_class &w) {
    mpz_class g;
    Hash_to_Group(g, vector<mpz_class>{u, w}, MODULUS);

    mpz_class z;
    mpz_powm(z.get_mpz_t(), g.get_mpz_t(), x.get_mpz_t(), MODULUS.get_mpz_t());

    mpz_class l;
    Hash_to_Prime(l, vector<mpz_class>{u, w, z});

    mpz_class alpha;
    Hash_to_Group(alpha, vector<mpz_class>{u, w, z, l}, MODULUS);

    mpz_class q, r;
    mpz_tdiv_qr(q.get_mpz_t(), r.get_mpz_t(), x.get_mpz_t(), l.get_mpz_t());

    mpz_class t1, t2, t3;
    mpz_powm(t1.get_mpz_t(), g.get_mpz_t(), alpha.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_mul(t2.get_mpz_t(), u.get_mpz_t(), t1.get_mpz_t());
    mpz_mod(t2.get_mpz_t(), t2.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_powm(t3.get_mpz_t(), t2.get_mpz_t(), q.get_mpz_t(), MODULUS.get_mpz_t());

    mpz_set(pi.z.get_mpz_t(), z.get_mpz_t());
    mpz_set(pi.Q.get_mpz_t(), t3.get_mpz_t());
    mpz_set(pi.r.get_mpz_t(), r.get_mpz_t());
}

bool Acc::NiPoKE2Verify(const NiPoKE2Pi &pi, const mpz_class &u, const mpz_class &w) {
    mpz_class g;
    Hash_to_Group(g, vector<mpz_class>{u, w}, MODULUS);

    mpz_class l;
    Hash_to_Prime(l, vector<mpz_class>{u, w, pi.z});

    mpz_class alpha;
    Hash_to_Group(alpha, vector<mpz_class>{u, w, pi.z, l}, MODULUS);

    mpz_class rhs;
    mpz_powm(rhs.get_mpz_t(), pi.z.get_mpz_t(), alpha.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_mul(rhs.get_mpz_t(), rhs.get_mpz_t(), w.get_mpz_t());
    mpz_mod(rhs.get_mpz_t(), rhs.get_mpz_t(), MODULUS.get_mpz_t());

    mpz_class lhs, lhs1, lhs2;
    mpz_powm(lhs1.get_mpz_t(), pi.Q.get_mpz_t(), l.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_powm(lhs2.get_mpz_t(), g.get_mpz_t(), alpha.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_mul(lhs2.get_mpz_t(), lhs2.get_mpz_t(), u.get_mpz_t());
    mpz_mod(lhs2.get_mpz_t(), lhs2.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_powm(lhs2.get_mpz_t(), lhs2.get_mpz_t(), pi.r.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_mul(lhs.get_mpz_t(), lhs1.get_mpz_t(), lhs2.get_mpz_t());
    mpz_mod(lhs.get_mpz_t(), lhs.get_mpz_t(), MODULUS.get_mpz_t());

    int status = mpz_cmp(lhs.get_mpz_t(), rhs.get_mpz_t());
    return status == 0;
}

void Acc::AggMemProvePoE(mpz_class &proof, mpz_class &poe_pi, mpz_class &I_prod, const mpz_class &digest, const vector<mpz_class> &I, const vector<mpz_class> &proofs) {
    AggMemProve(proof, I_prod, I, proofs);
    NiPoEProve(poe_pi, I_prod, proof, digest);
}

bool Acc::AggMemVerifyPoE(const mpz_class &digest, const vector<mpz_class> &I, const mpz_class &proof, const mpz_class &poe_pi) {
    mpz_class I_prod;
    MulMPZVec(I_prod, I);
    return NiPoEVerify(poe_pi, I_prod, proof, digest);
}

void Acc::AggNonMemProvePoE(NiNonMemPi &pi, mpz_class &I_prod, const mpz_class &digest, const vector<mpz_class> &I, const vector<pair<mpz_class, mpz_class> > &proofs) {
    AggNonMemProve(pi.nonmem_pi, I_prod, digest, I, proofs);
    mpz_powm(pi.V.get_mpz_t(), digest.get_mpz_t(), pi.nonmem_pi.first.get_mpz_t(), MODULUS.get_mpz_t());
    NiPoKE2Prove(pi.poke2_pi, pi.nonmem_pi.first, digest, pi.V);
    mpz_class VInv;
    mpz_invert(VInv.get_mpz_t(), pi.V.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_mul(VInv.get_mpz_t(), VInv.get_mpz_t(), BASE.get_mpz_t());
    mpz_mod(VInv.get_mpz_t(), VInv.get_mpz_t(), MODULUS.get_mpz_t());
    NiPoEProve(pi.poe_pi, I_prod, pi.nonmem_pi.second, VInv);
}

bool Acc::AggNonMemVerifyPoE(const NiNonMemPi pi, const mpz_class &digest, const vector<mpz_class> &I) {
    mpz_class I_prod;
    MulMPZVec(I_prod, I);

    mpz_class VInv;
    mpz_invert(VInv.get_mpz_t(), pi.V.get_mpz_t(), MODULUS.get_mpz_t());
    mpz_mul(VInv.get_mpz_t(), VInv.get_mpz_t(), BASE.get_mpz_t());
    mpz_mod(VInv.get_mpz_t(), VInv.get_mpz_t(), MODULUS.get_mpz_t());

    bool status;
    status = NiPoKE2Verify(pi.poke2_pi, digest, pi.V);
    status = status && NiPoEVerify(pi.poe_pi, I_prod, pi.nonmem_pi.second, VInv);
    return status;
}
#endif

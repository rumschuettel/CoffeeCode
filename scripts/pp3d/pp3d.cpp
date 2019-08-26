// postprocess-3dplot.py port to c++ for extended precision
// requires:
// conda install -c conda-forge nlohmann_json

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#ifdef F256
#include <boost/multiprecision/cpp_bin_float.hpp> 
using scalar = boost::multiprecision::number<
    boost::multiprecision::backends::cpp_bin_float<
        256,
        boost::multiprecision::backends::digit_base_2
    >
>;
#else
#include <boost/multiprecision/float128.hpp>
using scalar = boost::multiprecision::float128;
#endif
using integer = int64_t;

#include <iostream>
#include <assert.h>
#include <cmath>
#include <chrono>

using std::pow;
using std::log2;
using std::abs;

template<typename T>
scalar shannon_entropy_from_lambda(const T& lambda_list, const integer KSYS, const scalar fac, const scalar q1, const scalar q2, const scalar q3)
{
    scalar entropy = 0.;

    std::vector<scalar> pow_q1_lu(KSYS+1), pow_q2_lu(KSYS+1), pow_q3_lu(KSYS+1), pow_1_q1q2q3_lu(KSYS+1);
    for (size_t i = 0; i <= KSYS; i ++) {
        pow_q1_lu[i] = pow(q1, i);
        pow_q2_lu[i] = pow(q2, i);
        pow_q3_lu[i] = pow(q3, i);
        pow_1_q1q2q3_lu[i] = pow(static_cast<scalar>(1.) - q1 - q2 - q3, i);
    }

    for (const auto& _mult_poly : lambda_list) {
        const integer mult = _mult_poly[0];
        const auto& poly = _mult_poly[1];

        if (poly.size() == 0)
            continue;

        scalar term = 0.;

        for (const auto& _coeff_e1e2e3 : poly) {
            const integer coeff = _coeff_e1e2e3[0];
            const integer e1 = _coeff_e1e2e3[1][0];
            const integer e2 = _coeff_e1e2e3[1][1];
            const integer e3 = _coeff_e1e2e3[1][2];

            term += fac * coeff * pow_q1_lu[e1] * pow_q2_lu[e2] * pow_q3_lu[e3] * pow_1_q1q2q3_lu[KSYS - e1 - e2 - e3];
        }

        if (term != 0.)
            entropy -= mult * term * log2(term);
    }

    return entropy;
}

int main(int argc, char *argv[])
{
    // arguments
    assert(argc == 3);
    const int KSYS = atoi(argv[1]);
    const int KENV = atoi(argv[2]);

    std::cerr << "PP3D kSys=" << KSYS << " kEnv=" << KENV << " running... ";

    // deserialize from standard input
    json in;
    std::cin >> in;
    const size_t QCOUNT = in["qs"].size();

    // iterate over all qs; the double is not a typo since this just stores the output
    std::vector<double> S_a_minus_S(QCOUNT);

    const auto& lambda = in["lambda"];
    const auto& lambda_a = in["lambda_a"];
    const scalar fac_a = pow(static_cast<scalar>(2.0), -KENV);

    auto start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for (size_t i = 0; i < QCOUNT; i ++) {
        const scalar q1 = static_cast<double>(in["qs"][i][0]);
        const scalar q2 = static_cast<double>(in["qs"][i][1]);
        const scalar q3 = static_cast<double>(in["qs"][i][2]);
        
        const scalar S = shannon_entropy_from_lambda(in["lambda"], KSYS, 1.0, q1, q2, q3);
        const scalar S_a = shannon_entropy_from_lambda(in["lambda_a"], KSYS, fac_a, q1, q2, q3);

        S_a_minus_S[i] = static_cast<double>(S_a - S);
    }

    auto time = std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::high_resolution_clock::now() - start
    ).count();
    std::cerr << "done; T=" << time << " sec\n";

    // serialize to standard output
    json out;
    out["S_a-S"] = S_a_minus_S;
    std::cout << out << std::endl;

    return 0;
}
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// indval_cpp  –  single-partition IndVal (mean-abundance variant)
// Replaces R function indval(grp, spp, allscores=FALSE)
// grp: integer vector of 1s and 2s (length = nrow(spp))
// spp: numeric matrix (numUnit x numTxa)
// Returns: 2-column matrix [maxcls, ivmax]
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix indval_cpp(IntegerVector grp, NumericMatrix spp) {
  int numUnit = spp.nrow();
  int numTxa = spp.ncol();

  // count group sizes
  int n1 = 0, n2 = 0;
  for (int i = 0; i < numUnit; i++) {
    if (grp[i] == 1) n1++;
    else n2++;
  }

  NumericVector mnabnd1(numTxa, 0.0), mnabnd2(numTxa, 0.0);
  NumericVector frq1(numTxa, 0.0), frq2(numTxa, 0.0);

  for (int j = 0; j < numTxa; j++) {
    double sum1 = 0.0, sum2 = 0.0;
    double cnt1 = 0.0, cnt2 = 0.0;
    for (int i = 0; i < numUnit; i++) {
      double val = spp(i, j);
      if (grp[i] == 1) {
        sum1 += val;
        if (val > 0.0) cnt1 += 1.0;
      } else {
        sum2 += val;
        if (val > 0.0) cnt2 += 1.0;
      }
    }
    mnabnd1[j] = (n1 > 0) ? sum1 / n1 : 0.0;
    mnabnd2[j] = (n2 > 0) ? sum2 / n2 : 0.0;
    frq1[j] = (n1 > 0) ? cnt1 / n1 : 0.0;
    frq2[j] = (n2 > 0) ? cnt2 / n2 : 0.0;
  }

  NumericMatrix ivout(numTxa, 2);
  for (int j = 0; j < numTxa; j++) {
    double total = mnabnd1[j] + mnabnd2[j];
    double ra1 = (mnabnd1[j] == 0.0 || total == 0.0) ? 0.0 : mnabnd1[j] / total;
    double ra2 = (mnabnd2[j] == 0.0 || total == 0.0) ? 0.0 : mnabnd2[j] / total;
    double iv1 = ra1 * frq1[j] * 100.0;
    double iv2 = ra2 * frq2[j] * 100.0;
    double ivmax = (iv1 >= iv2) ? iv1 : iv2;
    int maxcls = 0;
    if (iv1 > iv2) maxcls = 1;
    else if (iv2 > iv1) maxcls = 2;
    ivout(j, 0) = maxcls;
    ivout(j, 1) = ivmax;
  }
  return ivout;
}


// ---------------------------------------------------------------------------
// indvals_cpp  –  single-partition IndVal (total-abundance variant)
// Replaces R function indvals(grp, spp, allscores=FALSE)
// BUG FIX: line 150 of R code sets frq1 instead of frq2 — corrected here
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix indvals_cpp(IntegerVector grp, NumericMatrix spp) {
  int numUnit = spp.nrow();
  int numTxa = spp.ncol();

  int n1 = 0, n2 = 0;
  for (int i = 0; i < numUnit; i++) {
    if (grp[i] == 1) n1++;
    else n2++;
  }

  NumericVector smabnd1(numTxa, 0.0), smabnd2(numTxa, 0.0);
  NumericVector frq1(numTxa, 0.0), frq2(numTxa, 0.0);

  for (int j = 0; j < numTxa; j++) {
    double sum1 = 0.0, sum2 = 0.0;
    double cnt1 = 0.0, cnt2 = 0.0;
    for (int i = 0; i < numUnit; i++) {
      double val = spp(i, j);
      if (grp[i] == 1) {
        sum1 += val;
        if (val > 0.0) cnt1 += 1.0;
      } else {
        sum2 += val;
        if (val > 0.0) cnt2 += 1.0;
      }
    }
    smabnd1[j] = sum1;
    smabnd2[j] = sum2;
    frq1[j] = (n1 > 0) ? cnt1 / n1 : 0.0;
    frq2[j] = (n2 > 0) ? cnt2 / n2 : 0.0;  // BUG FIX: R code had frq1 here
  }

  NumericMatrix ivout(numTxa, 2);
  for (int j = 0; j < numTxa; j++) {
    double total = smabnd1[j] + smabnd2[j];
    double ra1 = (smabnd1[j] == 0.0 || total == 0.0) ? 0.0 : smabnd1[j] / total;
    double ra2 = (smabnd2[j] == 0.0 || total == 0.0) ? 0.0 : smabnd2[j] / total;
    double iv1 = ra1 * frq1[j] * 100.0;
    double iv2 = ra2 * frq2[j] * 100.0;
    double ivmax = (iv1 >= iv2) ? iv1 : iv2;
    int maxcls = 0;
    if (iv1 > iv2) maxcls = 1;
    else if (iv2 > iv1) maxcls = 2;
    ivout(j, 0) = maxcls;
    ivout(j, 1) = ivmax;
  }
  return ivout;
}


// ---------------------------------------------------------------------------
// indvalp_cpp  –  matrix-based permuted IndVal (mean-abundance variant)
// Replaces R function indvalp(grpMatrix, spp)
// grpMatrix: numUnit x nPerm matrix of 1s and 2s
// spp: numUnit x numTxa matrix
// Returns: nPerm x numTxa matrix of IndVal maxima
// Uses Armadillo for the matrix multiplies
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix indvalp_cpp(IntegerMatrix grpMatrix, NumericMatrix spp) {
  int numUnit = spp.nrow();
  int numTxa = spp.ncol();
  int nPerm = grpMatrix.ncol();

  // Convert to arma matrices for fast multiply
  arma::mat Spp(spp.begin(), numUnit, numTxa, false);

  // Build group indicator matrices
  arma::mat G1(numUnit, nPerm);
  arma::mat G2(numUnit, nPerm);
  for (int p = 0; p < nPerm; p++) {
    for (int i = 0; i < numUnit; i++) {
      G1(i, p) = (grpMatrix(i, p) == 1) ? 1.0 : 0.0;
      G2(i, p) = (grpMatrix(i, p) == 2) ? 1.0 : 0.0;
    }
  }

  // group sizes (column sums)
  arma::rowvec n1 = arma::sum(G1, 0); // 1 x nPerm
  arma::rowvec n2 = arma::sum(G2, 0);

  // mean operator matrices: (1/n) * t(G) -> nPerm x numUnit
  arma::mat Op1 = G1.t();
  arma::mat Op2 = G2.t();
  for (int p = 0; p < nPerm; p++) {
    double inv1 = (n1(p) > 0) ? 1.0 / n1(p) : 0.0;
    double inv2 = (n2(p) > 0) ? 1.0 / n2(p) : 0.0;
    Op1.row(p) *= inv1;
    Op2.row(p) *= inv2;
  }

  // mean abundances: nPerm x numTxa
  arma::mat mnabnd1 = Op1 * Spp;
  arma::mat mnabnd2 = Op2 * Spp;

  // presence/absence
  arma::mat SppPA(numUnit, numTxa);
  for (int j = 0; j < numTxa; j++) {
    for (int i = 0; i < numUnit; i++) {
      SppPA(i, j) = (spp(i, j) > 0) ? 1.0 : 0.0;
    }
  }
  arma::mat frq1 = Op1 * SppPA;
  arma::mat frq2 = Op2 * SppPA;

  // relative abundance and IndVal
  arma::mat abndSum = mnabnd1 + mnabnd2;
  arma::mat ra1 = mnabnd1 / abndSum;
  arma::mat ra2 = mnabnd2 / abndSum;

  // zero out where abundance was zero (handle NaN from 0/0)
  ra1.elem(arma::find(mnabnd1 == 0.0)).zeros();
  ra2.elem(arma::find(mnabnd2 == 0.0)).zeros();
  // also handle NaN from division where both groups have 0
  ra1.replace(arma::datum::nan, 0.0);
  ra2.replace(arma::datum::nan, 0.0);

  arma::mat iv1 = ra1 % frq1 * 100.0;
  arma::mat iv2 = ra2 % frq2 * 100.0;
  arma::mat ivmax = arma::max(iv1, iv2);

  // Convert back to Rcpp
  NumericMatrix result(nPerm, numTxa);
  for (int p = 0; p < nPerm; p++) {
    for (int j = 0; j < numTxa; j++) {
      result(p, j) = ivmax(p, j);
    }
  }
  return result;
}


// ---------------------------------------------------------------------------
// indvalps_cpp  –  matrix-based permuted IndVal (total-abundance variant)
// Replaces R function indvalps(grpMatrix, spp)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix indvalps_cpp(IntegerMatrix grpMatrix, NumericMatrix spp) {
  int numUnit = spp.nrow();
  int numTxa = spp.ncol();
  int nPerm = grpMatrix.ncol();

  arma::mat Spp(spp.begin(), numUnit, numTxa, false);

  arma::mat G1(numUnit, nPerm);
  arma::mat G2(numUnit, nPerm);
  for (int p = 0; p < nPerm; p++) {
    for (int i = 0; i < numUnit; i++) {
      G1(i, p) = (grpMatrix(i, p) == 1) ? 1.0 : 0.0;
      G2(i, p) = (grpMatrix(i, p) == 2) ? 1.0 : 0.0;
    }
  }

  arma::rowvec n1 = arma::sum(G1, 0);
  arma::rowvec n2 = arma::sum(G2, 0);

  // Mean operator for frequency only
  arma::mat Op1 = G1.t();
  arma::mat Op2 = G2.t();
  for (int p = 0; p < nPerm; p++) {
    double inv1 = (n1(p) > 0) ? 1.0 / n1(p) : 0.0;
    double inv2 = (n2(p) > 0) ? 1.0 / n2(p) : 0.0;
    Op1.row(p) *= inv1;
    Op2.row(p) *= inv2;
  }

  // Total abundances (not mean): t(G) * spp
  arma::mat smabnd1 = G1.t() * Spp;
  arma::mat smabnd2 = G2.t() * Spp;

  // Frequencies (mean of presence)
  arma::mat SppPA(numUnit, numTxa);
  for (int j = 0; j < numTxa; j++) {
    for (int i = 0; i < numUnit; i++) {
      SppPA(i, j) = (spp(i, j) > 0) ? 1.0 : 0.0;
    }
  }
  arma::mat frq1 = Op1 * SppPA;
  arma::mat frq2 = Op2 * SppPA;

  arma::mat abndSum = smabnd1 + smabnd2;
  arma::mat ra1 = smabnd1 / abndSum;
  arma::mat ra2 = smabnd2 / abndSum;

  ra1.elem(arma::find(smabnd1 == 0.0)).zeros();
  ra2.elem(arma::find(smabnd2 == 0.0)).zeros();
  ra1.replace(arma::datum::nan, 0.0);
  ra2.replace(arma::datum::nan, 0.0);

  arma::mat iv1 = ra1 % frq1 * 100.0;
  arma::mat iv2 = ra2 % frq2 * 100.0;
  arma::mat ivmax = arma::max(iv1, iv2);

  NumericMatrix result(nPerm, numTxa);
  for (int p = 0; p < nPerm; p++) {
    for (int j = 0; j < numTxa; j++) {
      result(p, j) = ivmax(p, j);
    }
  }
  return result;
}


// ---------------------------------------------------------------------------
// Helper: compute permuted IndVals using BLAS matrix multiply (mean abundance)
// grpRandMatrix: numUnit x nPerm matrix of 1s and 2s (permuted group labels)
// Spp: numUnit x numTxa arma matrix
// SppPA: numUnit x numTxa presence/absence arma matrix
// Returns: nPerm x numTxa matrix of IndVal maxima
// ---------------------------------------------------------------------------
static arma::mat permiv_matmul_mean(const arma::mat& grpRandMatrix,
                                    const arma::mat& Spp,
                                    const arma::mat& SppPA,
                                    int numUnit, int nPerm, int numTxa) {
  // Build group indicator matrices
  arma::mat G1 = arma::conv_to<arma::mat>::from(grpRandMatrix == 1.0);
  arma::mat G2 = arma::conv_to<arma::mat>::from(grpRandMatrix == 2.0);

  arma::rowvec n1 = arma::sum(G1, 0);
  arma::rowvec n2 = arma::sum(G2, 0);

  // Mean operator: (1/n) * t(G) -> nPerm x numUnit
  arma::mat Op1 = G1.t();
  arma::mat Op2 = G2.t();
  for (int p = 0; p < nPerm; p++) {
    Op1.row(p) *= (n1(p) > 0) ? 1.0 / n1(p) : 0.0;
    Op2.row(p) *= (n2(p) > 0) ? 1.0 / n2(p) : 0.0;
  }

  // BLAS matrix multiply: nPerm x numTxa
  arma::mat mnabnd1 = Op1 * Spp;
  arma::mat mnabnd2 = Op2 * Spp;
  arma::mat frq1 = Op1 * SppPA;
  arma::mat frq2 = Op2 * SppPA;

  arma::mat abndSum = mnabnd1 + mnabnd2;
  arma::mat ra1 = mnabnd1 / abndSum;
  arma::mat ra2 = mnabnd2 / abndSum;

  ra1.elem(arma::find(mnabnd1 == 0.0)).zeros();
  ra2.elem(arma::find(mnabnd2 == 0.0)).zeros();
  ra1.replace(arma::datum::nan, 0.0);
  ra2.replace(arma::datum::nan, 0.0);

  return arma::max(ra1 % frq1, ra2 % frq2) * 100.0;
}

// ---------------------------------------------------------------------------
// Helper: compute permuted IndVals using BLAS matrix multiply (total abundance)
// ---------------------------------------------------------------------------
static arma::mat permiv_matmul_total(const arma::mat& grpRandMatrix,
                                     const arma::mat& Spp,
                                     const arma::mat& SppPA,
                                     int numUnit, int nPerm, int numTxa) {
  arma::mat G1 = arma::conv_to<arma::mat>::from(grpRandMatrix == 1.0);
  arma::mat G2 = arma::conv_to<arma::mat>::from(grpRandMatrix == 2.0);

  arma::rowvec n1 = arma::sum(G1, 0);
  arma::rowvec n2 = arma::sum(G2, 0);

  // Mean operator for frequency
  arma::mat Op1 = G1.t();
  arma::mat Op2 = G2.t();
  for (int p = 0; p < nPerm; p++) {
    Op1.row(p) *= (n1(p) > 0) ? 1.0 / n1(p) : 0.0;
    Op2.row(p) *= (n2(p) > 0) ? 1.0 / n2(p) : 0.0;
  }

  // Total abundance (not mean)
  arma::mat smabnd1 = G1.t() * Spp;
  arma::mat smabnd2 = G2.t() * Spp;
  arma::mat frq1 = Op1 * SppPA;
  arma::mat frq2 = Op2 * SppPA;

  arma::mat abndSum = smabnd1 + smabnd2;
  arma::mat ra1 = smabnd1 / abndSum;
  arma::mat ra2 = smabnd2 / abndSum;

  ra1.elem(arma::find(smabnd1 == 0.0)).zeros();
  ra2.elem(arma::find(smabnd2 == 0.0)).zeros();
  ra1.replace(arma::datum::nan, 0.0);
  ra2.replace(arma::datum::nan, 0.0);

  return arma::max(ra1 % frq1, ra2 % frq2) * 100.0;
}


// ---------------------------------------------------------------------------
// getivz_cpp  –  THE MAIN WIN
// Fuses observed IndVal + permutation + z-score + p-value into one C++ call.
// Uses BLAS-backed matrix multiply for the permutation step.
// Processes one partition at a time to avoid the full 3D array.
//
// clss: numUnit x numClass matrix of 1s and 2s
// spp: numUnit x numTxa matrix
// Returns: (4*numTxa) x numClass matrix = rbind(obsmx, ivz, obsiv, iv.p)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix getivz_cpp(IntegerMatrix clss, NumericMatrix spp,
                         bool ivTot, int nPerm, int numClass, bool imax) {
  int numUnit = spp.nrow();
  int numTxa = spp.ncol();

  // Convert spp to arma once
  arma::mat Spp(spp.begin(), numUnit, numTxa, false);
  arma::mat SppPA = arma::conv_to<arma::mat>::from(Spp > 0.0);

  // Output matrices
  NumericMatrix obsmx(numTxa, numClass);
  NumericMatrix obsiv(numTxa, numClass);
  NumericMatrix ivz(numTxa, numClass);
  NumericMatrix ivp(numTxa, numClass);

  // Process one partition at a time
  for (int c = 0; c < numClass; c++) {
    // Extract group vector for this partition
    IntegerVector grp(numUnit);
    for (int i = 0; i < numUnit; i++) {
      grp[i] = clss(i, c);
    }

    // --- Observed IndVal (scalar, single partition) ---
    int n1 = 0, n2 = 0;
    for (int i = 0; i < numUnit; i++) {
      if (grp[i] == 1) n1++; else n2++;
    }

    for (int j = 0; j < numTxa; j++) {
      double sum1 = 0.0, sum2 = 0.0, cnt1 = 0.0, cnt2 = 0.0;
      for (int i = 0; i < numUnit; i++) {
        double val = spp(i, j);
        if (grp[i] == 1) {
          sum1 += val;
          if (val > 0.0) cnt1 += 1.0;
        } else {
          sum2 += val;
          if (val > 0.0) cnt2 += 1.0;
        }
      }

      double abnd1, abnd2;
      if (ivTot) {
        abnd1 = sum1; abnd2 = sum2;
      } else {
        abnd1 = (n1 > 0) ? sum1 / n1 : 0.0;
        abnd2 = (n2 > 0) ? sum2 / n2 : 0.0;
      }
      double f1 = (n1 > 0) ? cnt1 / n1 : 0.0;
      double f2 = (n2 > 0) ? cnt2 / n2 : 0.0;

      double total = abnd1 + abnd2;
      double ra1 = (abnd1 == 0.0 || total == 0.0) ? 0.0 : abnd1 / total;
      double ra2 = (abnd2 == 0.0 || total == 0.0) ? 0.0 : abnd2 / total;

      double iv1 = ra1 * f1 * 100.0;
      double iv2 = ra2 * f2 * 100.0;
      obsiv(j, c) = std::max(iv1, iv2);

      int maxcls = 0;
      if (iv1 > iv2) maxcls = 1;
      else if (iv2 > iv1) maxcls = 2;
      obsmx(j, c) = maxcls;
    }

    // --- Permuted IndVal scores using BLAS matrix multiply ---
    // Generate nPerm permutations as a numUnit x nPerm matrix
    arma::mat grpRandMatrix(numUnit, nPerm);
    for (int p = 0; p < nPerm; p++) {
      IntegerVector pgrp = Rcpp::sample(grp, numUnit, false);
      for (int i = 0; i < numUnit; i++) {
        grpRandMatrix(i, p) = (double)pgrp[i];
      }
    }

    // Compute all permuted IndVals at once via matrix multiply
    arma::mat piv; // nPerm x numTxa
    if (ivTot) {
      piv = permiv_matmul_total(grpRandMatrix, Spp, SppPA, numUnit, nPerm, numTxa);
    } else {
      piv = permiv_matmul_mean(grpRandMatrix, Spp, SppPA, numUnit, nPerm, numTxa);
    }

    // --- Compute z-scores and p-values ---
    arma::rowvec mniv = arma::mean(piv, 0); // 1 x numTxa
    arma::rowvec sdiv = arma::stddev(piv, 0, 0); // 1 x numTxa, norm_type=0 (N-1)

    for (int j = 0; j < numTxa; j++) {
      double sd_j = sdiv(j);
      if (sd_j < 0.1) sd_j = 1.0;

      ivz(j, c) = (obsiv(j, c) - mniv(j)) / sd_j;

      // p-value: fraction of perm values >= observed, +1 as in R code
      double obs_val = obsiv(j, c);
      int count = 0;
      for (int p = 0; p < nPerm; p++) {
        if (piv(p, j) >= obs_val) count++;
      }
      ivp(j, c) = ((double)(count + 1)) / (double)nPerm;
    }
  }

  // Stack: rbind(obsmx, ivz, obsiv, iv.p)
  NumericMatrix result(4 * numTxa, numClass);
  for (int c = 0; c < numClass; c++) {
    for (int j = 0; j < numTxa; j++) {
      result(j, c) = obsmx(j, c);
      result(j + numTxa, c) = ivz(j, c);
      result(j + 2 * numTxa, c) = obsiv(j, c);
      result(j + 3 * numTxa, c) = ivp(j, c);
    }
  }
  return result;
}


// ---------------------------------------------------------------------------
// env_part_cpp  –  partition loop
// Replaces the for-loop in env.part() that creates the eclass matrix
// rankEnv: integer vector of ranks (1..numUnit)
// srtEnv2: sorted rankEnv
// numUnit: number of sites
// minSplt: minimum split size
// Returns: numUnit x numClass integer matrix of 1s and 2s
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
IntegerMatrix env_part_cpp(IntegerVector rankEnv, IntegerVector srtEnv2,
                           int numUnit, int minSplt) {
  int numClass = numUnit - 2 * minSplt + 1;
  IntegerMatrix eclass(numUnit, numClass);

  for (int c = 0; c < numClass; c++) {
    int threshold = srtEnv2[minSplt - 1 + c]; // R is 1-indexed, srtEnv2 passed as 0-indexed C++ vector
    for (int i = 0; i < numUnit; i++) {
      eclass(i, c) = (rankEnv[i] > threshold) ? 2 : 1;
    }
  }
  return eclass;
}

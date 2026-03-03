// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Forward declare getivz_cpp from titan_core.cpp
NumericMatrix getivz_cpp(IntegerMatrix clss, NumericMatrix spp,
                         bool ivTot, int nPerm, int numClass, bool imax);

// ---------------------------------------------------------------------------
// tboot_cpp  –  single bootstrap replicate
// Replaces the core of R function tboot()
// Performs: resample -> rank -> partition -> getivz -> NaN cleanup -> extract
//
// env: numeric vector of environmental values
// taxa: numUnit x numTxa matrix
// ivTot, minSplt, nPerm, imax: TITAN parameters
//
// Returns a List of 4: bt.metrics, ivzs, bSrti, rspdr
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List tboot_cpp(NumericVector env, NumericMatrix taxa,
               bool ivTot, int minSplt, int nPerm, bool imax) {

  int numUnit = env.size();
  int numTxa = taxa.ncol();

  // Resample with replacement using R's RNG
  IntegerVector row_idx(numUnit);
  for (int i = 0; i < numUnit; i++) {
    row_idx[i] = (int)(R::runif(0.0, (double)numUnit));
    if (row_idx[i] >= numUnit) row_idx[i] = numUnit - 1; // safety
  }

  // Build resampled taxa and env
  NumericMatrix booted_taxa(numUnit, numTxa);
  NumericVector permuted_env(numUnit);
  for (int i = 0; i < numUnit; i++) {
    permuted_env[i] = env[row_idx[i]];
    for (int j = 0; j < numTxa; j++) {
      booted_taxa(i, j) = taxa(row_idx[i], j);
    }
  }

  // Rank using R's rank function with ties.method="random"
  Function r_rank("rank");
  NumericVector bRankEnv = r_rank(permuted_env, Named("ties.method") = "random");

  // Sort env
  NumericVector bSrti = clone(permuted_env);
  std::sort(bSrti.begin(), bSrti.end());

  // Sort ranks
  IntegerVector bRankEnvInt(numUnit);
  for (int i = 0; i < numUnit; i++) {
    bRankEnvInt[i] = (int)bRankEnv[i];
  }
  IntegerVector bSrti2 = clone(bRankEnvInt);
  std::sort(bSrti2.begin(), bSrti2.end());

  // Partition
  int nClass = numUnit - 2 * minSplt + 1;
  IntegerMatrix boot_env(numUnit, nClass);
  for (int c = 0; c < nClass; c++) {
    int threshold = bSrti2[minSplt - 1 + c];
    for (int i = 0; i < numUnit; i++) {
      boot_env(i, c) = (bRankEnvInt[i] > threshold) ? 2 : 1;
    }
  }

  // Call getivz_cpp
  NumericMatrix ivzScores_bt = getivz_cpp(boot_env, booted_taxa, ivTot, nPerm, nClass, imax);

  // Replace NaNs with 0s in rows (numTxa+1) to (3*numTxa) (z-scores and IndVals)
  for (int k = numTxa; k < 3 * numTxa; k++) {
    for (int c = 0; c < nClass; c++) {
      if (R_IsNaN(ivzScores_bt(k, c))) {
        ivzScores_bt(k, c) = 0.0;
      }
    }
  }

  // Extract bt.metrics (numTxa x 4)
  NumericMatrix bt_metrics(numTxa, 4);
  for (int j = 0; j < numTxa; j++) {
    int maxSplt = 0;
    double maxVal = -1e30;
    if (!imax) {
      // z-score maxima
      for (int c = 0; c < nClass; c++) {
        double val = std::abs(ivzScores_bt(j + numTxa, c));
        if (val > maxVal) { maxVal = val; maxSplt = c; }
      }
    } else {
      // IndVal maxima
      for (int c = 0; c < nClass; c++) {
        double val = std::abs(ivzScores_bt(j + 2 * numTxa, c));
        if (val > maxVal) { maxVal = val; maxSplt = c; }
      }
    }
    bt_metrics(j, 0) = ivzScores_bt(j, maxSplt);           // maxgrp
    bt_metrics(j, 1) = (bSrti[minSplt - 1 + maxSplt] +
                         bSrti[minSplt + maxSplt]) / 2.0;   // env.cp
    bt_metrics(j, 2) = ivzScores_bt(j + numTxa, maxSplt);   // z.score
    bt_metrics(j, 3) = ivzScores_bt(j + 3 * numTxa, maxSplt); // iv.p
  }

  // Extract z-scores matrix (numTxa x nClass)
  NumericMatrix ivzs(numTxa, nClass);
  for (int j = 0; j < numTxa; j++) {
    for (int c = 0; c < nClass; c++) {
      ivzs(j, c) = ivzScores_bt(j + numTxa, c);
    }
  }

  // Extract response direction matrix (numTxa x nClass)
  NumericMatrix rspdr(numTxa, nClass);
  for (int j = 0; j < numTxa; j++) {
    for (int c = 0; c < nClass; c++) {
      rspdr(j, c) = ivzScores_bt(j, c);
    }
  }

  return List::create(bt_metrics, ivzs, bSrti, rspdr);
}


// ---------------------------------------------------------------------------
// small_boot_sumz_cpp  –  sum(z) maximization loop from small.boot()
// Replaces lines 440-461 of boot.functions.R
//
// zArray: numTxa x numClass x nBoot array (as 3D)
// rspdir: numTxa x numClass x nBoot array
// bEnvMatrix: nBoot x numUnit matrix
// sppSub1, sppSub2: integer vectors of 1-indexed taxon indices
// nBoot, minSplt: parameters
//
// Returns List of 4: maxSumz, maxFsumz (nBoot x 2 each)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List small_boot_sumz_cpp(NumericVector zArray_flat,
                         NumericVector rspdir_flat,
                         NumericMatrix bEnvMatrix,
                         IntegerVector sppSub1,
                         IntegerVector sppSub2,
                         int nTaxa, int numClass, int nBoot,
                         int minSplt) {

  // Access 3D arrays as flat vectors with manual indexing
  // Dimensions: nTaxa x numClass x nBoot (R column-major)
  // Index: [j, c, b] = j + nTaxa * c + nTaxa * numClass * b

  NumericMatrix maxSumz(nBoot, 2);
  NumericMatrix maxFsumz(nBoot, 2);

  int nSub1 = sppSub1.size();
  int nSub2 = sppSub2.size();

  for (int b = 0; b < nBoot; b++) {
    int base = nTaxa * numClass * b;

    // --- sumz (unfiltered) ---
    // direction 1
    {
      double best = -1e30;
      int bestc = 0;
      for (int c = 0; c < numClass; c++) {
        double s = 0.0;
        for (int j = 0; j < nTaxa; j++) {
          int idx = j + nTaxa * c + base;
          double zval = zArray_flat[idx];
          double rval = rspdir_flat[idx];
          if (!R_IsNA(zval) && rval == 1.0) {
            s += zval;
          }
        }
        if (s > best) { best = s; bestc = c; }
      }
      // R uses which.max (1-indexed), store env cp
      maxSumz(b, 0) = (bEnvMatrix(b, minSplt - 1 + bestc) +
                        bEnvMatrix(b, minSplt + bestc)) / 2.0;
    }

    // direction 2
    {
      double best = -1e30;
      int bestc = 0;
      for (int c = 0; c < numClass; c++) {
        double s = 0.0;
        for (int j = 0; j < nTaxa; j++) {
          int idx = j + nTaxa * c + base;
          double zval = zArray_flat[idx];
          double rval = rspdir_flat[idx];
          if (!R_IsNA(zval) && rval == 2.0) {
            s += zval;
          }
        }
        if (s > best) { best = s; bestc = c; }
      }
      maxSumz(b, 1) = (bEnvMatrix(b, minSplt - 1 + bestc) +
                        bEnvMatrix(b, minSplt + bestc)) / 2.0;
    }

    // --- filtered sumz ---
    // direction 1
    if (nSub1 > 1) {
      double best = -1e30;
      int bestc = 0;
      for (int c = 0; c < numClass; c++) {
        double s = 0.0;
        for (int k = 0; k < nSub1; k++) {
          int j = sppSub1[k] - 1; // convert to 0-indexed
          int idx = j + nTaxa * c + base;
          double zval = zArray_flat[idx];
          double rval = rspdir_flat[idx];
          if (!R_IsNA(zval) && rval == 1.0) {
            s += zval;
          }
        }
        if (s > best) { best = s; bestc = c; }
      }
      maxFsumz(b, 0) = (bEnvMatrix(b, minSplt - 1 + bestc) +
                         bEnvMatrix(b, minSplt + bestc)) / 2.0;
    } else if (nSub1 > 0) {
      int j = sppSub1[0] - 1;
      double best = -1e30;
      int bestc = 0;
      for (int c = 0; c < numClass; c++) {
        int idx = j + nTaxa * c + base;
        double zval = zArray_flat[idx];
        double rval = rspdir_flat[idx];
        double val = 0.0;
        if (!R_IsNA(zval) && rval == 1.0) val = zval;
        if (val > best) { best = val; bestc = c; }
      }
      maxFsumz(b, 0) = (bEnvMatrix(b, minSplt - 1 + bestc) +
                         bEnvMatrix(b, minSplt + bestc)) / 2.0;
    }

    // direction 2
    if (nSub2 > 1) {
      double best = -1e30;
      int bestc = 0;
      for (int c = 0; c < numClass; c++) {
        double s = 0.0;
        for (int k = 0; k < nSub2; k++) {
          int j = sppSub2[k] - 1;
          int idx = j + nTaxa * c + base;
          double zval = zArray_flat[idx];
          double rval = rspdir_flat[idx];
          if (!R_IsNA(zval) && rval == 2.0) {
            s += zval;
          }
        }
        if (s > best) { best = s; bestc = c; }
      }
      maxFsumz(b, 1) = (bEnvMatrix(b, minSplt - 1 + bestc) +
                         bEnvMatrix(b, minSplt + bestc)) / 2.0;
    } else if (nSub2 > 0) {
      int j = sppSub2[0] - 1;
      double best = -1e30;
      int bestc = 0;
      for (int c = 0; c < numClass; c++) {
        int idx = j + nTaxa * c + base;
        double zval = zArray_flat[idx];
        double rval = rspdir_flat[idx];
        double val = 0.0;
        if (!R_IsNA(zval) && rval == 2.0) val = zval;
        if (val > best) { best = val; bestc = c; }
      }
      maxFsumz(b, 1) = (bEnvMatrix(b, minSplt - 1 + bestc) +
                         bEnvMatrix(b, minSplt + bestc)) / 2.0;
    }
  }

  return List::create(Named("maxSumz") = maxSumz,
                      Named("maxFsumz") = maxFsumz);
}

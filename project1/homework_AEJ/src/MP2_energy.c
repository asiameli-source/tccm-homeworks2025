#include "MP2_energy.h"
// Eight permutation symmetry
double get_eri(int p, int q, int r, int s, int64_t nint, const int32_t* index, const double* value) {

  for (int64_t n = 0; n < nint; n++) {
    int i = index[4*n + 0];
    int j = index[4*n + 1];
    int k = index[4*n + 2];
    int l = index[4*n + 3];
    double v = value[n];

    if (
      (i==p && j==q && k==r && l==s) ||  // <ij|kl>
      (i==p && l==q && k==r && j==s) ||  // <il|kj>
      (k==p && l==q && i==r && j==s) ||  // <kl|ij>
      (k==p && j==q && i==r && l==s) ||  // <kj|il>
      (j==p && i==q && l==r && k==s) ||  // <ji|lk>
      (l==p && i==q && j==r && k==s) ||  // <li|jk>
      (l==p && k==q && j==r && i==s) ||  // <lk|ji>
      (j==p && k==q && l==r && i==s)     // <jk|li>
    ) {
      return v;
    }
  }
 return 0.0;
}

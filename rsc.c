#include "spherical.h"

static void check_for_ip_rights();

void compute_psidot_rsc(REAL* psidot, REAL* psi, param_list_t *param) {
  REAL a2[3][3],aa2[9];
  REAL da2[3][3];
  REAL collapsed_da2[3][3];
  REAL evec[9], eval[3];
  REAL inner;
  int i,j,k;

  check_for_ip_rights();

  tensor2(psi,a2);
  for (i=0;i<3;i++) for (j=0;j<3;j++)
    aa2[i*3+j] = a2[i][j];
  tensor2(psidot,da2);

/* collapsed_da2 = M:da2, where M = sum_{i=1}^3 e_i e_i e_i e_i, where
   e_i are the normalized eigenvectors of a2. */
  diagonalize_sym(3, aa2, eval, evec);
  memset(collapsed_da2,0,sizeof(collapsed_da2));
  for (i=0;i<3;i++) {
    inner = 0;
    for (j=0;j<3;j++) for (k=0;k<3;k++)
      inner += evec[i*3+j]*evec[i*3+k]*da2[k][j];
    for (j=0;j<3;j++) for (k=0;k<3;k++)
      collapsed_da2[j][k] += evec[i*3+j]*evec[i*3+k]*inner;
  }

/* Add to psidot: -(1-kappa) times collapsed_da2 converted to spherical
   harmonics. */
  for (i=0;i<3;i++) for (j=0;j<3;j++)
    collapsed_da2[i][j] *= -(1-param->kappa);
  reverse_tensor2(collapsed_da2,psidot);
}

static int first_time=1;
static void check_for_ip_rights() {
  char answer[10];

  if (first_time) {
    first_time = 0;
    if (getenv("MAY_USE_PATENT_7266469")==NULL) {
      printf("Do you understand that this program uses intellectual property protected by U.S. Patent No. 7,266,469? ");
      fgets(answer,10,stdin);
      if (strcmp(answer,"yes\n")!=0) {
        printf("This program will not run if given an answer other than \"yes\" (written in full).\n");
        exit(1);
      }
    }
  }
}

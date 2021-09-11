#include "spherical.h"

extern double kappa;
static void check_for_ip_rights();

void compute_psidot_rsc(COMPLEX* psidot, COMPLEX* psi) {
  double a2[3][3];
  double da2[3][3];
  double collapsed_da2[3][3];
  double evec[3][3], eval[3];
  double inner;
  int i,j,k;

  check_for_ip_rights();

  tensor2(psi,a2);
  tensor2(psidot,da2);

/* collapsed_da2 = M:da2, where M = sum_{i=1}^3 e_i e_i e_i e_i, where
   e_i are the normalized eigenvectors of a2. */
  diagonalize_sym(3, a2, eval, evec);
  memset(collapsed_da2,0,sizeof(collapsed_da2));
  for (i=0;i<3;i++) {
    inner = 0;
    for (j=0;j<3;j++) for (k=0;k<3;k++)
      inner += evec[i][j]*evec[i][k]*da2[k][j];
    for (j=0;j<3;j++) for (k=0;k<3;k++)
      collapsed_da2[j][k] += evec[i][j]*evec[i][k]*inner;
  }

/* Add to psidot: -(1-kappa) times collapsed_da2 converted to spherical
   harmonics. */
  for (i=0;i<3;i++) for (j=0;j<3;j++)
    collapsed_da2[i][j] *= -(1-kappa);
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

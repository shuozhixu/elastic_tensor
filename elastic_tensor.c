/*
 *  elastic_tensor.c
 *
 *  Created by Shuozhi Xu (shuozhixu@ucsb.edu) on 05/22/2020.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#undef I

int main()
{
  int i, j, k, l, aflag;
  int p, q, m, n, u;
  double muv, mur, lambdav, lambdar, kv, kr;
  double sum_temp, det_temp, delta_temp;
  double ctemp0, ctemp1, ctemp2, ctemp3, ctemp4, ctemp5;
  double CC[3][3][3][3], SS[3][3][3][3];
  double CI[6][6], K[6][6], KT[6][6], KC[6][6], KCK[6][6], CCI[6][6], SSI[6][6];
  double SI[6][6], CCV[6][6], CCR[6][6], CCH[6][6], SSIV[6][6], SSIR[6][6];
  double ed1[3], ed2[3], ed3[3];
  double e1[3], e2[3], e3[3];
  double omega[3][3];
  double bb[6][6], fac[6][6];

  double determinant(double [6][6], int);

  #define DELTA(i, j) ((i==j)?1:0)

  FILE *efp, *ofp;

  efp = fopen("elas.in", "r");

  if (efp == NULL) {
    fprintf(stderr, "Cannot open the input file %s\n", "elas.in");
    exit(1);
  }

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      for (k=0; k<3; k++) {
        for (l=0; l<3; l++) {
          CC[i][j][k][l] = 0.;
          SS[i][j][k][l] = 0.;
        }
      }
    }
  }

  // e_dagger

  for (i=0; i<3; i++) {
    ed1[i] = 0.;
    ed2[i] = 0.;
    ed3[i] = 0.;
  }

  fscanf(efp, "%d", &aflag);

  if(aflag < 0 || aflag > 3) {
    printf("Anisotropy type %d specified has not been developed yet \n", aflag);
    exit(1);
  }

  fscanf(efp, "%lf %lf %lf", &ctemp0, &ctemp1, &ctemp2);

  ed1[0] = ctemp0 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));
  ed1[1] = ctemp1 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));
  ed1[2] = ctemp2 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));

  fscanf(efp, "%lf %lf %lf", &ctemp0, &ctemp1, &ctemp2);

  ed2[0] = ctemp0 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));
  ed2[1] = ctemp1 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));
  ed2[2] = ctemp2 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));

  fscanf(efp, "%lf %lf %lf", &ctemp0, &ctemp1, &ctemp2);

  ed3[0] = ctemp0 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));
  ed3[1] = ctemp1 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));
  ed3[2] = ctemp2 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));

  // CI

  for (i=0; i<6; i++) {
    for (j=0; j<6; j++) {
      CI[i][j] = 0.;
    }
  }

  fscanf(efp, "%lf", &ctemp0);
  CI[0][0] = ctemp0;

  fscanf(efp, "%lf %lf", &ctemp0, &ctemp1);
  CI[1][0] = ctemp0;
  CI[1][1] = ctemp1;

  fscanf(efp, "%lf %lf %lf", &ctemp0, &ctemp1, &ctemp2);
  CI[2][0] = ctemp0;
  CI[2][1] = ctemp1;
  CI[2][2] = ctemp2;

  fscanf(efp, "%lf %lf %lf %lf", &ctemp0, &ctemp1, &ctemp2, &ctemp3);
  CI[3][0] = ctemp0;
  CI[3][1] = ctemp1;
  CI[3][2] = ctemp2;
  CI[3][3] = ctemp3;

  fscanf(efp, "%lf %lf %lf %lf %lf", &ctemp0, &ctemp1, &ctemp2, &ctemp3, &ctemp4);
  CI[4][0] = ctemp0;
  CI[4][1] = ctemp1;
  CI[4][2] = ctemp2;
  CI[4][3] = ctemp3;
  CI[4][4] = ctemp4;

  fscanf(efp, "%lf %lf %lf %lf %lf %lf", &ctemp0, &ctemp1, &ctemp2, &ctemp3, &ctemp4, &ctemp5);
  CI[5][0] = ctemp0;
  CI[5][1] = ctemp1;
  CI[5][2] = ctemp2;
  CI[5][3] = ctemp3;
  CI[5][4] = ctemp4;
  CI[5][5] = ctemp5;

  for (i=0; i<6; i++) {
    for (j=0; j<6; j++) {
      if(fabs(CI[i][j]) > 1.E-5){
        CI[j][i] = CI[i][j];
      }
    }
  }

  k = 0;
    
  for (i=0; i<6; i++) {
    for (j=0; j<6; j++) {
      if(fabs(CI[i][j]) < 1.E-5){
        k = k+1;
      }
    }
  }

  if(k == 36) {
    fprintf(stderr, "The input elastic tensor cannot be zero\n");
    exit(1);
  }

  for (i=0; i<4; i++) {
    e1[i] = 0.;
    e2[i] = 0.;
    e3[i] = 0.;
  }

  fscanf(efp, "%lf %lf %lf", &ctemp0, &ctemp1, &ctemp2);

  e1[0] = ctemp0 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));
  e1[1] = ctemp1 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));
  e1[2] = ctemp2 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));
  e1[3] = 1.;

  fscanf(efp, "%lf %lf %lf", &ctemp0, &ctemp1, &ctemp2);

  e2[0] = ctemp0 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));
  e2[1] = ctemp1 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));
  e2[2] = ctemp2 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));
  e2[3] = 1.;

  fscanf(efp, "%lf %lf %lf", &ctemp0, &ctemp1, &ctemp2);

  e3[0] = ctemp0 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));
  e3[1] = ctemp1 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));
  e3[2] = ctemp2 / sqrt(pow(ctemp0, 2.) + pow(ctemp1, 2.) + pow(ctemp2, 2.));
  e3[3] = 1.;

  fclose(efp);

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      omega[i][j] = 0.;
    }
  }

  // omega

  omega[0][0] = ed1[0]*e1[0] + ed1[1]*e1[1] + ed1[2]*e1[2];
  omega[0][1] = ed2[0]*e1[0] + ed2[1]*e1[1] + ed2[2]*e1[2];
  omega[0][2] = ed3[0]*e1[0] + ed3[1]*e1[1] + ed3[2]*e1[2];

  omega[1][0] = ed1[0]*e2[0] + ed1[1]*e2[1] + ed1[2]*e2[2];
  omega[1][1] = ed2[0]*e2[0] + ed2[1]*e2[1] + ed2[2]*e2[2];
  omega[1][2] = ed3[0]*e2[0] + ed3[1]*e2[1] + ed3[2]*e2[2];

  omega[2][0] = ed1[0]*e3[0] + ed1[1]*e3[1] + ed1[2]*e3[2];
  omega[2][1] = ed2[0]*e3[0] + ed2[1]*e3[1] + ed2[2]*e3[2];
  omega[2][2] = ed3[0]*e3[0] + ed3[1]*e3[1] + ed3[2]*e3[2];

  for (i=0; i<6; i++) {
    for (j=0; j<6; j++) {
      K[i][j] = 0.;
      KT[i][j] = 0.;
      KC[i][j] = 0.;
      KCK[i][j] = 0.;
    }
  }

  // K1

  K[0][0] = pow(omega[0][0], 2.);
  K[1][0] = pow(omega[1][0], 2.);
  K[2][0] = pow(omega[2][0], 2.);

  K[0][1] = pow(omega[0][1], 2.);
  K[1][1] = pow(omega[1][1], 2.);
  K[2][1] = pow(omega[2][1], 2.);

  K[0][2] = pow(omega[0][2], 2.);
  K[1][2] = pow(omega[1][2], 2.);
  K[2][2] = pow(omega[2][2], 2.);

  // K2

  K[0][3] = 2. * omega[0][1] * omega[0][2];
  K[1][3] = 2. * omega[1][1] * omega[1][2];
  K[2][3] = 2. * omega[2][1] * omega[2][2];

  K[0][4] = 2. * omega[0][2] * omega[0][0];
  K[1][4] = 2. * omega[1][2] * omega[1][0];
  K[2][4] = 2. * omega[2][2] * omega[2][0];

  K[0][5] = 2. * omega[0][0] * omega[0][1];
  K[1][5] = 2. * omega[1][0] * omega[1][1];
  K[2][5] = 2. * omega[2][0] * omega[2][1];

  // K3

  K[3][0] = omega[1][0] * omega[2][0];
  K[4][0] = omega[2][0] * omega[0][0];
  K[5][0] = omega[0][0] * omega[1][0];

  K[3][1] = omega[1][1] * omega[2][1];
  K[4][1] = omega[2][1] * omega[0][1];
  K[5][1] = omega[0][1] * omega[1][1];

  K[3][2] = omega[1][2] * omega[2][2];
  K[4][2] = omega[2][2] * omega[0][2];
  K[5][2] = omega[0][2] * omega[1][2];

  // K4

  K[3][3] = omega[1][1] * omega[2][2] + omega[1][2] * omega[2][1];
  K[4][3] = omega[2][1] * omega[0][2] + omega[2][2] * omega[0][1];
  K[5][3] = omega[0][1] * omega[1][2] + omega[0][2] * omega[1][1];

  K[3][4] = omega[1][2] * omega[2][0] + omega[1][0] * omega[2][2];
  K[4][4] = omega[2][2] * omega[0][0] + omega[2][0] * omega[0][2];
  K[5][4] = omega[0][2] * omega[1][0] + omega[0][0] * omega[1][2];

  K[3][5] = omega[1][0] * omega[2][1] + omega[1][1] * omega[2][0];
  K[4][5] = omega[2][0] * omega[0][1] + omega[2][1] * omega[0][0];
  K[5][5] = omega[0][0] * omega[1][1] + omega[0][1] * omega[1][0];

  // KT: transpose of K

  for (i=0; i<6; i++) {
    for (j=0; j<6; j++) {
      KT[j][i] = K[i][j];
    }
  }

  // KC = K * CI

  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {

      sum_temp = 0.;
      for (k = 0; k < 6; k++) {
        sum_temp = sum_temp + K[i][k] * CI[k][j];
      }
      KC[i][j] = sum_temp;

    }
  }

  // KCK = KC * KT

  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {

      sum_temp = 0.;
      for (k = 0; k < 6; k++) {
        sum_temp = sum_temp + KC[i][k] * KT[k][j];
      }
      KCK[i][j] = sum_temp;

    }
  }

  for (i=0; i<6; i++) { // check the symmetry of KCK
    for (j=0; j<6; j++) {
      if(fabs(KCK[j][i] - KCK[i][j]) > 1.){
        fprintf(stderr, "KCK does not pass the symmetry check %d %d \n", i, j);
        printf("KCK1 = %lf, KCK2 = %lf \n", KCK[j][i], KCK[i][j]);
        exit(1);
      }
    }
  }

  if (aflag == 0){
    for (i=0; i<6; i++) {
      for (j=0; j<6; j++) {
        CCI[i][j] = KCK[i][j];
      }
    }
  }

  for (i=0; i<6; i++) {
    for (j=0; j<6; j++) {
      CCV[i][j] = 0.;
      CCR[i][j] = 0.;
      CCH[i][j] = 0.;
    }
  }

  if(aflag == 1 || aflag == 3){ // Voigt or Hill isotropic

    muv = (1./15.) * (CI[0][0] + CI[1][1] + CI[2][2] - (CI[0][1] + CI[1][2] + CI[0][2]) + 3. * (CI[3][3] + CI[4][4] + CI[5][5]));
    kv = (1./9.) * (CI[0][0] + CI[1][1] + CI[2][2] + 2. * (CI[0][1] + CI[1][2] + CI[0][2]));
    lambdav = kv - (2./3.) * muv;

    CCV[0][0] = lambdav + 2. * muv;
    CCV[1][1] = lambdav + 2. * muv;
    CCV[2][2] = lambdav + 2. * muv;
    CCV[0][1] = lambdav;
    CCV[0][2] = lambdav;
    CCV[1][2] = lambdav;
    CCV[3][3] = muv;
    CCV[4][4] = muv;
    CCV[5][5] = muv;

    for (i=0; i<6; i++) {
      for (j=0; j<6; j++) {
        if(fabs(CCV[i][j]) > 1.E-5){
          CCV[j][i] = CCV[i][j];
        }
      }
    }

    for (i=0; i<6; i++) {
      for (j=0; j<6; j++) {
        CCI[i][j] = CCV[i][j];
      }
    }

  }

  if(aflag == 2 || aflag == 3){ // Reuss or Hill isotropic

    for (i=0; i<6; i++) {
      for (j = 0; j< 6; j++) {
        SI[i][j] = 0.;
        bb[i][j] = 0.;
        fac[i][j] = 0.;
      }
    }

    // calculate SI, which is the inverse of CI

    // the first step is to calculate the determinant of CI

    k = 6;

    det_temp = determinant(CI, k);

    // the second step is to calculate the cofactor of CI

    if (det_temp < 1.) {
      printf("CI is not invertible. \n");
      exit(1);
    }
    else {

      for (q = 0;q < k; q++) {
        for (p = 0;p < k; p++) {
          m = 0;
          n = 0;
          for (i = 0;i < k; i++) {
            for (j = 0;j < k; j++){
              if (i != q && j != p){
                bb[m][n] = CI[i][j];
                if (n < (k - 2)) {
                  n++;
                }
                else {
                  n = 0;
                  m++;
                }
              }
            }
          }
          fac[q][p] = pow(-1, q + p) * determinant(bb, k - 1);
        }
      }

      for (i = 0; i < k; i++){
        for (j = 0; j < k; j++){
          SI[i][j] = fac[j][i] / det_temp;
        }
      }
    }

    mur = 15. / (4. * (SI[0][0] + SI[1][1] + SI[2][2]) - 4. * (SI[0][1] + SI[1][2] + SI[0][2]) + 3. * (SI[3][3] + SI[4][4] + SI[5][5]));
    kr = 1. / (SI[0][0] + SI[1][1] + SI[2][2] + 2. * (SI[0][1] + SI[1][2] + SI[0][2]));
    lambdar = kr - (2./3.) * mur;

    CCR[0][0] = lambdar + 2. * mur;
    CCR[1][1] = lambdar + 2. * mur;
    CCR[2][2] = lambdar + 2. * mur;
    CCR[0][1] = lambdar;
    CCR[0][2] = lambdar;
    CCR[1][2] = lambdar;
    CCR[3][3] = mur;
    CCR[4][4] = mur;
    CCR[5][5] = mur;

    for (i=0; i<6; i++) {
      for (j=0; j<6; j++) {
        if(fabs(CCR[i][j]) > 1.E-5){
          CCR[j][i] = CCR[i][j];
        }
      }
    }

    for (i=0; i<6; i++) {
      for (j=0; j<6; j++) {
        CCI[i][j] = CCR[i][j];
      }
    }

  }

  if(aflag == 3) {

    for (i=0; i<6; i++) {
      for (j=0; j<6; j++) {
          CCH[i][j] = 0.5 * (CCV[i][j] + CCR[i][j]);
          CCI[i][j] = CCH[i][j];
      }
    }

  }

  // convert KCK, CCV, CCR, or CCH to CC

  ofp = fopen("elas.out", "w");

  if(aflag == 0){ // anisotropic

    fprintf(ofp, "Anisotropic \n");
    fprintf(ofp, "In units of GPa, the stiffness tensor is\n");
    for (i=0; i<6; i++) {
      fprintf(ofp, "%lf %lf %lf %lf %lf %lf \n", KCK[i][0], KCK[i][1], KCK[i][2], KCK[i][3], KCK[i][4], KCK[i][5]);
    }

    for (i = 0; i < 3; i++){
      CC[0][0][i][i] = KCK[0][i];
      CC[1][1][i][i] = KCK[1][i];
      CC[2][2][i][i] = KCK[2][i];
      CC[1][2][i][i] = KCK[3][i];
      CC[2][0][i][i] = KCK[4][i];
      CC[0][1][i][i] = KCK[5][i];
    }

    CC[0][0][1][2] = KCK[0][3];
    CC[1][1][1][2] = KCK[1][3];
    CC[2][2][1][2] = KCK[2][3];
    CC[1][2][1][2] = KCK[3][3];
    CC[2][0][1][2] = KCK[4][3];
    CC[0][1][1][2] = KCK[5][3];

    CC[0][0][2][0] = KCK[0][4];
    CC[1][1][2][0] = KCK[1][4];
    CC[2][2][2][0] = KCK[2][4];
    CC[1][2][2][0] = KCK[3][4];
    CC[2][0][2][0] = KCK[4][4];
    CC[0][1][2][0] = KCK[5][4];

    CC[0][0][0][1] = KCK[0][5];
    CC[1][1][0][1] = KCK[1][5];
    CC[2][2][0][1] = KCK[2][5];
    CC[1][2][0][1] = KCK[3][5];
    CC[2][0][0][1] = KCK[4][5];
    CC[0][1][0][1] = KCK[5][5];

  }
  else if(aflag == 1) {

    fprintf(ofp, "Isotropic in Voigt form \n");
    fprintf(ofp, "In units of GPa, the stiffness tensor is \n");
      for (i=0; i<6; i++) {
      fprintf(ofp, "%lf %lf %lf %lf %lf %lf \n", CCV[i][0], CCV[i][1], CCV[i][2], CCV[i][3], CCV[i][4], CCV[i][5]);
    }

    CC[0][0][0][0] = CCV[0][0];
    CC[1][1][1][1] = CCV[1][1];
    CC[2][2][2][2] = CCV[2][2];
    CC[0][0][1][1] = CCV[0][1];
    CC[0][0][2][2] = CCV[0][2];
    CC[1][1][2][2] = CCV[1][2];
    CC[0][1][0][1] = CCV[3][3];
    CC[1][2][1][2] = CCV[4][4];
    CC[0][2][0][2] = CCV[5][5];

  }
  else if(aflag == 2) {

    fprintf(ofp, "Isotropic in Reuss form \n");
    fprintf(ofp, "In units of GPa, the stiffness tensor is \n");
    for (i=0; i<6; i++) {
      fprintf(ofp, "%lf %lf %lf %lf %lf %lf \n", CCR[i][0], CCR[i][1], CCR[i][2], CCR[i][3], CCR[i][4], CCR[i][5]);
    }

    CC[0][0][0][0] = CCR[0][0];
    CC[1][1][1][1] = CCR[1][1];
    CC[2][2][2][2] = CCR[2][2];
    CC[0][0][1][1] = CCR[0][1];
    CC[0][0][2][2] = CCR[0][2];
    CC[1][1][2][2] = CCR[1][2];
    CC[0][1][0][1] = CCR[3][3];
    CC[1][2][1][2] = CCR[4][4];
    CC[0][2][0][2] = CCR[5][5];

  }
  else if(aflag == 3){

    fprintf(ofp, "Isotropic in Hill form \n");
    fprintf(ofp, "In units of GPa, the stiffness tensor is \n");
    for (i=0; i<6; i++) {
      fprintf(ofp, "%lf %lf %lf %lf %lf %lf \n", CCH[i][0], CCH[i][1], CCH[i][2], CCH[i][3], CCH[i][4], CCH[i][5]);
    }

    CC[0][0][0][0] = CCH[0][0];
    CC[1][1][1][1] = CCH[1][1];
    CC[2][2][2][2] = CCH[2][2];
    CC[0][0][1][1] = CCH[0][1];
    CC[0][0][2][2] = CCH[0][2];
    CC[1][1][2][2] = CCH[1][2];
    CC[0][1][0][1] = CCH[3][3];
    CC[1][2][1][2] = CCH[4][4];
    CC[0][2][0][2] = CCH[5][5];

  }

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      for (k=0; k<3; k++) {
        for (l=0; l<3; l++) {
          if(aflag == 0) { // check the symmetry of CC
            if(fabs(CC[k][l][i][j] - CC[i][j][k][l]) > 1.E-5){
              fprintf(stderr, "CC does not pass the first symmetry check %d %d %d %d \n", i, j, k, l);
              printf("CC1 = %lf, CC2 = %lf \n", CC[k][l][i][j], CC[i][j][k][l]);
              exit(1);
            }
          }
          else {
            if(fabs(CC[i][j][k][l]) > 1.E-5){
              CC[k][l][i][j] = CC[i][j][k][l];
            }
          }
        }
      }
    }
  }

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      for (k=0; k<3; k++) {
        for (l=0; l<3; l++) {
          if(fabs(CC[i][j][k][l]) > 1.E-5){
            CC[i][j][l][k] = CC[i][j][k][l];
            CC[j][i][k][l] = CC[i][j][k][l];
            CC[j][i][l][k] = CC[i][j][k][l];
          }
        }
      }
    }
  }

  // change units from GPa to Pa

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      for (k=0; k<3; k++) {
        for (l=0; l<3; l++) {
          CC[i][j][k][l] = CC[i][j][k][l] * 1.E9;
        }
      }
    }
  }

  for (i=0; i<6; i++) {
    for (j = 0; j< 6; j++) {
      SSI[i][j] = 0.;
    }
  }

  // calculate the inverse of KCK, CCV, or CCR

  if(aflag == 0) {

    for (i=0; i<6; i++) {
      for (j = 0; j< 6; j++) {
        bb[i][j] = 0.;
        fac[i][j] = 0.;
      }
    }

    // calculate SSI, which is the inverse of KCK

    // the first step is to calculate the determinant of KCK

    k = 6;

    det_temp = determinant(KCK, k);

    // the second step is to calculate the cofactor of KCK

    if (det_temp < 1.) {
      printf("KCK is not invertible. \n");
      exit(1);
    }
    else {

      for (q = 0;q < k; q++) {
        for (p = 0;p < k; p++) {
          m = 0;
          n = 0;
          for (i = 0;i < k; i++) {
            for (j = 0;j < k; j++){
              if (i != q && j != p){
                bb[m][n] = KCK[i][j];
                if (n < (k - 2)) {
                  n++;
                }
                else {
                  n = 0;
                  m++;
                }
              }
            }
          }
          fac[q][p] = pow(-1, q + p) * determinant(bb, k - 1);
        }
      }

      for (i = 0; i < k; i++){
        for (j = 0; j < k; j++){
          SSI[i][j] = fac[j][i] / det_temp;
        }
      }
    }

  }

  if(aflag == 1 || aflag == 3) {

    for (i=0; i<6; i++) {
      for (j = 0; j< 6; j++) {
        SSIV[i][j] = 0.;
        bb[i][j] = 0.;
        fac[i][j] = 0.;
      }
    }

    // calculate SSI, which is the inverse of CCV

    // the first step is to calculate the determinant of CCV

    k = 6;

    det_temp = determinant(CCV, k);

    // the second step is to calculate the cofactor of CCV

    if (det_temp < 1.) {
      printf("CCV is not invertible. \n");
      exit(1);
    }
    else {

      for (q = 0;q < k; q++) {
        for (p = 0;p < k; p++) {
          m = 0;
          n = 0;
          for (i = 0;i < k; i++) {
            for (j = 0;j < k; j++){
              if (i != q && j != p){
                bb[m][n] = CCV[i][j];
                if (n < (k - 2)) {
                  n++;
                }
                else {
                  n = 0;
                  m++;
                }
              }
            }
          }
          fac[q][p] = pow(-1, q + p) * determinant(bb, k - 1);
        }
      }

      for (i = 0; i < k; i++){
        for (j = 0; j < k; j++){
          SSI[i][j] = fac[j][i] / det_temp;
          SSIV[i][j] = SSI[i][j];
        }
      }
    }

  }

  if(aflag == 2 || aflag == 3) {

    for (i=0; i<6; i++) {
      for (j = 0; j< 6; j++) {
        SSIR[i][j] = 0.;
        bb[i][j] = 0.;
        fac[i][j] = 0.;
      }
    }

    // calculate SSI, which is the inverse of CCR

    // the first step is to calculate the determinant of CCR

    k = 6;

    det_temp = determinant(CCR, k);

    // the second step is to calculate the cofactor of CCR

    if (det_temp < 1.) {
      printf("CCR is not invertible. \n");
      exit(1);
    }
    else {

      for (q = 0;q < k; q++) {
        for (p = 0;p < k; p++) {
          m = 0;
          n = 0;
          for (i = 0;i < k; i++) {
            for (j = 0;j < k; j++){
              if (i != q && j != p){
                bb[m][n] = CCR[i][j];
                if (n < (k - 2)) {
                  n++;
                }
                else {
                  n = 0;
                  m++;
                }
              }
            }
          }
          fac[q][p] = pow(-1, q + p) * determinant(bb, k - 1);
        }
      }

      for (i = 0; i < k; i++){
        for (j = 0; j < k; j++){
          SSI[i][j] = fac[j][i] / det_temp;
          SSIR[i][j] = SSI[i][j];
        }
      }
    }

  }

  if(aflag == 3) {

    for (i=0; i<6; i++) {
      for (j = 0; j< 6; j++) {
        SSI[i][j] = 0.5 * (SSIV[i][j] + SSIR[i][j]);
      }
    }

  }

  // print the compliance tensor SSI

  fprintf(ofp, "In units of 1/GPa, the compliance tensor is \n");
  for (i=0; i<6; i++) {
    fprintf(ofp, "%lf %lf %lf %lf %lf %lf \n", SSI[i][0], SSI[i][1], SSI[i][2], SSI[i][3], SSI[i][4], SSI[i][5]);
  }

  fclose(ofp);

  // check the symmetry of SSI and (when aflag is not 3) whether SSI is the inverse of CCI

  for (i=0; i<6; i++) { // check the symmetry of SSI
    for (j=0; j<6; j++) {

      if(fabs(SSI[j][i] - SSI[i][j]) > 1./det_temp){
        fprintf(stderr, "SSI does not pass the symmetry check %d %d \n", i, j);
        printf("SSI1 = %lf, SSI2 = %lf \n", SSI[j][i], SSI[i][j]);
        exit(1);
      }

      if(aflag == 0 || aflag == 1 || aflag == 2) {

        sum_temp = 0.;
        for (k = 0; k < 6; k++) {
          sum_temp = sum_temp + CCI[i][k] * SSI[k][j];
        }

        if(fabs(sum_temp - DELTA(i,j)) > 1.E-5){
          fprintf(stderr, "SSI is not the inverse of CCI %d %d \n", i, j);
          printf("SSI = %lf, CCI = %lf \n", SSI[i][j], CCI[i][j]);
          exit(1);
        }
      }

    }
  }

  for (i = 0; i < 3; i++){
    SS[0][0][i][i] = SSI[0][i];
    SS[1][1][i][i] = SSI[1][i];
    SS[2][2][i][i] = SSI[2][i];
    SS[1][2][i][i] = SSI[3][i] / 2.;
    SS[2][0][i][i] = SSI[4][i] / 2.;
    SS[0][1][i][i] = SSI[5][i] / 2.;
  }

  SS[0][0][1][2] = SSI[0][3] / 2.;
  SS[1][1][1][2] = SSI[1][3] / 2.;
  SS[2][2][1][2] = SSI[2][3] / 2.;
  SS[1][2][1][2] = SSI[3][3] / 4.;
  SS[2][0][1][2] = SSI[4][3] / 4.;
  SS[0][1][1][2] = SSI[5][3] / 4.;

  SS[0][0][2][0] = SSI[0][4] / 2.;
  SS[1][1][2][0] = SSI[1][4] / 2.;
  SS[2][2][2][0] = SSI[2][4] / 2.;
  SS[1][2][2][0] = SSI[3][4] / 4.;
  SS[2][0][2][0] = SSI[4][4] / 4.;
  SS[0][1][2][0] = SSI[5][4] / 4.;

  SS[0][0][0][1] = SSI[0][5] / 2.;
  SS[1][1][0][1] = SSI[1][5] / 2.;
  SS[2][2][0][1] = SSI[2][5] / 2.;
  SS[1][2][0][1] = SSI[3][5] / 4.;
  SS[2][0][0][1] = SSI[4][5] / 4.;
  SS[0][1][0][1] = SSI[5][5] / 4.;

  for (i=0; i<3; i++) { // check the symmetry of SS
    for (j=0; j<3; j++) {
      for (k=0; k<3; k++) {
        for (l=0; l<3; l++) {
          if(fabs(SS[k][l][i][j] - SS[i][j][k][l]) > 1./det_temp){
            fprintf(stderr, "SS does not pass the first symmetry check %d %d %d %d \n", i, j, k, l);
            printf("SS1 = %lf, SS2 = %lf \n", SS[k][l][i][j], SS[i][j][k][l]);
            exit(1);
          }
        }
      }
    }
  }

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      for (k=0; k<3; k++) {
        for (l=0; l<3; l++) {
          if(fabs(SS[i][j][k][l]) > 1./det_temp){
            SS[i][j][l][k] = SS[i][j][k][l];
            SS[j][i][k][l] = SS[i][j][k][l];
            SS[j][i][l][k] = SS[i][j][k][l];
          }
        }
      }
    }
  }

  // change units from 1/GPa to 1/Pa

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      for (k=0; k<3; k++) {
        for (l=0; l<3; l++) {
          SS[i][j][k][l] = SS[i][j][k][l] / 1.E9;
        }
      }
    }
  }

  // when aflag is not 3, check whether SS is the inverse of CC

  if(aflag == 0 || aflag == 1 || aflag == 2) {

    for (m=0; m<3; m++) {
      for (n=0; n<3; n++) {
        for (k=0; k<3; k++) {
          for (l=0; l<3; l++) {

            sum_temp = 0.;

            for (i=0; i<3; i++) {
              for (j=0; j<3; j++) {

                sum_temp = sum_temp + CC[m][n][i][j] * SS[i][j][k][l];

              }
            }

            delta_temp = (1./2.) * (DELTA(m, k) * DELTA(n, l) + DELTA(m, l) * DELTA(n, k));

            if(fabs(sum_temp - delta_temp) > 1.E-5){
              fprintf(stderr, "SS is not the inverse of CC %d %d %d %d \n", m, n, k, l);
              printf("SS*CC = %lf, DELTA = %lf \n", sum_temp, delta_temp);
              exit(1);
            }
          }
        }
      }
    }
  }

  return 0;
}

// The following function is from https://www.sanfoundry.com/c-program-find-inverse-matrix/

double determinant(double a[6][6], int k)
{
  double s, det;
  double b[6][6];
  int i, j, m, n, c;

  s = 1.;
  det = 0.;

  if (k == 1) {
    return (a[0][0]);
  }
  else {
    det = 0;
    for (c = 0; c < k; c++) {
      m = 0;
      n = 0;
      for (i = 0;i < k; i++) {
        for (j = 0 ;j < k; j++) {
          b[i][j] = 0.;
          if (i != 0 && j != c){
            b[m][n] = a[i][j];
            if (n < (k - 2)) {
              n++;
            }
            else {
              n = 0;
              m++;
            }
          }
        }
      }
      det = det + s * (a[0][c] * determinant(b, k - 1));
      s = -1. * s;
    }
  }

  return (det);
}

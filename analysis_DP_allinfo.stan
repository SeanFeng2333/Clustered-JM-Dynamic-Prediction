functions{
// ------------------------------------------------------
//     LINEAR PREDICTOR FOR THE LONGITUDINAL SUBMODEL                
// ------------------------------------------------------ 
//  vector linear_predictor(vector X, vector visits, int<lower=1, upper=n> ID, vector beta, vector bi){
//         int N = num_elements(visits);
//         vector[N] out;

//         out = beta[1] + X[ID]*beta[2]  + beta[3]*visits + bi[ID];

//         return out;
      real linear_predictor1(real X1, real X2, real X3, real X4, real X5, real visits, int IDtooth, vector beta1, vector bi, matrix bij){
        real out;
        out = beta1[1] + X1*beta1[2] + X2*beta1[3] + X3*beta1[4] + X4*beta1[5] + X5*beta1[6] + beta1[7]*visits + bi[1] + bij[IDtooth,1] + bij[IDtooth,2]*visits;
        
        return out;
    } 
    
      real linear_predictor2(real X1, real X2, real X3, real X4, real X5, real visits, int IDtooth, vector beta2, vector bi, matrix bij){
        real out;
        out = beta2[1] + X1*beta2[2] + X2*beta2[3] + X3*beta2[4] + X4*beta2[5] + X5*beta2[6] + beta2[7]*visits + bi[2] + bij[IDtooth,3] + bij[IDtooth,4]*visits;
        return out;
    } 
// ------------------------------------------------------ 
}


data{
  int N; // total number of longitudinal outcomes 
  int nteeth; // total number of teeth
  vector[N] y1; // longitudinal outcomes
  array[N] int<lower=0, upper=1> y2;  // longitudinal outcomes
  vector[nteeth] X1; // covariate X1 
  vector[nteeth] X2; // covariate X2
  vector[nteeth] X3; // covariate X3
  vector[nteeth] X4; // covariate X4
  vector[nteeth] X5; // covariate X5

  vector[N] visits; // visit times for repeated observations
  vector[nteeth] times; // times to event
  vector[nteeth] status; // status indicator
//  int<lower=1,upper=n> ID[N];
  array[N] int IDtooth;
  array[nteeth] int IDtooth_s;
  int length_GK;
  vector[length_GK] sk;
  vector[length_GK] wk;
  
  vector[7] beta1;
  vector[7] beta2;
  vector[6] gamma;
  real<lower=0> xi;
  
  real alpha1;
  real alpha2;

  real<lower=0> sigma_e;  
  real<lower=0> sigma_ri;


  vector<lower=0>[2] sigma_bi;        // Standard deviations for bi
  vector<lower=0>[4] sigma_bij;      // Standard deviations for bijk

  matrix[2,2] L_bi;       // Cholesky factor of correlation matrix for bi
  matrix[4,4] L_bij;     // Cholesky factor of correlation matrix for bijk
  
 
}

parameters{
//  vector[2] z_bi_mat;
//  matrix[4,nteeth] z_bij_mat;
  vector[2] bi;
  matrix[nteeth,4] bij;
  real ri;
}

//transformed parameters{
//  vector[2] bi;
//  matrix[nteeth,4] bij;
//  bi = (diag_pre_multiply(sigma_bi, L_bi) * z_bi_mat);
//  bij = (diag_pre_multiply(sigma_bij, L_bij) * z_bij_mat)';
//}



model{
// ------------------------------------------------------
//        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL                
// ------------------------------------------------------
{
   vector[N] linpred1; 
   vector[N] linpred2; 
   vector[N] p2; 
   for (i in 1:N) {
     int IDtoothi = IDtooth[i];
     real X1_i = X1[IDtoothi];
     real X2_i = X2[IDtoothi];
     real X3_i = X3[IDtoothi];
     real X4_i = X4[IDtoothi];
     real X5_i = X5[IDtoothi];
     real visit_i = visits[i];
     int y2_i = y2[i];
     linpred1[i] = linear_predictor1(X1_i, X2_i, X3_i, X4_i, X5_i, visit_i, IDtoothi, beta1, bi, bij);
     linpred2[i] = linear_predictor2(X1_i, X2_i, X3_i, X4_i, X5_i, visit_i, IDtoothi, beta2, bi, bij);
     p2[i] = inv_logit(linpred2[i]);
     target += bernoulli_lpmf(y2_i | p2[i]);
   }
   
   // Longitudinal Normal log-likelihood
   target += normal_lpdf(y1 | linpred1, sigma_e);
}  
// ------------------------------------------------------
//          LOG-LIKELIHOOD FOR SURVIVAL SUBMODEL                
// ------------------------------------------------------
{
   vector[nteeth] linpred1_at_event; 
   vector[nteeth] linpred2_at_event; 
   vector[nteeth] loghazard;
   vector[nteeth] logsurv;
   for (i in 1:nteeth){
      int IDtoothi = IDtooth_s[i];
      real X1_i = X1[IDtoothi];
      real X2_i = X2[IDtoothi];
      real X3_i = X3[IDtoothi];
      real X4_i = X4[IDtoothi];
      real X5_i = X5[IDtoothi];

      real time_i = times[i];
      real status_i = status[i];
      real linpred1_at_event_i = beta1[1] + beta1[2]*X1_i + beta1[3]*X2_i + beta1[4]*X3_i + beta1[5]*X4_i + beta1[6]*X5_i + beta1[7]*time_i + bi[1] + bij[IDtoothi,1] + bij[IDtoothi,2]*time_i;
      real linpred2_at_event_i = beta2[1] + beta2[2]*X1_i + beta2[3]*X2_i + beta2[4]*X3_i + beta2[5]*X4_i + beta2[6]*X5_i + beta2[7]*time_i + bi[2] + bij[IDtoothi,3] + bij[IDtoothi,4]*time_i;

      linpred1_at_event[i] = linpred1_at_event_i;
      linpred2_at_event[i] = linpred2_at_event_i;
      loghazard[i] = status_i * (log(xi) + (xi-1)*log(time_i) +  
       ri + gamma[1] + gamma[2]*X1_i + gamma[3]*X2_i + gamma[4]*X3_i + gamma[5]*X4_i + gamma[6]*X5_i + alpha1 * linpred1_at_event_i + alpha2 * linpred2_at_event_i);
      
      vector[length_GK] logsurv1;
      for (j in 1:length_GK){
        real sk_scaled_j = (sk[j]+1)*time_i/2;
        real wk_scaled_j = wk[j]*time_i/2;
        logsurv1[j] = wk_scaled_j*xi*sk_scaled_j^(xi-1)*exp(alpha1*(beta1[7]+bij[IDtoothi,2])*sk_scaled_j + 
          alpha2*(beta2[7]+bij[IDtoothi,4])*sk_scaled_j);
      }
      real approx_logsurv1 = sum(logsurv1);
      
      real logsurv2 = exp(ri + gamma[1] + gamma[2]*X1_i + gamma[3]*X2_i + 
        alpha1*(beta1[1]+beta1[2]*X1_i+beta1[3]*X2_i+beta1[4]*X3_i+beta1[5]*X4_i+ beta1[6]*X5_i+ bi[1] + bij[IDtoothi,1]) + 
        alpha2*(beta2[1]+beta2[2]*X1_i+beta2[3]*X2_i+beta2[4]*X3_i+beta2[5]*X4_i+ beta2[6]*X5_i+ bi[2] + bij[IDtoothi,3]));
        
      logsurv[i] = approx_logsurv1 * logsurv2;
   }
//   linpred_at_event = beta[1] + beta[2]*X1 + beta[3]*times + bi;
   
//   hazard = exp(gamma[1] + gamma[2]*X1 + alpha * linpred_at_event);
   // Log-hazard function
//   loghazard = rows_dot_product(status, (gamma[1] + gamma[2]*X1 + alpha * linpred_at_event));
   // Log-survival function
//   logsurv = -rows_dot_product(times, hazard);
    
   // Survival log-likelihood
   target += sum(loghazard) - sum(logsurv); 
} 
// ------------------------------------------------------
//                       LOG-PRIORS                       
// ------------------------------------------------------
   // Random-effects
   target += normal_lpdf(ri | 0, sigma_ri);
   target += multi_normal_cholesky_lpdf(bi | rep_vector(0,2), L_bi);
   for (t in 1:nteeth)
   target += multi_normal_cholesky_lpdf(bij[t]' | rep_vector(0, 4), L_bij);

//   target += normal_lpdf(to_vector(z_bi_mat) | 0, 1); 
//   target += normal_lpdf(to_vector(z_bij_mat) | 0, 1);    
}


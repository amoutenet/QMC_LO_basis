#pragma once

#include <cmath>
#include <triqs/mc_tools/random_generator.hpp>
#include "lo_vertex.hpp"
#include "configuration.hpp"

struct change_order_haar_scattering_alpha_sym {
 configuration* config;
 triqs::mc_tools::random_generator* RNG;

 int order_max;
 int order_min;

 double t_max;
 double lambda;
 double c_t;
 double c_a;
 dcomplex diag;

 std::complex<double> i{0,1};
 std::array<std::complex<double>, 4> coeffs;
 double norm_cauchy;

 bool moving_up;

 std::complex<double> det;
 std::vector<lo_vertex> V;

 std::vector<triqs::arrays::matrix<int>> perm_vecs;
 //std::vector<std::vector<int>> min_perm_vec;

 change_order_haar_scattering_alpha_sym(configuration *config_, triqs::mc_tools::random_generator* RNG_, int order_max_, int order_min_, double t_max_, double lambda_, double cauchy_t, double cauchy_a, dcomplex diag_, std::vector<triqs::arrays::matrix<int>> p_):
   config(config_), RNG(RNG_), order_max(order_max_), order_min(order_min_), t_max(t_max_), lambda(lambda_), c_t(cauchy_t), c_a(cauchy_a), diag(diag_), coeffs{i,{-1,0},-i,{1,0}}, perm_vecs(p_) {
     
  // normalisation factor for the Cauchy law
  norm_cauchy = std::atan((t_max_ - cauchy_t)/cauchy_a) - std::atan(- cauchy_t/cauchy_a);
 }


// ----------------------------------------------------------------
// Attempt() function
// ----------------------------------------------------------------

 std::complex<double> attempt(){
  //std::cout << "Attempting move change order" << std::endl;
  //std::cout << "We are at order " << config->order << std::endl;

  double proba_mv_up;
  if (config->order == order_max) proba_mv_up = 0;
  else if (config->order > order_min) proba_mv_up = (*RNG)(1.);
  else proba_mv_up = 1;

  moving_up = (proba_mv_up > 0.5);

  if (moving_up) {
   V.resize(config->order + 1);

   for (int i = 0; i < V.size() - 1; i++) V[i] = lo_vertex{config->vertices[i].t, (*RNG)(2), (*RNG)(2)};

   V[V.size() - 1] = get_random_lo_vertex(t_max, RNG, c_t, c_a);
   auto ratio = (V[V.size() - 1].t - c_t) / c_a;
   auto proba = c_a * norm_cauchy * (1 + ratio* ratio);

   std::sort(V.begin(), V.end());

   // quick exit
   if (V[V.size() - 1].itau == 0) return 0;

   if (config->order < 3){
     bool b = (config->order > 0);
     for (int i = 0; (i < V.size() - 1) and b; i++) b = (V[i].itau == (config->order %2));
     if (b) return 0;
   }

   // Compute first det
   det = 1;
   for (int sigma = 0; sigma <2; sigma++) det *= config->compute_det(sigma, V); 

   int pos_v = 0;
   for (int i = 0; i < V.size(); i++) {
     if (V[i].itau) pos_v += (1 << (2*V.size() - 1 - 2*i));
     if (V[i].l) pos_v += (1 << (2*V.size() - 2 - 2*i));
   }
   
   std::vector<lo_vertex> W(V.size());

   int j = 0;
   while (j < first_dim(perm_vecs[config->order + 1])) {
     int pos_w = perm_vecs[config->order + 1](j, pos_v);
     if (pos_w == -1) break;
     else {
       for (int i = 0; i < W.size(); i++) {
         W[i].t = V[i].t;
         W[i].itau = (pos_w >> (2*W.size() - 1 - 2*i)) & 1;
         W[i].l = (pos_w >> (2*W.size() - 2 - 2*i)) & 1;
       }
     }

     // Compute new det
     std::complex<double> det_temp = 1;
     for (int sigma = 0; sigma <2; sigma++) det_temp *= config->compute_det(sigma, W); 

     det += det_temp;
     j += 1;
   }
   det /= (j + 1);

   // get accurate 0.5 factors in the detailed balance
   double c_to_cp = (proba_mv_up < 1) ? 0.5 : 1;
   double cp_to_c = ((config->order + 1) < order_max) ? 0.5 : 1;

   auto weight = real(coeffs[(config->order + 1) % 4] * det) / real(coeffs[config->order % 4] * config->det) * 
     2 * lambda * proba * cp_to_c / ((config->order + 1) * c_to_cp);

   //std::cout << "MC weight is " << weight << std::endl;
   return isnan(weight) ? 0 : weight;
   //return det_ratio * 2 * i * t_max * lambda / (config->order + 1); //imaginary weight
  }

  else {
   V.resize(config->order - 1);

   int p = (*RNG)(config->order);
   auto ratio = (config->vertices[p].t - c_t) / c_a;
   auto proba = c_a * norm_cauchy * (1 + ratio* ratio);

   for (int i = 0; i < p; i++) V[i] = lo_vertex{config->vertices[i].t, (*RNG)(2), (*RNG)(2)};
   for (int i = p; i < V.size(); i++) V[i] = lo_vertex{config->vertices[i + 1].t, (*RNG)(2), (*RNG)(2)};

   std::sort(V.begin(), V.end());

   // quick exit
   if ((config->order > 1) and (V[V.size() - 1].itau == 0)) return 0;

   if (config->order == 3 ) {
     bool b = (config->order > 2);
     for (int i = 0; (i < V.size() - 1) and b; i++) b = (V[i].itau == (config->order %2));
     if (b) return 0;
   }
   
   if (config->order == 1) det = 2*diag;
   else {
    det = 1;
    for (int sigma = 0; sigma <2; sigma++) det *= config->compute_det(sigma, V); 

    int pos_v = 0;
    for (int i = 0; i < V.size(); i++) {
      if (V[i].itau) pos_v += (1 << (2*V.size() - 1 - 2*i));
      if (V[i].l) pos_v += (1 << (2*V.size() - 2 - 2*i));
    }

    std::vector<lo_vertex> W(V.size());

    int j = 0;
    while (j < first_dim(perm_vecs[config->order - 1])) {
      int pos_w = perm_vecs[config->order - 1](j, pos_v);
      if (pos_w == -1) break;
      else {
        for (int i = 0; i < W.size(); i++) {
          W[i].t = V[i].t;
          W[i].itau = (pos_w >> (2*W.size() - 1 - 2*i)) & 1;
          W[i].l = (pos_w >> (2*W.size() - 2 - 2*i)) & 1;
        }
      }

      // Compute new det
      std::complex<double> det_temp = 1;
      for (int sigma = 0; sigma <2; sigma++) det_temp *= config->compute_det(sigma, W); 

      det += det_temp;
      j += 1;
    }
    det /= (j + 1);

   }

   // get accurate 0.5 factors in the detailed balance
   double c_to_cp = (proba_mv_up > 0) ? 0.5 : 1;
   double cp_to_c = ((config->order - 1) > order_min) ? 0.5 : 1;

   //std::cout << "sum of new dets is :" << det_new + det_group_new << std::endl;

   auto weight =  real(coeffs[(config->order - 1)%4] * det) / real(coeffs[config->order %4] * config->det) * 
     config->order * cp_to_c / (2 * lambda * proba * c_to_cp);
   //std::cout << "MC weights is " << weight << std::endl;

   return isnan(weight) ? 0 : weight;
   //return det_ratio * config->order / (2 * i * t_max * lambda);
  }

}

// ----------------------------------------------------------------
// Accept() function
// ----------------------------------------------------------------

std::complex<double> accept(){
  //std::cout << "Move change order accepted" << std::endl;
  moving_up ? config->order++ : config->order--;

  config->det = det;
  config->vertices = V;

  return 1.0;
 }

// ----------------------------------------------------------------
// Reject() function
// ----------------------------------------------------------------

void reject() {}

};

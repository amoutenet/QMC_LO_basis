#pragma once

#include <cmath>
#include <triqs/mc_tools/random_generator.hpp>
#include "lo_vertex.hpp"
#include "configuration.hpp"

struct change_order {
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

 std::complex<double> new_det;
 std::vector<lo_vertex> V;

 change_order(configuration *config_, triqs::mc_tools::random_generator* RNG_, int order_max_, int order_min_, double t_max_, double lambda_, double cauchy_t, double cauchy_a, dcomplex diag_):
   config(config_), RNG(RNG_), order_max(order_max_), order_min(order_min_), t_max(t_max_), lambda(lambda_), c_t(cauchy_t), c_a(cauchy_a), diag(diag_), coeffs{i,{-1,0},-i,{1,0}} {

  // normalisation factor for the Cauchy law
  norm_cauchy = std::atan((t_max_ - cauchy_t)/cauchy_a) - std::atan(- cauchy_t/cauchy_a);
 }

 std::complex<double> attempt(){

  //std::cout << "Order is " << config->order << std::endl;
  //std::cout << "Current det is " << config->det << " and current proba is " << config->proba << std::endl;

  new_det = 1;

  double proba_mv_up;
  if (config->order == order_max) proba_mv_up = 0;
  else if (config->order > order_min) proba_mv_up = (*RNG)(1.);
  else proba_mv_up = 1;

  moving_up = (proba_mv_up > 0.5);

  if (moving_up) {
   //std::cout << "Attempting to add a vertex" << std::endl;
   V.resize(config->order + 1);

   for (int i =0; i < V.size() - 1; i++) {
     V[i] = lo_vertex{config->vertices[i].t, (*RNG)(2), (*RNG)(2)};
   }
   
   V[V.size() - 1] = get_random_lo_vertex(t_max, RNG, c_t, c_a);
   auto ratio = (V[V.size() - 1].t - c_t) / c_a;
   auto proba = c_a * norm_cauchy * (1 + ratio * ratio);

   std::sort(V.begin(), V.end());

   //for (int i = 0; i < V.size(); i++) std::cout << "time " << V[i].t << " - itau " << V[i].itau << " - l " << V[i].l <<std::endl;
   //for (int i = 0; i < config->order; i++) std::cout << "time " << config->vertices[i].t << " - itau " << config->vertices[i].itau << " - l " << config->vertices[i].l <<std::endl;

   // quick exit
   if (V[V.size() - 1].itau == 0) return 0;

   // Construt det up and det down
   for (int sigma = 0; sigma < 2; sigma++) new_det *= config->compute_det(sigma, V);

   //get accurate 0.5 factors in the detailed balance
   double c_to_cp = (proba_mv_up < 1) ? 0.5 : 1;
   double cp_to_c = ((config->order + 1) < order_max) ? 0.5 : 1;
   //TRIQS_PRINT(c_to_cp);
   //TRIQS_PRINT(cp_to_c);

   //TRIQS_PRINT(new_det);
   //TRIQS_PRINT(config->det);
   //TRIQS_PRINT(new_proba);
   //TRIQS_PRINT(config->proba);

   auto det_ratio = real(coeffs[(config->order + 1) % 4] * new_det) / real(coeffs[config->order % 4] * config->det); 
   //TRIQS_PRINT(det_ratio);
   auto proba_ratio = 2 * lambda * proba * cp_to_c / (c_to_cp * (config->order + 1));
   //TRIQS_PRINT(proba_ratio);

   //auto weight = real(coeffs[(config->order + 1) % 4] * new_det) / real(coeffs[config->order % 4] * config->det) *
   //  2 * lambda * new_proba * cp_to_c / (c_to_cp * config->proba);
   auto weight = det_ratio * proba_ratio;

   //TRIQS_PRINT(weight);

   return isnan(weight) ? 0 : weight;
   //return det_ratio * 2 * i * t_max * lambda / (config->order + 1); //imaginary weight
  }

  else{
   //std::cout << "Attempting to remove a vertex" << std::endl;
   V.resize(config->order - 1);

   int p = (*RNG)(config->order); 
   //TRIQS_PRINT(p);
   auto ratio = (config->vertices[p].t - c_t) / c_a;
   auto proba = c_a * norm_cauchy * (1 + ratio * ratio);

   for (int i =0; i < p; i++) V[i] = lo_vertex{config->vertices[i].t, (*RNG)(2), (*RNG)(2)};
   for (int i =p; i < V.size(); i++) V[i] = lo_vertex{config->vertices[i + 1].t, (*RNG)(2), (*RNG)(2)};

   std::sort(V.begin(), V.end());

   // quick exit
   if ((config->order > 1) and (V[V.size() - 1].itau == 0)) return 0;

   // Construt det up and det down
   if (config->order == 1) new_det = 2*diag; 
   else {
     for (int sigma = 0; sigma < 2; sigma++) new_det *= config->compute_det(sigma, V);
   }

   //for (int i = 0; i < V.size(); i++) std::cout << "time " << V[i].t << " - itau " << V[i].itau << " - l " << V[i].l <<std::endl;
   //for (int i = 0; i < config->order; i++) std::cout << "time " << config->vertices[i].t << " - itau " << config->vertices[i].itau << " - l " << config->vertices[i].l <<std::endl;

   double c_to_cp = (proba_mv_up > 0) ? 0.5 : 1;
   double cp_to_c = ((config->order - 1) > order_min) ? 0.5 : 1;
   //TRIQS_PRINT(c_to_cp);
   //TRIQS_PRINT(cp_to_c);

   //TRIQS_PRINT(new_det);
   //TRIQS_PRINT(config->det);
   //TRIQS_PRINT(new_proba);
   //TRIQS_PRINT(config->proba);

   auto det_ratio = real(coeffs[(config->order - 1)%4] * new_det) / real(coeffs[config->order %4] * config->det); 
   //TRIQS_PRINT(det_ratio);
   auto proba_ratio = cp_to_c * config->order / (2 * lambda * proba * c_to_cp);
   //TRIQS_PRINT(proba_ratio);

   //auto weight = real(coeffs[(config->order - 1)%4] * new_det) / real(coeffs[config->order %4] * config->det) * 
     //cp_to_c * config->proba / (2 * lambda * new_proba * c_to_cp);
   auto weight = det_ratio*proba_ratio;
   //TRIQS_PRINT(weight);

   return isnan(weight) ? 0 : weight;

   //return det_ratio * config->order / (2 * i * t_max * lambda);
  }

 }

 std::complex<double> accept(){
  config->det = new_det;
  config->vertices = V;

  (moving_up) ? config->order++ : config->order--;
  return 1.0;
 }

 void reject() {}
};

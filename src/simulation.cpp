#include <cpp11.hpp>
#include <iostream>
#include <random>
#include <map>
#include <string>

//   https://cpp11.r-lib.org/index.html

auto lm_fit = cpp11::package("stats")["lm.fit"];

void update_cases(std::map<std::string, int> &cases, const int x, const int v) {
  switch(x) {
  case 1:
    v == 1 ? cases["MM"]++ : v == 2 ? cases["MW"]++ : cases["MX"]++;
    break;
  case 2:
    break;
  case 3:
    break;
  }
}

[[cpp11::register]]
cpp11::writable::list make_x_v_cpp(const int n, cpp11::doubles probs, 
                            cpp11::doubles pmat_m, cpp11::doubles pmat_w, 
                            cpp11::doubles pmat_x) {
  std::random_device rd;
  std::mt19937 gen(rd());
  
  std::discrete_distribution<> x_dist(probs.begin(), probs.end());
  std::discrete_distribution<> vm_dist(pmat_m.begin(), pmat_m.end());
  std::discrete_distribution<> vw_dist(pmat_w.begin(), pmat_w.end());
  std::discrete_distribution<> vx_dist(pmat_x.begin(), pmat_x.end());
  
  std::map<std::string, int> cases {
    {"MM", 0}, {"MW", 0}, {"MX", 0},
    {"WM", 0}, {"WW", 0}, {"WX", 0},
    {"XM", 0}, {"XW", 0}, {"XX", 0}
  };
  
  using namespace cpp11::literals;
  cpp11::writable::integers x(n);
  cpp11::writable::integers v(n);
  cpp11::writable::list out;
  auto v_i = v.begin();
  for (auto i : x) {
    i = x_dist(gen) + 1;
    *v_i = i == 1 ? vm_dist(gen) : i == 2 ? vw_dist(gen) : vx_dist(gen);
    update_cases(cases, i, *v_i);
    ++v_i;
  }
  out.push_back({"X"_nm = x});
  out.push_back({"V"_nm = v});
  cpp11::writable::data_frame counts(
      {"MM"_nm = cases["MM"], "MW"_nm = cases["MW"]}
  );
  out.push_back({"counts"_nm = counts});
  return(out);
}

[[cpp11::register]]
cpp11::list fit_models_cpp(cpp11::doubles_matrix<> mat_act, 
                           cpp11::doubles_matrix<> mat_obs, cpp11::doubles y) {
  using namespace cpp11::literals;
  cpp11::writable::list models;
  models.push_back({"actual"_nm = lm_fit(mat_act, y)});
  models.push_back({"observed"_nm = lm_fit(mat_obs, y)});
  return(models);
}

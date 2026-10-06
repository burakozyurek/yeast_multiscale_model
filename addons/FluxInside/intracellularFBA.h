#ifndef intracellularFBA_H
#define intracellularFBA_H
#include <vector>

void buildModelFromMat(const char *file_path);
void changeBound(int col, double lb, double ub);
int optimizeModel(std::string obj_fun = "biomass", int obj_dir = 1, bool warm_start = true, bool print_error = false);
double getFVAByName(const std::string &rxn_name, const std::string &obj_name, const double flux);
std::vector<double> getFluxes();
#endif

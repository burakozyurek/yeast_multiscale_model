#include "intracellularFBA.h"
// #include "buildModel.h"
#include <iostream>
#include <string>
#include <chrono>
#include <vector>
#include <string>
#include <glpk.h>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <cstdlib>
#include <cstdarg>
#include <csetjmp>
#include <cstdio>
#include <cstdarg>
#include <unordered_map>

extern "C"
{
#include <matio.h>

}

using namespace std;

// static void my_glpk_error_hook(void *info, const char *msg);


class IntracellularFBA
{
private:
public:
    ~IntracellularFBA()
    {
        // if (lp)
        //     glp_delete_prob(lp);
    }

    glp_prob *lp;
    glp_smcp parm;
    // store reaction names and a fast lookup map for name->column index (zero-based)
    std::vector<std::string> reaction_names;
    std::unordered_map<std::string, int> reaction_index_map;
    // Initialize GLPK parameters in the constructor
    IntracellularFBA() : lp(nullptr)
    {
        glp_init_smcp(&parm);
        parm.msg_lev = GLP_MSG_OFF;
        parm.meth = GLP_DUALP;
        // parm.it_lim = 2000;
    }
    // IntracellularFBA(glp_prob *lp)
    // {
    //     this->lp = lp;
    //

    void buildModelFromMat(const char *file_path)
    {
        std::cout << "Reading MAT file...\n";
        mat_t *matfp = Mat_Open(file_path, MAT_ACC_RDONLY);
        if (!matfp)
        {
            std::cerr << "Error: Could not open MAT file: " << file_path << std::endl;
        }

        // Print MAT file version
        if (Mat_GetVersion(matfp) == MAT_FT_MAT4)
            std::cout << "MAT file version: v4\n";
        else if (Mat_GetVersion(matfp) == MAT_FT_MAT5)
            std::cout << "MAT file version: v5\n";
        else if (Mat_GetVersion(matfp) == MAT_FT_MAT73)
            std::cout << "MAT file version: v7.3\n";
        else
            std::cout << "MAT file version: unknown\n";
        // Mat_Close(matfp);
        // return nullptr;

        std::cout << "Reading variables from: " << file_path << "\n";
        matvar_t *var;
        var = Mat_VarReadNext(matfp); // Read the first variable (assumed struct)
        std::cout << "Variable name: " << var->name << "\n";

        // List all field names in the struct
        char *const *field_names = Mat_VarGetStructFieldnames(var);
        int num_fields = Mat_VarGetNumberOfFields(var);
        std::cout << "Number of fields: " << num_fields << "\n\n";
        for (int i = 0; i < num_fields; i++)
            std::cout << "Field " << i << ": " << field_names[i] << "\n";

        // std::cout << "Freeing variable: " << var->name << "...\n";

        matvar_t *S = Mat_VarGetStructFieldByName(var, "S", 0);
        matvar_t *lb = Mat_VarGetStructFieldByName(var, "lb", 0);
        matvar_t *ub = Mat_VarGetStructFieldByName(var, "ub", 0);
        matvar_t *rxn = Mat_VarGetStructFieldByName(var, "rxns", 0);
        matvar_t *objective = Mat_VarGetStructFieldByName(var, "c", 0);
        // matvar_t *ver = Mat_VarGetStructFieldByName(var, "version", 0);
        if (!var)
        {
            Mat_Close(matfp);
            throw runtime_error("MAT variable read failed");
        }
        if (!S || !lb || !ub)
        {
            Mat_Close(matfp);
            throw runtime_error("Missing required fields S/lb/ub");
        }

        // char *version_str = (char *)ver->data;
        //  printf("Human-GEM version: %s\n", version_str);
        //  // Print model version to console
        //  if (ver && ver->data_type == MAT_T_DOUBLE && ver->data != nullptr)
        //  {
        //      double *ver_data = static_cast<double *>(ver->data);
        //      std::cout << "Model version: " << ver_data[0] << "\n";
        //  }
        //  else
        //  {
        //      std::cerr << "Version field not found or not of type double.\n";
        //  }

        // store reaction names
        std::vector<std::string> reaction_names;
        if (rxn->class_type == MAT_C_CELL)
        {
            size_t num_rxns = rxn->dims[0];
            matvar_t **rxn_cells = static_cast<matvar_t **>(rxn->data);
            for (size_t i = 0; i < num_rxns; ++i)
            {
                matvar_t *rxn_str = rxn_cells[i];
                if (rxn_str && (rxn_str->data_type == MAT_T_UTF8 ||
                                rxn_str->data_type == MAT_T_UINT8))
                {
                    std::string reaction_name(static_cast<char *>(rxn_str->data), rxn_str->nbytes);
                    reaction_names.push_back(reaction_name);
                    // std::cout << "Reaction " << i << ": " << reaction_name << "\n";
                }
                else if (rxn_str && rxn_str->data_type == MAT_T_DOUBLE)
                    reaction_names.push_back("[numeric data]");
                // std::cout << "Reaction " << i << ": [numeric data]" << "\n";
                else
                    reaction_names.push_back("[unknown type]");
                // std::cout << "Reaction " << i << ": [unknown type]" << "\n";
            }
        }
        else
            std::cout << "rxnNames is not a cell array.\n";

        // write reaction names and bounds to csv file
        std::ofstream csv_file("c:\\Users\\burak\\My Drive\\ecoli-multiscale\\yeast-GEM-9.0.2\\reaction_bounds.csv");
        if (csv_file.is_open())
        {
            csv_file << "Reaction,Lower Bound,Upper Bound,Objective Coefficient\n";
            for (size_t i = 0; i < reaction_names.size(); ++i)
            {
                double lb_value = (lb->data_type == MAT_T_DOUBLE && lb->data != nullptr && i < (size_t)lb->dims[0]) ? static_cast<double *>(lb->data)[i] : 0.0;
                double ub_value = (ub->data_type == MAT_T_DOUBLE && ub->data != nullptr && i < (size_t)ub->dims[0]) ? static_cast<double *>(ub->data)[i] : 0.0;
                double obj_coef = (objective && objective->data_type == MAT_T_DOUBLE && objective->data != nullptr && i < (size_t)objective->dims[0]) ? static_cast<double *>(objective->data)[i] : 0.0;
                csv_file << reaction_names[i] << "," << lb_value << "," << ub_value << "," << obj_coef << "\n";
            }
            csv_file.close();
            std::cout << "Reaction bounds written to reaction_bounds.csv\n";
        }
        else
        {
            std::cerr << "Unable to open file for writing: reaction_bounds.csv\n";
        }

        // if (rxn->class_type == MAT_C_CELL)

        // Persist reaction names to the class and build a fast lookup map
        this->reaction_names = reaction_names;
        this->reaction_index_map.clear();
        for (size_t i = 0; i < this->reaction_names.size(); ++i)
        {
            this->reaction_index_map[this->reaction_names[i]] = (int)i; // zero-based
        }


        // Build the stoichiometric matrix S_matrix from the sparse matrix S
        std::vector<std::vector<double>> S_matrix;

        // read matrix from mat file with matio and handle sparse and dense conditions
        if (S)
        {
            if (S->class_type == MAT_C_SPARSE)
            {
                mat_sparse_t *sp = static_cast<mat_sparse_t *>(S->data);
                size_t n_rows = S->dims[0];
                size_t n_cols = S->dims[1];
                size_t n_nz = sp->ndata;
                std::cout << "Sparse matrix detected.\n";
                std::cout << "\nDense dimensions of matrix: " << n_rows << " x " << n_cols << "\n";
                std::cout << "Number of non-zero elements: " << n_nz << "\n";
                std::cout << "1" << std::endl;
                double *values = static_cast<double *>(sp->data);
                std::cout << "2" << std::endl;
                mat_uint32_t *ir = sp->ir;
                mat_uint32_t *jc = sp->jc;
                std::cout << "3" << std::endl;
                S_matrix.clear(); // Clear any existing data in S_matrix
                std::cout << "444" << std::endl;
                // Convert sparse CSC (jc, ir, values) into triplet list (row, col, value)
                for (size_t col = 0; col < n_cols; ++col)
                {
                    size_t start = jc[col];
                    size_t end = jc[col + 1];
                    for (size_t idx = start; idx < end; ++idx)
                    {
                        size_t row = ir[idx];
                        double val = values[idx];
                        // store as (row, col, value)
                        S_matrix.push_back({static_cast<double>(row), static_cast<double>(col), val});
                    }
                }
            }
            else if ((S->class_type == MAT_C_DOUBLE || S->class_type == MAT_C_DOUBLE) && S->data != nullptr)
            {
                std::cout << "Dense matrix detected.\n";
                // handle dense matrix
                double *dense_data = static_cast<double *>(S->data);
                size_t n_rows = S->dims[0];
                size_t n_cols = S->dims[1];
                S_matrix.clear();
                for (size_t i = 0; i < n_rows; ++i)
                    for (size_t j = 0; j < n_cols; ++j)
                    {
                        double value = dense_data[i + j * n_rows]; // MATLAB column-major order
                        if (value != 0.0)
                            S_matrix.push_back({static_cast<double>(i), static_cast<double>(j), value});
                    }
            }
            else
                std::cout << "\nField is not a double or sparse matrix.\n";
        }
        else
        {
            std::cout << "Field 'S' not found in the variable.\n";
            // return 1;
        }
        std::cout << "6" << std::endl;
        // Determine the number of columns and rows in S_matrix
        int real_col = 0;
        for (const auto &col : S_matrix)
            if (!col.empty() && col[1] > real_col)
                real_col = static_cast<int>(col[1]);

        real_col += 1; // Adjust for 1-based indexing in GLPK
        std::cout << "7" << std::endl;
        int real_row = 0;
        for (const auto &col : S_matrix)
            if (!col.empty() && col[0] > real_row)
                real_row = static_cast<int>(col[0]);

        real_row += 1; // Adjust for 1-based indexing in GLPK
        std::cout << "8" << std::endl;
        // Extract lower and upper bounds for each reaction
        double *lb_data = nullptr;
        double *ub_data = nullptr;
        size_t lb_size = 0;
        size_t ub_size = 0;

        // compute sizes (product of dims) and get data pointers safely
        if (lb && lb->data_type == MAT_T_DOUBLE && lb->data != nullptr)
        {
            lb_data = static_cast<double *>(lb->data);
            lb_size = 1;
            for (size_t i = 0; i < lb->rank; ++i)
                lb_size *= lb->dims[i];
        }
        else
        {
            std::cerr << "Lower bounds not found or not of type double.\n";
        }
        std::cout << "9" << std::endl;
        if (ub && ub->data_type == MAT_T_DOUBLE && ub->data != nullptr)
        {
            ub_data = static_cast<double *>(ub->data);
            ub_size = 1;
            for (size_t i = 0; i < ub->rank; ++i)
                ub_size *= ub->dims[i];
        }
        else
        {
            std::cerr << "Upper bounds not found or not of type double.\n";
        }

        std::cout << "10" << std::endl;

        // If sizes don't match, handle a common case: a scalar upper bound should be broadcast
        if (lb_size != (size_t)real_col)
        {
            std::cerr << "Warning: lower-bounds length (" << lb_size << ") does not match number of reactions (" << real_col << ").\n";
        }

        std::vector<double> ub_broadcast;
        double *ub_ptr = ub_data;
        if (ub_size == 1 && ub_data != nullptr)
        {
            // broadcast scalar upper bound to all reactions
            ub_broadcast.assign(real_col, ub_data[0]);
            ub_ptr = ub_broadcast.data();
        }
        else if (ub_size == (size_t)real_col)
        {
            ub_ptr = ub_data;
        }
        else
        {
            std::cerr << "Warning: upper-bounds length (" << ub_size << ") does not match number of reactions (" << real_col << ").\n";
            // fall back to using provided ub_data where possible (may cause out-of-range access if incorrect)
            ub_ptr = ub_data;
        }

        std::cout << "11" << std::endl;
        Mat_Close(matfp);

        std::cout << "13" << std::endl;

        this->lp = glp_create_prob();
        if (!this->lp)
            throw std::runtime_error("glp_create_prob failed");

        // glp_set_prob_name(lp, "FBA");

        // dimensions of the matrix
        std::cout << "S_matrix size: " << S_matrix.size() << " x " << S_matrix[0].size() << std::endl;
        std::cout << "14" << std::endl;
        std::cout << "Dense size: " << real_row << " x " << real_col << std::endl;
        glp_add_cols(lp, real_col); // +1 for the objective variable

        // use rxn names as glp col name
        std::cout << "15" << std::endl;
        double *obj_data = (objective && objective->data_type == MAT_T_DOUBLE && objective->data != nullptr) ? static_cast<double *>(objective->data) : nullptr;
        for (int i = 1; i <= real_col; i++)
        {
            // glp_set_col_name(lp, i, ("v" + std::to_string(i)).c_str()); // Name each reaction variable
            glp_set_col_name(lp, i, reaction_names[i - 1].c_str()); // Name each reaction variable
            glp_set_col_kind(lp, i, GLP_CV);                        // Continuous variable

            double coef = obj_data ? obj_data[i - 1] : 0.0;
            glp_set_obj_coef(lp, i, coef); // Default objective coefficient
            if ((obj_data[i - 1] == 1))
            {
                cout << reaction_names[i - 1] << ", obj coef: " << obj_data[i - 1] << endl;
            }
        }
        std::cout << "5" << std::endl;
        glp_set_obj_dir(lp, GLP_MAX);


        for (int i = 1; i <= real_col; i++)
        {
            // When lower bound equals upper bound, use GLP_FX instead of GLP_DB
            if (lb_data[i - 1] == ub_data[i - 1])
            {
                glp_set_col_bnds(lp, i, GLP_FX, lb_data[i - 1], lb_data[i - 1]);
            }
            else if (lb_data[i - 1] > ub_data[i - 1])
            {
                // If lower bound > upper bound, that's an error
                std::cerr << "Error: Column " << i << " has lb (" << lb_data[i - 1]
                          << ") > ub (" << ub_data[i - 1] << ")" << std::endl;
                // return 1;
            }
            else
            {
                // Normal case: lb < ub
                glp_set_col_bnds(lp, i, GLP_DB, lb_data[i - 1], ub_data[i - 1]);
            }
        }
        std::cout << "16" << std::endl;
        // glp_set_col_bnds(lp, 47, GLP_DB, -10, 1000);
        // glp_set_col_bnds(lp, 52, GLP_DB, -1, 0);
        // glp_set_col_bnds(lp, 55, GLP_DB, -1000, 1000);
        // glp_set_col_bnds(lp, 56, GLP_DB, -1000, 1000);
        // // glp_set_col_bnds(lp, 7735, GLP_FX, glutamate, glutamate);
        // glp_set_col_bnds(lp, 94, GLP_DB, -20, 0); // O2 difusion

        glp_add_rows(lp, real_row);
        std::cout << "5" << std::endl;
        for (int i = 1; i <= real_row; i++)
            glp_set_row_bnds(lp, i, GLP_FX, 0.0, 0.0);
        int k = 0;

        // Declare arrays with appropriate sizes
        std::vector<int> ia(S_matrix.size() + 1, 0);
        std::vector<int> ja(S_matrix.size() + 1, 0);
        std::vector<double> ar(S_matrix.size() + 1, 0.0);

        for (size_t i = 0; i < S_matrix.size(); i++)
        {
            k++;
            ia[k] = static_cast<int>(S_matrix[i][0] + 1); // Convert to 1-based index
            ja[k] = static_cast<int>(S_matrix[i][1] + 1); // Convert to 1-based index
            ar[k] = S_matrix[i][2];
        }
        std::cout << "17" << std::endl;
        glp_load_matrix(lp, k, ia.data(), ja.data(), ar.data());
        std::cout << "18" << std::endl;

        auto start_opt = std::chrono::high_resolution_clock::now();
        glp_simplex(lp, &parm);
        auto end_opt = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_opt - start_opt;
        std::cout << "Optimization took " << elapsed.count() << " seconds.\n";
        // model_initiated = true;
        //  return lp; // Return the pointer as is
        //   Mat_VarFree(var);
        //   std::cout << "Variable freed.\n";
    }

    void changeBound(int col, double lb, double ub)
    {
        // Validate bounds before setting
        if (col + 1 < 0 || col + 1 >= glp_get_num_cols(lp))
        {
            throw runtime_error("Invalid column index in changeBound: " + std::to_string(col + 1));
        }

        // Handle NaN/Inf
        if (std::isnan(lb) || std::isnan(ub) || std::isinf(lb) || std::isinf(ub))
        {
            throw runtime_error("NaN or Inf detected in changeBound bounds");
        }

        if (lb == ub)
            glp_set_col_bnds(lp, col + 1, GLP_FX, lb, lb);
        else if (lb < ub)
            glp_set_col_bnds(lp, col + 1, GLP_DB, lb, ub);
        else
            throw runtime_error("lb > ub in changeBound: lb=" + std::to_string(lb) + " ub=" + std::to_string(ub));
    }

    int optimizeModel(std::string obj_fun = "biomass", int obj_dir = 1, bool warm_start = true, bool print_error = false)
    {
        auto start_opt = std::chrono::high_resolution_clock::now();

        glp_init_smcp(&parm);
        if (print_error)
            parm.msg_lev = GLP_MSG_ERR;
        else
            parm.msg_lev = GLP_MSG_OFF;
        parm.meth = GLP_DUALP;
        


        if (obj_dir == 1)
        {
            glp_set_obj_coef(lp, getReactionIndex(obj_fun) + 1, 1.0);
        }
        if (obj_dir == -1)
        {
            glp_set_obj_coef(lp, getReactionIndex(obj_fun) + 1, -1.0);
        }
        
        if (!warm_start)
        {
            glp_std_basis(lp); // Reset the basis so each solve starts cold
        }

        // glp_iptcp iparm;
        // glp_init_iptcp(&iparm);
        // iparm.msg_lev = GLP_MSG_ERR;

        
        glp_simplex(lp, &parm);
        // glp_interior(lp, &iparm);
        auto end_opt = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed_opt = end_opt - start_opt;
        //std::cout << "Model optimized in " << elapsed_opt.count() << " ms." << endl;
        
        int status = glp_get_status(lp);
        return status;
    }

    const std::vector<double> getFluxes() const
    {
        int real_col = glp_get_num_cols(lp);
        // print all flux values with their colnames
        for (int i = 1; i <= real_col; i++)
        {
            const char *reaction_name = glp_get_col_name(lp, i);
            // std::cout << "flx_" << i << ": " << round(glp_get_col_prim(lp, i) * 100.0) / 100 << "\t" << reaction_name << endl;
        }

        std::vector<double> fluxes;
        for (int i = 1; i <= real_col; ++i)
            fluxes.push_back(glp_get_col_prim(lp, i));

        // glp_delete_prob(lp);
        // std::cout << "GLPK problem deleted.\n";
        return fluxes;
    }

    void print_reactions()
    {
        int real_col = glp_get_num_cols(lp);
        for (int i = 1; i <= real_col; i++)
        {
            const char *reaction_name = glp_get_col_name(lp, i);
            std::cout << "Reaction " << i << ": " << glp_get_col_lb(lp, i) << "\t- " << glp_get_col_ub(lp, i) << ", \tname: " << reaction_name << "\n";
        }
    }

    // Return zero-based column index for a given reaction name, or -1 if not found
    int getReactionIndex(const std::string &name) const
    {
        auto it = reaction_index_map.find(name);
        if (it == reaction_index_map.end())
            return -1;
        return it->second; // Convert to 1-based index for GLPK
    }

    // Convenience: change bounds by reaction name
    void changeBoundByName(const std::string &name, double lb, double ub)
    {
        int idx = getReactionIndex(name);
        if (idx < 0)
            throw std::runtime_error("Unknown reaction name in changeBoundByName: " + name);
        changeBound(idx, lb, ub);
    }

    // Get flux by reaction name (returns primal value)
    double getFluxByName(const std::string &name) const
    {
        int idx = -1;
        auto it = reaction_index_map.find(name);
        if (it == reaction_index_map.end())
            throw std::runtime_error("Unknown reaction name in getFluxByName: " + name);
        idx = it->second;

            return glp_get_col_prim(lp, idx + 1);
        
    }

    // Set objective coefficient by reaction name
    void setObjectiveByName(const std::string &name, double coef)
    {
        int idx = getReactionIndex(name);
        if (idx < 0)
            throw std::runtime_error("Unknown reaction name in setObjectiveByName: " + name);
        glp_set_obj_coef(lp, idx + 1, coef);
    }

    double getFVAByName(const std::string &rxn_name, const std::string &obj_name, const double flux)
    {
        // Implement FVA logic
        // change objective
        int idx = getReactionIndex(rxn_name);

        //store original objective and value    
        //const string obj_name = glp_get_obj_name(lp);
        int obj_dir = glp_get_obj_dir(lp);

        changeBoundByName(obj_name, flux, flux);
        setObjectiveByName(rxn_name, -1);

        glp_simplex(lp, &parm);

        //set back to original objective
        setObjectiveByName(obj_name, obj_dir);
        changeBoundByName(obj_name, 0, 1000);
        return getFluxByName(rxn_name);
    }
};
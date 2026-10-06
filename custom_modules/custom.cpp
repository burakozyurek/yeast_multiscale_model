/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
#include "../addons/FluxInside/fba_timing_logger.h"
#include "../addons/FluxInside/intracellularFBA.hpp"
#include <fstream>
#include <unordered_map>
#include <mutex>
#include <algorithm>

static std::mutex fba_mutex;
static IntracellularFBA yeast_fba;
static bool fba_initialized = false;

void initialize_fba(void)
{
	if (!fba_initialized)
	{
		std::lock_guard<std::mutex> lock(fba_mutex);

		// Double-check before initializing (in case another thread won the race)
		if (!fba_initialized)
		{
			// get path from parameters
			std::string model_path = parameters.strings("GEM_path");
			yeast_fba.buildModelFromMat(model_path.c_str());
			fba_initialized = true;
			yeast_fba.print_reactions();
		}
	}
}

void create_cell_types(void)
{
	// Initialize FBA model ONCE at startup (before OpenMP parallelization begins)
	// This avoids thread-safety issues during GLPK/matio initialization
	initialize_fba();

	// set the random seed
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}

	/*
	   Put any modifications to default cell definition here if you
	   want to have "inherited" by other cell types.

	   This is a good place to set default functions.
	*/

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment(&microenvironment);

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL;
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based;
	cell_defaults.functions.custom_cell_rule = NULL;
	cell_defaults.functions.contact_function = NULL;

	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;

	/*
	   This parses the cell definitions in the XML config file.
	*/

	initialize_cell_definitions_from_pugixml();

	/*
	   This builds the map of cell definitions and summarizes the setup.
	*/

	build_cell_definitions_maps();

	/*
	   This intializes cell signal and response dictionaries
	*/

	setup_signal_behavior_dictionaries();

	/*
	   Cell rule definitions
	*/

	setup_cell_rules();

	/*
	   Put any modifications to individual cell definitions here.

	   This is a good place to set custom functions.
	*/

	cell_defaults.functions.update_phenotype = phenotype_function;
	cell_defaults.functions.custom_cell_rule = custom_function;
	cell_defaults.functions.contact_function = contact_function;

	/*
	   This builds the map of cell definitions and summarizes the setup.
	*/

	display_cell_definitions(std::cout);

	return;
}

void setup_microenvironment(void)
{
	// set domain parameters

	// put any custom code to set non-homogeneous initial conditions or
	// extra Dirichlet nodes here.

	// initialize BioFVM

	initialize_microenvironment();

	return;
}

void setup_tissue(void)
{
	double Xmin = microenvironment.mesh.bounding_box[0];
	double Ymin = microenvironment.mesh.bounding_box[1];
	double Zmin = microenvironment.mesh.bounding_box[2];

	double Xmax = microenvironment.mesh.bounding_box[3];
	double Ymax = microenvironment.mesh.bounding_box[4];
	double Zmax = microenvironment.mesh.bounding_box[5];

	if (default_microenvironment_options.simulate_2D == true)
	{
		Zmin = 0.0;
		Zmax = 0.0;
	}

	double Xrange = Xmax - Xmin;
	double Yrange = Ymax - Ymin;
	double Zrange = Zmax - Zmin;

	// create some of each type of cell

	Cell *pC;

	for (int k = 0; k < cell_definitions_by_index.size(); k++)
	{
		Cell_Definition *pCD = cell_definitions_by_index[k];
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
		for (int n = 0; n < parameters.ints("number_of_cells"); n++)
		{
			std::vector<double> position = {0, 0, 0};
			position[0] = Xmin + UniformRandom() * Xrange;
			position[1] = Ymin + UniformRandom() * Yrange;
			position[2] = Zmin + UniformRandom() * Zrange;

			pC = create_cell(*pCD);
			pC->assign_position(position);
		}
	}
	std::cout << std::endl;

	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	set_parameters_from_distributions();

	return;
}

void phenotype_function(Cell *pCell, Phenotype &phenotype, double dt)
{
	// skip entire function if cell is dead
	if (phenotype.death.dead)
		return;

	// parameter initialization------------------------------------------------
	double t = PhysiCell_globals.current_time;

	double Vbase = parameters.ints("base_volume");
	double Vcurr = pCell->phenotype.volume.total;
	double biomass_scale = Vcurr / Vbase;
	double metabolic_stress_limit = parameters.ints("metabolic_stress_limit");
	
	bool warm_start = parameters.bools("warm_start");
	bool use_Fva = parameters.bools("use_Fva");

	// O2 Kinetic Limit (Michaelis-Menten)
	double o2_vmax = parameters.doubles("o2_vmax"); // mmol/gDW/h
	double o2_km = parameters.doubles("o2_km");		// mmol/L

	// GLC Kinetic Limit (Michaelis-Menten)
	double glc_vmax = parameters.doubles("glc_vmax"); // mmol/gDW/h
	double glc_km = parameters.doubles("glc_km");	  // mmol/L

	//=========================================================================
	// A. Environment Sensing (BioFVM → PhysiCell)
	//=========================================================================
	static int Noxygen = microenvironment.find_density_index("oxygen");
	double oxygen = pCell->nearest_density_vector()[Noxygen];

	static int Nglucose = microenvironment.find_density_index("glucose");
	double glucose = pCell->nearest_density_vector()[Nglucose];

	// static int Nferm = microenvironment.find_density_index("fermentation_activity");
	// double fermentation_product = pCell->nearest_density_vector()[Nferm];

	//=========================================================================
	// B. Constraint Setting (PhysiCell → FluxInside)
	//=========================================================================
	// calculate uptake bounds based on local oxygen concentration and cell volume
	double o2_uptake_bound = biomass_scale * o2_vmax * (oxygen / (o2_km + oxygen));		 // makaleler biomass scalingini fbaden sonra kullanıyor
	double glc_uptake_bound = biomass_scale * glc_vmax * (glucose / (glc_km + glucose)); // convert from 1/h to 1/min
	// double lac_secretion_bound = biomass_scale * vmax_lac * Km_lac / (Km_lac + lactate);

	bool optimization_success = false;
	std::vector<double> fluxes;

	std::lock_guard<std::mutex> lock(fba_mutex);

	yeast_fba.changeBoundByName("r_1992", std::min(0.0, -o2_uptake_bound), 0);		 // o2
	yeast_fba.changeBoundByName("r_1714", std::min(0.0, -glc_uptake_bound), 0); // Glucose

	//=====================================================================
	// C. FBA (FluxInside)
	//=====================================================================

	auto solve_start = std::chrono::high_resolution_clock::now();
	double working_set_before_mb = fba_timing_logger.get_process_working_set_mb();
	double private_bytes_before_mb = fba_timing_logger.get_process_private_bytes_mb();

	// Optimize the FBA model
	int optimization_status = yeast_fba.optimizeModel("r_2111", 1, warm_start, false);
	
	auto solve_end = std::chrono::high_resolution_clock::now();
	double solve_ms = std::chrono::duration<double, std::milli>(solve_end - solve_start).count();
	double working_set_after_mb = fba_timing_logger.get_process_working_set_mb();
	double private_bytes_after_mb = fba_timing_logger.get_process_private_bytes_mb();
	double working_set_delta_mb = working_set_after_mb - working_set_before_mb;
	double private_bytes_delta_mb = private_bytes_after_mb - private_bytes_before_mb;
	auto record_fba_timing = [&]() {
		int dead_cell_count = 0;
		for (auto cell : *all_cells)
		{
			if (cell->phenotype.death.dead)
				dead_cell_count++;
		}
		fba_timing_logger.record(t, solve_ms, working_set_delta_mb, private_bytes_delta_mb, dead_cell_count);
	};

	if (optimization_status == GLP_OPT || optimization_status == GLP_FEAS || optimization_status == GLP_NOFEAS)
	{
		// cout << fluxes << endl;
		// ---- Biomass growth
		double mu = yeast_fba.getFluxByName("r_2111"); // biomass flux

		if (!std::isfinite(mu))
		{
			std::cerr << "Non-finite biomass flux for cell " << pCell->ID << std::endl;
			pCell->start_death(1);
			record_fba_timing();
			return;
		}

		//cout << "mu:" << mu << endl;
		double metabolic_stress = pCell->custom_data["metabolic_stress"];

		if (mu < 0.1)
			metabolic_stress += dt;
		else
			metabolic_stress -= dt;

		if (metabolic_stress >= metabolic_stress_limit)
		{
			pCell->start_death(1);
			record_fba_timing();
		
			return;
		}
		metabolic_stress = std::max(0.0, metabolic_stress);
		pCell->custom_data["metabolic_stress"] = metabolic_stress;

		double growth_factor = exp(mu * dt / 60.0); // convert from 1/h to 1/min
		if (!std::isfinite(growth_factor))
		{
			std::cerr << "Non-finite growth factor for cell " << pCell->ID << std::endl;
			pCell->start_death(1);
			record_fba_timing();
			return;
		}

		//=====================================================================
		// 	D.	Phenotype Update (FluxInside → PhysiCell)
		//=====================================================================

		double Vnew = Vcurr * growth_factor;

		pCell->set_total_volume(Vnew);

		if (Vnew >= 2.0 * Vbase)
		{
			pCell->divide();
		}

		//=====================================================================
		// 	E. Secretion/Uptake (FluxInside → BioFVM)
		//=====================================================================

		double o2_flux;
		double glc_flux;
		if (use_Fva)
		{
			o2_flux = std::max(0.0,  -yeast_fba.getFVAByName("r_1992", "r_2111", mu)); 
			glc_flux = std::max(0.0, -yeast_fba.getFVAByName("r_1714", "r_2111", mu));
		}
		else
		{	
			o2_flux = std::max(0.0, -yeast_fba.getFluxByName("r_1992")); // mmol/gDW/h	 
			glc_flux = std::max(0.0, -yeast_fba.getFluxByName("r_1714")); // mmol/gDW/h
		}

		//cout << "Cell " << pCell->ID << " Oxygen: " << oxygen << " mmol/L, uptake rate: " << o2_flux << " mmol/gDW/h" << endl;
		//cout << "Cell " << pCell->ID << " Glucose: " << glucose << " mmol/L, uptake rate: " << glc_flux << " mmol/gDW/h" << endl;
		
		if (oxygen > 1e-8)
		{
			double uptake_rate_1_per_min = (o2_flux / oxygen) / 60.0; // convert from 1/h to 1/min
			set_single_behavior(pCell, "oxygen uptake", uptake_rate_1_per_min);
			// cout << "Oxygen: " << oxygen << " mmol/L, uptake rate: " << uptake_rate_1_per_min << " 1/min" << endl;
		}
		else
		{
			set_single_behavior(pCell, "oxygen uptake", 0.0);
		}

		if (glucose > 1e-8)
		{
			double uptake_rate_1_per_min = (glc_flux / glucose) / 60.0; // convert from 1/h to 1/min
			set_single_behavior(pCell, "glucose uptake", uptake_rate_1_per_min);
			// cout << "Glucose: " << glucose << " mmol/L, uptake rate: " << uptake_rate_1_per_min << " 1/min" << endl;
		}

		else
		{
			set_single_behavior(pCell, "glucose uptake", 0.0);
		}
	}
 
	else
	{
		std::cout << "Optimization failed for cell " << pCell->ID << " with status " << optimization_status << std::endl;
		// set_single_behavior(pCell, "necrosis", 0.1);
		// pCell->start_death(1); // start necrosis
		record_fba_timing();
		return;
	}
	record_fba_timing();
}

void custom_function(Cell *pCell, Phenotype &phenotype, double dt)
{
	return;
}

void contact_function(Cell *pMe, Phenotype &phenoMe, Cell *pOther, Phenotype &phenoOther, double dt)
{
	return;
}

std::vector<std::string> my_coloring_function(Cell *pCell)
{
	std::vector<std::string> output(4, "black");
	if ((bool)get_single_signal(pCell, "dead"))
	{
		output[0] = "rgb(128,70,35)"; // brown cytoplasm fill
		output[1] = "rgb(128,35,35)"; // brown cytoplasm outline
		output[2] = "rgb(64,35,35)";  // dark brown nucleus fill
		output[3] = "rgb(64,35,35)";  // dark brown nucleus outline
		return output;
	}

	if (pCell->type_name == "yeast")
	{
		output[0] = "rgb(123, 255, 0)"; // cytoplasm fill
		output[2] = "rgb(63, 116, 29)"; // nuclear fill
		// outlines will be set based on stress below
		double stress = pCell->custom_data["metabolic_stress"]; // accumulated stress
		// cout << "Cell " << pCell->ID << " metabolic stress: " << stress << "\n";
		std::string stress_color = "rgb(" + std::to_string(200) + "," + std::to_string(20) + "," + std::to_string(20) + ")";
		if (stress > 0)
		{
			output[1] = stress_color; // cytoplasm outline
			output[3] = stress_color; // nuclear outline
		}

		output[1] = "rgb(185, 255, 119)";
		output[3] = "rgb(85, 139, 52)";
	}

	return output;
}

std::vector<std::string> (*cell_coloring_function)(Cell *) = my_coloring_function;
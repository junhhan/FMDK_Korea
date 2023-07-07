#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mtwister_single.h"

// VARIABLES
#define _simperiod 20 // 9 days of maximum period between E and Ic + 11 days of observation
#define _pigs_per_pen 50

#define _ncol_disease 8
#define _s 0
#define _e 1
#define _is 2
#define _ic 3
#define _sm 4
#define _em 5
#define _ism 6
#define _r 7

#define _ncol_pen 2
#define _col_s 0
#define _col_m 1

#define _ncol_transition 3
#define _col_eis 0
#define _col_isic 1 
#define _col_ismr 2 

// E to Is: Binomial (97, 0.02) [1 ~ 7] - Kinsley et al. 2016
#define _n_eis 97
#define _p_eis 0.02
#define _min_eis 1
#define _max_eis 7

// Is to Ic: Binomial (66, 0.02) - Kinsley et al. 2016
#define _n_isic 66
#define _p_isic 0.02

// E to Ic: Binomial (35, 0.06) [1 ~ 9] - Kinsley et al. 2016
#define _n_eic 35
#define _p_eic 0.06
#define _min_eic 1
#define _max_eic 9

// E(m) to Is(m) or Ic(m): Uniform (1, 2) - Parida et al. 2007
// Is(m) to R: Uniform (3, 7) - Parida et al. 2007
#define _min_ismr 3
#define _max_ismr 7

// FUNCTIONS
int _Random_binomial(int n, double p);
double _Random_uniform(double a, double b);
void _Infection(int **animal_transition, int d_eis, int d_isic, int d_ismr, double alpha, int s_or_m, int i);
void _Intro(int *animal_status, int **animal_transition, int *compartment_size, int d_eis, int d_isic, int d_ismr, double alpha, int s_or_m);
void _Calculate_lambda(int **pen_inf, int n_pig, int n_pen, double beta_w, double beta_b, double psi, double *lambda);
void _Spread(int *animal_status, int **animal_transition, int **pen_inf, int *compartment_size, int n_pig, int n_pen, int mod_id,
	double beta_w, double beta_b, double psi, double alpha, int d_eis, int d_isic, int d_ismr);

// MAIN FUNCTION
void _FMD_sim(int *_seed, int *_n_pig, int *_mod_id,
	double *_rho, double *_beta, double *_omega, double *_psi, double *_alpha, int *_d_eis, int *_d_isic, int *_d_ismr, int *_tau, int *_res_ts) {
	int i;
	int seed= (*_seed);
	_MTinit(seed);
	
	int n_pig= (*_n_pig);
	int n_pen= (int)ceil((double)n_pig / _pigs_per_pen);
	double rho= (*_rho);
	double beta_w= (*_beta);
	double omega= (*_omega);
	double beta_b= beta_w * omega;
	double psi= (*_psi);
	double alpha= (*_alpha);
	int d_eis= (*_d_eis);
	int d_isic= (*_d_isic);
	int d_ismr= (*_d_ismr);
	int tau= (*_tau);

	int mod_id= (*_mod_id);
	if (mod_id == 1) {
		rho= 0.0;
	} else if (mod_id == 2) {
		psi= 1.0;
	} else if (mod_id == 3) {
		alpha= 0.0;
	}
	
	// An array of animal infection status
	int *animal_status;
	animal_status= calloc(n_pig, sizeof(int));
	
	// An array of days-left for transition between infection status
	int **animal_transition;
	animal_transition= calloc(n_pig, sizeof(int *));

	// An array of the number of animals in each infection status
	int *compartment_size;
	compartment_size= calloc(_ncol_disease, sizeof(int));

	for (i= 0; i < n_pig; ++i) {
		if (_Random_binomial(1, rho) == 1 && i > 0) {
			animal_status[i]= _sm;
			++compartment_size[_sm];
		} else {
			animal_status[i]= _s;
			++compartment_size[_s];
		}
		animal_transition[i]= calloc(_ncol_transition, sizeof(int));
	}

	// An array of infectious pigs per pen
	int **pen_inf;
	pen_inf= calloc(n_pen, sizeof(int *));
	for (i= 0; i < n_pen; ++i) {
		pen_inf[i]= calloc(_ncol_pen, sizeof(int));
	}

	// Simulation
	for (i= tau; i < _simperiod; ++i) {
		if (i == tau) {
			_Intro(animal_status, animal_transition, compartment_size, d_eis, d_isic, d_ismr, alpha, _col_s);
		}
		_Spread(animal_status, animal_transition, pen_inf, compartment_size, n_pig, n_pen, mod_id, beta_w, beta_b, psi, alpha, d_eis, d_isic, d_ismr);
		_res_ts[i]= compartment_size[_ic];
		compartment_size[_ic]= 0;
	}

	// Free memories
	for (i= 0; i < n_pig; ++i) {
		free(animal_transition[i]);
	}
	for (i= 0; i < n_pen; ++i) {
		free(pen_inf[i]);
	}
	free(animal_status);
	free(animal_transition);
	free(pen_inf);
	free(compartment_size);
}

// FUNCTIONS IN DETAIL
int _Random_binomial(int n, double p) {
    int index;
	int n_pos= 0;
	int i= 0;
	if (p > 1) {p= 1.0;}
	if (p < 0) {p= 0.0;}

	double x;

	if (p == 0.0 || n == 0) {
		return 0;
	} else if (p == 1.0) {
		return n;
	} else {
	    do {
		    do {
				x= _Rand();
			} while (x > 1.0 || x < 0.0);

			if(p >= x) {index= 1;} else {index= 0;}
			n_pos += index;
			i++;
	    } while (i < n);
	
		return n_pos;
	}
}

void _Infection(int **animal_transition, int d_eis, int d_isic, int d_ismr, double alpha, int s_or_m, int i) {
	animal_transition[i][_col_eis]= d_eis;
	animal_transition[i][_col_isic]= d_isic;
	animal_transition[i][_col_ismr]= d_ismr;
	if (s_or_m == _col_s) {
		animal_transition[i][_col_ismr]= -1;
	} else {
		if (_Random_binomial(1, alpha) == 1) {
			animal_transition[i][_col_isic]= -1;
		} else {
			animal_transition[i][_col_ismr]= -1;
		}
	}
}

void _Intro(int *animal_status, int **animal_transition, int *compartment_size, int d_eis, int d_isic, int d_ismr, double alpha, int s_or_m) {
	int i= 0;
	animal_status[i]= _e;
	--compartment_size[_s];
	++compartment_size[_e];
	_Infection(animal_transition, d_eis, d_isic, d_ismr, alpha, _col_s, i);
}

void _Calculate_lambda(int **pen_inf, int n_pig, int n_pen, double beta_w, double beta_b, double psi, double *lambda) {
	int i, j, penid, pencrs;
	double sum;

	if (n_pen % 2 == 0) {
		pencrs= n_pen / 2;
	} else {
		pencrs= (n_pen + 1) / 2;
	}
	
	for (i= 0; i < n_pig; ++i) {
		penid= (int)floor((double)i / _pigs_per_pen);

		// Infectious pressure
		sum= 0.0;
		for (j= 0; j < n_pen; ++j) {
			if (abs(penid - j) == 1 || abs(penid - j) == (n_pen - 1) || abs(penid - j) == pencrs) {
				sum += ((double)pen_inf[j][_col_s] + (psi * (double)pen_inf[j][_col_m]));
			} 
		}
		lambda[i]= (beta_w * ((double)pen_inf[penid][_col_s] + (psi * (double)pen_inf[penid][_col_m]))) + (beta_b * sum);
	}
}

void _Spread(int *animal_status, int **animal_transition, int **pen_inf, int *compartment_size, int n_pig, int n_pen, int mod_id,
	double beta_w, double beta_b, double psi, double alpha, int d_eis, int d_isic, int d_ismr) {
	int i, penid;
	double *lambda, p;
	lambda= calloc(n_pig, sizeof(double));
	_Calculate_lambda(pen_inf, n_pig, n_pen, beta_w, beta_b, psi, lambda);
	
	for (i= 0; i < n_pig; ++i) {
		penid= (int)floor((double)i / _pigs_per_pen);
		// Susceptible
		if (animal_status[i] == _s) {
			p= 1.0 - exp(-1.0 * lambda[i]);
			if (p > 0.0 && _Random_binomial(1, p) == 1) {
				animal_status[i]= _e;
				--compartment_size[_s];
				++compartment_size[_e];
				_Infection(animal_transition, d_eis, d_isic, d_ismr, alpha, _col_s, i);
			}
		// Exposed
		} else if (animal_status[i] == _e) {
			if (animal_transition[i][_col_eis] > 0) {
				--animal_transition[i][_col_eis];
			} else {
				animal_status[i]= _is;
				--compartment_size[_e];
				++compartment_size[_is];
				++pen_inf[penid][_col_s];
			}
		// Subclinically infectious
		} else if (animal_status[i] == _is) {
			if (animal_transition[i][_col_isic] > 0) {
				--animal_transition[i][_col_isic];
			} else {
				animal_status[i]= _ic;
				--compartment_size[_is];
				++compartment_size[_ic];
				--pen_inf[penid][_col_s];
			}
		// Immunised susceptible
		} else if (animal_status[i] == _sm) {
			if (mod_id != 5) {
				p= 1.0 - exp(-1.0 * lambda[i]);
				if (p > 0.0 && _Random_binomial(1, p) == 1) {
					animal_status[i]= _em;
					--compartment_size[_sm];
					++compartment_size[_em];
					_Infection(animal_transition, d_eis, d_isic, d_ismr, alpha, _col_m, i);
				}
			}
		// Immunised exposed
		} else if (animal_status[i] == _em) {
			if (animal_transition[i][_col_eis] > 0) {
				--animal_transition[i][_col_eis];
			} else {
				animal_status[i]= _ism;
				--compartment_size[_em];
				++compartment_size[_ism];
				++pen_inf[penid][_col_m];
			}
		// Immunised subclinically infectious
		} else if (animal_status[i] == _ism) {
			if (animal_transition[i][_col_isic] > 0) {
				--animal_transition[i][_col_isic];
			} else if (animal_transition[i][_col_isic] == 0) {
				animal_status[i]= _ic;
				--compartment_size[_ism];
				++compartment_size[_ic];
				--pen_inf[penid][_col_m];
			}
			
			if (animal_transition[i][_col_ismr] > 0) {
				--animal_transition[i][_col_ismr];
			} else if (animal_transition[i][_col_ismr] == 0) {
				animal_status[i]= _r;
				--compartment_size[_ism];
				++compartment_size[_r];
				--pen_inf[penid][_col_m];
			}
		}
	}
	free(lambda);
}

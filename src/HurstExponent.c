#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_statistics.h> 
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_matrix.h>


void getHurstExponent(gsl_matrix *DFASpace, gsl_vector *c, gsl_matrix *cov, double *chisq)
{
	/*
	 * fit the space of DFA and s and get the Hurst Exponent
	 */
	gsl_vector *logDFA = gsl_vector_alloc(DFASpace->size1);
	gsl_vector *weights = gsl_vector_alloc(DFASpace->size1);
	gsl_matrix *logS = gsl_matrix_alloc(DFASpace->size1, 2);
	gsl_matrix_set_all(logS, 1.0f);	
	size_t i;
	for (i=0; i<DFASpace->size1 - 1; i++)
	{
		gsl_vector_set(logDFA, i, log(gsl_matrix_get(DFASpace, i, 1)));
		gsl_matrix_set(logS, i, 1, log(gsl_matrix_get(DFASpace, i, 0)));
		gsl_vector_set(weights, i, abs(gsl_matrix_get(logS, i+1, 1) - gsl_matrix_get(logS, i, 1)));
	}

	gsl_vector_set(logDFA, i, log(gsl_matrix_get(DFASpace, i, 1)));
	gsl_matrix_set(logS, i, 1, log(gsl_matrix_get(DFASpace, i, 0)));
	gsl_vector_set(weights, i, abs(gsl_matrix_get(logS, i, 1) - gsl_matrix_get(logS, i-1, 1)));

	gsl_multifit_linear_workspace *work  = gsl_multifit_linear_alloc(DFASpace->size1, 2);  
  	gsl_multifit_wlinear(logS, weights,logDFA, c, cov, chisq, work);
	
	gsl_vector_free(logDFA);
	gsl_matrix_free(logS);
	gsl_vector_free(weights);
	gsl_multifit_linear_free(work);

	gsl_vector_set(c, 0, exp(gsl_vector_get(c, 0)));

	*chisq /= DFASpace->size1;
}

int main( int argc, char *argv[] )
{
	// Get terminal arguments
	if ( argc != 6 )
	{
		return 1;
	}
	
	size_t rows = atoi(argv[1]);
	size_t offset = atoi(argv[2]);
	size_t sublen = atoi(argv[3]);
	char *inputFileName = argv[4];
	char *paramsFileName = argv[5];


	gsl_matrix *DFAMatrix = gsl_matrix_alloc( rows, 2 );
    	FILE *file;
	if ( ( file = fopen( inputFileName, "r" ) ) == 0 ) return 1;
	gsl_matrix_fscanf( file, DFAMatrix );
	fclose( file );
	
	gsl_matrix_view subDFAMatrix = gsl_matrix_submatrix(DFAMatrix, offset, 0, sublen, DFAMatrix->size2);

	// Fit DFA for calculate Hurst Exponent
	gsl_vector *params = gsl_vector_alloc(2);
	gsl_matrix *covarianceMatrix = gsl_matrix_alloc(2, 2);
	double chiSquare;
	getHurstExponent(&subDFAMatrix.matrix, params, covarianceMatrix, &chiSquare);
	
	gsl_matrix_free( DFAMatrix );
	gsl_matrix_free(covarianceMatrix);
	
	if ( ( file = fopen( paramsFileName, "w" ) ) == 0 ) return 1;
	gsl_vector_fprintf(file, params, "%f");
	gsl_vector_free(params);
	fprintf(file, "%f", chiSquare);
	fclose(file);

	return 0;
}


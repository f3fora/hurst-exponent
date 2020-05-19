#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_statistics.h> 
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_matrix.h>

#define TIME_COLUMN 0
#define X_COLUMN 1
#define profileFileName "profile.dat"
#define DFAFileName "DFA.dat"
#define covFileName "covariance.dat"
#define paramsFileName "params.dat"

     
int print_matrix(FILE *f, const gsl_matrix *m)
{
	/*
	 * https://gist.github.com/jmbr/668067
	 */
        int status, n = 0;
	size_t i, j;
        for (i = 0; i < m->size1; i++) {
                for (j = 0; j < m->size2; j++) {
                        if ((status = fprintf(f, "%g ", gsl_matrix_get(m, i, j))) < 0)
                                return -1;
                        n += status;
                }

                if ((status = fprintf(f, "\n")) < 0)
                        return -1;
                n += status;
        }

        return n;
}


gsl_vector *getSubVector(gsl_vector *vector, size_t offset, size_t len)
{
	/*
	 * create a copy of a vector v' = v[offset: offset+len]]
	 */
	gsl_vector *subVector = gsl_vector_alloc(len);
	size_t i;
	for (i=offset; (i<(len+offset) && (i<vector->size)); i++)
	{
		gsl_vector_set(subVector, i-offset, gsl_vector_get(vector, i));
	}

	return subVector;
}


gsl_vector *determinateProfile(gsl_matrix *data, size_t toProfileColumn)
{
	/*
	 * determinate the profile (step1)
	 */
	gsl_vector *profiledData = gsl_vector_alloc( data->size1 );
	gsl_vector_view column = gsl_matrix_column(data, toProfileColumn);
	double mean = gsl_stats_mean(column.vector.data, column.vector.stride, column.vector.size);

        size_t i;
        double sum = 0;
        for ( i=0; i<column.vector.size; i++ )
        {
                sum += ( gsl_vector_get(&column.vector, i) - mean );
                gsl_vector_set (profiledData, i, sum);
        }

        return profiledData;
}


double getVarianceOfSegment(gsl_vector *partialProfile, gsl_matrix *partialTime)
{
	/*
	 * determinate variance of ONE segment (actually multifit_linear return chisq aka 'The sum of squares of the residuals from the best-fit', Its's exactly what we are looking for!!!)
	 */
	double chisq;
	gsl_matrix  *cov;
	gsl_vector *c;

	size_t n, m;
	n = partialTime->size1;
	m = partialTime->size2;


	c = gsl_vector_alloc(m);
	cov = gsl_matrix_alloc(m, m);
	gsl_multifit_linear_workspace *work  = gsl_multifit_linear_alloc(n, m);  
  	gsl_multifit_linear(partialTime, partialProfile, c, cov, &chisq, work);

    	gsl_multifit_linear_free(work);
	gsl_matrix_free(cov);
	gsl_vector_free(c);
	return chisq;
}


double getLocalFit(gsl_vector *time, gsl_vector* profile, size_t i, size_t s, size_t orderOfDetrend)
{
	/*
	 * fit a local part of the series and return the chisq correctly normalized
	 */
	if (s <= orderOfDetrend)	// if the number of points is minor or less of grade of freedom the chisq is zero
		return 0.0f;

	gsl_vector *partialProfile;
	gsl_vector *ausiliarTime;
	gsl_matrix *partialTime;
	double k;

	partialProfile = getSubVector(profile, i, s);	// No need a copy, just to be sure
	ausiliarTime = getSubVector(time, i, s);	// I need a copy cause it will modify
	partialTime = gsl_matrix_alloc(s, orderOfDetrend+1);
	gsl_matrix_set_all(partialTime, 1.0f);
		
	for (k=1; k<=orderOfDetrend; k++)
	{
		gsl_matrix_set_col(partialTime, k, ausiliarTime);
		gsl_vector_mul(ausiliarTime, ausiliarTime);
	}
	
	gsl_vector_free(ausiliarTime);
	double chisq = getVarianceOfSegment(partialProfile, partialTime) / (s - orderOfDetrend - 1); 
	gsl_vector_free(partialProfile);
	gsl_matrix_free(partialTime);

	return chisq;
}


gsl_vector *getVarianceOfEachSegment(gsl_vector *time, gsl_vector *profile, size_t orderOfDetrend, size_t numberOfSegment, size_t lenghtOfSegment, size_t numberOfWindows)
{
	/*
	 * Split the dataset into small segments and store the chisq of each in a vector for each windows
	 */
	double begin = gsl_vector_get(time, 0);

	double step = 1.0f/numberOfWindows;
	gsl_vector *reducedChiSquared = gsl_vector_alloc((numberOfSegment - 2)  * numberOfWindows + 1);

	size_t i, j;
	size_t n, s, a;
	double chisq;
	
	for  (a=0; a<numberOfWindows;a++)
	{
		j=0;
		for (n=0; n<numberOfSegment; n++)
		{
			for (i=j; ((lenghtOfSegment * ((double)n - a*step)  + begin) <= gsl_vector_get(time, j)) 
					&& (gsl_vector_get(time, j) < (lenghtOfSegment * ((double)n+ 1.0f - a*step) + begin) 
					&& (j < time->size) ); j++);

			if ((a==0) || ((a!=0) && (n>0) && (n<(numberOfSegment-1)))) 
			{
				s = j-i;
				chisq = getLocalFit(time, profile, i, s, orderOfDetrend);
				gsl_vector_set( reducedChiSquared, n , chisq); 
			}
		}

	}
	return reducedChiSquared;	
}


double getVarianceOfSeries(gsl_vector *time, gsl_vector *profile, size_t orderOfDetrend, size_t orderOfFluctuation, size_t numberOfSegment, size_t lenghtOfSegment, size_t numberOfWindows)
{
	/*
	 * determinate the variance of a series for a fixed number of segments
	 */
	gsl_vector *reducedChiSquared = getVarianceOfEachSegment(time, profile, orderOfDetrend, numberOfSegment, lenghtOfSegment, numberOfWindows);
	
	size_t i;
	double variance = 0;
	double aus;
	for (i=0; i<numberOfSegment; i++)
	{
		aus = gsl_vector_get(reducedChiSquared, i);
		variance += pow(aus, (double)orderOfFluctuation / 2.0f);
	}

	variance = pow(variance / (reducedChiSquared->size), 1.0f/orderOfFluctuation);
	gsl_vector_free(reducedChiSquared);
	return variance;
}


gsl_matrix *getDFASpace(gsl_vector *time, gsl_vector *profile, size_t orderOfDetrend, size_t orderOfFluctuation, size_t minNumberOfSegment, size_t maxNumberOfSegment, size_t numberOfWindows)
{
	/*
	 * calculate the space of DFA and s;
	 */
	gsl_matrix *DFASpace = gsl_matrix_alloc( maxNumberOfSegment - minNumberOfSegment, 2 );

	double begin = gsl_vector_get(time, 0);
	double end = gsl_vector_get(time, time->size - 1);
	double lenghtOfSegment;

	size_t i;
	double aus;
	for (i=0; i<maxNumberOfSegment-minNumberOfSegment; i++) 
	{
		lenghtOfSegment = (end - begin) / (double)(i + minNumberOfSegment);
		aus = getVarianceOfSeries(time, profile, orderOfDetrend, orderOfFluctuation, i + minNumberOfSegment, lenghtOfSegment, numberOfWindows);
		gsl_matrix_set(DFASpace, i, 0, lenghtOfSegment);
		gsl_matrix_set(DFASpace, i, 1, aus);
	}

	return DFASpace;
}


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


int main ( int argc, char *argv[] )
{
	// Get terminal arguments
	if ( argc != 9 )
	{
		return 1;
	}

	size_t rows = atoi(argv[1]); // First arg is number of rows of input file
	size_t cols = atoi(argv[2]); // Second arg is number of cols
	size_t detrend = atoi(argv[3]);
	size_t fluctuation = atoi(argv[4]);
	size_t minSeg = atoi(argv[5]);
	size_t maxSeg = atoi(argv[6]);
	size_t windows = atoi(argv[7]);
	char *inputFileName = argv[8];

	// Import initial matrix from the input file
	gsl_matrix *inputData = gsl_matrix_alloc( rows,cols );
    	FILE *file;
	if ( ( file = fopen( inputFileName, "r" ) ) == 0 ) return 1;
	gsl_matrix_fscanf( file, inputData );
	fclose( file );
	
	// Determinate profile of the time series
	gsl_vector *profile = determinateProfile(inputData, X_COLUMN) ;
	gsl_vector *time = gsl_vector_alloc( inputData->size1);
	gsl_matrix_get_col(time, inputData, TIME_COLUMN);
	gsl_matrix_free( inputData );	

	// Export profile data
	if ( ( file = fopen( profileFileName, "w" ) ) == 0 ) return 1;
	gsl_matrix *profileMatrix = gsl_matrix_alloc(time->size, 2);
	gsl_matrix_set_col(profileMatrix, 0, time);
	gsl_matrix_set_col(profileMatrix, 1, profile);
	print_matrix(file, profileMatrix);
	fclose( file );
	gsl_matrix_free( profileMatrix );

	// Calculate DFA
	gsl_matrix *DFAMatrix = getDFASpace(time, profile, detrend, fluctuation, minSeg, maxSeg, windows);

	gsl_vector_free( time );
	gsl_vector_free( profile );

	// Export dfa data
	if ( ( file = fopen( DFAFileName, "w" ) ) == 0 ) return 1;
	print_matrix(file, DFAMatrix);
	fclose( file );
	
	// Fit DFA for calculate Hurst Exponent
	gsl_vector *params = gsl_vector_alloc(2);
	gsl_matrix *covarianceMatrix = gsl_matrix_alloc(2, 2);
	double chiSquare;
	getHurstExponent(DFAMatrix, params, covarianceMatrix, &chiSquare);
	gsl_matrix_free( DFAMatrix );
	
	// Export fit result
	if ( ( file = fopen( covFileName, "w" ) ) == 0 ) return 1;
	print_matrix(file, covarianceMatrix);
	fclose(file);
	gsl_matrix_free(covarianceMatrix);
	
	if ( ( file = fopen( paramsFileName, "w" ) ) == 0 ) return 1;
	gsl_vector_fprintf(file, params, "%f");
	fprintf(file, "%f\n", gsl_vector_get(params, 1) - 1.0f); 
	gsl_vector_free(params);
	fprintf(file, "%f", chiSquare);
	fclose(file);

	return 0;
}

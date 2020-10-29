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
#define PROFILE_NAME "profile.dat"
#define DFA_NAME "DFA.dat"

#define INVALID -42.0f
     

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


int main ( int argc, char *argv[] )
{
	// Get terminal arguments
	if ( argc != 5) return 1;
		
	size_t rows = atoi(argv[1]); 
	size_t cols = atoi(argv[2]); 
	char *inputFileName = argv[3];
	char *outputPath = argv[4];

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

	char *profileFileName = malloc(strlen(outputPath) + strlen(PROFILE_NAME) + 1);
	strcpy(profileFileName, outputPath);
	strcat(profileFileName, PROFILE_NAME);

	// Export profile data
	if ( ( file = fopen( profileFileName, "w" ) ) == 0 ) return 1;
	gsl_matrix *profileMatrix = gsl_matrix_alloc(time->size, 2);
	gsl_matrix_set_col(profileMatrix, 0, time);
	gsl_matrix_set_col(profileMatrix, 1, profile);
	print_matrix(file, profileMatrix);
	fclose( file );
	gsl_matrix_free( profileMatrix );

	free(profileFileName);

	gsl_vector_free( time );
	gsl_vector_free( profile );

	return 0;
}

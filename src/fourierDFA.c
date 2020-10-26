#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>


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

void fourierDFA(gsl_matrix *data, size_t removingTerms)
{
	gsl_fft_real_wavetable * real;
	gsl_fft_halfcomplex_wavetable * hc;
	gsl_fft_real_workspace * work;
	work = gsl_fft_real_workspace_alloc(data->size1);
	real = gsl_fft_real_wavetable_alloc(data->size1);
	gsl_vector_view column = gsl_matrix_column(data, 1);

	gsl_fft_real_transform(column.vector.data, column.vector.stride, column.vector.size, real, work);
	gsl_fft_real_wavetable_free (real);

	size_t i;
	for (i = 0; i < removingTerms; i++)
	{
	       	gsl_vector_set(&column.vector, i, 0);
	}

	hc = gsl_fft_halfcomplex_wavetable_alloc(data->size1);
	gsl_fft_halfcomplex_inverse(column.vector.data, column.vector.stride, column.vector.size, hc, work);
	gsl_fft_halfcomplex_wavetable_free(hc);
	gsl_fft_real_workspace_free(work);
}


int main( int argc, char *argv[] )
{
	// Get terminal arguments
	if ( argc != 6 )
	{
		return 1;
	}
	
	size_t rows = atoi(argv[1]);
	size_t cols = atoi(argv[2]);
	size_t removingTerms = atoi(argv[3]);
	char *inputFileName = argv[4];
	char *outputFileName = argv[5];


	gsl_matrix *inputData = gsl_matrix_alloc( rows, cols );
    	FILE *file;
	if ( ( file = fopen( inputFileName, "r" ) ) == 0 ) return 1;
	gsl_matrix_fscanf( file, inputData );
	fclose( file );

	fourierDFA( inputData, removingTerms);
	
	if ( ( file = fopen( outputFileName, "w" ) ) == 0 ) return 1;
	print_matrix(file, inputData);
	fclose(file);
	gsl_matrix_free( inputData );
	
	return 0;
}


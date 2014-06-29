#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 1


void arguments( int argc, char *argv[] )
{
    int c;

	nblosum = 62;
	fmodel = 0;
	calledByXced = 0;
	scoremtx = 0;
	dorp = NOTSPECIFIED;
	ppenalty = NOTSPECIFIED;
	ppenalty_ex = NOTSPECIFIED;
	poffset = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					fprintf( stderr, "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
					fprintf( stderr, "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = atoi( *++argv );
					fprintf( stderr, "kimuraR = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = atoi( *++argv );
					scoremtx = 1;
					fprintf( stderr, "blosum %d\n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = atoi( *++argv );
					scoremtx = 0;
					fprintf( stderr, "jtt %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					fmodel = 1;
					break;
				case 'r':
					fmodel = -1;
					break;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
                default:
                    fprintf( stderr, "illegal option %c\n", c );
                    argc = 0;
                    break;
            }
		}
		nextoption:
			;
	}
    if( argc == 1 )
    {
        cut = atof( (*argv) );
        argc--;
    }
    if( argc != 0 ) 
    {
        fprintf( stderr, "options: Check source file !\n" );
        exit( 1 );
    }
	if( tbitr == 1 && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : o, m or u\n" );
		exit( 1 );
	}
	if( alg == 'C' && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : C, o\n" );
		exit( 1 );
	}
}

int main( int argc, char *argv[] )
{
	static int  nlen[M];	
	static char name[M][B], **seq;
	int i, j, alloclen, c;
	double **mtx;
	double *self;
	double tmpdouble;
	FILE *fp;

	arguments( argc, argv );

	getnumlen( stdin );
	rewind( stdin );

	if( njob < 2 )
	{
		fprintf( stderr, "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}

	seq = AllocateCharMtx( njob, nlenmax*9+1 );
	mtx = AllocateDoubleMtx( njob, njob );
	self = AllocateDoubleVec( njob );
	alloclen = nlenmax*9;

	readData( stdin, name, nlen, seq );
	constants( njob, seq );




	c = seqcheck( seq );
	if( c )
	{
		fprintf( stderr, "Illeagal character %c\n", c );
		exit( 1 );
	}

	for( i=0; i<njob; i++ ) 
	{
		self[i] = (double)substitution_nid( seq[i], seq[i] );
//		fprintf( stdout, "self[%d] = %f\n", i, self[i] );
	}

	for( i=0; i<njob-1; i++ ) 
		for( j=i+1; j<njob; j++ ) 
		{
			tmpdouble = (double)substitution_score( seq[i], seq[j] );
//			fprintf( stdout, "tmpdouble = %f\n", tmpdouble );
			mtx[i][j] = ( 1.0 - tmpdouble / MIN( self[i], self[j] ) );
			if( mtx[i][j] < 0.95 )
				mtx[i][j] = - log( 1.0 - mtx[i][j] );
			else
				mtx[i][j] = 3.0;
		}
	
#if TEST
	for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ ) 
		fprintf( stdout, "i=%d, j=%d, mtx[][] = %f\n", i, j, mtx[i][j] );
#endif

	fp = fopen( "hat2", "w" );
	WriteHat2( fp, njob, name, mtx );
	fclose( fp );
	exit( 0 );

	return( 0 );
}

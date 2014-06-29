#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0

#define NODIST -9999

static char *whereispairalign;
static char *laraparams;
static char foldalignopt[1000];
static int stdout_align;
static int stdout_dist;
static int store_localhom;
static int store_dist;

typedef struct _jobtable
{
	int i;
	int j;
} Jobtable;

#ifdef enablemultithread
typedef struct _thread_arg
{
	int thread_no;
	int njob;
	Jobtable *jobpospt;
	char **name;
	char **seq;
	LocalHom **localhomtable;
	double **distancemtx;
	double *selfscore;
	char ***bpp;
	int alloclen;
	pthread_mutex_t *mutex_counter;
	pthread_mutex_t *mutex_stdout;
} thread_arg_t;
#endif

static void t2u( char *seq )
{
	while( *seq )
	{
		if     ( *seq == 'A' ) *seq = 'a';
		else if( *seq == 'a' ) *seq = 'a';
		else if( *seq == 'T' ) *seq = 'u';
		else if( *seq == 't' ) *seq = 'u';
		else if( *seq == 'U' ) *seq = 'u';
		else if( *seq == 'u' ) *seq = 'u';
		else if( *seq == 'G' ) *seq = 'g';
		else if( *seq == 'g' ) *seq = 'g';
		else if( *seq == 'C' ) *seq = 'c';
		else if( *seq == 'c' ) *seq = 'c';
		else *seq = 'n';
		seq++;
	}
}

static float recallpairfoldalign( char **mseq1, char **mseq2, int m1, int m2, int *of1pt, int *of2pt, int alloclen )
{
	static FILE *fp = NULL;
	float value;
	char *aln1;
	char *aln2;
	int of1tmp, of2tmp;

	if( fp == NULL )
	{
		fp = fopen( "_foldalignout", "r" );
		if( fp == NULL )
		{
			fprintf( stderr, "Cannot open _foldalignout\n" );
			exit( 1 );
		}
	}

	aln1 = calloc( alloclen, sizeof( char ) );
	aln2 = calloc( alloclen, sizeof( char ) );

	readpairfoldalign( fp, *mseq1, *mseq2, aln1, aln2, m1, m2, &of1tmp, &of2tmp, alloclen );

	if( strstr( foldalignopt, "-global") )
	{
		fprintf( stderr, "Calling G__align11\n" );
		value = G__align11( mseq1, mseq2, alloclen, outgap, outgap );
		*of1pt = 0;
		*of2pt = 0;
	}
	else
	{
		fprintf( stderr, "Calling L__align11\n" );
		value = L__align11( mseq1, mseq2, alloclen, of1pt, of2pt );
	}

//	value = (float)naivepairscore11( *mseq1, *mseq2, penalty ); // nennnotame

	if( aln1[0] == 0 )
	{
		fprintf( stderr, "FOLDALIGN returned no alignment between %d and %d.  Sequence alignment is used instead.\n", m1+1, m2+1 );
	}
	else
	{
		strcpy( *mseq1, aln1 );
		strcpy( *mseq2, aln2 );
		*of1pt = of1tmp;
		*of2pt = of2tmp;
	}

//	value = naivepairscore11( *mseq1, *mseq2, penalty ); // v6.511 ha kore wo tsukau, global nomi dakara.

//	fclose( fp ); // saigo dake yatta houga yoi.

//	fprintf( stderr, "*mseq1 = %s\n", *mseq1 );
//	fprintf( stderr, "*mseq2 = %s\n", *mseq2 );


	free( aln1 );
	free( aln2 );

	return( value );
}

static void callfoldalign( int nseq, char **mseq )
{
	FILE *fp;
	int i;
	int res;
	static char com[10000];

	for( i=0; i<nseq; i++ )
		t2u( mseq[i] );

	fp = fopen( "_foldalignin", "w" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open _foldalignin\n" );
		exit( 1 );
	}
	for( i=0; i<nseq; i++ )
	{
		fprintf( fp, ">%d\n", i+1 );
		fprintf( fp, "%s\n", mseq[i] );
	}
	fclose( fp );

	sprintf( com, "env PATH=%s  foldalign210 %s _foldalignin > _foldalignout ", whereispairalign, foldalignopt );
	res = system( com );
	if( res )
	{
		fprintf( stderr, "Error in foldalign\n" );
		exit( 1 );
	}

}

static void calllara( int nseq, char **mseq, char *laraarg )
{
	FILE *fp;
	int i;
	int res;
	static char com[10000];

//	for( i=0; i<nseq; i++ )

	fp = fopen( "_larain", "w" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open _larain\n" );
		exit( 1 );
	}
	for( i=0; i<nseq; i++ )
	{
		fprintf( fp, ">%d\n", i+1 );
		fprintf( fp, "%s\n", mseq[i] );
	}
	fclose( fp );


//	fprintf( stderr, "calling LaRA\n" );
	sprintf( com, "env PATH=%s:/bin:/usr/bin mafft_lara -i _larain -w _laraout -o _lara.params %s", whereispairalign, laraarg );
	res = system( com );
	if( res )
	{
		fprintf( stderr, "Error in lara\n" );
		exit( 1 );
	}
}

static float recalllara( char **mseq1, char **mseq2, int alloclen )
{
	static FILE *fp = NULL;
	static char *ungap1;
	static char *ungap2;
	static char *ori1;
	static char *ori2;
//	int res;
	static char com[10000];
	float value;


	if( fp == NULL )
	{
		fp = fopen( "_laraout", "r" );
		if( fp == NULL )
		{
			fprintf( stderr, "Cannot open _laraout\n" );
			exit( 1 );
		}
		ungap1 = AllocateCharVec( alloclen );
		ungap2 = AllocateCharVec( alloclen );
		ori1 = AllocateCharVec( alloclen );
		ori2 = AllocateCharVec( alloclen );
	}


	strcpy( ori1, *mseq1 );
	strcpy( ori2, *mseq2 );

	fgets( com, 999, fp );
	myfgets( com, 9999, fp );
	strcpy( *mseq1, com );
	myfgets( com, 9999, fp );
	strcpy( *mseq2, com );

	gappick0( ungap1, *mseq1 );
	gappick0( ungap2, *mseq2 );
	t2u( ungap1 );
	t2u( ungap2 );
	t2u( ori1 );
	t2u( ori2 );

	if( strcmp( ungap1, ori1 ) || strcmp( ungap2, ori2 ) )
	{
		fprintf( stderr, "SEQUENCE CHANGED!!\n" );
		fprintf( stderr, "*mseq1  = %s\n", *mseq1 );
		fprintf( stderr, "ungap1  = %s\n", ungap1 );
		fprintf( stderr, "ori1    = %s\n", ori1 );
		fprintf( stderr, "*mseq2  = %s\n", *mseq2 );
		fprintf( stderr, "ungap2  = %s\n", ungap2 );
		fprintf( stderr, "ori2    = %s\n", ori2 );
		exit( 1 );
	}

	value = (float)naivepairscore11( *mseq1, *mseq2, penalty );

//	fclose( fp ); // saigo dake yatta houga yoi.

	return( value );
}


static float callmxscarna_giving_bpp( char **mseq1, char **mseq2, char **bpp1, char **bpp2, int alloclen, int i, int j )
{
	FILE *fp;
	int res;
	char *com;
	float value;
	char *dirname;


	dirname = calloc( 100, sizeof( char ) );
	com = calloc( 1000, sizeof( char ) );
	sprintf( dirname, "_%d-%d", i, j );
	sprintf( com, "rm -rf %s", dirname );
	system( com );
	sprintf( com, "mkdir %s", dirname );
	system( com );


	sprintf( com, "%s/_bpporg", dirname );
	fp = fopen( com, "w" );
	if( !fp )
	{
		fprintf( stderr, "Cannot write to %s/_bpporg\n", dirname );
		exit( 1 );
	}
	fprintf( fp, ">a\n" );
	while( *bpp1 )
		fprintf( fp, "%s", *bpp1++ );

	fprintf( fp, ">b\n" );
	while( *bpp2 )
		fprintf( fp, "%s", *bpp2++ );
	fclose( fp );

	sprintf( com, "tr -d '\\r' < %s/_bpporg > %s/_bpp", dirname, dirname );
	system( com ); // for cygwin, wakaran

	t2u( *mseq1 );
	t2u( *mseq2 );

	sprintf( com, "%s/_mxscarnainorg", dirname );
	fp = fopen( com, "w" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open %s/_mxscarnainorg\n", dirname );
		exit( 1 );
	}
	fprintf( fp, ">1\n" );
//	fprintf( fp, "%s\n", *mseq1 );
	write1seq( fp, *mseq1 );
	fprintf( fp, ">2\n" );
//	fprintf( fp, "%s\n", *mseq2 );
	write1seq( fp, *mseq2 );
	fclose( fp );

	sprintf( com, "tr -d '\\r' < %s/_mxscarnainorg > %s/_mxscarnain", dirname, dirname );
	system( com ); // for cygwin, wakaran

#if 0
	sprintf( com, "cd %s; %s/mxscarnamod -readbpp _mxscarnain > _mxscarnaout 2>_dum", dirname, whereispairalign );
#else
	sprintf( com, "_mxscarnash%s", dirname );
	fp = fopen( com, "w" );
	fprintf( fp, "cd %s\n", dirname );
	fprintf( fp, "%s/mxscarnamod -readbpp _mxscarnain > _mxscarnaout 2>_dum\n", whereispairalign );
	fprintf( fp, "exit $tatus\n" );
	fclose( fp );

	sprintf( com, "tr -d '\\r' < _mxscarnash%s > _mxscarnash%s.unix", dirname, dirname );
	system( com ); // for cygwin, wakaran

	sprintf( com, "sh _mxscarnash%s.unix 2>_dum%s", dirname, dirname );
#endif
	res = system( com );
	if( res )
	{
		fprintf( stderr, "Error in mxscarna\n" );
		exit( 1 );
	}

	sprintf( com, "%s/_mxscarnaout", dirname );

	fp = fopen( com, "r" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open %s/_mxscarnaout\n", dirname );
		exit( 1 );
	}

	fgets( com, 999, fp );
	load1SeqWithoutName_new( fp, *mseq1 );
	fgets( com, 999, fp );
	load1SeqWithoutName_new( fp, *mseq2 );

	fclose( fp );

//	fprintf( stderr, "*mseq1 = %s\n", *mseq1 );
//	fprintf( stderr, "*mseq2 = %s\n", *mseq2 );

	value = (float)naivepairscore11( *mseq1, *mseq2, penalty );

#if 0
	sprintf( com, "rm -rf %s > /dev/null 2>&1", dirname );
	if( system( com ) )
	{
		fprintf( stderr, "retrying to rmdir\n" );
		usleep( 2000 );
		system( com );
	}
#endif

	free( dirname );
	free( com );


	return( value );
}

#if 0
static float callmxscarna_slow( char **mseq1, char **mseq2, int alloclen )
{
	FILE *fp;
	int res;
	static char com[10000];
	float value;


	t2u( *mseq1 );
	t2u( *mseq2 );
	fp = fopen( "_mxscarnain", "w" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open _mxscarnain\n" );
		exit( 1 );
	}
	fprintf( fp, ">1\n" );
	fprintf( fp, "%s\n", *mseq1 );
	fprintf( fp, ">2\n" );
	fprintf( fp, "%s\n", *mseq2 );
	fclose( fp );

	sprintf( com, "env PATH=%s mxscarnamod _mxscarnain > _mxscarnaout 2>_dum", whereispairalign );
	res = system( com );
	if( res )
	{
		fprintf( stderr, "Error in mxscarna\n" );
		exit( 1 );
	}

	fp = fopen( "_mxscarnaout", "r" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open _mxscarnaout\n" );
		exit( 1 );
	}

	fgets( com, 999, fp );
	load1SeqWithoutName_new( fp, *mseq1 );
	fgets( com, 999, fp );
	load1SeqWithoutName_new( fp, *mseq2 );

	fclose( fp );

//	fprintf( stderr, "*mseq1 = %s\n", *mseq1 );
//	fprintf( stderr, "*mseq2 = %s\n", *mseq2 );

	value = (float)naivepairscore11( *mseq1, *mseq2, penalty );

	return( value );
}
#endif

static void readhat4( FILE *fp, char ***bpp )
{
	char oneline[1000];
	int bppsize;
	int onechar;
//	double prob;
//	int posi, posj;

	bppsize = 0;
//	fprintf( stderr, "reading hat4\n" );
	onechar = getc(fp);
//	fprintf( stderr, "onechar = %c\n", onechar );
	if( onechar != '>' )
	{
		fprintf( stderr, "Format error\n" );
		exit( 1 );
	}
	ungetc( onechar, fp );
	fgets( oneline, 999, fp );
	while( 1 )
	{
		onechar = getc(fp);
		ungetc( onechar, fp );
		if( onechar == '>' || onechar == EOF )
		{
//			fprintf( stderr, "Next\n" );
			*bpp = realloc( *bpp, (bppsize+2) * sizeof( char * ) );
			(*bpp)[bppsize] = NULL;
			break;
		}
		fgets( oneline, 999, fp );
//		fprintf( stderr, "oneline=%s\n", oneline );
//		sscanf( oneline, "%d %d %f", &posi, &posj, &prob );
//		fprintf( stderr, "%d %d -> %f\n", posi, posj, prob );
		*bpp = realloc( *bpp, (bppsize+2) * sizeof( char * ) );
		(*bpp)[bppsize] = calloc( 100, sizeof( char ) );
		strcpy( (*bpp)[bppsize], oneline );
		bppsize++;
	}
}

static void preparebpp( int nseq, char ***bpp )
{
	FILE *fp;
	int i;

	fp = fopen( "hat4", "r" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open hat4\n" );
		exit( 1 );
	}
	for( i=0; i<nseq; i++ )
		readhat4( fp, bpp+i );
	fclose( fp );
}

void arguments( int argc, char *argv[] )
{
    int c;

	nthread = 1;
	foldalignopt[0] = 0;
	laraparams = NULL;
	inputfile = NULL;
	fftkeika = 0;
	pslocal = -1000.0;
	constraint = 0;
	nblosum = 62;
	fmodel = 0;
	calledByXced = 0;
	devide = 0;
	use_fft = 0;
	fftscore = 1;
	fftRepeatStop = 0;
	fftNoAnchStop = 0;
    weight = 3;
    utree = 1;
	tbutree = 1;
    refine = 0;
    check = 1;
    cut = 0.0;
    disp = 0;
    outgap = 1;
    alg = 'A';
    mix = 0;
	tbitr = 0;
	scmtd = 5;
	tbweight = 0;
	tbrweight = 3;
	checkC = 0;
	treemethod = 'x';
	contin = 0;
	scoremtx = 1;
	kobetsubunkatsu = 0;
	divpairscore = 0;
	stdout_align = 0;
	stdout_dist = 0;
	store_dist = 1;
	store_localhom = 1;
	dorp = NOTSPECIFIED;
	ppenalty = NOTSPECIFIED;
	ppenalty_OP = NOTSPECIFIED;
	ppenalty_ex = NOTSPECIFIED;
	ppenalty_EX = NOTSPECIFIED;
	poffset = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	RNAppenalty = NOTSPECIFIED;
	RNApthr = NOTSPECIFIED;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( ( c = *++argv[0] ) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'O':
					ppenalty_OP = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'E':
					ppenalty_EX = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = atoi( *++argv );
//					fprintf( stderr, "kimuraR = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = atoi( *++argv );
					scoremtx = 1;
//					fprintf( stderr, "blosum %d\n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = atoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
					fprintf( stderr, "jtt %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = atoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
					fprintf( stderr, "TM %d\n", pamN );
					--argc;
					goto nextoption;
				case 'l':
					ppslocal = (int)( atof( *++argv ) * 1000 + 0.5 );
					pslocal = (int)( 600.0 / 1000.0 * ppslocal + 0.5);
//					fprintf( stderr, "ppslocal = %d\n", ppslocal );
//					fprintf( stderr, "pslocal = %d\n", pslocal );
					--argc;
					goto nextoption;
				case 'd':
					whereispairalign = *++argv;
					fprintf( stderr, "whereispairalign = %s\n", whereispairalign );
					--argc; 
					goto nextoption;
				case 'p':
					laraparams = *++argv;
					fprintf( stderr, "laraparams = %s\n", laraparams );
					--argc; 
					goto nextoption;
				case 'C':
					nthread = atoi( *++argv );
					fprintf( stderr, "nthread = %d\n", nthread );
					--argc; 
					goto nextoption;
				case 'c':
					stdout_dist = 1;
					break;
				case 'n':
					stdout_align = 1;
					break;
				case 'x':
					store_localhom = 0;
					store_dist = 0;
					break;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'r':
					fmodel = -1;
					break;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
				case 'e':
					fftscore = 0;
					break;
#if 0
				case 'O':
					fftNoAnchStop = 1;
					break;
#endif
				case 'Q':
					calledByXced = 1;
					break;
#if 0
				case 'x':
					disp = 1;
					break;
				case 'a':
					alg = 'a';
					break;
#endif
				case 'S':
					alg = 'S';
					break;
				case 't':
					alg = 't';
					store_localhom = 0;
					break;
				case 'L':
					alg = 'L';
					break;
				case 's':
					alg = 's';
					break;
				case 'B':
					alg = 'B';
					break;
				case 'T':
					alg = 'T';
					break;
				case 'H':
					alg = 'H';
					break;
				case 'M':
					alg = 'M';
					break;
				case 'R':
					alg = 'R';
					break;
				case 'N':
					alg = 'N';
					break;
				case 'K':
					alg = 'K';
					break;
				case 'A':
					alg = 'A';
					break;
				case 'V':
					alg = 'V';
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'v':
					tbrweight = 3;
					break;
				case 'y':
					divpairscore = 1;
					break;
/* Modified 01/08/27, default: user tree */
				case 'J':
					tbutree = 0;
					break;
/* modification end. */
				case 'o':
//					foldalignopt = *++argv;
					strcat( foldalignopt, " " );
					strcat( foldalignopt, *++argv );
					fprintf( stderr, "foldalignopt = %s\n", foldalignopt );
					--argc; 
					goto nextoption;
				case 'z':
					fftThreshold = atoi( *++argv );
					--argc; 
					goto nextoption;
				case 'w':
					fftWinSize = atoi( *++argv );
					--argc;
					goto nextoption;
				case 'Z':
					checkC = 1;
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
}

int countamino( char *s, int end )
{
	int val = 0;
	while( end-- )
		if( *s++ != '-' ) val++;
	return( val );
}

#if enablemultithread
static void *athread( void *arg )
{
	thread_arg_t *targ = (thread_arg_t *)arg;
	int i, j;
	int clus1, clus2;
	int off1, off2;
	int intdum;
	double bunbo;
	float pscore = 0.0; // by D.Mathog
	double *effarr;
	double *effarr1;
	double *effarr2;
	char **pair;
	char *indication1, *indication2;
	char **mseq1, **mseq2;
	char **aseq;

// thread_arg
	int thread_no = targ->thread_no;
	int njob = targ->njob;
	Jobtable *jobpospt = targ->jobpospt;
	char **name = targ->name;
	char **seq = targ->seq;
	LocalHom **localhomtable = targ->localhomtable;
	double **distancemtx = targ->distancemtx;
	double *selfscore = targ->selfscore;
	char ***bpp = targ->bpp;
	int alloclen = targ->alloclen;

//	fprintf( stderr, "thread %d start!\n", thread_no );

	effarr = AllocateDoubleVec( njob );
	for( i=0; i<njob; i++ ) effarr[i] = 1.0;
	effarr1 = AllocateDoubleVec( njob );
	effarr2 = AllocateDoubleVec( njob );
	indication1 = AllocateCharVec( 150 );
	indication2 = AllocateCharVec( 150 );
	pair = AllocateCharMtx( njob, njob );
	for( i=0; i<njob; i++ ) for( j=0; j<njob; j++ ) pair[i][j] = 0;
	for( i=0; i<njob; i++ ) pair[i][i] = 1;
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );
	aseq = AllocateCharMtx( njob, alloclen+10 );

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex_counter );
		j = jobpospt->j;
		i = jobpospt->i;
		j++;
		if( j == njob )
		{
			i++;
			j = i + 1;
			if( i == njob-1 )
			{
//				fprintf( stderr, "thread %d end!\n", thread_no );
				pthread_mutex_unlock( targ->mutex_counter );

				if( commonIP ) FreeIntMtx( commonIP );
				commonIP = NULL;
				if( commonJP ) FreeIntMtx( commonJP );
				commonJP = NULL;
				Falign( NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
				G__align11_noalign( NULL, 0, 0, NULL, NULL, 0 );
				L__align11( NULL, NULL, 0, NULL, NULL );
				genL__align11( NULL, NULL, 0, NULL, NULL );
				free( effarr );
				free( effarr1 );
				free( effarr2 );
				free( indication1 );
				free( indication2 );
				FreeCharMtx( pair  );
				free( mseq1 );
				free( mseq2 );
				FreeCharMtx( aseq  );
				return( NULL );
			}
		}
		jobpospt->j = j;
		jobpospt->i = i;
		pthread_mutex_unlock( targ->mutex_counter );


		if( j == i+1 || j % 100 == 0 ) 
		{
			fprintf( stderr, "% 5d / %d (by thread %3d) \r", i, njob, thread_no );
//			fprintf( stderr, "% 5d - %5d / %d (thread %d)\n", i, j, njob, thread_no );
		}


		if( strlen( seq[i] ) == 0 || strlen( seq[j] ) == 0 )
		{
			if( store_dist ) distancemtx[i][j] = 2.0;
			if( stdout_dist) 
			{
				pthread_mutex_lock( targ->mutex_stdout );
				fprintf( stdout, "%d %d d=%.3f\n", i+1, j+1, 2.0 );
				pthread_mutex_unlock( targ->mutex_stdout );
			}
			continue;
		}

		strcpy( aseq[i], seq[i] );
		strcpy( aseq[j], seq[j] );
		clus1 = conjuctionfortbfast( pair, i, aseq, mseq1, effarr1, effarr, indication1 );
		clus2 = conjuctionfortbfast( pair, j, aseq, mseq2, effarr2, effarr, indication2 );
//		fprintf( stderr, "mseq1 = %s\n", mseq1[0] );
//		fprintf( stderr, "mseq2 = %s\n", mseq2[0] );
	
#if 0
		fprintf( stderr, "group1 = %.66s", indication1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "group2 = %.66s", indication2 );
		fprintf( stderr, "\n" );
#endif
//		for( l=0; l<clus1; l++ ) fprintf( stderr, "## STEP-eff for mseq1-%d %f\n", l, effarr1[l] );

		if( use_fft )
		{
			pscore = Falign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, &intdum, NULL, 0, NULL );
//			fprintf( stderr, "pscore (fft) = %f\n", pscore );
			off1 = off2 = 0;
		}
		else
		{
			switch( alg )
			{
				case( 'L' ):
					pscore = G__align11_noalign( amino_dis, penalty, penalty_ex, mseq1, mseq2, alloclen );
					L__align11( mseq1, mseq2, alloclen, &off1, &off2 );
					break;
				case( 'A' ):
					pscore = G__align11( mseq1, mseq2, alloclen, outgap, outgap );
					off1 = off2 = 0;
					break;
				case( 'N' ):
					pscore = G__align11_noalign( amino_dis, penalty, penalty_ex, mseq1, mseq2, alloclen );
					genL__align11( mseq1, mseq2, alloclen, &off1, &off2 );
					break;
				case( 't' ):
					pscore = G__align11_noalign( amino_dis, penalty, penalty_ex, mseq1, mseq2, alloclen );
					off1 = off2 = 0;
					break;
				case( 's' ):
					pscore = callmxscarna_giving_bpp( mseq1, mseq2, bpp[i], bpp[j], alloclen, i, j );
					off1 = off2 = 0;
					break;
#if 0 
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
					off1 = off2 = 0;
					break;
				case( 'K' ):
					pscore = genG__align11( mseq1, mseq2, alloclen );
					off1 = off2 = 0;
					break;
				case( 'H' ):
					pscore = recallpairfoldalign( mseq1, mseq2, i, j, &off1, &off2, alloclen );
					break;
				case( 'B' ):
				case( 'T' ):
					pscore = recalllara( mseq1, mseq2, alloclen );
					off1 = off2 = 0;
					break;
				case( 'M' ):
					pscore = MSalign11( mseq1, mseq2, alloclen );
					break;
#endif
				default:
					ErrorExit( "\n\nERROR IN SOURCE FILE\n\n" );
			}
		}

		if( alg == 't' || ( mseq1[0][0] != 0 && mseq2[0][0] != 0  ) ) // 't' no jouken ha iranai to omou. if( ( mseq1[0][0] != 0 && mseq2[0][0] != 0  ) )
		{
#if SCOREOUT
			fprintf( stderr, "score = %10.2f (%d,%d)\n", pscore, i, j );
#endif
			if( !store_localhom )
				;
			else if( alg == 'H' )
				putlocalhom_ext( mseq1[0], mseq2[0], localhomtable[i]+j, off1, off2, (int)pscore, strlen( mseq1[0] ) );
			else if( alg != 'S' && alg != 'V' )
			{
				putlocalhom2( mseq1[0], mseq2[0], localhomtable[i]+j, off1, off2, (int)pscore, strlen( mseq1[0] ) );
			}

			if( (bunbo=MIN( selfscore[i], selfscore[j] )) == 0.0 || bunbo < pscore )
				pscore = 2.0;
			else
				pscore = ( 1.0 - pscore / bunbo ) * 2.0;
		}
		else
		{
			pscore = 2.0;
		}

#if 1 // mutex
		if( stdout_align )
		{
			pthread_mutex_lock( targ->mutex_stdout );
			if( alg != 't' )
			{
				fprintf( stdout, "sequence %d - sequence %d, pairwise distance = %10.5f\n", i+1, j+1, pscore );
				fprintf( stdout, ">%s\n", name[i] );
				write1seq( stdout, mseq1[0] );
				fprintf( stdout, ">%s\n", name[j] );
				write1seq( stdout, mseq2[0] );
				fprintf( stdout, "\n" );
			}
			pthread_mutex_unlock( targ->mutex_stdout );
		}
		if( stdout_dist )
		{
			pthread_mutex_lock( targ->mutex_stdout );
			if( j == i+1 ) fprintf( stdout, "%d %d d=%.3f\n", i+1, i+1, 0.0 );
			fprintf( stdout, "%d %d d=%.3f\n", i+1, j+1, pscore );
			pthread_mutex_unlock( targ->mutex_stdout );
		}
#endif // mutex
		if( store_dist) distancemtx[i][j] = pscore;
	}
}
#endif

static void pairalign( char **name, int nlen[M], char **seq, char **aseq, char **mseq1, char **mseq2, double *effarr, int alloclen )
{
	int i, j, ilim;
	int clus1, clus2;
	int off1, off2;
	float pscore = 0.0; // by D.Mathog
	static char *indication1, *indication2;
	FILE *hat2p, *hat3p;
	double **distancemtx;
	double *selfscore;
	double *effarr1;
	double *effarr2;
	char *pt;
	char *hat2file = "hat2";
	LocalHom **localhomtable = NULL, *tmpptr;
	static char **pair;
	int intdum;
	double bunbo;
	char ***bpp = NULL; // mxscarna no toki dake

	if( store_localhom )
	{
		localhomtable = (LocalHom **)calloc( njob, sizeof( LocalHom *) );
		for( i=0; i<njob; i++)
		{
			localhomtable[i] = (LocalHom *)calloc( njob, sizeof( LocalHom ) );
			for( j=0; j<njob; j++)
			{
				localhomtable[i][j].start1 = -1;
				localhomtable[i][j].end1 = -1;
				localhomtable[i][j].start2 = -1; 
				localhomtable[i][j].end2 = -1; 
				localhomtable[i][j].opt = -1.0;
				localhomtable[i][j].next = NULL;
				localhomtable[i][j].nokori = 0;
			}
		}
	}

	if( store_dist ) distancemtx = AllocateDoubleMtx( njob, njob );
	else distancemtx = NULL;
	selfscore = AllocateDoubleVec( njob );
	effarr1 = AllocateDoubleVec( njob );
	effarr2 = AllocateDoubleVec( njob );
	indication1 = AllocateCharVec( 150 );
	indication2 = AllocateCharVec( 150 );
#if 0
#else
	pair = AllocateCharMtx( njob, njob );
#endif

#if 0
	fprintf( stderr, "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		fprintf( stderr, "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif


//	writePre( njob, name, nlen, aseq, 0 );

	for( i=0; i<njob; i++ ) for( j=0; j<njob; j++ ) pair[i][j] = 0;
	for( i=0; i<njob; i++ ) pair[i][i] = 1;

	if( alg == 'H' )
	{
		fprintf( stderr, "Calling FOLDALIGN with option '%s'\n", foldalignopt );
		callfoldalign( njob, seq );
		fprintf( stderr, "done.\n" );
	}
	if( alg == 'B' )
	{
		fprintf( stderr, "Running LARA (Bauer et al. http://www.planet-lisa.net/)\n" );
		calllara( njob, seq, "" );
		fprintf( stderr, "done.\n" );
	}
	if( alg == 'T' )
	{
		fprintf( stderr, "Running SLARA (Bauer et al. http://www.planet-lisa.net/)\n" );
		calllara( njob, seq, "-s" );
		fprintf( stderr, "done.\n" );
	}
	if( alg == 's' )
	{
		fprintf( stderr, "Preparing bpp\n" );
//		bpp = AllocateCharCub( njob, nlenmax, 0 );
		bpp = calloc( njob, sizeof( char ** ) );
		preparebpp( njob, bpp );
		fprintf( stderr, "done.\n" );
		fprintf( stderr, "Running MXSCARNA (Tabei et al. http://www.ncrna.org/software/mxscarna)\n" );
	}

	for( i=0; i<njob; i++ )
	{
		pscore = 0.0;
		for( pt=seq[i]; *pt; pt++ )
			pscore += amino_dis[(int)*pt][(int)*pt];
		selfscore[i] = pscore;

	}

#if enablemultithread
	if( nthread > 0 )
	{
		Jobtable jobpos;
		pthread_t *handle;
		pthread_mutex_t mutex_counter;
		pthread_mutex_t mutex_stdout;
		thread_arg_t *targ;

		jobpos.i = 0;
		jobpos.j = 0;

		targ = calloc( nthread, sizeof( thread_arg_t ) );
		handle = calloc( nthread, sizeof( pthread_t ) );
		pthread_mutex_init( &mutex_counter, NULL );
		pthread_mutex_init( &mutex_stdout, NULL );

		for( i=0; i<nthread; i++ )
		{
			targ[i].thread_no = i;
			targ[i].njob = njob;
			targ[i].jobpospt = &jobpos;
			targ[i].name = name;
			targ[i].seq = seq;
			targ[i].localhomtable = localhomtable;
			targ[i].distancemtx = distancemtx;
			targ[i].selfscore = selfscore;
			targ[i].bpp = bpp; 
			targ[i].alloclen = alloclen;
			targ[i].mutex_counter = &mutex_counter;
			targ[i].mutex_stdout = &mutex_stdout;

//			athread( (void *)targ );
			pthread_create( handle+i, NULL, athread, (void *)(targ+i) );
//			pthread_create( handle+i, NULL, bthread, (void *)(targ+i) );
		}


		for( i=0; i<nthread; i++ )
		{
			pthread_join( handle[i], NULL );
		}
		pthread_mutex_destroy( &mutex_counter );
		pthread_mutex_destroy( &mutex_stdout );
		free( handle );
		free( targ );
	}
	else
#endif
	{
		ilim = njob - 1;
		for( i=0; i<ilim; i++ ) 
		{
			if( stdout_dist) fprintf( stdout, "%d %d d=%.3f\n", i+1, i+1, 0.0 );
			fprintf( stderr, "% 5d / %d\r", i, njob );
			fflush( stderr );
			for( j=i+1; j<njob; j++ )
			{
	
				if( strlen( seq[i] ) == 0 || strlen( seq[j] ) == 0 )
				{
					if( store_dist ) distancemtx[i][j] = 2.0;
					if( stdout_dist) fprintf( stdout, "%d %d d=%.3f\n", i+1, j+1, 2.0 );
					continue;
				}
	
				strcpy( aseq[i], seq[i] );
				strcpy( aseq[j], seq[j] );
				clus1 = conjuctionfortbfast( pair, i, aseq, mseq1, effarr1, effarr, indication1 );
				clus2 = conjuctionfortbfast( pair, j, aseq, mseq2, effarr2, effarr, indication2 );
	//			fprintf( stderr, "mseq1 = %s\n", mseq1[0] );
	//			fprintf( stderr, "mseq2 = %s\n", mseq2[0] );
		
#if 0
				fprintf( stderr, "group1 = %.66s", indication1 );
				fprintf( stderr, "\n" );
				fprintf( stderr, "group2 = %.66s", indication2 );
				fprintf( stderr, "\n" );
#endif
	//			for( l=0; l<clus1; l++ ) fprintf( stderr, "## STEP-eff for mseq1-%d %f\n", l, effarr1[l] );
	
				if( use_fft )
				{
					pscore = Falign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, &intdum, NULL, 0, NULL );
//					fprintf( stderr, "pscore (fft) = %f\n", pscore );
					off1 = off2 = 0;
				}
				else
				{
					switch( alg )
					{
						case( 'a' ):
							pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
							off1 = off2 = 0;
							break;
						case( 't' ):
							pscore = G__align11_noalign( amino_dis, penalty, penalty_ex, mseq1, mseq2, alloclen );
							off1 = off2 = 0;
							break;
						case( 'A' ):
							pscore = G__align11( mseq1, mseq2, alloclen, outgap, outgap );
							off1 = off2 = 0;
							break;
						case( 'N' ):
							pscore = G__align11_noalign( amino_dis, penalty, penalty_ex, mseq1, mseq2, alloclen );
							genL__align11( mseq1, mseq2, alloclen, &off1, &off2 );
							break;
						case( 'K' ):
							pscore = genG__align11( mseq1, mseq2, alloclen );
							off1 = off2 = 0;
							break;
						case( 'L' ):
							pscore = G__align11_noalign( amino_dis, penalty, penalty_ex, mseq1, mseq2, alloclen );
							L__align11( mseq1, mseq2, alloclen, &off1, &off2 );
							break;
						case( 'H' ):
							pscore = recallpairfoldalign( mseq1, mseq2, i, j, &off1, &off2, alloclen );
							break;
						case( 'B' ):
						case( 'T' ):
							pscore = recalllara( mseq1, mseq2, alloclen );
							off1 = off2 = 0;
							break;
						case( 's' ):
							pscore = callmxscarna_giving_bpp( mseq1, mseq2, bpp[i], bpp[j], alloclen, i, j );
							off1 = off2 = 0;
							break;
						case( 'M' ):
							pscore = MSalign11( mseq1, mseq2, alloclen );
							break;
						default:
							ErrorExit( "ERROR IN SOURCE FILE" );
					}
				}
	
				if( alg == 't' || ( mseq1[0][0] != 0 && mseq2[0][0] != 0  ) ) // 't' no jouken ha iranai to omou. if( ( mseq1[0][0] != 0 && mseq2[0][0] != 0  ) )
				{
#if SCOREOUT
					fprintf( stderr, "score = %10.2f (%d,%d)\n", pscore, i, j );
#endif
					if( !store_localhom )
						;
					else if( alg == 'H' )
						putlocalhom_ext( mseq1[0], mseq2[0], localhomtable[i]+j, off1, off2, (int)pscore, strlen( mseq1[0] ) );
					else if( alg != 'S' && alg != 'V' )
						putlocalhom2( mseq1[0], mseq2[0], localhomtable[i]+j, off1, off2, (int)pscore, strlen( mseq1[0] ) );
	
	
					if( (bunbo=MIN( selfscore[i], selfscore[j] )) == 0.0 || bunbo < pscore )
						pscore = 2.0;
					else
						pscore = ( 1.0 - pscore / bunbo ) * 2.0;
				}
				else
				{
					pscore = 2.0;
				}
	
				if( stdout_align )
				{
					if( alg != 't' )
					{
						fprintf( stdout, "sequence %d - sequence %d, pairwise distance = %10.5f\n", i+1, j+1, pscore );
						fprintf( stdout, ">%s\n", name[i] );
						write1seq( stdout, mseq1[0] );
						fprintf( stdout, ">%s\n", name[j] );
						write1seq( stdout, mseq2[0] );
						fprintf( stdout, "\n" );
					}
				}
				if( stdout_dist ) fprintf( stdout, "%d %d d=%.3f\n", i+1, j+1, pscore );
				if( store_dist) distancemtx[i][j] = pscore;
			}
		}
	}


	if( store_dist )
	{
		hat2p = fopen( hat2file, "w" );
		if( !hat2p ) ErrorExit( "Cannot open hat2." );
		WriteHat2_pointer( hat2p, njob, name, distancemtx );
		fclose( hat2p );
	}

	hat3p = fopen( "hat3", "w" );
	if( !hat3p ) ErrorExit( "Cannot open hat3." );
	if( store_localhom )
	{
		fprintf( stderr, "\n\n##### writing hat3\n" );
		ilim = njob-1;	
		for( i=0; i<ilim; i++ ) 
		{
			for( j=i+1; j<njob; j++ )
			{
				for( tmpptr=localhomtable[i]+j; tmpptr; tmpptr=tmpptr->next )
				{
					if( tmpptr->opt == -1.0 ) continue;
// tmptmptmptmptmp
//					if( alg == 'B' || alg == 'T' )
//						fprintf( hat3p, "%d %d %d %7.5f %d %d %d %d %p\n", i, j, tmpptr->overlapaa, 1.0, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, (void *)tmpptr->next ); 
//					else
						fprintf( hat3p, "%d %d %d %7.5f %d %d %d %d h\n", i, j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2 ); 
				}
			}
		}
#if DEBUG
		fprintf( stderr, "calling FreeLocalHomTable\n" );
#endif
		FreeLocalHomTable( localhomtable, njob );
#if DEBUG
		fprintf( stderr, "done. FreeLocalHomTable\n" );
#endif
	}
	fclose( hat3p );

	if( alg == 's' )
	{
		char **ptpt;
		for( i=0; i<njob; i++ )
		{
			ptpt = bpp[i];
			while( 1 )
			{
				if( *ptpt ) free( *ptpt );
				else break;
				ptpt++;
			}
			free( bpp[i] );
		}
		free( bpp );
	}
	free( selfscore );
	free( effarr1 );
	free( effarr2 );
	free( indication1 );
	free( indication2 );
	if( store_dist ) FreeDoubleMtx( distancemtx );
}

static void WriteOptions( FILE *fp )
{

	if( dorp == 'd' ) fprintf( fp, "DNA\n" );
	else
	{
		if     ( scoremtx ==  0 ) fprintf( fp, "JTT %dPAM\n", pamN );
		else if( scoremtx ==  1 ) fprintf( fp, "BLOSUM %d\n", nblosum );
		else if( scoremtx ==  2 ) fprintf( fp, "M-Y\n" );
	}
    fprintf( stderr, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );
    if( use_fft ) fprintf( fp, "FFT on\n" );

	fprintf( fp, "tree-base method\n" );
	if( tbrweight == 0 ) fprintf( fp, "unweighted\n" );
	else if( tbrweight == 3 ) fprintf( fp, "clustalw-like weighting\n" );
	if( tbitr || tbweight ) 
	{
		fprintf( fp, "iterate at each step\n" );
		if( tbitr && tbrweight == 0 ) fprintf( fp, "  unweighted\n" ); 
		if( tbitr && tbrweight == 3 ) fprintf( fp, "  reversely weighted\n" ); 
		if( tbweight ) fprintf( fp, "  weighted\n" ); 
		fprintf( fp, "\n" );
	}

   	 fprintf( fp, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );

	if( alg == 'a' )
		fprintf( fp, "Algorithm A\n" );
	else if( alg == 'A' ) 
		fprintf( fp, "Algorithm A+\n" );
	else if( alg == 'S' ) 
		fprintf( fp, "Apgorithm S\n" );
	else
		fprintf( fp, "Unknown algorithm\n" );

    if( use_fft )
    {
        fprintf( fp, "FFT on\n" );
        if( dorp == 'd' )
            fprintf( fp, "Basis : 4 nucleotides\n" );
        else
        {
            if( fftscore )
                fprintf( fp, "Basis : Polarity and Volume\n" );
            else
                fprintf( fp, "Basis : 20 amino acids\n" );
        }
        fprintf( fp, "Threshold   of anchors = %d%%\n", fftThreshold );
        fprintf( fp, "window size of anchors = %dsites\n", fftWinSize );
    }
	else
        fprintf( fp, "FFT off\n" );
	fflush( fp );
}
	 

int main( int argc, char *argv[] )
{
	int  nlen[M];	
	char **name, **seq;
	char **mseq1, **mseq2;
	char **aseq;
	char **bseq;
	double *eff;
	int i;
	FILE *infp;
	char c;
	int alloclen;

	arguments( argc, argv );
#ifndef enablemultithread
	nthread = 0;
#endif

	if( inputfile )
	{
		infp = fopen( inputfile, "r" );
		if( !infp )
		{
			fprintf( stderr, "Cannot open %s\n", inputfile );
			exit( 1 );
		}
	}
	else
		infp = stdin;

	getnumlen( infp );
	rewind( infp );

	if( njob < 2 )
	{
		fprintf( stderr, "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}
	if( njob > M )
	{
		fprintf( stderr, "The number of sequences must be < %d\n", M );
		fprintf( stderr, "Please try the splittbfast program for such large data.\n" );
		exit( 1 );
	}

	alloclen = nlenmax*2;
	seq = AllocateCharMtx( njob, alloclen+10 );
	aseq = AllocateCharMtx( njob, alloclen+10 );
	bseq = AllocateCharMtx( njob, alloclen+10 );
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );
	name = AllocateCharMtx( njob, B );

	eff = AllocateDoubleVec( njob );

#if 0
	Read( name, nlen, seq );
#else
	readData_pointer( infp, name, nlen, seq );
#endif
	fclose( infp );

	constants( njob, seq );

#if 0
	fprintf( stderr, "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM();

	initFiles();

	WriteOptions( trap_g );

	c = seqcheck( seq );
	if( c )
	{
		fprintf( stderr, "Illegal character %c\n", c );
		exit( 1 );
	}

//	writePre( njob, name, nlen, seq, 0 );

	for( i=0; i<njob; i++ ) eff[i] = 1.0;


	for( i=0; i<njob; i++ ) gappick0( bseq[i], seq[i] );

	pairalign( name, nlen, bseq, aseq, mseq1, mseq2, eff, alloclen );

	fprintf( trap_g, "done.\n" );
#if DEBUG
	fprintf( stderr, "closing trap_g\n" );
#endif
	fclose( trap_g );

//	writePre( njob, name, nlen, aseq, !contin );
#if 0
	writeData( stdout, njob, name, nlen, aseq );
#endif
#if IODEBUG
	fprintf( stderr, "OSHIMAI\n" );
#endif
	SHOWVERSION;

	if( stdout_dist && nthread > 1 )
	{
		fprintf( stderr, "\nThe order of distances is not identical to that in the input file, because of the parallel calculation.  Reorder them by yourself, using sort -n -k 2 | sort -n -k 1 -s\n" );
	}
	if( stdout_align && nthread > 1 )
	{
		fprintf( stderr, "\nThe order of pairwise alignments is not identical to that in the input file, because of the parallel calculation.  Reorder them by yourself.\n" );
	}
	FreeCharMtx( seq );
	FreeCharMtx( aseq );
	FreeCharMtx( bseq );
	FreeCharMtx( name );
	free( mseq1 );
	free( mseq2 );
	free( eff );

	return( 0 );
}

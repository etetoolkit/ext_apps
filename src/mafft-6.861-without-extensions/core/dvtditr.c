 /* Tree-dependent-iteration */
 /* Devide to segments       */ 

#include "mltaln.h"

extern char **seq_g;
extern char **res_g;

static int intop;
static int intree;

void arguments( int argc, char *argv[] )
{
	int c;
	char *argkey;

	outnumber = 0;
	nthread = 1;
	randomseed = 0;
	scoreout = 0;
	parallelizationstrategy = BAATARI1;
	intop = 0;
	intree = 0;
	inputfile = NULL;
	rnakozo = 0;
	rnaprediction = 'm';
	nevermemsave = 0;
	score_check = 1;
	fftkeika = 1;
	constraint = 0;
	fmodel = 0;
	kobetsubunkatsu = 1;
	bunkatsu = 1;
	nblosum = 62;
	niter = 100;
	calledByXced = 0;
	devide = 1;
	divWinSize = 20; /* 70 */
	divThreshold = 65;
	fftscore = 1;
	fftRepeatStop = 0;
	fftNoAnchStop = 0;
    scmtd = 5;
	cooling = 1;
    weight = 4;
    utree = 1;
    refine = 1;
    check = 1;
    cut = 0.0;
	disp = 0;
	outgap = 1;
	use_fft = 0; // CHUUI dochira demo mafft.tmpl deha F
	force_fft = 0;
	alg = 'A';  /* chuui */
	mix = 0;
	checkC = 0;
	tbitr = 0;
	treemethod = 'X';
	scoremtx = 1;
	dorp = NOTSPECIFIED;
	ppenalty = NOTSPECIFIED;
	ppenalty_ex = NOTSPECIFIED;
	poffset = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	RNAppenalty = NOTSPECIFIED;
	RNAppenalty_ex = NOTSPECIFIED;
	RNApthr = NOTSPECIFIED;
	TMorJTT = JTT;
	consweight_multi = 1.0;
	consweight_rna = 0.0;

	while( --argc > 0 && (*++argv)[0] == '-' )
	{
		while ( (c = *++argv[0]) )
		{
			switch( c )
			{
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'I':
					niter = atoi( *++argv );
					fprintf( stderr, "niter = %d\n", niter );
					--argc;
					goto nextoption;
				case 'e':
					RNApthr = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'o':
					RNAppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
					fprintf( stderr, "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = atoi( *++argv );
					fprintf( stderr, "kappa = %d\n", kimuraR );
					--argc; 
					goto nextoption;
				case 'b':
					nblosum = atoi( *++argv );
					scoremtx = 1;
					fprintf( stderr, "blosum %d / kimura 200\n", nblosum );
					--argc; 
					goto nextoption;
				case 'j':
					pamN = atoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
					fprintf( stderr, "jtt/kimura %d\n", pamN );
					--argc; 
					goto nextoption;
				case 'm':
					pamN = atoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
					fprintf( stderr, "tm %d\n", pamN );
					--argc; 
					goto nextoption;
				case 'l':
					fastathreshold = atof( *++argv );
					constraint = 2;
					--argc;
					goto nextoption;
				case 'r':
					consweight_rna = atof( *++argv );
					rnakozo = 1;
					--argc;
					goto nextoption;
				case 'c':
					consweight_multi = atof( *++argv );
					--argc;
					goto nextoption;
				case 'C':
					nthread = atoi( *++argv );
					fprintf( stderr, "nthread = %d\n", nthread );
					--argc; 
					goto nextoption;
				case 't':
					randomseed = atoi( *++argv );
					fprintf( stderr, "randomseed = %d\n", randomseed );
					--argc; 
					goto nextoption;
				case 'p':
					argkey = *++argv;
					if( !strcmp( argkey, "BESTFIRST" ) ) parallelizationstrategy = BESTFIRST;
					else if( !strcmp( argkey, "BAATARI0" ) ) parallelizationstrategy = BAATARI0;
					else if( !strcmp( argkey, "BAATARI1" ) ) parallelizationstrategy = BAATARI1;
					else if( !strcmp( argkey, "BAATARI2" ) ) parallelizationstrategy = BAATARI2;
					else
					{
						fprintf( stderr, "Unknown parallelization strategy, %s\n", argkey );
						exit( 1 );
					}
//					exit( 1 );
					--argc; 
					goto nextoption;
				case 'S' :
					scoreout = 1;
					break;
				case 's' :
					RNAscoremtx = 'r';
					break;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'N':
					nevermemsave = 1;
					break;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
				case 'Q':
					alg = 'Q';
					break;
				case 'R':
					rnaprediction = 'r';
					break;
				case 'O':
					fftNoAnchStop = 1;
					break;
#if 0
				case 'e':
					fftscore = 0;
					break;
				case 'r':
					fmodel = -1;
					break;
				case 'R':
					fftRepeatStop = 1;
					break;
#endif
				case 'T':
					kobetsubunkatsu = 0;
					break;
				case 'B':
					bunkatsu = 0;
					break;
#if 0
				case 'c':
					cooling = 1;
					break;
				case 'a':
					alg = 'a';
					break;
				case 's' :
					treemethod = 's';
					break;
#endif
				case 'H':
					alg = 'H';
					break;
				case 'A':
					alg = 'A';
					break;
				case 'M':
					alg = 'M';
					break;
				case 'F':
					use_fft = 1;
					break;
#if 0
				case 't':
					weight = 4;
					break;
#endif
				case 'u':
					weight = 0;
					break;
				case 'U':
					intree = 1;
					break;
				case 'V':
					intop = 1;
					break;
				case 'J':
					utree = 0;
					break;
				case 'd':
					disp = 1;
					break;
				case 'Z':
					score_check = 0;
					break;
				case 'Y':
					score_check = 2;
					break;
#if 0
				case 'n' :
					treemethod = 'n';
					break;
#endif
				case 'n' :
					outnumber = 1;
					break;
				case 'X' :
					treemethod = 'X';
					break;
				case 'E' :
					treemethod = 'E';
					break;
				case 'q' :
					treemethod = 'q';
					break;
				case 'z':
					fftThreshold = atoi( *++argv );
					--argc;
					goto nextoption;
				case 'w':
					fftWinSize = atoi( *++argv );
					--argc;
					goto nextoption;
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
		fprintf( stderr, "options : Check source file!\n" );
		exit( 1 );
	}
#if 0
	if( alg == 'A' && weight == 0 ) 
		ErrorExit( "ERROR : Algorithm A+ and un-weighted\n" ); 
#endif
}


int main( int argc, char *argv[] )
{
    int identity;
	static int nlen[M];
	static char **name, **seq, **aseq, **bseq;
	static Segment *segment = NULL;
	static int anchors[MAXSEG];
	int i, j;
	int iseg, nseg;
	int ***topol;
	double **len;
	double **eff;
	FILE *prep;
	FILE *infp;
	int alloclen;
	int returnvalue;
	char c;
	int ocut;
	char **seq_g_bk;
	LocalHom **localhomtable = NULL; // by D.Mathog
	RNApair ***singlerna;
	int nogaplen;
	static char **nogap1seq;
	static char *kozoarivec;
	int nkozo;

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

#if 0
	PreRead( stdin, &njob, &nlenmax );
#else
	getnumlen( infp );
#endif
	rewind( infp );

	nkozo = 0;

	if( njob < 2 ) 
	{
		fprintf( stderr, "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}

	ocut = cut;

	segment = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
//	if( treemethod == 'X' || treemethod == 'E' || treemethod == 'q' )
	topol = AllocateIntCub( njob, 2, njob );
	len = AllocateDoubleMtx( njob, 2 );
	eff = AllocateDoubleMtx( njob, njob );
	seq = AllocateCharMtx( njob, nlenmax*9+1 );
	name = AllocateCharMtx( njob, B+1 );
	seq_g = AllocateCharMtx( njob, nlenmax*9+1 );
	res_g = AllocateCharMtx( njob, nlenmax*9+1 );
	aseq = AllocateCharMtx( njob, nlenmax*9+1 );
	bseq = AllocateCharMtx( njob, nlenmax*9+1 );
	nogap1seq = AllocateCharMtx( 1, nlenmax*1+1 );
	alloclen = nlenmax * 9;
	seq_g_bk = AllocateCharMtx( njob, 0 );
	for( i=0; i<njob; i++ ) seq_g_bk[i] = seq_g[i];
	kozoarivec = AllocateCharVec( njob );

	if( constraint )
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
				localhomtable[i][j].overlapaa = -1.0; 
				localhomtable[i][j].opt = -1.0;
				localhomtable[i][j].importance = -1.0; 
				localhomtable[i][j].next = NULL; 
				localhomtable[i][j].nokori = 0;
				localhomtable[i][j].extended = -1;
				localhomtable[i][j].last = localhomtable[i]+j;
				localhomtable[i][j].korh = 'h';
			}
		}
		fprintf( stderr, "Loading 'hat3' ... " );
		fflush( stderr );
		prep = fopen( "hat3", "r" );
		if( prep == NULL ) ErrorExit( "Make hat3." );
		readlocalhomtable2( prep, njob, localhomtable, kozoarivec );
		fclose( prep ); 
//		for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ )
//			fprintf( stdout, "%d %d %d %d %d %d %d\n", i, j, localhomtable[i][j].opt, localhomtable[i][j].start1, localhomtable[i][j].end1, localhomtable[i][j].start2, localhomtable[i][j].end2 );
		fprintf( stderr, "done.\n" );
		fflush( stderr );
#if 0
		fprintf( stderr, "Extending localhom ... " );
		extendlocalhom( njob, localhomtable );
		fprintf( stderr, "done.\n" );
#endif
		nkozo = 0;
		for( i=0; i<njob; i++ ) if( kozoarivec[i] ) nkozo++;
	}


//		if( nlenmax > 30000 )
		if( nlenmax > 50000 ) // version >= 6.823
		{
#if 0
			if( constraint )
			{
				fprintf( stderr, "\nnlenmax=%d, nagasugi!\n", nlenmax ); 
				exit( 1 );
			}
			if( nevermemsave )
			{
				fprintf( stderr, "\nnevermemsave=1, nlenmax=%d, nagasugi!\n", nlenmax ); 
				exit( 1 );
			}
#endif
			if( !constraint && !nevermemsave && alg != 'M' )
			{
				fprintf( stderr, "\nnlenmax=%d, Switching to the memsave mode\n", nlenmax ); 
				alg = 'M';
			}
		}

#if 0
	Read( name, nlen, seq_g );
#else
	readData_pointer( infp, name, nlen, seq_g );
#endif
	fclose( infp );


	for( i=0; i<njob; i++ )
	{
		res_g[i][0] = 0;
	}

	identity = 1;
	for( i=0; i<njob; i++ ) 
	{
		identity *= ( nlen[i] == nlen[0] );
	}
	if( !identity ) 
	{
		fprintf( stderr, "Input pre-aligned data\n" );
		exit( 1 );
	}
	constants( njob, seq_g );

#if 0
	fprintf( stderr, "penalty = %d\n", penalty ); 
	fprintf( stderr, "penalty_ex = %d\n", penalty_ex ); 
	fprintf( stderr, "offset = %d\n", offset ); 
#endif

	initSignalSM();

	initFiles();

#if 0
	if( njob == 2 )
	{
		writePre( njob, name, nlen, seq_g, 1 );
		SHOWVERSION;
		return( 0 );
	}
#else
	if( njob == 2 )
	{
		weight = 0;
		niter = 1;
	}
#endif

	c = seqcheck( seq_g );
	if( c )
	{
		fprintf( stderr, "Illegal character %c\n", c );
		exit( 1 );
	}
	commongappick( njob, seq_g );

	if( rnakozo && rnaprediction == 'm' )
	{
		singlerna = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
		prep = fopen( "hat4", "r" );
		if( prep == NULL ) ErrorExit( "Make hat4 using mccaskill." );
		fprintf( stderr, "Loading 'hat4' ... " );
		fflush( stderr );
		for( i=0; i<njob; i++ )
		{
			gappick0( nogap1seq[0], seq_g[i] );
			nogaplen = strlen( nogap1seq[0] );
			singlerna[i] = (RNApair **)calloc( nogaplen, sizeof( RNApair * ) );
			for( j=0; j<nogaplen; j++ )
			{
				singlerna[i][j] = (RNApair *)calloc( 1, sizeof( RNApair ) );
				singlerna[i][j][0].bestpos = -1;
				singlerna[i][j][0].bestscore = -1.0;
			}
			readmccaskill( prep, singlerna[i], nogaplen );
		}
		fclose( prep );
		fprintf( stderr, "\ndone.\n" );
		fflush( stderr );
	}
	else
		singlerna = NULL;



	if( utree )
	{
		prep = fopen( "hat2", "r" );
		if( !prep ) ErrorExit( "Make hat2." );
		readhat2_pointer( prep, njob, name, eff );
		fclose( prep );
#if DEBUG
		for( i=0; i<njob-1; i++ ) 
		{
			for( j=i+1; j<njob; j++ ) 
			{
				printf( " %f", eff[i][j] );
			}
			printf( "\n" );
		}
#endif
		if( intree )
			veryfastsupg_double_loadtree( njob, eff, topol, len );
		else if( intop ) // v6.528 deha if( intop ) dattanode intree ga mukou datta.
			veryfastsupg_double_loadtop( njob, eff, topol, len );
		else if( treemethod == 'X' || treemethod == 'E' || treemethod == 'q' ) 
//			veryfastsupg_double_outtree( njob, eff, topol, len, name );
			fixed_musclesupg_double_treeout( njob, eff, topol, len, name );
		else if( treemethod == 'n' ) 
			nj( njob, eff, topol, len );
		else if( treemethod == 's' )
			spg( njob, eff, topol, len );
		else if( treemethod == 'p' )
			upg2( njob, eff, topol, len );
		else ErrorExit( "Incorrect treemethod.\n" );
	}
#if DEBUG
	printf( "utree = %d\n", utree );
#endif


	fprintf( stderr, "\n" );
	if( ( !utree && kobetsubunkatsu ) || constraint || !bunkatsu )
	{
		nseg = 0;
		anchors[0] =0;
		anchors[1] =strlen( seq_g[0] );
		nseg += 2;
	}
	else
	{
		nseg = searchAnchors( njob, seq_g, segment );
#if 0
		fprintf( stderr, "### nseg = %d\n", nseg );
		fprintf( stderr, "### seq_g[0] = %s\n", seq_g[0] );
		fprintf( stderr, "### nlenmax = %d\n", nlenmax );
		fprintf( stderr, "### len = %d\n", strlen( seq_g[0] ) );

#endif

		anchors[0] = 0;
		for( i=0; i<nseg; i++ ) anchors[i+1] = segment[i].center;
		anchors[nseg+1] = strlen( seq_g[0] );
		nseg +=2;

#if 0
		for( i=0; i<nseg; i++ )
			fprintf( stderr, "anchor[%d] = %d\n", i, anchors[i] );
#endif
	}

	for( i=0; i<njob; i++ ) res_g[i][0] = 0;

	for( iseg=0; iseg<nseg-1; iseg++ )
	{
		int tmplen = anchors[iseg+1]-anchors[iseg];
		int pos = strlen( res_g[0] );
		for( j=0; j<njob; j++ )
		{
			strncpy( seq[j], seq_g[j], tmplen );
			seq[j][tmplen]= 0;
			seq_g[j] += tmplen;	

		}
		fprintf( stderr, "Segment %3d/%3d %4d-%4d\n", iseg+1, nseg-1, pos+1, pos+1+tmplen );
		fflush( stderr );
		fprintf( trap_g, "Segment %3d/%3d %4d-%4d\n", iseg+1, nseg-1, pos+1, pos+1+tmplen );
	
		cut = ocut;
		returnvalue = TreeDependentIteration( njob, name, nlen, seq, bseq, topol, len, alloclen, localhomtable, singlerna, nkozo, kozoarivec );

		for( i=0; i<njob; i++ )
			strcat( res_g[i], bseq[i] );
	}
	FreeCharMtx( seq_g_bk );
	FreeIntCub( topol );
	FreeDoubleMtx( len );
	FreeDoubleMtx( eff );
	if( constraint ) FreeLocalHomTable( localhomtable, njob );

#if 0
	Write( stdout, njob, name, nlen, bseq );
#endif

	fprintf( stderr, "done\n" );
	fprintf( trap_g, "done\n" );
	fclose( trap_g );


	devide = 0; 
	writePre( njob, name, nlen, res_g, 1 );
#if 0
	writeData( stdout, njob, name, nlen, res_g, 1 );
#endif


	SHOWVERSION;
	return( 0 );
}

#if 0
signed int main( int argc, char *argv[] )
{
	int i, nlen[M];
	char b[B];
	char a[] = "=";
	int value;

	gets( b ); njob = atoi( b );

/*
	scoremtx = 0;
	if( strstr( b, "ayhoff" ) ) scoremtx = 1;
	else if( strstr( b, "dna" ) || strstr( b, "DNA" ) ) scoremtx = -1;
	else if( strstr( b, "M-Y" ) || strstr( b, "iyata" ) ) scoremtx = 2;
	else scoremtx = 0;
*/
	if( strstr( b, "constraint" ) ) cnst = 1;

	nlenmax = 0;
	i = 0;
	while( i<njob )
	{
		gets( b );
		if( !strncmp( b, a, 1 ) ) 
		{
			gets( b ); nlen[i] = atoi( b );
			if( nlen[i] > nlenmax ) nlenmax = nlen[i];
			i++;
		}
	}
	if( nlenmax > N || njob > M ) 
	{
		fprintf( stderr, "ERROR in main\n" );
		exit( 1 );
	}
	/*
	nlenmax = Na;
	*/
	rewind( stdin );
	value = main1( nlen, argc, argv );
	exit( 0 );
}
#endif

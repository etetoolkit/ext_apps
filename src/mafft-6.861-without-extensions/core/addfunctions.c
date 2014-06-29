#include "mltaln.h"

void profilealignment( int n0, int n1, int n2, char **aln0, char **aln1, char **aln2, int alloclen, char alg ) // n1 ha allgap
{
	int i, newlen;
	double *effarr0, *effarr2;
	float dumfl;
	double eff;
	effarr0 = AllocateDoubleVec( n0 );
	effarr2 = AllocateDoubleVec( n2 );

	commongappick( n0, aln0 );
	commongappick( n2, aln2 );


	eff = 1.0 / (double)n0; for( i=0; i<n0; i++ ) effarr0[i] = eff;
	eff = 1.0 / (double)n2; for( i=0; i<n2; i++ ) effarr2[i] = eff;

	if( alg == 'M' )
		MSalignmm( aln0, aln2, effarr0, effarr2, n0, n2, alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 ); //outgap=1??
	else
		A__align( aln0, aln2, effarr0, effarr2, n0, n2, alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 ); //outgap=1??

	newlen = strlen( aln0[0] );

	for( i=0; i<newlen; i++ ) aln1[0][i] = '-';
	aln1[0][i] = 0;
	for( i=1; i<n1; i++ ) strcpy( aln1[i], aln1[0] );


	free( effarr0 );
	free( effarr2 );
}

void eq2dash( char *s )
{
	while( *s )
	{
		if( *s == '=' ) *s = '-';
		s++;
	}
}

void findnewgaps( int n, char **seq, int *gaplen )
{
	int i, pos, len;
	len = strlen( seq[0] );	

//	fprintf( stderr, "seq[0] = %s\n", seq[0] );
	for( i=0; i<len; i++ ) gaplen[i] = 0;
	
	pos = 0;
	for( i=0; i<len; i++ )
	{
		if( seq[0][i] == '=' ) 
		{
//			fprintf( stderr, "Newgap! pos = %d\n", pos );
			gaplen[pos]++;
		}
		else
			pos++;
	}
}

void findcommongaps( int n, char **seq, int *gaplen )
{
	int i, j, pos, len, len1;
	len = strlen( seq[0] );	
	len1 = len+1;

//	fprintf( stderr, "seq[0] = %s\n", seq[0] );
	for( i=0; i<len1; i++ ) gaplen[i] = 0;
	
	pos = 0;
	for( i=0; i<len; i++ )
	{
		for( j=0; j<n; j++ )
			if( seq[j][i] != '-' ) break;

		if( j == n ) gaplen[pos]++;
		else
			pos++;
	}
#if 0
	for( i=0; i<pos; i++ )
	{
		fprintf( stderr, "vec[%d] = %d\n", i, gaplen[i] );
	}
#endif
}

void adjustgapmap( int newlen, int *gapmap, char *seq )
{
	int j;
	int pos;
	int newlen1 = newlen+1;
	int *tmpmap;

	tmpmap = AllocateIntVec( newlen+2 );
	j = 0;
	pos = 0;
	while( *seq )
	{
//		fprintf( stderr, "j=%d *seq = %c\n", j, *seq );
		if( *seq++ == '=' )
			tmpmap[j++] = 0;
		else
		{
			tmpmap[j++] = gapmap[pos++];
		}
	}
	tmpmap[j++] = gapmap[pos];

	for(j=0; j<newlen1; j++)
		gapmap[j] = tmpmap[j];

	free( tmpmap );
}

void insertnewgaps( int njob, int *alreadyaligned, char **seq, int *ex1, int *ex2, int *gaplen, int *gapmap, int alloclen, char alg )
{
	int *mar;
	char *gaps;
	char *cptr;
	int i, j, k, len, rep, len0;
	char **mseq2, **mseq0, **mseq1;
	char **aseq;
	int ngroup2, ngroup0, ngroup1;
	int *list0, *list1, *list2;
	int posin12, gapshift, newpos;

	mar = calloc( njob, sizeof( int ) );
	list0 = calloc( njob, sizeof( int ) );
	list1 = calloc( njob, sizeof( int ) );
	list2 = calloc( njob, sizeof( int ) );

	for( i=0; i<njob; i++ ) mar[i] = 0;
	for( i=0; i<njob; i++ ) 
	{
		if( alreadyaligned[i]==0 ) mar[i] = 3;
	}
	for( i=0; (k=ex1[i])>-1; i++ ) 
	{
		mar[k] = 1;
//		fprintf( stderr, "excluding %d\n", ex1[i] );
	}
	for( i=0; (k=ex2[i])>-1; i++ ) 
	{
		mar[k] = 2;
//		fprintf( stderr, "excluding %d\n", ex2[i] );
	}

	ngroup2 = ngroup1 = ngroup0 = 0;
	for( i=0; i<njob; i++ )
	{
		if( mar[i] == 2 ) 
		{
			list2[ngroup2] = i;
			ngroup2++;
		}
		if( mar[i] == 1 ) 
		{
			list1[ngroup1] = i;
			ngroup1++;
		}
		if( mar[i] == 0 ) 
		{
			list0[ngroup0] = i;
			ngroup0++;
		}
	}
	list0[ngroup0] = list1[ngroup1] = list2[ngroup2] = -1;
	if( ngroup0 == 0 )
	{
		fprintf( stderr, "Nothing to do\n" );
		free( mar );
		free( list0 );
		free( list1 );
		free( list2 );
		return;
	}

	for( i=0; i<njob; i++ ) if( mar[i] == 0 ) break;
	rep = i;
	len = strlen( seq[rep] );
	len0 = len+1;

//
//	if( i == njob )
//	{
////		fprintf( stderr, "Nothing to do\n" );
//		free( mar );
//		return;
//	}

	mseq2 = AllocateCharMtx( ngroup2, alloclen );
	mseq1 = AllocateCharMtx( ngroup1, alloclen );
	mseq0 = AllocateCharMtx( ngroup0, alloclen );
	aseq = AllocateCharMtx( njob, alloclen );
	gaps = calloc( alloclen, sizeof( char ) );

	for( i=0; i<njob; i++ ) aseq[i][0] = 0;
	posin12 = 0;
	for( j=0; j<len0; j++ )
	{
		if( gaplen[j] )
		{
			for( i=0; i<ngroup0; i++ ) mseq0[i][0] = 0;
			for( i=0; i<ngroup1; i++ ) mseq1[i][0] = 0;
			for( i=0; i<ngroup2; i++ ) mseq2[i][0] = 0;

			gapshift = gaplen[j];
			cptr = gaps;
			while( gapshift-- ) *cptr++ = '-';
			*cptr = 0;
			gapshift = gaplen[j];

			for( i=0; i<ngroup0; i++ ) strncat( mseq0[i], gaps, gapshift );
			for( i=0; i<ngroup1; i++ ) strncat( mseq1[i], seq[list1[i]]+posin12, gapshift );
			for( i=0; i<ngroup2; i++ ) strncat( mseq2[i], seq[list2[i]]+posin12, gapshift );
			posin12 += gapshift;

			gapshift = gapmap[posin12];
//			fprintf( stderr, "gapmap[%d] kouho = %d\n", posin12, gapmap[posin12] );


			for( i=0; i<ngroup0; i++ ) strncat( mseq0[i], seq[list0[i]]+j, gapshift );
			for( i=0; i<ngroup1; i++ ) strncat( mseq1[i], seq[list1[i]]+posin12, gapshift );
			for( i=0; i<ngroup2; i++ ) strncat( mseq2[i], seq[list2[i]]+posin12, gapshift );
#if 0
			for( i=0; i<ngroup0; i++ ) fprintf( stderr, "### mseq0[%d] = %s\n", i, mseq0[i] );
			for( i=0; i<ngroup1; i++ ) fprintf( stderr, "### mseq1[%d] = %s\n", i, mseq1[i] );
			for( i=0; i<ngroup2; i++ ) fprintf( stderr, "### mseq2[%d] = %s\n", i, mseq2[i] );
#endif

			if( gapshift ) profilealignment( ngroup0, ngroup1, ngroup2, mseq0, mseq1, mseq2, alloclen, alg );

			j += gapshift;
			posin12 += gapshift;

			for( i=0; i<ngroup0; i++ ) strcat( aseq[list0[i]], mseq0[i] );
			for( i=0; i<ngroup1; i++ ) strcat( aseq[list1[i]], mseq1[i] );
			for( i=0; i<ngroup2; i++ ) strcat( aseq[list2[i]], mseq2[i] );
		}

		newpos = strlen( aseq[rep] );
		for( i=0; i<ngroup0; i++ ) aseq[list0[i]][newpos] = seq[list0[i]][j];
		for( i=0; i<ngroup1; i++ ) aseq[list1[i]][newpos] = seq[list1[i]][posin12];
		for( i=0; i<ngroup2; i++ ) aseq[list2[i]][newpos] = seq[list2[i]][posin12];
		for( i=0; i<ngroup0; i++ ) aseq[list0[i]][newpos+1] = 0;
		for( i=0; i<ngroup1; i++ ) aseq[list1[i]][newpos+1] = 0;
		for( i=0; i<ngroup2; i++ ) aseq[list2[i]][newpos+1] = 0;

		posin12++;
	}

//	for( i=0; i<njob; i++ ) if( mar[i] != 3 ) strcpy( seq[i], aseq[i] );
	for( i=0; i<ngroup0; i++ ) strcpy( seq[list0[i]], aseq[list0[i]] );
	for( i=0; i<ngroup1; i++ ) strcpy( seq[list1[i]], aseq[list1[i]] );
	for( i=0; i<ngroup2; i++ ) strcpy( seq[list2[i]], aseq[list2[i]] );

	free( mar );
	free( gaps );
	free( list0 );
	free( list1 );
	free( list2 );
	FreeCharMtx( mseq2 );
	FreeCharMtx( mseq0 );
}


void restorecommongaps( int njob, char **seq, int *ex1, int *ex2, int *gaplen, int alloclen )
{
	int *mar;
	char *tmpseq;
	char *cptr;
	int *iptr;
	int *tmpgaplen;
	int i, j, k, len, rep, len1;

	mar = calloc( njob, sizeof( int ) );
	tmpseq = calloc( alloclen, sizeof( char ) );
	tmpgaplen = calloc( alloclen, sizeof( int ) );
//	tmpseq = calloc( alloclen+2, sizeof( char ) );
//	tmpgaplen = calloc( alloclen+2, sizeof( int ) );


	for( i=0; i<njob; i++ ) mar[i] = 0;
	for( i=0; (k=ex1[i])>-1; i++ ) 
	{
		mar[k] = 1;
//		fprintf( stderr, "excluding %d\n", ex1[i] );
	}
	for( i=0; (k=ex2[i])>-1; i++ ) 
	{
		mar[k] = 1;
//		fprintf( stderr, "excluding %d\n", ex2[i] );
	}

	for( i=0; i<njob; i++ )
		if( mar[i] ) break;

	if( i == njob )
	{
//		fprintf( stderr, "Nothing to do\n" );
		free( mar );
		free( tmpseq );
		free( tmpgaplen );
		return;
	}
	rep = i;
	len = strlen( seq[rep] );
	len1 = len+1;


	for( i=0; i<njob; i++ )
	{
		if( mar[i] == 0 ) continue;
		cptr = tmpseq;
		for( j=0; j<len1; j++ )
		{
			for( k=0; k<gaplen[j]; k++ )
				*(cptr++) = '-';
			*(cptr++) = seq[i][j];
		}
		*cptr = 0;
		strcpy( seq[i], tmpseq );
	}

	iptr = tmpgaplen;
	for( j=0; j<len1; j++ )
	{
		*(iptr++) = gaplen[j];
		for( k=0; k<gaplen[j]; k++ )
			*(iptr++) = 0;
	}
	*iptr = -1;

	iptr = tmpgaplen;
	while( *iptr != -1 ) *gaplen++ = *iptr++;

	free( mar );
	free( tmpseq );
	free( tmpgaplen );
}

#if 0
int samemember( int *mem, int *cand )
{
	int i, j;

#if 0
	fprintf( stderr, "mem = " );
	for( i=0; mem[i]>-1; i++ )	fprintf( stderr, "%d ", mem[i] );
	fprintf( stderr, "\n" );

	fprintf( stderr, "cand = " );
	for( i=0; cand[i]>-1; i++ )	fprintf( stderr, "%d ", cand[i] );
	fprintf( stderr, "\n" );
#endif

	for( i=0, j=0; mem[i]>-1; )	
	{
		if( mem[i++] != cand[j++] ) return( 0 );
	}

	if( cand[j] == -1 )
	{
		return( 1 );
	}
	else
	{
		return( 0 );
	}
}
#else
int samemember( int *mem, int *cand )
{
	int i, j;
	int nm, nc;
#if 0
	fprintf( stderr, "mem = " );
	for( i=0; mem[i]>-1; i++ )	fprintf( stderr, "%d ", mem[i] );
	fprintf( stderr, "\n" );

	fprintf( stderr, "cand = " );
	for( i=0; cand[i]>-1; i++ )	fprintf( stderr, "%d ", cand[i] );
	fprintf( stderr, "\n" );
#endif

	nm = 0; for( i=0; mem[i]>-1; i++ ) nm++;
	nc = 0; for( i=0; cand[i]>-1; i++ ) nc++;

	if( nm != nc ) return( 0 );

	for( i=0; mem[i]>-1; i++ )	
	{
		for( j=0; cand[j]>-1; j++ )
			if( mem[i] == cand[j] ) break;
		if( cand[j] == -1 ) return( 0 );
	}

	if( mem[i] == -1 )
	{
		return( 1 );
	}
	else
	{
		return( 0 );
	}
}
#endif


int includemember( int *mem, int *cand ) // mem in cand 
{
	int i, j;

#if 0
	fprintf( stderr, "mem = " );
	for( i=0; mem[i]>-1; i++ )	fprintf( stderr, "%d ", mem[i] );
	fprintf( stderr, "\n" );

	fprintf( stderr, "cand = " );
	for( i=0; cand[i]>-1; i++ )	fprintf( stderr, "%d ", cand[i] );
	fprintf( stderr, "\n" );
#endif

	for( i=0; mem[i]>-1; i++ )
	{
		for( j=0; cand[j]>-1; j++ )
			if( mem[i] == cand[j] ) break;
		if( cand[j] == -1 ) return( 0 );
	}
//	fprintf( stderr, "INCLUDED! mem[0]=%d\n", mem[0] );
	return( 1 );
}

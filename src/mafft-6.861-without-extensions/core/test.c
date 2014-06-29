#include <stdio.h>
#include <unistd.h>
main()
{
	fprintf( stdout, "nprocess = %d\n", sysconf( _SC_NPROCESSORS_ONLN ) );
	fprintf( stdout, "nprocess = %d\n", sysconf( _SC_NPROCESSORS_CONF ) );
}


$infile = $ARGV[0];

open(INFILE,$infile);

@a = <INFILE>;

$flag = 0;
foreach $line (@a)	{
	if ($flag)	{
		@b = split(/\s/,$line);
		if ((@b) != 2)	{
			die "error : (@b)\n";
		}
		chomp $b[1];
		print "$b[0]\t", 100 * ($b[1] =~ tr/\?//) / length($b[1]), "\n";
	}
	if ($line =~ /matrix/)	{
		$flag = 1;
	}
}
		
		

	

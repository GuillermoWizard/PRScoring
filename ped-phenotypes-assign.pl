#! /bin/perl


open(PED,'<',$ARGV[0])or die $!;


#@test = (1..10);
#print "@test\n" ;
#shift @test for 0..3;
#print "@test\n" ;


while($a=<PED>){
	chomp($a);
	@ped_f = split /\t/,$a;
	@ped_a = @ped_f ; 
	$cols = scalar(@ped_a); 
	#print "Number of columns: $cols \n";
	shift @ped_a for 0..6;

	##print("$ped_f[0]\t$ped_f[1]\t$ped_f[2]\t$ped_f[3]\t$ped_f[4]\n");

	open(PHEN,'<',$ARGV[1]) or die $!;
	while($line = <PHEN>){
		chomp($line);
		@pheno_f = split /\t/,$line;
		if($ped_f[1] eq $pheno_f[0]){
			
			#print("found...");
			@ped_info=("$ped_f[1]","$ped_f[1]","$ped_f[2]","$ped_f[3]","$ped_f[4]","$pheno_f[1]");
			print "@ped_info\t";
			print "@ped_a\n";
			#$i=1;
			#foreach $j (@ped_f) {
			#	print "$j\t";
			#	$i++;
			#	if($i==6){
			#		print "\n";
			#		next;
			#	}
			#	next if $i == 6;
			#}

		}

	}
	close(PHEN);
}
close(PED);
use strict;
use threads;

if (@ARGV < 2){
	die "\n\tUsage: perl $0 <FILE_TAB_RUTES> <THREADS_NUMBER>
			Example: perl $0 z_scrit1_argv0_paths_htseq.list 20\n\n";
}

#### ------------------------------------------------ OPTIONS
my $num_of_threads = "10";
if (defined $ARGV[1]){
	$num_of_threads = $ARGV[1];
}




#### ------------------------------------------------ RECOVER FILES SPLITTED

my @lines_rute = `cat $ARGV[0]`;
chomp @lines_rute;

my @threads = initThreads();



#### ------------------------------------------------ DO PROCESS
my $number_xfiles = scalar (@lines_rute);

my $i2 =  $number_xfiles + $num_of_threads;
my $index = 0;

for (my $j=1; $j<=$i2; $j+=$num_of_threads){
	foreach my $each_thread (@threads){
		# Tell each thread to perform our 'doOperation()' subroutine.
		$each_thread = threads->create(\&doOperation);
		$index++;
	}

	# This tells the main program to keep running until all threads have finished.
	foreach(@threads){
		$_->join();
	}
}





#### ------------------------------------------------ SUB RUTINAS ------------------------------------------------ ###

sub initThreads{
	my @initThreads;
	for(my $i=1; $i<=$num_of_threads; $i++){
		push(@initThreads,$i);
	}
	return @initThreads;
}


sub doOperation{
	# Get the thread id. Allows each thread to be identified.
	my $id = threads->tid();

	my ($gtf, $sam, $out_folder,$t, $i) = split ("\t", $lines_rute[$index]);
	my @sam_name = split ("\/", $sam);
	my @name_file = split ("\_", $sam_name[-1]);

	if (defined $lines_rute[$index]){
		system "htseq-count -s no -m intersection-nonempty -t $t -i $i $sam $gtf \> $out_folder/$name_file[0]\.htseq\n";
		print "Work $index is done!\n";
	}

	# Exit the thread
	threads->exit();

}


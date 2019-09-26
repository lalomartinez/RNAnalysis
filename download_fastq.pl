open (IN,"<",$ARGV[0]);
#$rute_SRA = "/media/lalo/9938e3a6-9b4c-487f-9a33-6c4e9448ad24/Lmajor_rnaseq/SRA_accession_1.txt";
$rute_SRA2 = (IN,"<",$ARGV[0]);  

while (<IN>){
	$rute_file_1 = "";
	$rute_file_2 = "";
	chomp $_;
	($rute_file_1, $rute_file_2) = split ("\t", $_);
	#print "$rute_file_1  $rute_file_2$rute_file_1  \n";

	#system "nohup parallel-fastq-dump --sra-id  $rute_file_1 --threads 10 --split-files  --gzip   --outdir $rute_file_2$rute_file_1  &"
	 system "nohup fastq-dump -A $rute_file_1  --split-files  --gzip   --outdir $rute_file_2$rute_file_1  &"
	}


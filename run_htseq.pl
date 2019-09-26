open (IN,"<",$ARGV[0]);

while (<IN>){
	chomp $_; 
	($gtf, $sam, $out_folder,$t, $i)= split("\t", $_);
	@sam_name=split("\/",$sam);
	@name_file=split("\_", $sam_name[-1]);

	system "htseq-count -s no -m intersection-nonempty -t $t -i $i $sam $gtf \> $out_folder$name_file[0]\.htseq\n";

}

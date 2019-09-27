
open FILE,'experiments.txt';
%ref=();
while ($line=<FILE>){
	chomp $line;
	@t=split '\.',$line;
	$ref{$t[0]}=0;
}


# download tissue information
@exps=keys %ref;

foreach $exp (@exps){
	if (-e "tissue_html/$exp.html"){}else{
		print "$exp\n";
		system "wget -O $exp.html http://www.ncbi.nlm.nih.gov/sra/?term=$exp";
		system "mv $exp.html tissue_html/";
	}
}


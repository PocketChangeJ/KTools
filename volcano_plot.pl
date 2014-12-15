#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

# integrade it to DE tools later

my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
DE_plot_volcano(\%options, \@ARGV);

sub DE_plot_volcano
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t DE_plot_volcano -p [column] -f [column]  input_file 

	-p column number of pvalue
	-f column number of fold change

';

	print $usage and exit unless defined $$files[0];
	my $input_file = $$files[0];
	die "[ERR]input file not exist\n" unless -s $input_file;
	
	print $usage and exit unless (defined $$options{'p'} && defined $$options{'f'});
	my $fcol = $$options{'f'};
	my $pcol = $$options{'p'};

	my $output_file = $input_file.".pdf";
	$output_file = $$options{'o'} if defined $$options{'o'};
	# check the column number and input files

	# generate R code
	my $r_code = qq'
data<-read.delim(\"$input_file\", header=T)
data<-data[,c(1,$fcol,$pcol)]
data<-data[!is.na(data[,3]),]
x<-data[,2]
y<-data[,3]
data1<-data[-log(data[,3],10)<5,]
x1<-data1[,2]
y1<-data1[,3]
pdf(\"$output_file\", width=8, height=10)
plot(log(x),  -log(y,10),  main="Volcano Plot",          pch=16, col=ifelse(((x>2 & y<0.05) | (x<0.5 & y<0.05)), "red", 200), cex=0.5, xlab="log2fc", ylab="-log10pvalue")
plot(log(x1), -log(y1,10), main="Filtered Volcano Plot", pch=16, col=ifelse(((x1>2 & y1<0.05) | (x1<0.5 & y1<0.05)), "red", 200), cex=0.5, xlab="log2fc", ylab="-log10pvalue")
dev.off()

';	

	# run R code and
	my $tmp = IO::File->new(">temp.R") || die "Can not open temp.R file $!\n";
	print $tmp $r_code;
	$tmp->close;
	system("R --no-save < temp.R") && die "Error at cmd R --no-save < temp.R\n";
}


#!/usr/bin/perl

open FILE, "<".$ARGV[0];

$flg=0;
while(<FILE>){
    if($_ =~ /\/\//){
	$_ =~ s/\/\/#/#/g;
    }
    for($i=1;$i<$#ARGV+1;$i++){
	if($_ =~/$ARGV[$i]/ && $_ =~/#define/ ){
	    $flg=1;
	}
    }
    if($flg==1){
	print "//", $_;	
    }else{
	print $_;
    }
    $flg=0;
}

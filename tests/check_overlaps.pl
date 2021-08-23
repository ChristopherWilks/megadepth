#!/usr/bin/env perl
use strict;
use warnings;

my $target_bg=shift;
my $coordsF=shift;
#e.g. wgEncodeHaibTfbsH1hescSp1Pcr1xRawRep1.bg.bz
open(IN,"<$coordsF");

my $total_sum = 0;
while(my $line1=<IN>) 
{
    chomp($line1) ; 
    my ($c,$s,$e,$n)=split(/\t/,$line1); 
    my $output=`tabix $target_bg $c:$s-$e`;
    my @output = split(/\n/,$output);
    if(scalar(@output) == 0)
    {
        print STDERR "no matches for $line1, skipping\n";
        next;
    }
    my $local_sum = 0;
    while(my $line2=shift(@output))
    { 
        chomp($line2); 
        my ($c2,$s2,$e2,$v)=split(/\t/,$line2);

        #left side of annotation overhangs
        if($s < $s2 && $e <= $e2)
        {
            $local_sum+= ($e - $s2)*$v;
        }
        #right side of annotation overhangs
        elsif($s >= $s2 && $e > $e2)
        {
            $local_sum+= ($e2 - $s)*$v;
        }
        #target either exactly matches or is strictly contained in annotation
        elsif($s <= $s2 && $e >= $e2)
        {
            $local_sum+= ($e2 - $s2)*$v;
        }
        #annotation either exactly matches or is strictly contained in target
        elsif($s >= $s2 && $e <= $e2)
        {
            $local_sum+= ($e - $s)*$v;
        }
        else
        {
            print STDERR "ran out of matching possibilities $c:$s-$e (annotation) vs. $c2:$s2-$e2 (target)\n";
            next;
        }
    }
    print "$line1\t$local_sum\n";
    $total_sum += $local_sum;
}
close(IN);
print "total_sum\t$total_sum\n";
            

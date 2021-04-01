#!/usr/bin/perl

use POSIX;
use List::Util qw(max);
use List::UtilsBy qw(max_by);
use List::Util 'sum';

$inputfile = $ARGV[0];
$treshold = $ARGV[1];
#$per = $ARGV[2];

open (BLAST,  "<$inputfile");
open (VIRUS, ">$inputfile.$treshold.contigs_virus.txt");
open (VIRUS_MAX, ">$inputfile.$treshold.contigs_virus_max.txt");
open (BACTERIA_MAX, ">$inputfile.$treshold.contigs_bacteria_max.txt");
open (FUNGI, ">$inputfile.$treshold.contigs_fungi_max.txt");
open (EUKA, ">$inputfile.$treshold.contigs_eukaryota_max.txt");
open (VIRDPLA, ">$inputfile.$treshold.contigs_virdplantae_max.txt");
open (UNKW, ">$inputfile.$treshold.contigs_unknown_max.txt");
open (METAZ, ">$inputfile.$treshold.contigs_metazoo_max.txt");
open (ARCHAEA, ">$inputfile.$treshold.contigs_archaea_max.txt");
open (ALLTAX, ">$inputfile.$treshold.contigs_all_tax.txt");

my @ncontigs;
my %hash_contigs;

my @tax_c;
my %hash_taxonomy;
my %hash_order;
my %hash_fam;
my %hash_genus;
my %hash_sp;

while(<BLAST>){
        chomp $_;
        @result_line = split ("\t", $_);
        $contign = $result_line[0];
        $eval = $result_line[6] + 0;
        #$treshold = 1e-1 + 0;
        next if ($treshold < $eval);
        #print $result_line[6]."\n";

        $hash_contigs{$contign} = {} if !$hash_contigs{$contign};
        $hash_taxonomy{$contign} = {} if !$hash_taxonomy{$contign};
        $hash_order{$contign} = {} if !$hash_order{$contign};
        $hash_fam{$contign} = {} if !$hash_fam{$contign};
        $hash_genus{$contign} = {} if !$hash_genus{$contign};
        $hash_sp{$contign} = {} if !$hash_sp{$contign};


        @$contign = () if !@$contign;
        @$taxtemp = () if !@$contign;
        @$ordertemp = () if !@$contign;
        @$famtemp = () if !@$contign;
        @$genustemp = () if !@$contign;
        @$sptemp = () if !@$contign;


        push @ncontigs, $contign if !@$contign;


        @tax_tree = split (" /\ ", $result_line[8]);

                ##checking spaces in the taxa
                $mycountarray = 0;
                foreach $a (@tax_tree){
                        $a =~ s/\s/-/g;
                        $tax_tree[$mycountarray] = $a;
                        $mycountarray++;
                }

        my $last3 = $tax_tree[3]."/".$tax_tree[4]."/".$tax_tree[5]."/".$tax_tree[6];

        push @$ordertemp, $tax_tree[3];
        push @$famtemp, $tax_tree[4];
        push @$genustemp, $tax_tree[5];
        push @$sptemp, $tax_tree[6];
        push @$taxtemp, $last3;
        push @$contign, $tax_tree[0];


        $hash_contigs{$contign} = "@$contign";
        $hash_taxonomy{$contign} = "@$taxtemp";
        $hash_order{$contign} = "@$ordertemp";
        $hash_fam{$contign} = "@$famtemp";
        $hash_genus{$contign} = "@$genustemp";
        $hash_sp{$contign} = "@$sptemp";
        $old_id = $contign;
}

while ( ($k,$v) = each %hash_contigs ) {
        print $k.":\t";
        print ALLTAX $k."\t";

        @array = split (' ', $v);
        @tax_sp  = split(' ', $hash_taxonomy{$k}); ### the hash that contain the tax ids  are entereed with the key that is the same key used, so the contig ID.
        my %count_tax;                     ### hash to bulk the counts
        $count_tax{$_}++ foreach @tax_sp;  ## counts of every tax level on the hash %count_tax
        my $highest_tax = max_by { $count_tax{$_} } keys %count_tax; ## selecting the most abundant element per hash, that is teh number of species per contig
        my $len_tax = keys %count_tax;   ## counting the number of taxas per contig
        my @sorted_tax = sort { $count_tax{$b} <=> $count_tax{$a} } keys %count_tax;

        @order_dump = split (' ', $hash_order{$k});
        @fam_dump = split (' ', $hash_fam{$k});
        @genus_dump = split (' ', $hash_genus{$k});
        @sp_dump = split (' ', $hash_sp{$k});

        my %count_order;
        $count_order{$_}++ foreach @order_dump;
        my $highest_order = max_by { $count_order{$_} } keys %count_order;

        my %count_fam;
        $count_fam{$_}++ foreach @fam_dump;
        my $highest_fam = max_by { $count_fam{$_} } keys %count_fam;


        my %count_genus;
        $count_genus{$_}++ foreach @genus_dump;
        my $highest_genus = max_by { $count_genus{$_} } keys %count_genus;

        my %count_sp;
        $count_sp{$_}++ foreach @sp_dump;
        my $highest_sp = max_by { $count_sp{$_} } keys %count_sp;

        $len = @array;
        my %count;
        $count{$_}++ foreach @array;
        print $len."\t";

        my %classif =(
                "Metazoa" => 0,
                "Archaea" => 0,
                "Bacteria" => 0,
                "Eukaryota" => 0,
                "Viridiplantae" => 0,
                "Fungi" => 0,
                "Viruses" => 0,
                "unknown" => 0,
                "kingdom" => 0,
        );

                while (my ($key, $value) = each(%count)) {
                $per = sprintf("%.2f", $value/$len);
                #print "$key:$value:$per\t";
                print VIRUS "$k\t$key\t$value\t$per\n" if $key eq "Viruses" && $per >= 0.6;
                $classif{$key} = $value;
                }

        my $highest = max_by { $classif{$_} } keys %classif;
        my $value_count = sum values %classif;
        #print "the highest is: $highest\n";

                foreach my $name (keys %classif) {
                #printf "%-8s %s\t", $name, $classif{$name};
                $per = sprintf("%.2f", $classif{$name}/$value_count);
                #print "$per\t";
                printf "%-8s %s\t", $name, $classif{$name};
                print $per."\t";

                printf ALLTAX "%-8s %s\t", $name, $classif{$name};
                print ALLTAX $per."\t";




                print VIRUS_MAX "$k\t$name\t$classif{$name}\t$per\t$highest_order\t$highest_fam\t$highest_genus\t$highest_sp\ttax_cat:$len_tax\t@sorted_tax[0]:$count_tax{$sorted_tax[0]}\t@sorted_tax[1]:$count_tax{$sorted_tax[1]}\n" if $name eq $highest && $name eq "Viruses";
                print BACTERIA_MAX "$k\t$name\t$classif{$name}\t$per\t$highest_order\t$highest_fam\t$highest_genus\t$highest_sp\ttax_cat:$len_tax\t@sorted_tax[0]:$count_tax{$sorted_tax[0]}\t@sorted_tax[1]:$count_tax{$sorted_tax[1]}\n" if $name eq $highest && $name eq "Bacteria";
                print FUNGI  "$k\t$name\t$classif{$name}\t$per\t$highest_order\t$highest_fam\t$highest_genus\t$highest_sp\ttax_cat:$len_tax\t@sorted_tax[0]:$count_tax{$sorted_tax[0]}\t@sorted_tax[1]:$count_tax{$sorted_tax[1]}\n" if $name eq $highest && $name eq "Fungi";
                print EUKA "$k\t$name\t$classif{$name}\t$per\t$highest_order\t$highest_fam\t$highest_genus\t$highest_sp\ttax_cat:$len_tax\t@sorted_tax[0]:$count_tax{$sorted_tax[0]}\t@sorted_tax[1]:$count_tax{$sorted_tax[1]}\n" if $name eq $highest && $name eq "Eukaryota";
                print UNKW "$k\t$name\t$classif{$name}\t$per\t$highest_order\t$highest_fam\t$highest_genus\t$highest_sp\ttax_cat:$len_tax\t@sorted_tax[0]:$count_tax{$sorted_tax[0]}\t@sorted_tax[1]:$count_tax{$sorted_tax[1]}\n" if $name eq $highest && $name eq "unknown";
                print METAZ "$k\t$name\t$classif{$name}\t$per\t$highest_order\t$highest_fam\t$highest_genus\t$highest_sp\ttax_cat:$len_tax\t@sorted_tax[0]:$count_tax{$sorted_tax[0]}\t@sorted_tax[1]:$count_tax{$sorted_tax[1]}\n" if $name eq $highest && $name eq "Metazoa";
                print VIRDPLA "$k\t$name\t$classif{$name}\t$per\t$highest_order\t$highest_fam\t$highest_genus\t$highest_sp\ttax_cat:$len_tax\t@sorted_tax[0]:$count_tax{$sorted_tax[0]}\t@sorted_tax[1]:$count_tax{$sorted_tax[1]}\n" if $name eq $highest && $name eq "Viridiplantae";
                print ARCHAEA "$k\t$name\t$classif{$name}\t$per\t$highest_order\t$highest_fam\t$highest_genus\t$highest_sp\ttax_cat:$len_tax\t@sorted_tax[0]:$count_tax{$sorted_tax[0]}\t@sorted_tax[1]:$count_tax{$sorted_tax[1]}\n" if $name eq $highest && $name eq "Archaea";


                }
        print "$highest_order\t$highest_fam\t$highest_genus\t$highest_sp\ttax_cat:$len_tax\t@sorted_tax[0]:$count_tax{$sorted_tax[0]}\t@sorted_tax[1]:$count_tax{$sorted_tax[1]}\n";

        print ALLTAX "$highest_order\t$highest_fam\t$highest_genus\t$highest_sp\n";

}

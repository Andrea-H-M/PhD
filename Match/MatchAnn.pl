#!/usr/bin/perl -w
use strict;

my $reporte = "Trinity.isoform.counts.matrix";
open (FILE1, "$reporte");
my @Count1 = <FILE1>;
chomp @Count1;
my $header = shift @Count1;
my $c1 = scalar @Count1;
my @segmentos4 = split ("\t", $header);
        chomp  @segmentos4;
        my $header1corregido = join(",",@segmentos4);

my $reporte2 = "TablaVanilla_vFinal.csv";
open (FILE2, "$reporte2");
my @Count2 = <FILE2>;
chomp @Count2;
my $header2 = shift @Count2;

my $c2 = scalar @Count2; 

my $datafinal = "Anotada.".$reporte.".csv";
open (OUTPUT, ">$datafinal");
my $headerF= "ID".",".$header1corregido.",".$header2;
print OUTPUT "$headerF\n";
my $Info = "";

for (my $i=0; $i < $c1; $i++){ 
        my $renglon = $Count1[$i];
        chomp $renglon;
        my @segmentos = split ("\t", $renglon);
        chomp  @segmentos;
        my $mew = join(",",@segmentos);
        my $gene = $segmentos[0];
        $gene =~ s/"//g;
        $gene =~ s/\s//g;
     
for (my $ii=0; $ii < $c2; $ii++){ 
        my $renglon2 = $Count2[$ii];
        chomp $renglon2;
        my @segmentos1 = split (",", $renglon2);
        chomp  @segmentos1;
        my $gene2 = $segmentos1[0];
        $gene2 =~ s/\s//g;
        #print "$gene  eq  $gene2\n";
        if($gene  eq  $gene2){
         #print "$gene    $gene2\n";
              $Info = $renglon2;
         last;
            }

        }
 my $linea = $mew.",".$Info;
              print OUTPUT "$linea\n";
            $Info = "";
              
}
close OUTPUT;
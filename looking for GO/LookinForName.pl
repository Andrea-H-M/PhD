#!/usr/bin/perl -w
use strict;
my $archivoAnotada="TablaMaestraVanilla.csv";
open(INPUT1,"$archivoAnotada"); 
my @contenidoAnotada = <INPUT1>; 
chomp @contenidoAnotada;
my $longitudarreglo1=scalar@contenidoAnotada;
print "Imprimiendo la longitud del arreglo de la tabla $longitudarreglo1\n";

my $archivoSalida="Auxinasv2.csv";
open(OUTPUT1,">$archivoSalida"); 


for(my $i=0;$i<$longitudarreglo1;$i++){
    my $FilaAnotada =$contenidoAnotada[$i];
    my @segmentos =split(",",$FilaAnotada);
        chomp@segmentos;
        my$BP= $segmentos[25];
        my$MF= $segmentos[27];
        my$CC= $segmentos[29];
        
        if ($BP =~ "auxin" or $BP =~ "auxin" or $BP =~ "auxin"){
            print "$FilaAnotada\n";
            print OUTPUT1 "$FilaAnotada\n";
        }
}
#!/usr/bin/perl -w
#
#
# 
#  script: calculate_number_hits_per_bin_EUseq_v1.pl
#  
#  published in: Macheret M., and Halazonetis, T.D.
#                Intragenic origins due to short G1 phases underlie oncogene-induced DNA replication stress.
#                Nature 2018.
#
#  disclaimers: no claims are made that this script performs any specific task
#               not responsible for any damage/loss caused by using this script
#               script is to be used only for research purposes
#               
#  prerequisites: good knowledge of linux and perl are required to run this script
#                 no support can be provided to troubleshoot this script
#                 if you have difficulties, please contact your local linux/perl expert
#
#  IMPORTANT: works with the masked human genome assembly: GRCh37/hg19
#             the sam files must be aligned to the above assembly 
#
#  copyright (2017): Thanos D. Halazonetis, University of Geneva
#
#  
#
#
#  usage: calculate_number_hits_per_bin_EUseq_v1.pl    bin_size   quality   read_length   sam_file   
#         calculate_number_hits_per_bin_EUseq_v1.pl    10000      60        100           EUseq_file.sam
#
#
#
print "\n";
##
$bin_size=$ARGV[0];
$quality_score=$ARGV[1];
$read_length=join("",$ARGV[2],"M");
$sam_repl_File=$ARGV[3];
#
#
#
# select appropriate adjust file, by commenting out the inappropriate file names
#
#$sam_adjust_File="RPE1_adjust__bin-size_10000_chr1-X_adjbin_0b.csv";
#$sam_adjust_File="HELA_adjust__bin-size_10000_chr1-X_adjbin_0b.csv";
$sam_adjust_File="U2OS_adjust__bin-size_10000_chr1-X_adjbin_0b.csv";
#
# checks that input bin_size matches bin_size of adjust file
#
@temp_bin_size_1=split("bin-size_",$sam_adjust_File);
@temp_bin_size_2=split("_",$temp_bin_size_1[1]);
if ($bin_size ne $temp_bin_size_2[0]) {
    print "\n\nusage:  calculate_number_hits_per_bin_EUseq_v1.pl    bin_size   quality   read_length   sam_file\n";
    print "usage:  calculate_number_hits_per_bin_EUseq_v1.pl    10000      60        100           EUseq_file.sam\n";
    print "input bin-size and bin-size of adjust file do not match\n\n";
    exit;
}
#
#
#
@file_name=split(".sam",$sam_repl_File);
print "\nSam file:    $sam_repl_File\n";
print "Adjust file: $sam_adjust_File\n";
#
# exclude values for adjust hit counts
#
if ($bin_size==10000) {
   $exclude_low=25;     #these values are good for bin_size of 10000
   $exclude_high=10000; #these values are good for bin_size of 10000
} elsif ($bin_size==1000) {
   $exclude_low=10;     #these values are good for bin_size of 1000
   $exclude_high=2000;  #these values are good for bin_size of 1000
} else {
   $exclude_low=1;         #limit 0
   $exclude_high=1000000;  #no limit
}
#
#
#
#
if ($#ARGV<3) {
    print "\n\nusage:  calculate_number_hits_per_bin_EUseq_v1.pl    bin_size   quality   read_length   sam_file\n";
    print "usage:  calculate_number_hits_per_bin_EUseq_v1.pl    10000      60        100           EUseq_file.sam\n";
    print "        incorrect number of variables in command line\n\n";
    exit;
}
#
if (!-r "$sam_repl_File" ) {
    print "\n\nusage:  calculate_number_hits_per_bin_EUseq_v1.pl    bin_size   quality   read_length   sam_file\n";
    print "usage:  calculate_number_hits_per_bin_EUseq_v1.pl    10000      60        100           EUseq_file.sam\n";
    print "        input sam file non-existent or not readable\n\n";
    exit;
}
#
if (!-r "$sam_adjust_File" ) {
    print "\n\nusage:  calculate_number_hits_per_bin_EUseq_v1.pl    bin_size   quality   read_length   sam_file\n";
    print "usage:  calculate_number_hits_per_bin_EUseq_v1.pl    10000      60        100           EUseq_file.sam\n";
    print "        input adjust file non-existent or not readable\n\n";
    exit;
}
#
# read adjust hits
#
$chr_size=250000000;
$number_bins=$chr_size/$bin_size;
$loopb=0;     #### zero base
$number_OK_adjust_lines=0;
$number_notOK_adjust_lines=0;
#
open FILE, $sam_adjust_File or die "Can't open input file: $:\n";
while ( <FILE> ) {
       chomp; #remove newline
       @line = split(/,/,$_); # split each line
       if ($#line==22) {
          for ($loopc=1; $loopc<=23; $loopc++){
             $adjust_hits[$loopc][$loopb]=$line[$loopc-1];
          }
          $loopb++;
          $number_OK_adjust_lines++;
       } else {
          $number_notOK_adjust_lines++;
       }
}
close FILE;
#
print "\nBin size:      $bin_size\n";
print "Quality score: $quality_score\n";
print "Read length:   $read_length\n\n";
#
print "Number of bins: $number_bins\n";
print "Number of OK adjust lines (should be equal to number of bins or one more): $number_OK_adjust_lines\n";
print "Number of not OK adjust lines (should be zero or one):                     $number_notOK_adjust_lines\n\n";
#
# create array for bin hits
#
$chr_size=250000000;
$number_bins=$chr_size/$bin_size;
##
for ($loopc=1; $loopc<=23; $loopc++){
   for ($loopb=0; $loopb<=$number_bins; $loopb++){   
       $bin_hits_f[$loopc][$loopb]=0;
       $bin_hits_r[$loopc][$loopb]=0;
   }
}
#
#
# read bin hits
#
#
$qualOK = 0;
$qualOK_f = 0;
$qualOK_r = 0;
$qualnotOK = 0;
#
open FILE, $sam_repl_File or die "Can't open input file: $:\n";
while ( <FILE> ) {
    chomp; #remove newline
    @line = split(/\t/,$_); # split each line
    if ($#line>10) {
        $chromosome=$line[2];
        @temp_chr=split("chr",$chromosome);
        if ($#temp_chr==0) { $chr_number=26 } else { $chr_number=$temp_chr[1] }
        if ($chr_number eq "1") {$chr_number=1}
        elsif ($chr_number eq "2") {$chr_number=2}
        elsif ($chr_number eq "3") {$chr_number=3}
        elsif ($chr_number eq "4") {$chr_number=4}
        elsif ($chr_number eq "5") {$chr_number=5}
        elsif ($chr_number eq "6") {$chr_number=6}
        elsif ($chr_number eq "7") {$chr_number=7}
        elsif ($chr_number eq "8") {$chr_number=8}
        elsif ($chr_number eq "9") {$chr_number=9}
        elsif ($chr_number eq "10") {$chr_number=10}
        elsif ($chr_number eq "11") {$chr_number=11}
        elsif ($chr_number eq "12") {$chr_number=12}
        elsif ($chr_number eq "13") {$chr_number=13}
        elsif ($chr_number eq "14") {$chr_number=14}
        elsif ($chr_number eq "15") {$chr_number=15}
        elsif ($chr_number eq "16") {$chr_number=16}
        elsif ($chr_number eq "17") {$chr_number=17}
        elsif ($chr_number eq "18") {$chr_number=18}
        elsif ($chr_number eq "19") {$chr_number=19}
        elsif ($chr_number eq "20") {$chr_number=20}
        elsif ($chr_number eq "21") {$chr_number=21}
        elsif ($chr_number eq "22") {$chr_number=22}
        elsif ($chr_number eq "X") {$chr_number=23}
        elsif ($chr_number eq "Y") {$chr_number=24}
        elsif ($chr_number eq "M") {$chr_number=25}
        else {$chr_number=26}
        #if ($line[4]>=$quality_score && ($line[5] eq "99M" || $line[5] eq "100M") && $chr_number<24 && $chromosome!~/mchr/) {  
              #use above if sam file is merge of sam files with different read lengths
        if ($line[4]>=$quality_score && $line[5] eq $read_length && $chr_number<24 && $chromosome!~/mchr/) {  #some sam files may contain mouse sequences for calibration
             $qualOK ++;
             $bin_allocate=$line[3]/$bin_size;
             $bin_allocate_int=int($bin_allocate);
             if ($line[1]==16) { 
                 $bin_hits_f[$chr_number][$bin_allocate_int] ++;
                 $qualOK_f ++;
             } elsif ($line[1]==0) {
                 $bin_hits_r[$chr_number][$bin_allocate_int] ++;
                 $qualOK_r ++;
             } 
        } else {
             $qualnotOK ++;
        }
    }
}
close FILE;
$qualOK_fr=$qualOK_f+$qualOK_r;
print "Sam file $sam_repl_File\n";
print "quality forward OK $qualOK_f\n";
print "quality reverse OK $qualOK_r\n";
print "quality OK $qualOK\n";
print "validation: f+r: $qualOK_fr\n";
print "quality not OK $qualnotOK\n\n";
#
#
# eliminate bin and adjust hits outside low and high exclude range
#
#
for ($loopc=1; $loopc<=23; $loopc++){
   for ($loopb=0; $loopb<=$number_bins; $loopb++){
       if ($adjust_hits[$loopc][$loopb]<$exclude_low || $adjust_hits[$loopc][$loopb]>$exclude_high) {
           $bin_hits_f[$loopc][$loopb]=0;
           $bin_hits_r[$loopc][$loopb]=0;
           $adjust_hits[$loopc][$loopb]=0;
       }
   }
}
#
#
# calculate total hits
#
#
$bin_total_hits_f=0;
$bin_total_hits_r=0;
$adjust_total_hits=0;
for ($loopc=1; $loopc<=23; $loopc++){
   for ($loopb=0; $loopb<=$number_bins; $loopb++){
       $bin_total_hits_f=$bin_total_hits_f+$bin_hits_f[$loopc][$loopb];
       $bin_total_hits_r=$bin_total_hits_r+$bin_hits_r[$loopc][$loopb];
       $adjust_total_hits=$adjust_total_hits+$adjust_hits[$loopc][$loopb];
   }
}
$Acorr_factor=$adjust_total_hits*2/($bin_total_hits_f+$bin_total_hits_r);
#
#
print "Total sample hits forward:  $bin_total_hits_f\n";
print "Total sample hits reverse:  $bin_total_hits_r\n";
print "Adjust total hits:          $adjust_total_hits\n";
print "Correction factor:          $Acorr_factor\n";
#
#
# adjust bin hits and calculate total adjbin hits
#
#
$adjbin_total_hits_f=0;
$adjbin_total_hits_r=0;
for ($loopc=1; $loopc<=23; $loopc++){
   for ($loopb=0; $loopb<=$number_bins; $loopb++){
       if ($adjust_hits[$loopc][$loopb]>0) {
           $adjbin_hits_f[$loopc][$loopb]=int($bin_hits_f[$loopc][$loopb]*$Acorr_factor*1000/$adjust_hits[$loopc][$loopb]+0.5);
           $adjbin_total_hits_f=$adjbin_total_hits_f+$adjbin_hits_f[$loopc][$loopb];
           $adjbin_hits_r[$loopc][$loopb]=int($bin_hits_r[$loopc][$loopb]*$Acorr_factor*1000/$adjust_hits[$loopc][$loopb]+0.5);
           $adjbin_total_hits_r=$adjbin_total_hits_r+$adjbin_hits_r[$loopc][$loopb];
       } else {
           $adjbin_hits_f[$loopc][$loopb]=0;
           $adjbin_hits_r[$loopc][$loopb]=0;
       }
   }
}
print "Total adjusted sample hits for: $adjbin_total_hits_f\n";
print "Total adjusted sample hits rev: $adjbin_total_hits_r\n\n";
#
#
# write files
#
$prepare_file_name=join("",$file_name[0],"--EUfor");
$result_File = join("_",$prepare_file_name,"bin-size",$bin_size,"quality",$quality_score,"chr1-X_adjbin_qual_counts_EU_0b.txt");
#
open RESULT, "> $result_File";
#
printf RESULT "\n";
printf RESULT "file $sam_repl_File\n";
printf RESULT "quality OK $qualOK_f\n";
printf RESULT "quality not OK $qualnotOK\n\n";
#
printf RESULT "\n";
printf RESULT "Total sample hits: $bin_total_hits_f\n";
printf RESULT "Adjust total hits: $adjust_total_hits\n";
printf RESULT "Correction factor: $Acorr_factor\n";
printf RESULT "Total adjusted sample hits: $adjbin_total_hits_f\n\n";
#
close RESULT;
#
#
$prepare_file_name=join("",$file_name[0],"--EUrev");
$result_File = join("_",$prepare_file_name,"bin-size",$bin_size,"quality",$quality_score,"chr1-X_adjbin_qual_counts_EU_0b.txt");
#
open RESULT, "> $result_File";
#
printf RESULT "\n";
printf RESULT "file $sam_repl_File\n";
printf RESULT "quality OK $qualOK_r\n";
printf RESULT "quality not OK $qualnotOK\n\n";
#
printf RESULT "\n";
printf RESULT "Total sample hits: $bin_total_hits_r\n";
printf RESULT "Adjust total hits: $adjust_total_hits\n";
printf RESULT "Correction factor: $Acorr_factor\n";
printf RESULT "Total adjusted sample hits: $adjbin_total_hits_r\n\n";
#
close RESULT;
#
#
#
$prepare_file_name=join("",$file_name[0],"--EUfor");
$result_File = join("_",$prepare_file_name,"bin-size",$bin_size,"quality",$quality_score,"chr1-X_adjbin_EU_0b.csv");
#
open RESULT, "> $result_File";
#
for ($loopb=0; $loopb<=$number_bins; $loopb++){       #### zero base
   for ($loopc=1; $loopc<=23; $loopc++){
       printf RESULT "$bin_hits_f[$loopc][$loopb],$adjbin_hits_f[$loopc][$loopb],$adjust_hits[$loopc][$loopb],,";
   }
   printf RESULT "\n";
}
#
close RESULT;
#
#
$prepare_file_name=join("",$file_name[0],"--EUrev");
$result_File = join("_",$prepare_file_name,"bin-size",$bin_size,"quality",$quality_score,"chr1-X_adjbin_EU_0b.csv");
#
open RESULT, "> $result_File";
#
for ($loopb=0; $loopb<=$number_bins; $loopb++){       #### zero base
   for ($loopc=1; $loopc<=23; $loopc++){
       printf RESULT "$bin_hits_r[$loopc][$loopb],$adjbin_hits_r[$loopc][$loopb],$adjust_hits[$loopc][$loopb],,";
   }
   printf RESULT "\n";
}
#
close RESULT;
#
#
#
exit;
#
#

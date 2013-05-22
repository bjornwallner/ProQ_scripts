#!/usr/bin/perl
#
#       Author:         Arlo Randall
#       File:           blast2index.pl
#       Description:    parse blast results for best match, make index file of the form:
#               
#       First line:     (basic info on query)
#                       query name,query length
#       Additional lines:       (detailed info on one segment of alignment)
#
#       $aligned,$gaps,$eval,$ident,$posi,$subj,$subj_len,$subj_start,$subj_end,$query_start,$query_end
#       0       $aligned        # number of amino acids aligned
#       1       $gaps
#       2       $eval
#       3       $ident
#       4       $posi
#       5       $subj
#       6       $subj_len
#       7       $subj_start
#       8       $subj_end
#       9       $query_start
#       10      $query_end
#*************************************************************************    
    
use strict;

my $infile = $ARGV[0];
#my $infile = "testout.txt";
my @text;

my @unsorted_index;
my @sorted_index;
my @foo;

my @one_result;
my $query; # these vars correspond to DB fields
my $subj;
my $query_len;
my $subj_len;
my $eval;
my $subj_start;
my $subj_end;
my $query_start;
my $query_end;
my $ident;
my $posi = 0;
my $gaps = 0;
my $p_align = 0;
my $aligned;    # number of residues aligned in blast alignment
my $align_length;

my $temp;
my @row_tokens;
my $done;
my $shorter;

open(INFILE, $infile);
@text = <INFILE>;
close(INFILE);

# parse query
@row_tokens = split(/\s+/,$text[8]);
$query = $row_tokens[1];

# parse query length
@row_tokens = split(/\s+|\(/,$text[9]);
$query_len = $row_tokens[2];


my $i = 0;
my $j = 0;
while($i < $#text)
{
        # new results start with '>' or 'Score' if multiple records
        if(($text[$i] =~ /^>/) or ($text[$i] =~ /Score\s+=/))
        {       
                if($text[$i] =~ /^>/)
                {
                        @row_tokens = split(/\s+|>/,$text[$i]);    
                        $subj = $row_tokens[1];
                
                        # parse subject length
                        if($text[$i+1] =~ /Length/)
                        {       
                                @row_tokens = split(/\s+/,$text[$i+1]); 
                                $subj_len = $row_tokens[3];                     
                        }
                        elsif($text[$i+2] =~ /Length/)
                        {
                                @row_tokens = split(/\s+/,$text[$i+2]); 
                                $subj_len = $row_tokens[3];
                                $i++;
                        }
                        $i++;
                        $i++;
                        $i++;
                }
                
                # parse e-value
                @row_tokens = split(/\s+/,$text[$i]);
                $eval = $row_tokens[8];

                # parse identities
                @row_tokens = split(/\(|%/,$text[$i+1]);
                $ident = $row_tokens[1];
                
                # parse positives
                $posi = $row_tokens[3];
                if ($text[$i+1] =~ /Gaps/)
                {       $gaps = $row_tokens[5];         }
                else
                {       $gaps = 0;                      }      
                if ($gaps !~ /[0-9]+/)
                {       $gaps = 0;      }
              
                # set p_align
                @row_tokens = split(/\s+|\//,$text[$i+1]);
                $aligned = $row_tokens[4];
                
                if ($query_len < $subj_len)
                {       $shorter = $query_len;       }
                else
                {       $shorter = $subj_len;        }

                # DEBUG
                if ($shorter > 0)
                {
                        $p_align = $aligned / $shorter;
                        $p_align = substr($p_align, 0, 5);
                }       
                else
                {       $p_align = 0;                }
                
                
                
                if($p_align > 1)
                {       $p_align = 1;   }
                                
                # get start position of query
                @row_tokens = split(/\s+/, $text[$i+3]);
                $query_start = $row_tokens[1];
                $query_end = $row_tokens[3];

                # get start position of subject
                @row_tokens = split(/\s+/, $text[$i+5]);
                $subj_start = $row_tokens[1];
                $subj_end = $row_tokens[3];

                # different length alignments require looking for last rows for end pos.
                $done = 0; # flag
                $j = $i + 2;      # start just before first query line

                while(!$done)
                {
                        if($text[$j] =~ /Query/)
                        {
                                @row_tokens = split(/\s+/, $text[$j]);
                                $query_end = $row_tokens[3];
                        }
                        elsif($text[$j] =~ /Sbjct/)
                        {
                                @row_tokens = split(/\s+/, $text[$j]);
                                $subj_end = $row_tokens[3];
                        }

                        # finished if 
                        if((($text[$j] =~ /Lambda/) or ($text[$j] =~/^>/)) or ($text[$j] =~ /Score\s+=/ ))
                        {       $done = 1;        }
        
                        $j++;
                }
                
        # put in top of file basic info in query: $query, $query_len
        @one_result = 1;
        shift(@one_result);
        push (@one_result, $aligned);
        push (@one_result, $gaps);
        push (@one_result, $eval);
        push (@one_result, $ident);
        push (@one_result, $posi);
        push (@one_result, $subj);
        push (@one_result, $subj_len);
        push (@one_result, $subj_start);
        push (@one_result, $subj_end);
        push (@one_result, $query_start);
        push (@one_result, $query_end);
        
        # build data structure and sort by aligned, then display when finished!
        push (@unsorted_index,[@one_result]); 
        
        }

$i++;
}


# sort records by length of alignment, then output
@sorted_index = sort { $b->[0] <=> $a->[0] } @unsorted_index;
    
# display first line    
print "$query $query_len \n";

# display lines w/ aligment indices    
for $i (0..$#sorted_index)
{
        print $sorted_index[$i][0]."\t";
        print $sorted_index[$i][1]."\t";
        print $sorted_index[$i][2]."\t";
        print $sorted_index[$i][3]."\t";
        print $sorted_index[$i][4]."\t";
        print $sorted_index[$i][5]."\t";
        print $sorted_index[$i][6]."\t";
        print $sorted_index[$i][7]."\t";
        print $sorted_index[$i][8]."\t";
        print $sorted_index[$i][9]."\t";
        print $sorted_index[$i][10]."\n";
}

#!/usr/bin/perl
use strict;
use warnings;
use LWP::UserAgent;
use Tk;
use Bio::SeqIO;

my ($aligned_seq1, $aligned_seq2, $score);
# Create main window
my $mw = MainWindow->new;

$mw->title("Main Window");

$mw->overrideredirect(1);
$mw->geometry($mw->screenwidth . "x" . $mw->screenheight . "+0+0");
$mw->configure(-borderwidth => 15, -relief => 'groove',-background => 'light blue');

my $right_frame = $mw->Frame(-borderwidth => 10, -relief => 'solid',-background => 'cyan');
$right_frame->place(-x => 890, -y => 15, -width => 600, -height => 800);

my $canvasrf = $right_frame->Canvas(
    -background => 'white',-highlightthickness => 0);
$canvasrf->place(-x => 0, -y => 0, -width => 580, -height => 780);

my $textshow = $canvasrf->Scrolled('Text',
    -background =>"powder blue",
    -foreground => "black",
    -scrollbars => 'e',
    -wrap => 'word',
    -highlightthickness => 1
)->pack(-expand => 1, -fill => 'both');

my $text1 = $mw->Scrolled('Text',
    -background =>"powder blue",
    -foreground => "black",
    -scrollbars => 'e',
    -wrap => 'word',
    -highlightthickness => 1
)->place(-x => 15, -y => 15, -width => 850, -height => 300);
$text1->insert('end',"sequence 1 : paste the sequence in FASTA format or paste the sequence id");

my $text2 = $mw->Scrolled('Text',
    -background =>"powder blue",
    -foreground => "black",
    -scrollbars => 'e',
    -wrap => 'word',
    -highlightthickness => 1
)->place(-x => 15, -y => 330, -width =>850, -height => 300);
$text2->insert('end',"sequence 2 : paste the sequence in FASTA format or paste the sequence id");

my $localbutton = $mw->Button(
    -text    => "RUN LOCAL ALIGNMENT",
    -command => sub { runlocalalignment() },
    -foreground => "black",
    -background => "cyan",
    -font => '-family {Times New Roman} -weight bold -size 14',
    -height => 2.8,
    -width => 15);
$localbutton->place(-x => 18, -y => 690, -width => 270, -height => 80);

my $globalbutton = $mw->Button(
    -text    => "RUN GLOBAL ALIGNMENT",
    -command => sub { runglobalalignment() },
    -foreground => "black",
    -background => "cyan",
    -font => '-family {Times New Roman} -weight bold -size 14',
    -height => 2.8,
    -width => 15);
$globalbutton->place(-x => 303, -y => 690, -width => 270, -height => 80);

my $exitbutton = $mw->Button(
    -text    => "EXIT",
    -command => sub { exit_confirmation() },
    -foreground => "black",
    -background => "cyan",
    -font => '-family {Times New Roman} -weight bold -size 14',
    -height => 2.8,
    -width => 15);
$exitbutton->place(-x => 588, -y => 690, -width => 270, -height => 80);

sub runglobalalignment
{
$textshow->delete('1.0', 'end');
my $seq1 = $text1->get("1.0", "end-1c");
my $seq2 = $text2->get("1.0", "end-1c");
$seq1=~ s/^>.*\n//mg;
$seq2=~ s/^>.*\n//mg;

if (contains_digit_underscore_or_dot($seq1)) {
    $seq1=retrieve_ncbi_sequence($seq1);
    $seq1=~ s/^>.*\n//mg;

}

if (contains_digit_underscore_or_dot($seq2)) {
    $seq2=retrieve_ncbi_sequence($seq2);
    $seq2=~ s/^>.*\n//mg;
}

my $gap_penalty = -2;
my ($aligned_seq1, $aligned_seq2, $alignment_score) = perform_sequence_alignment($seq1, $seq2, $gap_penalty);
print_sequences_with_line_length($aligned_seq1, $aligned_seq2, $alignment_score, 60);
}

sub perform_sequence_alignment {
    my ($seq1, $seq2, $gap_penalty) = @_;
    my @matrix;
    $matrix[0][0] = 0;
    for my $i (1 .. length($seq1)) {
        $matrix[$i][0] = $matrix[$i-1][0] + $gap_penalty;
    }
    for my $j (1 .. length($seq2)) {
        $matrix[0][$j] = $matrix[0][$j-1] + $gap_penalty;
    }
    for my $i (1 .. length($seq1)) {
        for my $j (1 .. length($seq2)) {
            my $match = substr($seq1, $i-1, 1) eq substr($seq2, $j-1, 1) ? 1 : 0;
            my $diagonal_score = $matrix[$i-1][$j-1] + $match;
            my $left_score = $matrix[$i][$j-1] + $gap_penalty;
            my $up_score = $matrix[$i-1][$j] + $gap_penalty;
            $matrix[$i][$j] = maximum($diagonal_score, $left_score, $up_score);
        }
    }

    return traceback($seq1, $seq2, $gap_penalty, @matrix);
}

sub maximum {
    my $max = shift;
    for (@_) {
        $max = $_ if $_ > $max;
    }
    return $max;
}

sub traceback {
    my ($seq1, $seq2, $gap_penalty, @matrix) = @_;

    my $i = length($seq1);
    my $j = length($seq2);

    my $aligned_seq1 = "";
    my $aligned_seq2 = "";

    my $alignment_score = $matrix[$i][$j];

    while ($i > 0 || $j > 0) {
        my $current_score = $matrix[$i][$j];
        my $diagonal_score = $i > 0 && $j > 0 ? $matrix[$i-1][$j-1] : 0;
        my $left_score = $j > 0 ? $matrix[$i][$j-1] : 0;
        my $up_score = $i > 0 ? $matrix[$i-1][$j] : 0;

        if ($current_score == $diagonal_score + (substr($seq1, $i-1, 1) eq substr($seq2, $j-1, 1) ? 1 : 0)) {
            $aligned_seq1 = substr($seq1, $i-1, 1) . $aligned_seq1;
            $aligned_seq2 = substr($seq2, $j-1, 1) . $aligned_seq2;
            $i--;
            $j--;
        } elsif ($current_score == $left_score + $gap_penalty) {
            $aligned_seq1 = "-" . $aligned_seq1;
            $aligned_seq2 = substr($seq2, $j-1, 1) . $aligned_seq2;
            $j--;
        } else {
            $aligned_seq1 = substr($seq1, $i-1, 1) . $aligned_seq1;
            $aligned_seq2 = "-" . $aligned_seq2;
            $i--;
        }
    }

    return ($aligned_seq1, $aligned_seq2, $alignment_score);
}

sub print_sequences_with_line_length {
    my ($seq1, $seq2, $alignment_score, $line_length) = @_;

    $textshow->insert('end', "GLOAL ALIGNMENT SCORE : $alignment_score\n\n");
    my $alignment_line = "";

    while (length($seq1) > 0) {
        my $print_length = length($seq1) > $line_length ? $line_length : length($seq1);
        my $substr_seq1 = substr($seq1, 0, $print_length);
        my $substr_seq2 = substr($seq2, 0, $print_length);

        $alignment_line = join("", map { substr($substr_seq1, $_, 1) eq substr($substr_seq2, $_, 1) ? "|" : " " } (0 .. $print_length-1));

        $textshow->insert('end', "Sequence 1: $substr_seq1\n");
        $textshow->insert('end', "            $alignment_line\n");
        $textshow->insert('end', "Sequence 2: $substr_seq2\n\n");
        $textshow->insert('end', "\n");

        $seq1 = substr($seq1, $print_length);
        $seq2 = substr($seq2, $print_length);
    }
}

sub runlocalalignment{
$textshow->delete('1.0', 'end');
my $sequence1 = $text1->get("1.0", "end-1c");
my $sequence2 = $text2->get("1.0", "end-1c");


$sequence1=~ s/^>.*\n//mg;
$sequence1=~ s/^>.*\n//mg;

if (contains_digit_underscore_or_dot($sequence1)) {
    $sequence1=retrieve_ncbi_sequence($sequence1);
    $sequence1=~ s/^>.*\n//mg;

}

if (contains_digit_underscore_or_dot($sequence2)) {
    $sequence2=retrieve_ncbi_sequence($sequence2);
    $sequence2=~ s/^>.*\n//mg;

}
    my ($aligned_seq1, $aligned_seq2, $score) = smith_waterman($sequence1, $sequence2, 1, -1, -1);

    $textshow->insert('end', "LOCAL ALIGNMENT SCORE : $score\n\n");
    print_aligned_sequences($aligned_seq1, $aligned_seq2);
    return ($aligned_seq1, $aligned_seq2, $score);
}

sub print_aligned_sequences {
    my ($aligned_seq1, $aligned_seq2) = @_;

    my $len = length($aligned_seq1 // '');  # Use the defined operator to handle undefined values

    for (my $i = 0; $i < $len; $i += 60) {
        my $end_index = $i + 59;
        $end_index = $len - 1 if $end_index >= $len;

        my $substr_seq1 = substr($aligned_seq1, $i, $end_index - $i + 1) // '';  
        my $substr_seq2 = substr($aligned_seq2, $i, $end_index - $i + 1) // ''; 

        my @matches = map {
            substr($aligned_seq1, $_, 1) eq substr($aligned_seq2, $_, 1) ? "|" : " "
        } $i .. $end_index;

        $textshow->insert('end', "Sequence1: $substr_seq1\n");
        $textshow->insert('end', "           " . join("", @matches) . "\n");
        $textshow->insert('end', "Sequence2: $substr_seq2\n\n");
        $textshow->insert('end', "\n");
    }
}

sub smith_waterman {
    my ($seq1, $seq2, $match_score, $mismatch_penalty, $gap_penalty) = @_;

    my $len1 = length($seq1);
    my $len2 = length($seq2);

    my @matrix;
    my $max_score = 0;
    my ($max_i, $max_j) = (0, 0);

    my ($aligned_seq1, $aligned_seq2)=("","");
    for my $i (0 .. $len1) {
        for my $j (0 .. $len2) {
            if ($i == 0 || $j == 0) {
                $matrix[$i][$j] = 0;
            } else {
                my $match = substr($seq1, $i - 1, 1) eq substr($seq2, $j - 1, 1) ? $matrix[$i - 1][$j - 1] + $match_score : 0;
                my $delete = $matrix[$i - 1][$j] + $gap_penalty;
                my $insert = $matrix[$i][$j - 1] + $gap_penalty;
                my $score = $match > $delete ? ($match > $insert ? $match : $insert) : ($delete > $insert ? $delete : $insert);
                $matrix[$i][$j] = $score > 0 ? $score : 0;

                if ($matrix[$i][$j] > $max_score) {
                    $max_score = $matrix[$i][$j];
                    ($max_i, $max_j) = ($i, $j);
                }
            }
        }
    }
    while ($matrix[$max_i][$max_j] != 0) {
        if ($matrix[$max_i][$max_j] == $matrix[$max_i - 1][$max_j - 1] + $match_score &&
            substr($seq1, $max_i - 1, 1) eq substr($seq2, $max_j - 1, 1)) {
            $aligned_seq1 = substr($seq1, $max_i - 1, 1) . $aligned_seq1;
            $aligned_seq2 = substr($seq2, $max_j - 1, 1) . $aligned_seq2;
            $max_i--;
            $max_j--;
        } elsif ($matrix[$max_i][$max_j] == $matrix[$max_i - 1][$max_j] + $gap_penalty) {
            $aligned_seq1 = substr($seq1, $max_i - 1, 1) . $aligned_seq1;
            $aligned_seq2 = "-" . $aligned_seq2;
            $max_i--;
        } else {
            $aligned_seq1 = "-" . $aligned_seq1;
            $aligned_seq2 = substr($seq2, $max_j - 1, 1) . $aligned_seq2;
            $max_j--;
        }
    }

    return ($aligned_seq1, $aligned_seq2, $max_score);
}

sub exit_confirmation {
    my $answer = $mw->messageBox(
        -title   => 'Exit Confirmation',
        -message => 'Are you sure you want to exit?',
        -type    => 'yesno',
        -icon    => 'question',);

    if ($answer eq 'Yes') {
        $mw->destroy;
    }
}

sub contains_digit_underscore_or_dot {
    my ($string) = @_;
    return $string =~ /[0-9_.]/;
}

sub retrieve_ncbi_sequence {
    my ($accession_number) = @_;
    my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$accession_number&rettype=fasta&retmode=text";
    my $ua = LWP::UserAgent->new;
    my $response = $ua->get($url);
    if ($response->is_success) {
        my $fasta_data = $response->decoded_content;
        my $seqio = Bio::SeqIO->new(-string => $fasta_data, -format => 'fasta');
        my $sequence = "";
        while (my $seq = $seqio->next_seq) {
            $sequence .= ">$accession_number\n";
            $sequence .= $seq->seq . "\n";
        }

        return $sequence;
    } else {
        die "Error fetching sequence for $accession_number: ", $response->status_line;
    }
}

MainLoop;


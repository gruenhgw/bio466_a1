#! /usr/bin/perl
# Author: George Gruenhagen
# Assignment No.1
$file     = $ARGV[0];
$fileType = $ARGV[1];
open( INPUT, "<$file" ) or die "Can't open file: $file\n";

# Count the number of sequences and store them
$seqCounter = 0;
@sequences  = ();
if ( $fileType eq "fastq" ) {
	# for fastq file type also store the fragment count
	$fragmentCounter = 0;
	$lineCounter     = 1;

	while ( $data = <INPUT> ) {
		chomp($data);
		if ( $lineCounter % 4 == 1 && substr( $data, 4, 1 ) eq $seqCounter + 1 ) {
			$seqCounter++;
			$fragmentCounter = 0;
		}
		if ( $lineCounter % 4 == 1 && substr( $data, 6, 1 ) eq $fragmentCounter + 1 ) {
			$fragmentCounter++;
		}
		if ( $lineCounter % 4 == 2 ) {
			$sequences[ $seqCounter - 1 ][ $fragmentCounter - 1 ] = $data;
		}

		$lineCounter++;
	}
}
else {
	# fasta case
	while ( $data = <INPUT> ) {
		if ( substr( $data, 0, 1 ) eq ">" ) {
			$seqCounter++;
		}
		else {
			chomp($data);
			$sequences[ $seqCounter - 1 ][0] .= $data;
		}
	}
}
print "There are $seqCounter sequences in the file: $file\n";

# ===========================================================================Q1
# Ask the user if they would like to view the file in fasta format and process
# their input.
do {
	print "FASTA format or not (Y/N)? \n";
	$isFasta = <STDIN>;
	$isFasta = lc($isFasta);
	chomp($isFasta);
	$validResponse = 1;
	if ( $isFasta ne "y" && $isFasta ne "yes" && $isFasta ne "n" && $isFasta ne "no" ) {
		print "Invalid Input: Enter Valid Input Y or N.\n";
		$validResponse = 0;
	}
} while ( !$validResponse );

if ( $isFasta eq "y" || $isFasta eq "yes" ) {
	# Print the file in a fasta format
	$isFirstLine = 1;
	open( INPUT, "<$file" );
	while ( $data = <INPUT> ) {
		$firstChar = substr( $data, 0, 1 );

		# If the file is not a fasta file, autogenerate a header, including
		# a name and a description to mimic a fasta file
		if ( $isFirstLine && $firstChar ne ">" ) {
			print ">defSeq|Auto-generated sequence name and description\n";
		}
		elsif ( $firstChar eq "A" || $firstChar eq "T" || $firstChar eq "G" ||
			$firstChar eq "C" || $firstChar eq "N" || $firstChar eq ">" ) {
			print $data;
		}
		$isFirstLine = 0;
	}
}
else {
	# Else the user doesn't want to see the file in fasta format

	#========================================================================Q2
	# Ask for how many nucleotides per line
	do {
		print "How many nucleotides per line [50|100]?\n";
		$nts           = <STDIN>;
		$validResponse = 1;
		if ( $nts != 50 && $nts != 100 ) {
			print "Invalid Input: Enter Valid Input 50 or 100.\n";
			$validResponse = 0;
		}
	} while ( !$validResponse );
	$cols = $nts / 10;

	#========================================================================Q3
	# Ask for spacer
	do {
		print "Do you need a spacer (Y/N)? \n";
		$spacer = <STDIN>;
		$spacer = lc($spacer);
		chomp($spacer);
		$validResponse = 1;
		if ( $spacer eq "y" || $spacer eq "yes" ) {
			$spacer = 1;
		}
		elsif ( $spacer eq "n" || $spacer eq "no" ) {
			$spacer = 0;
		}
		else {
			print "Invalid Input: Enter Valid Input Y or N.\n";
			$validResponse = 0;
		}
	} while ( !$validResponse );

	#========================================================================Q4
	# Ask what sequence should be examined
	do {
		print "Which sequence do you want to examine [";
		for ( $i = 1 ; $i < $seqCounter ; $i++ ) {
			print "$i|";
		}
		print "$seqCounter]\n";
		$seqNum = <STDIN>;
		chomp($seqNum);
		$seqNum--;
		$validResponse = 1;
		if ( $seqNum > $seqCounter || $seqCounter < 1 ) {
			print "Invalid Input: Enter Valid Number 1-$seqCounter.\n";
			$validResponse = 0;
		}
	} while ( !$validResponse );
	
	# If the the file is fastq then assemble the sequences
	# Assuming comment line is @SEQ1_1
	if ( $fileType eq "fastq" ) {
		# Store the overlaps between the fragments into a hash. Where the 
		# sequence that overlaps is the key and the value is an array where
		# the first element is the index of the fragment that has the overlap at
		# its end and the last element is the index of the fragment that has the 
		# overlap at its beginning.
		%overlapFrag = ();
		$maxNt = 15;
		
		# Loop through all the fragments and compare them to the others for overlaps.
		for ( $frag1 = 0 ; $frag1 < $fragmentCounter ; $frag1++ ) {
			for ( $frag2 = $frag1 + 1 ; $frag2 < $fragmentCounter ; $frag2++ ) {

				$overlapBegOfFrag1 = "";
				$overlapEndOfFrag1 = "";
				for ( $i = $maxNt ; $i > 0 ; $i-- ) {
					$begOfFrag1 = substr( $sequences[$seqNum][$frag1], 0, ($maxNt+1) - $i );
					$endOfFrag1 = substr( $sequences[$seqNum][$frag1], length( $sequences[$seqNum][$frag1] ) - ($maxNt+1) + $i, (15 - $i)+1 );
					$begOfFrag2 = substr( $sequences[$seqNum][$frag2], 0, ($maxNt+1) - $i );
					$endOfFrag2 = substr( $sequences[$seqNum][$frag2], length( $sequences[$seqNum][$frag2] ) - ($maxNt+1) + $i, (15 - $i)+1 );
					if ( $begOfFrag1 eq $endOfFrag2 ) {
			  			# checks for overlap at the start of frag1 and at the end of frag2
						$overlapBegOfFrag1 = $begOfFrag1;
					}
					if ( $begOfFrag2 eq $endOfFrag1 ) {
						# checks for overlap at the end of frag1 and at the start of frag2
						$overlapEndOfFrag1 = $begOfFrag2;
					}
				}
				
				# If the overlap is 5 nt or longer store it in the hash.
				if ( length($overlapBegOfFrag1) >= 5 ) {
					$overlapFrag{$overlapBegOfFrag1} = [ $frag2, $frag1 ];
				}
				if ( length($overlapEndOfFrag1) >= 5 ) {
					$overlapFrag{$overlapEndOfFrag1} = [ $frag1, $frag2 ];
				}
			}
		}
	
	# Find start and end fragments by looping through the hash and checking to
	# see if the fragment has overlaps at its beginning and end. If it has an
	# overlap only at its end, it is the start fragment. If it has an overlap
	# only at its beginning its the end fragment.
	@keys = keys(%overlapFrag);
	$startFrag = -1;
	$endFrag = -1;
	for ($frag = 0; $frag < $fragmentCounter; $frag++) {
		$hasBegOverlap = 0;
		$hasEndOverlap = 0;
		foreach $key (@keys) {
			if ($overlapFrag{$key}[0] == $frag) {
				$hasEndOverlap = 1;	
			}
			if ($overlapFrag{$key}[1] == $frag) {
				$hasBegOverlap = 1;	
			}
		}
		if ($hasBegOverlap && !$hasEndOverlap) {
			$endFrag = $frag
		}
		if (!$hasBegOverlap && $hasEndOverlap) {
			$startFrag = $frag
		}
	}
	
	# String the fragments together and delete the extra overlap. Do this by 
	# looping through the hash starting with the start fragment and then going
	# to the fragment it overlaps with until all fragments have been added.
	$sequence = $sequences[$seqNum][$startFrag];
	$visitedFrags = 0;
	$currFrag = $startFrag;
	while ($visitedFrags < $fragmentCounter) {
		foreach $key (@keys) {
			if ($overlapFrag{$key}[0] == $currFrag) {
				$sequence .= substr($sequences[$seqNum][$overlapFrag{$key}[1]], length($key), length($sequences[$seqNum][$overlapFrag{$key}[1]])-length($key));
				$currFrag = $overlapFrag{$key}[1];
				last;
			}
		}
		$visitedFrags++;
	}
	
	
	} # end fastq if
	else {
		# else the file is fasta
		$sequence = $sequences[ $seqNum ][0];
	}
	

	# print the top row: the column numbers
	print "    " . ( $spacer ? "" : " " );
	for ( $col = 1 ; $col <= $cols ; $col++ ) {
		for ( $space = 0 ; $space < 8 + $spacer ; $space++ ) {
			print " ";
		}
		print $col == 10 ? $col : " $col";
	}

	# print the 2nd to the top row: the sub-column numbers
	print "\nLine";
	for ( $col = 1 ; $col <= $cols ; $col++ ) {
		print( !$spacer && $col != 1 ? "" : " " );
		for ( $subCol = 1 ; $subCol <= 10 ; $subCol++ ) {
			print $subCol % 10;
		}
	}

	# Set the counter variables for the nucleotides to 0
	$line_number = 0;
	$a           = 0;
	$t           = 0;
	$g           = 0;
	$c           = 0;
	$n           = 0;
	$o           = 0;
	
	# Print the sequence properly spaced
	for ( $i = 0 ; $i < length($sequence) ; $i++ ) {

		# count nucleotides
		$nt = substr( $sequence, $i, 1 );
		if    ( $nt eq "A" ) { $a++; }
		elsif ( $nt eq "T" ) { $t++; }
		elsif ( $nt eq "G" ) { $g++; }
		elsif ( $nt eq "C" ) { $c++; }
		elsif ( $nt eq "N" ) { $n++; }
		else                 { $o++; }

		# print the line number, but properly spaced for multiple digits
		if ( $i % $nts == 0 ) {
			print "\n";
			for ( $space = 0 ; $space < 3 - ($line_number+1) / 10 ; $space++ ) {
				print " ";
			}
			print ( ($line_number+1). " " );
			$line_number++;
		}
		# Print the spaces if needed
		elsif ( $spacer && $i % 10 == 0 && $i % $nts != 0 ) {
			print " ";
		}
		
		# Print the actual nucleotide
		print $nt;

	}
	# Print the nucleotide counts
	print "\n\nThe sequence in the file $file has nucleotide counts:";
	print "\nA=$a";
	print "\nT=$t";
	print "\nG=$g";
	print "\nC=$c";
	print "\nN=$n";
	print "\nOther=$o";
}
# 🧬 Pairwise Sequence Alignment Tool

![Perl](https://img.shields.io/badge/Language-Perl-blue?style=for-the-badge&logo=perl)
![License](https://img.shields.io/badge/License-MIT-green?style=for-the-badge)
![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-orange?style=for-the-badge)
![GUI](https://img.shields.io/badge/Interface-Tkinter-purple?style=for-the-badge)

A powerful Perl-based graphical application for performing both local and global sequence alignments using dynamic programming algorithms! 🚀

## 🔬 Overview

This tool provides a user-friendly GUI interface for comparing biological sequences using two different alignment algorithms:

- 🌍 **Global Alignment**: Needleman-Wunsch algorithm for end-to-end sequence comparison
- 📍 **Local Alignment**: Smith-Waterman algorithm for finding optimal local alignments

## ✨ Features

- 📁 **Dual Input Methods**: Accept sequences in FASTA format or NCBI accession numbers
- 🌐 **Automatic NCBI Retrieval**: Fetch sequences directly from NCBI database using accession numbers
- 🖥️ **Interactive GUI**: User-friendly Tkinter-based interface with scrollable text areas
- 👁️ **Visual Alignment Display**: Shows aligned sequences with match indicators (|) and gaps (-)
- 📊 **Scoring System**: Displays alignment scores for both local and global alignments
- ⚡ **Real-time Results**: View alignment results immediately in the application

## 📁 Files

- 📜 `seq align.pl` - Main Perl script containing the GUI application and alignment algorithms
- 📋 `idname.txt` - Sample data file containing NCBI accession numbers and example sequences including:
  - 🧬 BRCA1 gene sequences from various species (Macropus rufus, Didelphis virginiana)
  - 💉 Insulin mRNA sequences (Octodon degus, Aplysia californica)

## 🛠️ Prerequisites

### 📦 Required Perl Modules

```perl
use strict;
use warnings;
use LWP::UserAgent;    # For NCBI sequence retrieval
use Bio::SeqIO;        # For FASTA sequence parsing
use Tk;                # For GUI interface
```

### 💻 Installation

Install required modules using CPAN:

```bash
cpan install LWP::UserAgent
cpan install Bio::SeqIO
cpan install Tk
```

## 🚀 Usage

### 🏃‍♂️ Running the Application

```bash
perl "seq align.pl"
```

### 🖼️ Interface Overview

The application window contains:

- 📝 **Two text input areas**: For entering sequences or accession numbers
- 🔘 **Three main buttons**:
  - `RUN LOCAL ALIGNMENT` - Performs Smith-Waterman local alignment
  - `RUN GLOBAL ALIGNMENT` - Performs Needleman-Wunsch global alignment
  - `EXIT` - Close the application
- 📊 **Results panel**: Displays alignment results with scores and visual representation

### 📥 Input Formats

#### 🧬 FASTA Format

```fasta
>Sequence_Name
ATCGATCGATCGATCG
```

#### 🔢 NCBI Accession Numbers

Simply enter the accession number (e.g., `AY211956.1`, `NM_001204686.1`)

### 🧪 Example Sequences

The `idname.txt` file contains several example sequences you can use for testing:

- 🦘 `AY211956.1` - Macropus rufus BRCA1 gene
- 🐭 `AY211955.1` - Didelphis virginiana BRCA1 gene
- 🐹 `M57671.1` - Octodon degus insulin mRNA
- 🐌 `NM_001204686.1` - Aplysia californica insulin precursor

## ⚙️ Algorithm Details

### 🌍 Global Alignment (Needleman-Wunsch)

- ✅ **Match Score**: +1
- ❌ **Mismatch Score**: 0
- 🕳️ **Gap Penalty**: -2
- 🎯 Aligns sequences end-to-end for overall similarity

### 📍 Local Alignment (Smith-Waterman)

- ✅ **Match Score**: +1
- ❌ **Mismatch Penalty**: -1
- 🕳️ **Gap Penalty**: -1
- 🎯 Finds the best local alignment within sequences

## 📋 Output Format

Results are displayed showing:

- 📊 Alignment score
- 🧬 Aligned sequences with gaps represented by dashes (-)
- ✅ Match indicators (|) showing identical positions
- 📏 60 characters per line for easy reading

## 🔧 Technical Implementation

- 🧮 **Dynamic Programming**: Both algorithms use matrix-based dynamic programming
- 🔄 **Traceback**: Reconstructs optimal alignment path
- 🖥️ **GUI Framework**: Built with Perl/Tk for cross-platform compatibility
- 🌐 **Web Integration**: Uses LWP::UserAgent for NCBI database queries
- 🧬 **Bioinformatics Support**: Leverages BioPerl for sequence handling

## 🐛 Troubleshooting

- ✅ Ensure all required Perl modules are installed
- 🌐 Check internet connectivity for NCBI sequence retrieval
- 📝 Verify that input sequences are in valid FASTA format
- ⏱️ For large sequences, processing may take additional time

## 📄 License

MIT License - This project is licensed under the MIT License, making it free to use, modify, and distribute for both educational and commercial purposes.

## 🤝 Contributing

Feel free to enhance the tool by:

- 🧪 Adding support for protein sequences
- 📊 Implementing additional scoring matrices (BLOSUM, PAM)
- 🔄 Adding batch processing capabilities
- 🎨 Improving the GUI interface

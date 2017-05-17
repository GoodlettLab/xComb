# xComb  
## What is xComb?  
xComb is a program that computes all possible cross-links based on proteins of
interest(Max. 30 proteins) for further interrogation using standard search
engine (e.g. Phenyx, Mascot, Sequest, OMSSA, X!Tandem). If you use xComb, please
cite our original publication:

Panchaud A, Singh P, Shaffer S a, Goodlett DR. xComb: A Cross-Linked Peptide
Database Approach to Protein-Protein Interaction Analysis. J Proteome Res
2010;9:2508–15. [[Link]](http://pubs.acs.org/doi/abs/10.1021/pr9011816)

## Licensing Information  
xComb is available as open source software here on GitHub, under the GNU
General Public License Version 3.0 (See the `LICENSE.md` file for details).
Though we offer xComb free for non-profit and academic use, commercial interests
should be directed through the Unversity of Washington Express Licensing Program
at:  

UW CoMotion
4545 Roosevelt Way NE, Suite 400
Seattle, WA 98105
Phone: (206) 543-3970
Fax: (206) 616-3322
Email: license@uw.edu


## How do I use it? (2 options)

### 1. On the web  
The easiest way to use xComb is through our [website](https://goodlett.umaryland.edu/xcomb/index.php).

### 2. Download and use the software locally
First, you will need to install Perl. Which can be downloaded here: [perl.org](https://www.perl.org/)  

Once you have Perl installed proceed with the `protein2digest.pl` script which
takes in your fasta file(s) and produces a directory with .digest files in it.

Below is an example of the usage for `protein2digest.pl`:

```
protein2digest.pl (xComb v1.4)

FUNCTION:
Transforms DAT or FASTA protein files into DIGEST format for further use with
compatible softwares.

USAGE:
protein2digest.pl -in=[value] -ez=[value] -mc=[value] [file_name]

        -in=[value]: Specify input format

          0= UniProt DAT
          1= UniProt FASTA
          2= Any FASTA

        -ez=[value]: Specify enzyme

          0= Trypsin
          1= Lys-C
          2= Arg-C
          3= Glu-C (E)
          4= Glu-C (DE)
          5= Asp-N
          6= Trypsin, no P rule
          7= ARG-C, no P rule
          8= Chymotrypsin, no P rule
                9= Thermolysin

        -mc=[value]: Specify nb of missed cleavage allowed (max. = 3)

        [file_name]: Specify file to be converted

Ex: protein2digest.pl -in=1 -ez=0 -mc=2 filename
```
The resulting directory can now be used by `digest2cxdb.pl` to create your new
cross-linked database.

Below is an example of the usage for `digest2cxdb.pl`:

```
digest2cxdb.pl (xComb v1.4)

FUNCTION:
Creates inter & intra cross-linking databases for search with standard algorithm
such as Sequest, Mascot or Phenyx. Uses pre-digested protein sequences in DIGEST
format. If a single digest is present, performs only intra database. The maximum
protein sequences allowed per directory in digest format is 50!!!

USAGE:
digest2cxdb.pl -i=[value] -x=[value] -mc=[value] -pl=[value] -f=[value]
-t=[value] [digest_dir]

        -i=[value]: Inter/Intra type mode

          0= Both intra and inter crosslinks are generated
          1= Only inter crosslinks are generated
          2= Only intra crosslinks are generated

        -x=[value]: Specify type of cross-linker

          0= Amine/Amine (e.g. BS2, DST, DFDNB)
          1= Sulfhydryl/Sulfhydryl (e.g. BMOE, BMDB or S-S bridge)
          2= Amine/carboxyl (e.g. EDC)
          3= Amine/Non-selective (e.g. SDA)
          4= Amine/Sulfhydryl (e.g. AMAS, SPDP)
          5= Tyrosine/Tyrosine (e.g. dityrosine)
          6= Lysine/Lysine

        -mc=[value]: Specify whether amine cross-linking needs 1 trypsin missed
cleavage at least. on=1 / off=0

        -pl=[value]: Specify the minimum peptide length for building a
        permutation [1-6]

        -f=[value]: Phenyx FASTA formatted output. on=1 / off=0

        -t=[value]: Test mode which adds "|" between the two peptides for easier
        reading. on=1 / off=0

        [digest_dir]: Specifies directory containing files in DIGEST format

Ex: digest2cxdb.pl -i=1 -x=1 -mc=1 -pl=4 -f=0 -t=0 digest_dir
```

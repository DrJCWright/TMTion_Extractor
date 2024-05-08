TMT Ion Extractor and PLink2 Mapper

Dr James C Wright
Institute of Cancer Research - London
2024

This tools can quickly parse all MGF files in a directory extracting TMT ion peaks present in spectra without the need to identify the spectrum first.
A table delimited table of spectra and ions is generated for further analyses.
The tool can optionally parse Protein and CrossLinked Peptide results from PLink2 and map TMT ions into the results files.

usage: 

Basic:
python ExtractTMTions_V1.py --mgfdir MGFDIR

PLink Ion Mapping:
ExtractTMTions_V1.py --mgfdir MGFDIR --proteins PROTEINFILE --crosslinks CROSSLINKFILE


Full cmd line parameters:
usage: ExtractTMTions_V1.py [-h] --mgfdir MGFDIR [--tol ITOL] [--proteins PROTEINFILE] [--crosslinks CROSSLINKFILE] [--outprefix OUTPRE] [--norm]

Tool for extracting TMT ion intensities and optionally mapping to plink2 output. Dr James Wright - 2024 - Institute of Cancer Research

optional arguments:
  
  --mgfdir MGFDIR, -i MGFDIR    Directory conatining MGF spectrum files matching the plink2 input (required)

  --tol ITOL, -t ITOL   TMT ion matching tolerance inppm, default = 15ppm (optional)

  --proteins PROTEINFILE, -p PROTEINFILE  csv plink2 output File Containing proteins and peptides (optional)
  --crosslinks CROSSLINKFILE, -c CROSSLINKFILE  csv plink2 output File Containing crosslinked peptides (optional)

  --outprefix OUTPRE, -o OUTPRE prefix for output files (optional)

  --norm, -n  Run normalisation and scaling of TMT channels using total TMT channel intensity (beta method). Default=False (optional)

  -h, --help  show help message and exit
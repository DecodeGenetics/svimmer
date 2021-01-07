## svimmer - SV merging tool

Merges similar SVs from multiple single sample VCF files. The tool was written for merging SVs discovered using Manta calls, but should support (almost) any SV VCFs. The output is a VCF file containing all merged SV sites (with no calls). The output can be given as input into [GraphTyper](https://github.com/DecodeGenetics/graphtyper) to genotype the sites.

### Requirements

* Python 3.4+
* pysam

### Usage

```sh
python3 svimmer input_vcfs chrA chrB chrC ...
```

where input is a list of tabix indexed+bgzipped VCF files and chromosomes are the chromosomes to merge. For further details see the help page:

```sh
python3 svimmer -h
```

## Test data example

```sh
python3 svimmer test_vcfs chr20 > test/actual_output.vcf
diff test/actual_output.vcf test/expected_output.vcf
```

## License
GNU GPLv3

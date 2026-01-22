import sys

def parse_info(info_str):
    info = {}
    for entry in info_str.split(';'):
        if '=' in entry:
            key, value = entry.split('=', 1)
            info[key] = value
        else:
            info[entry] = True
    return info

def process_vcf(input_vcf, output_txt):
    with open(input_vcf, 'r') as fin, open(output_txt, 'w') as fout:
        fout.write('Chr\tstart\tEnd\tSVTYPE\tID\tBreakpoint\n')
        for line in fin:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            chrom = fields[0].replace('chr', '')  # Remove 'chr' prefix if present
            start = fields[1]
            var_id = fields[2]
            info = parse_info(fields[7])
            svtype = info.get('SVTYPE', '')
            end = info.get('END', start)
            svlen = info.get('SVLEN', None)
            # Write start breakpoint
            fout.write(f"{chrom}\t{start}\t{end}\t{svtype}\t{var_id}\tstart\n")
            # Write end breakpoint if SVLEN is present and numeric
            if svlen is not None:
                try:
                    svlen_int = int(svlen)
                    end_pos = int(start) + svlen_int
                    fout.write(f"{chrom}\t{start}\t{end_pos}\t{svtype}\t{var_id}\tend\n")
                except ValueError:
                    pass

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Extract breakpoints from VCF")
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--out', required=True, help='Output file')
    args = parser.parse_args()
    process_vcf(args.vcf, args.out)

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

def process_vcf(input_vcf, output_bed):
    with open(input_vcf, 'r') as fin, open(output_bed, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            chrom = fields[0]
            start = int(fields[1])
            var_id = fields[2]
            info = parse_info(fields[7])
            svtype = info.get('SVTYPE', '')
            end = int(info.get('END', start))
            svlen = info.get('SVLEN', None)
            # Write start breakpoint (BED: 0-based start, exclusive end)
            bed_start = start - 1
            bed_end = start
            name = f"{var_id}|{svtype}|start"
            fout.write(f"{chrom}\t{bed_start}\t{bed_end}\t{name}\t.\t.\n")
            # Write end breakpoint if SVLEN is present and numeric
            if svlen is not None:
                try:
                    svlen_int = int(svlen)
                    end_pos = start + svlen_int
                    bed_end_start = end_pos - 1
                    bed_end_end = end_pos
                    name_end = f"{var_id}|{svtype}|end"
                    fout.write(f"{chrom}\t{bed_end_start}\t{bed_end_end}\t{name_end}\t.\t.\n")
                except ValueError:
                    pass

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Extract breakpoints from VCF to BED")
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--out', required=True, help='Output BED file')
    args = parser.parse_args()
    process_vcf(args.vcf, args.out)

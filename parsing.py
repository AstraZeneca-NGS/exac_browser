"""
Utils for reading flat files that are loaded into database
"""
import re
import traceback
from utils import *
import copy


def get_base_coverage_from_file(base_coverage_file, canonical_transcripts):
    """
    Read a base coverage file and return iter of dicts that look like:
    {
        'xpos': 1e9+1,
        'mean': 0.0,
        'median': 0.0,
        '1': 0.0,
        '5': 0.0,
        '10': 0.0,
        '15': 0.0,
        '20': 0.0,
        '25': 0.0,
        '30': 0.0,
        '50': 0.0,
        '100': 0.0,
    }
    """

    float_header_fields = ['mean', 'median', '1', '5', '10', '15', '20', '25', '30', '50', '100']
    for line in base_coverage_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')
        d = {
            'xpos': get_xpos(fields[0], int(fields[1])),
            'pos': int(fields[1]),
        }
        for i, k in enumerate(float_header_fields):
            d[k] = float(fields[i+2])
        yield d


def get_variants_from_sites_vcf(sites_vcf, canonical_transcripts):
    """
    Parse exac sites VCF file and return iter of variant dicts
    sites_vcf is a file (gzipped), not file path
    """
    ann_field_names = None
    samples = []

    for line in sites_vcf:
        try:
            line = line.strip('\n')
            if line.startswith('##INFO=<ID=ANN'):
                ann_field_names = line.split(': ')[-1].strip('">').strip("'").split('|')
                ann_field_names = [f.strip() for f in ann_field_names]
            if line.startswith('##INFO=<ID=DP_HIST'):
                dp_mids = map(float, line.split('Mids: ')[-1].strip('">').split('|'))
            if line.startswith('##INFO=<ID=GQ_HIST'):
                gq_mids = map(float, line.split('Mids: ')[-1].strip('">').split('|'))
            if line.startswith('#CHROM'):
                fs = line.split('\t')
                samples = fs[fs.index('FORMAT') + 1:]
            if line.startswith('#'):
                continue

            # If we get here, it's a variant line
            if ann_field_names is None:
                raise Exception("ANN_field_names is None. Make sure VCF header is present.")
            # This elegant parsing code below is copied from https://github.com/konradjk/loftee
            fields = line.split('\t')
            info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', fields[7])])
            annotation_array = info_field['ANN'].split(',') if 'ANN' in info_field else []
            all_annotations = [dict(zip(ann_field_names, x.split('|'))) for x in annotation_array if len(ann_field_names) == len(x.split('|'))]
            coding_annotations = [ann for ann in all_annotations if ann['Feature_ID'].startswith('NM')]

            alt_alleles = fields[4].split(',')

            # different variant for each alt allele
            for i, alt_allele in enumerate(alt_alleles):

                annotations = [ann for ann in coding_annotations if (ann['Allele']) == alt_allele]

                # Variant is just a dict
                # Make a copy of the info_field dict - so all the original data remains
                # Add some new keys that are allele-specific
                pos, ref, alt = get_minimal_representation(fields[1], fields[3], alt_allele)

                variant = {}
                variant['chrom'] = fields[0]
                variant['pos'] = pos
                rs_ids = re.findall(r'rs\d+', fields[2])
                if rs_ids:
                    variant['rsid'] = rs_ids[0]
                variant['xpos'] = get_xpos(variant['chrom'], variant['pos'])
                variant['ref'] = ref
                variant['alt'] = alt
                variant['xstart'] = variant['xpos']
                variant['xstop'] = variant['xpos'] + len(variant['alt']) - len(variant['ref'])
                variant['variant_id'] = '{}-{}-{}-{}'.format(variant['chrom'], variant['pos'], variant['ref'], variant['alt'])
                variant['orig_alt_alleles'] = [
                    '{}-{}-{}-{}'.format(variant['chrom'], *get_minimal_representation(fields[1], fields[3], x))
                    for x in alt_alleles
                ]
                variant['site_quality'] = float(fields[5])
                variant['filter'] = fields[6]
                variant['vep_annotations'] = [dict((k.replace('.', '_'), v) for k, v in annotation.iteritems()) for annotation in annotations]

                variant['allele_freq'] = float(info_field['AF'].split(',')[i])
                variant['allele_count'] = int(info_field['AC'].split(',')[i])
                variant['allele_num'] = int(info_field['AN'])

                if 'AC_MALE' in info_field:
                    variant['ac_male'] = info_field['AC_MALE']
                if 'AC_FEMALE' in info_field:
                    variant['ac_female'] = info_field['AC_FEMALE']
                if 'AN_MALE' in info_field:
                    variant['an_male'] = info_field['AN_MALE']
                if 'AN_FEMALE' in info_field:
                    variant['an_female'] = info_field['AN_FEMALE']
                variant['sample_names'] = samples
                samples_data = fields[-len(samples):]
                samples_data_info = fields[-len(samples) - 1].split(':')
                variant['sample_data'] = []
                for idx, sample in enumerate(samples):
                    fs = samples_data[idx].split(':')
                    if len(fs) < 6:
                        variant['sample_data'].append('')
                        continue
                    freq_col = samples_data_info.index('AF')
                    depth_col = samples_data_info.index('DP')
                    freq, depth = fs[freq_col], fs[depth_col]
                    if freq.replace('.', '', 1).isdigit():
                        freq = str(float(freq) * 100) + '%'
                    variant['sample_data'].append('AF:' + freq + ',DP:' + depth)

                if variant['chrom'] in ('X', 'Y'):
                    # variant['pop_hemis'] = dict([(POPS[x], int(info_field['Hemi_%s' % x].split(',')[i])) for x in POPS])
                    variant['hemi_count'] = sum(variant['pop_hemis'].values())
                variant['quality_metrics'] = dict([(x.replace('.', '_'), info_field[x]) for x in METRICS if x in info_field])

                variant['genes'] = set()
                variant['transcripts'] = set()
                for vep_annotation, annotation in zip(variant['vep_annotations'], annotations):
                    gene = annotation['Gene_Name']
                    transcript = annotation['Feature_ID'].split('.')[0]
                    if gene in canonical_transcripts and canonical_transcripts[gene] == transcript:
                        vep_annotation['CANONICAL'] = 'YES'
                    else:
                        vep_annotation['CANONICAL'] = 'NO'
                    if 'LoF' not in annotation:
                        vep_annotation['LoF'] = ''
                    variant['genes'].add(gene)
                    variant['transcripts'].add(transcript)
                    vep_annotation['Feature_ID'] = transcript

                variant['genes'] = list(variant['genes'])
                variant['transcripts'] = list(variant['transcripts'])
                if 'DP_HIST' in info_field:
                    hists_all = [info_field['DP_HIST'].split(',')[0], info_field['DP_HIST'].split(',')[i+1]]
                    variant['genotype_depths'] = [zip(dp_mids, map(int, x.split('|'))) for x in hists_all]
                if 'GQ_HIST' in info_field:
                    hists_all = [info_field['GQ_HIST'].split(',')[0], info_field['GQ_HIST'].split(',')[i+1]]
                    variant['genotype_qualities'] = [zip(gq_mids, map(int, x.split('|'))) for x in hists_all]

                yield variant
        except Exception:
            print("Error parsing vcf line: " + line)
            traceback.print_exc()
            break


def get_regions(regions_file, canonical_transcripts):
    float_header_fields = ['size', 'gene', 'depth', 'samples', 'annotation']
    depth_threshold = 0
    for line in regions_file:
        if line.startswith('#'):
            if 'Coverage threshold Nx is ' in line:
                threshold_pattern = 'Coverage threshold Nx is (\d+)x'
                depth_threshold = re.findall(threshold_pattern, line)[0]
            continue
        fields = line.strip('\n').split('\t')
        chrom = fields[0]
        if chrom not in CHROMOSOME_TO_CODE:
            continue
        start = int(fields[1])
        stop = int(fields[2])
        gene = fields[4]
        if not gene or gene == '.' or gene == 'None':
            continue
        d = {
            'chrom': chrom,
            'start': format_value(start, is_html=True, human_readable=True),
            'stop': format_value(start, is_html=True, human_readable=True),
            'depth_threshold': depth_threshold
        }
        for i, k in enumerate(float_header_fields):
            if i + 3 < len(float_header_fields):
                d[k] = float(fields[i+3])
        yield d


def get_mnp_data(mnp_file):
    header = mnp_file.readline().strip().split('\t')
    for line in mnp_file:
        data = dict(zip(header, line.split('\t')))
        if any(map(lambda x: x == 'True', data['QUESTIONABLE_PHASING'])): continue
        chroms = data['CHROM'].split(',')
        chrom = chroms[0]
        sites = data['SITES'].split(',')
        refs = data['REF'].split(',')
        alts = data['ALT'].split(',')
        for i, site in enumerate(sites):
            all_sites = zip(chroms, sites, refs, alts)
            all_sites.remove(all_sites[i])
            mnp = {}
            mnp['xpos'] = get_xpos(chrom, site)
            mnp['ref'] = refs[i]
            mnp['alt'] = alts[i]
            mnp['site2'] = '-'.join(all_sites[0])
            if len(all_sites) > 1:
                mnp['site3'] = all_sites[1]
            mnp['combined_codon_change'] = data['COMBINED_CODON_CHANGE']
            mnp['category'] = data['CATEGORY']
            mnp['number_samples'] = data['NSAMPS']
            yield mnp


def get_constraint_information(constraint_file):
    _, _, _, header = constraint_file.readline().strip().split(None, 3)
    header = header.split()
    for line in constraint_file:
        transcript, gene, chrom, info = line.strip().split(None, 3)
        transcript_info = dict(zip(header, map(float, info.split())))
        transcript_info['transcript'] = transcript.split('.')[0]
        yield transcript_info


def get_canonical_transcripts(canonical_transcript_file):
    for line in canonical_transcript_file:
        gene, transcript = line.strip().split()
        yield gene, transcript


def get_omim_associations(omim_file):
    for line in omim_file:
        fields = line.strip().split('\t')
        if len(fields) == 4:
            yield fields
        else:
            yield None


def get_genes_from_features(features_file):
    """
    Parse features bed file;
    Returns iter of gene dicts
    """
    gene_ids = set()
    for line in features_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')
        feature_type = fields[6]

        if feature_type != 'Gene':
            continue

        chrom = fields[0]
        if chrom not in CHROMOSOME_TO_CODE:
            continue
        start = int(fields[1])
        stop = int(fields[2])
        gene_id = fields[3]
        strand = fields[5]
        if (gene_id, chrom) in gene_ids:
            continue
        gene_ids.add((gene_id, chrom))
        gene = {
            'gene_id': gene_id,
            'gene_name': gene_id,
            'gene_name_upper': gene_id.upper(),
            'chrom': chrom[3:],
            'start': start,
            'stop': stop,
            'strand': strand,
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield gene


def get_transcripts_from_features(features_file):
    """
    Parse features bed file;
    Returns iter of transcript dicts
    """
    for line in features_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')
        feature_type = fields[6]

        if feature_type != 'Transcript':
            continue

        chrom = fields[0]
        if chrom not in CHROMOSOME_TO_CODE:
            continue
        start = int(fields[1])
        stop = int(fields[2])
        gene_id = fields[3]
        strand = fields[5]
        transcript_id = fields[8]

        gene = {
            'transcript_id': transcript_id,
            'gene_id': gene_id,
            'chrom': chrom[3:],
            'start': start,
            'stop': stop,
            'strand': strand,
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield gene


def get_exons_from_features(features_file):
    """
    Parse features bed file;
    Returns iter of transcript dicts
    """
    for line in features_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')
        feature_type = fields[6]
        if fields[7] == 'UTR':
            feature_type = 'UTR'

        if feature_type not in ['Exon', 'CDS', 'UTR']:
            continue

        chrom = fields[0]
        if chrom not in CHROMOSOME_TO_CODE:
            continue
        start = int(fields[1])
        stop = int(fields[2])
        gene_id = fields[3]
        strand = fields[5]
        transcript_id = fields[8]

        exon = {
            'feature_type': feature_type,
            'transcript_id': transcript_id,
            'gene_id': gene_id,
            'chrom': chrom[3:],
            'start': start,
            'stop': stop,
            'strand': strand,
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield exon


def get_dbnsfp_info(dbnsfp_file):
    """
    Parse dbNSFP_gene file;
    Returns iter of transcript dicts
    """
    header = dbnsfp_file.next().split('\t')
    fields = dict(zip(header, range(len(header))))
    for line in dbnsfp_file:
        line = line.split('\t')
        other_names = line[fields["Gene_old_names"]].split(';') if line[fields["Gene_old_names"]] != '.' else []
        if line[fields["Gene_other_names"]] != '.':
            other_names.extend(line[fields["Gene_other_names"]].split(';'))
        gene_info = {
            'gene_name': line[fields["Gene_name"]],
            'ensembl_gene': line[fields["Ensembl_gene"]],
            'gene_full_name': line[fields["Gene_full_name"]],
            'gene_other_names': other_names
        }
        yield gene_info


def get_snp_from_dbsnp_file(dbsnp_file, canonical_transcripts):
    for line in dbsnp_file:
        fields = line.split('\t')
        if len(fields) < 3: continue
        rsid = int(fields[0])
        chrom = fields[1].rstrip('T')
        if chrom == 'PAR': continue
        start = int(fields[2]) + 1
        snp = {
            'xpos': get_xpos(chrom, start),
            'rsid': rsid
        }
        yield snp

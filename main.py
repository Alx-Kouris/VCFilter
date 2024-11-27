import pysam
from multiprocessing import Pool
from functools import partial
from scipy.stats import fisher_exact
import os
import argparse
import pandas as pd

clinical_settings = {
    "DP": 10,
    "AF": 0.5
}

research_settings = {
    "DP": 5,
    "AF": 0.4
}

# Mapping between BAMs and VCFs
file_mapping = {}

class UserProfile:
    def __init__(self):
        self.settings_dictionary = {}
        self.profile_name = ""

    @staticmethod
    def get_clinical_profile():
        result = UserProfile()
        result.settings_dictionary = clinical_settings
        result.profile_name = "Clinical"
        return result

    @staticmethod
    def get_research_profile():
        result = UserProfile()
        result.settings_dictionary = research_settings
        result.profile_name = "Research"
        return result

    def has_strand_bias_set(self):
        return 'strand_bias_threshold' in self.settings_dictionary.keys()

    def strand_bias_threshold(self):
        return self.settings_dictionary['strand_bias_threshold']


def setup_output_filename(filename):
    # Generate output filename with "filtered" inserted
    base, ext = os.path.splitext(filename)
    if ext == ".gz":  # Handle .vcf.gz files
        base, ext2 = os.path.splitext(base)
        output_vcf = f"{base}.filtered{ext2}.gz"
    else:
        output_vcf = f"{base}.filtered{ext}"

    return output_vcf

def load_variant_list(variant_file):
    # Load variants from a CSV file with columns: chromosome, position, ref, alt
    variants = pd.read_csv(variant_file, sep="\\t")
    # Convert to a set of tuples for quick lookup (chromosome, position, ref, alt)
    variant_set = set(zip(variants['chromosome'], variants['position'], variants['ref'], variants['alt']))
    return variant_set

def calculate_strand_bias(bam_file, region):
    bam = pysam.AlignmentFile(bam_file, "rb")

    strand_counts = {"forward": 0, "reverse": 0}

    # Iterate through reads in the specified region
    for read in bam.fetch(region=region):
        if read.is_unmapped:
            continue
        if read.is_reverse:
            strand_counts["reverse"] += 1
        else:
            strand_counts["forward"] += 1

    bam.close()

    # Perform Fisher's exact test
    forward = strand_counts["forward"]
    reverse = strand_counts["reverse"]
    table = [[forward, reverse], [reverse, forward]]
    _, p_value = fisher_exact(table)

    return p_value

def associate_bams_with_vcfs(bams, vcfs):
    for bam in bams:
        bam_base = os.path.basename(bam)
        sample_name = bam_base.split('.')[0]
        for vcf in vcfs:
            vcf_base = os.path.basename(vcf)
            if vcf_base.startswith(sample_name):
                print(f"BAM file={bam_base} with sample_name={sample_name} was associated with VCF={vcf_base}")
                file_mapping[vcf] = bam

def process_vcf(vcf_file_path, profile, debug = False, variant_file = None):
    print("Opening VCF: ", vcf_file_path)
    vcf_in = pysam.VariantFile(vcf_file_path)

    output_vcf = setup_output_filename(vcf_file_path)
    vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)

    filter_out_variant_list = load_variant_list(variant_file) if variant_file else set()

    for record in vcf_in:
        if debug:
            print(f"\nVariant at {record.chrom}:{record.pos}")

        ################################################################################################################
        ## Filter out variants that exist in custom provided variant list
        ################################################################################################################
        if len(filter_out_variant_list):
            found = False
            for alt in record.alts:
                variant_tuple = (record.chrom, record.pos, record.ref, alt)
                if variant_tuple in filter_out_variant_list:
                    print(f"Filtering out variant {record.chrom}:{record.pos}:{record.ref}:{alt}")
                    found = True

            if found:
                continue

        ################################################################################################################
        ## Filter out variants that overlap with region that stand bias may exist
        ## BAM file must be provided in order for this filtering to work
        ################################################################################################################
        if len(file_mapping) and user_profile.has_strand_bias_set():
            region = f"{record.chrom}:{record.pos}-{record.pos + len(record.ref)}"
            strand_bias = calculate_strand_bias(file_mapping[vcf_file_path], region)
            if strand_bias < profile.strand_bias_threshold:
                print(f"Strand bias confirmed for region={region}")
                continue

        if len(record.alts) > 1:
            for alt_index, alt_allele in enumerate(record.alts):
                split_record = record.copy()
                split_record.alts = (alt_allele,)

                passes_info_filters = True
                for field, threshold in profile.settings_dictionary.items():
                    if field in record.info:
                        if isinstance(record.info[field], tuple):
                            split_record.info[field] = (record.info[field][alt_index],)
                        else:
                            split_record.info[field] = record.info[field]

                        if (isinstance(split_record.info[field], tuple) and split_record.info[field][0] <= threshold) or \
                            (not isinstance(split_record.info[field], tuple) and split_record.info[field] <= threshold):
                            passes_info_filters = False
                            break

                passes_format_filters = True
                for sample_name, sample in split_record.samples.items():
                    for field, threshold in profile.settings_dictionary.items():
                        if field in sample:
                            if isinstance(sample[field], tuple):
                                sample[field] = (sample[field][alt_index],)
                            if (isinstance(sample[field], tuple) and sample[field][0] <= threshold) or \
                                    (not isinstance(sample[field], tuple) and sample[field] <= threshold):
                                passes_format_filters = False
                                break
                    if not passes_format_filters:
                        break

                if passes_info_filters or passes_format_filters:
                    vcf_out.write(split_record)

        else:
            passes_info_filters = True
            for field, threshold in profile.settings_dictionary.items():
                if field in record.info:
                    if (isinstance(record.info[field], tuple) and record.info[field][0] <= threshold) or \
                            (not isinstance(record.info[field], tuple) and record.info[field] <= threshold):
                        passes_info_filters = False
                        break

            passes_format_filters = True
            for sample_name, sample in record.samples.items():
                for field, threshold in profile.settings_dictionary.items():
                    if field in sample:
                        if (isinstance(sample[field], tuple) and sample[field][0] <= threshold) or \
                                (not isinstance(sample[field], tuple) and sample[field] <= threshold):
                            passes_format_filters = False
                            break
                if not passes_format_filters:
                    break

            # Write the record only if it passes both INFO and FORMAT filters
            if passes_info_filters or passes_format_filters:
                vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()

parser = argparse.ArgumentParser(description="Process some integers.")

# Add arguments
parser.add_argument("--files_path", type=str, help="The directory containing the VCFs to process")
parser.add_argument("--user_profile", choices=["clinical", "research"], type=str, help="User profile specifying the settings, clinical or research")
parser.add_argument("--debug", type=bool, default=False, help="Print debug messages")
parser.add_argument("--variant_list", type=str, help="Path to file containing variants to exclude (columns: chromosome, position, ref, alt)")
parser.add_argument("--bam_files_path", type=str, help="The directory containing the BAMs of the VCFs. VCF file and BAM should have the same sample name in order to be associated")
parser.add_argument("--strand_bias_threshold", default=0.05, type=float, help="The p-value from Fisher's exact test that would be used to determine strand bias")

args = parser.parse_args()

print("Files path: ", args.files_path)
print("User profile: ", args.user_profile)
print("Variant list file: ", args.variant_list)
print("BAM files path: ", args.bam_files_path)
print("Strand Bias threshold: ", args.strand_bias_threshold)
print("Debug: ", args.debug)

vcf_files = [
    os.path.join(args.files_path, f)
    for f in os.listdir(args.files_path)
    if f.endswith(".vcf") or f.endswith(".vcf.gz")
]

bam_files = [
    os.path.join(args.bam_files_path, f)
    for f in os.listdir(args.bam_files_path)
    if f.endswith(".bam")
]

if len(bam_files):
    associate_bams_with_vcfs(bam_files, vcf_files)

if args.user_profile.lower() == "clinical":
    user_profile = UserProfile.get_clinical_profile()
else:
    user_profile = UserProfile.get_research_profile()

if args.strand_bias_threshold > 0:
    user_profile.settings_dictionary['strand_bias_threshold'] = args.strand_bias_threshold

with Pool(processes=4) as pool:
    process_vcf_with_profile = partial(process_vcf, profile=user_profile)
    results = pool.map(process_vcf_with_profile, vcf_files)

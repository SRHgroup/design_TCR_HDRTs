import io
import os
import re
import subprocess
import argparse

import zipfile
import tempfile

import logging
from datetime import datetime

import pandas as pd
import numpy as np

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

# =========================
# Logging setup
# =========================
def setup_logger(log_path):
    logger = logging.getLogger("TCR_pipeline")
    logger.setLevel(logging.INFO)

    if logger.handlers:
        return logger   # prevents duplicate handlers

    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s"
    )

    fh = logging.FileHandler(log_path)
    fh.setFormatter(formatter)

    ch = logging.StreamHandler()
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger

# =========================
# Codon utilities
# =========================

# Human codon usage obtaind from http://www.kazusa.or.jp/codon/
human_codon_usage = {
    "F": {"TTT": 0.46, "TTC": 0.54},
    "S": {"TCT": 0.19, "TCC": 0.22, "TCA": 0.15, "TCG": 0.05, "AGT": 0.15, "AGC": 0.24},
    "Y": {"TAT": 0.44, "TAC": 0.56},
    "C": {"TGT": 0.46, "TGC": 0.54},
    "L": {"TTA": 0.08, "TTG": 0.13, "CTT": 0.13, "CTC": 0.20, "CTA": 0.07, "CTG": 0.40},
    "*": {"TAA": 0.30, "TGA": 0.47, "TAG": 0.24},
    "W": {"TGG": 1.00},
    "P": {"CCT": 0.29, "CCC": 0.32, "CCA": 0.28, "CCG": 0.11},
    "H": {"CAT": 0.42, "CAC": 0.58},
    "Q": {"CAA": 0.27, "CAG": 0.73},
    "R": {"CGT": 0.08, "CGC": 0.18, "CGA": 0.11, "CGG": 0.20, "AGA": 0.21, "AGG": 0.21},
    "I": {"ATT": 0.36, "ATC": 0.47, "ATA": 0.17},
    "T": {"ACT": 0.25, "ACC": 0.36, "ACA": 0.28, "ACG": 0.11},
    "N": {"AAT": 0.47, "AAC": 0.53},
    "K": {"AAA": 0.43, "AAG": 0.57},
    "M": {"ATG": 1.00},
    "V": {"GTT": 0.18, "GTC": 0.24, "GTA": 0.12, "GTG": 0.46},
    "A": {"GCT": 0.27, "GCC": 0.40, "GCA": 0.23, "GCG": 0.11},
    "D": {"GAT": 0.46, "GAC": 0.54},
    "E": {"GAA": 0.42, "GAG": 0.58},
    "G": {"GGT": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25}
}

def select_codon(aa, current_codon):
    # Get codons sorted by usage frequency (highest first)
    sorted_codons = sorted(
        human_codon_usage[aa].items(),
        key=lambda x: x[1],
        reverse=True
    )
    for codon, _ in sorted_codons:
        if codon != current_codon.upper():
            return codon
    # fallback: return current codon if no alternative (rare)
    return current_codon.upper()

# =========================
# restriction sites
# =========================

def has_bsmBI_site(seq: str, motifs=["CGTCTC", "GAGACG"]) -> bool:
    seq = seq.upper()
    return any(motif in seq for motif in motifs)


def replace_restriction_sites(var_seq_series, var_aa_series, motifs=["CGTCTC","GAGACG"], max_iter=50):
    res = var_seq_series.copy()
    pattern = "|".join(motifs)

    for _ in range(max_iter):
        matches = res.str.contains(pattern, case=False)
        if not matches.any():
            break

        for idx in res[matches].index:
            dna_seq = res[idx]
            aa_seq = var_aa_series[idx]

            for m in re.finditer(pattern, dna_seq, flags=re.IGNORECASE):
                start, end = m.start(), m.end()
                codon_start = (start // 3) * 3
                codon_end = ((end-1) // 3) * 3 + 3

                for pos in range(codon_start, codon_end, 3):
                    aa_idx = pos // 3
                    if aa_idx >= len(aa_seq):
                        continue
                    current_codon = dna_seq[pos:pos+3]
                    replacement_codon = select_codon(aa_seq[aa_idx], current_codon)
                    dna_seq = dna_seq[:pos] + replacement_codon + dna_seq[pos+3:]

            res[idx] = dna_seq.upper()

    return res

# =========================
# genbank record
# =========================

def add_cdr3_feature(features, parent_start, nt_sequence, aa_sequence, label):
    """
    Adds a CDR3 feature by locating the AA sequence inside translated NT sequence.
    
    parent_start = start position of the variable region in full construct
    nt_sequence = nucleotide sequence of TRB or TRA variable region
    aa_sequence = CDR3 amino acid sequence
    """

    if not aa_sequence:
        return

    translated = str(Seq(nt_sequence).translate(to_stop=False))

    aa_index = translated.find(aa_sequence)
    if aa_index == -1:
        return  # CDR3 not found in translation

    nt_start = parent_start + aa_index * 3
    nt_end = nt_start + len(aa_sequence) * 3

    features.append(
        SeqFeature(
            FeatureLocation(nt_start, nt_end),
            type="misc_feature",
            qualifiers={
                "label": label,
                "translation": aa_sequence
            }
        )
    )

def make_genbank_records(TCR_df, constants, logger):

    gb_files = {
        "Clonal_Gene": {},
        "BsmBI_Golden_Gate": {}
    }

    def build_record(seq, record_id, description, feature_blocks, row=None):
        record = SeqRecord(
            Seq(seq),
            id=record_id,
            name=record_id,
            description=description
        )
        record.annotations["molecule_type"] = "DNA"
        record.annotations["topology"] = "linear"
        record.annotations["date"] = datetime.today().strftime("%d-%b-%Y").upper()

        features = []
        pos = 0

        for fragment, label in feature_blocks:
            if not fragment:
                continue

            length = len(fragment)
            start_pos = pos
            end_pos = pos + length

            features.append(
                SeqFeature(
                    FeatureLocation(start_pos, end_pos),
                    type="misc_feature",
                    qualifiers={"label": label}
                )
            )

            # --- Add CDR3 annotation if variable region ---
            if label == "TRB_variable" and row is not None:
                add_cdr3_feature(
                    features,
                    parent_start=start_pos,
                    nt_sequence=fragment,
                    aa_sequence=row.get("TRB_CDR3"),
                    label="TRB_CDR3"
                )

            if label == "TRA_variable" and row is not None:
                add_cdr3_feature(
                    features,
                    parent_start=start_pos,
                    nt_sequence=fragment,
                    aa_sequence=row.get("TRA_CDR3"),
                    label="TRA_CDR3"
                )

            pos = end_pos

        record.features = features
        return record

    for _, row in TCR_df.iterrows():

        # -------------------------
        # Clonal construct
        # -------------------------
        standard_blocks = [
            (constants["hom_5"], "hom_5"),
            (constants["Linker_1"], "Linker_1"),
            (constants["T2A"], "T2A"),
            (row["TRB_nt_trim"], "TRB_variable"),
            (row["TRBC_nt"], "TRBC"),
            (constants["Furin"], "Furin"),
            (constants["Linker_2"], "Linker_2"),
            (constants["P2A"], "P2A"),
            (row["TRA_nt_trim"], "TRA_variable"),
            (constants["TRAC_const"], "TRAC"),
            (constants["hom_3"], "hom_3"),
        ]

        record_std = build_record(
            seq=row["complete_insert_no_BsmBI"],
            record_id=row["TCR_name"],
            description="HDRT construct",
            feature_blocks=standard_blocks,
            row=row
        )

        buffer = io.StringIO()
        SeqIO.write(record_std, buffer, "genbank")
        gb_files["Clonal_Gene"][f"{row['TCR_name']}.gb"] = buffer.getvalue().encode()

        # -------------------------
        # BsmBI Golden Gate construct
        # -------------------------
        bsmBI_blocks = [
            (constants["BsmBI_5"], "BsmBI_5"),
            (row["TRB_nt_trim"], "TRB_variable"),
            (row["TRBC_nt"], "TRBC"),
            (constants["Furin"], "Furin"),
            (constants["Linker_2"], "Linker_2"),
            (constants["P2A"], "P2A"),
            (row["TRA_nt_trim"], "TRA_variable"),
            (constants["BsmBI_3"], "BsmBI_3"),
        ]

        record_bsmBI = build_record(
            seq=row["BsmBI_fragment"],
            record_id=row["TCR_name"] + "_BsmBI",
            description="HDRT construct with BsmBI sites",
            feature_blocks=bsmBI_blocks,
            row=row
        )

        buffer_bsmBI = io.StringIO()
        SeqIO.write(record_bsmBI, buffer_bsmBI, "genbank")
        gb_files["BsmBI_Golden_Gate"][f"{row['TCR_name']}_BsmBI.gb"] = buffer_bsmBI.getvalue().encode()

    return gb_files

# =========================
# HDRT Assembly
# =========================

def assemble_tcr_hdrt(
    TCR_df,
    lookup_table,
    logger,
    use_D112K=True
):

    logger.info("Running HDRT assembly")
    
    logger.info("Loading lookup table")
    lookup = pd.read_excel(lookup_table)
    seq_dict = pd.Series(lookup.sequence.values, index=lookup.name).to_dict()

    def get_seq(name):
        if name not in seq_dict:
            raise ValueError(f"Sequence lookup failed for: {name}")
        return seq_dict[name].upper()

    # -------------------------
    # Constants
    # -------------------------
    
    logger.info("Preparing constants")

    labels = [
        "hom_5","Linker_1","T2A","Furin",
        "Linker_2","P2A","TRAC","hom_3",
        "BsmBI_5","BsmBI_3","tra_c_terminal"
    ]

    constants = {label: get_seq(label) for label in labels}
    constants["TRAC_const"] = constants.pop("TRAC")

    # -------------------------
    # TRBC mapping
    # -------------------------
    logger.info("Mapping TRBC regions")

    TCR_df["TRBC_nt"] = np.select(
        [
            TCR_df["TRBC"] == "TRBC1*01",
            TCR_df["TRBC"] == "TRBC2*01",
            TCR_df["TRBC"] == "No TRBC"
        ],
        [
            get_seq("TRBC1_D112K" if use_D112K else "TRBC1"),
            get_seq("TRBC2_D112K" if use_D112K else "TRBC2"),
            ""
        ],
        default=""
    )

    # -------------------------
    # Trim sequences
    # -------------------------
    logger.info("Trimming TRA/TRB sequences")

    TCR_df["TRA_nt_trim"] = (
        TCR_df["TRA_nt"]
        .str.replace(r"ATATCCAGAACC.*", "", regex=True)
        .str.upper()
    )

    TCR_df["TRB_nt_trim"] = (
        TCR_df["TRB_nt"]
        .str.replace(r"GAGGACCTGAA.*", "", regex=True)
        .str.upper()
    )

    # -------------------------
    # Assemble constructs
    # -------------------------
    logger.info("Assembling HDRT constructs")

    TCR_df["synthesis_5_arm"] = (
        constants["hom_5"] +
        constants["Linker_1"] +
        constants["T2A"]
    ).upper()

    TCR_df["synthesis_3_arm"] = (
        constants["TRAC_const"] +
        constants["hom_3"]
    ).upper()

    TCR_df["variable_synthesis"] = (
        TCR_df["TRB_nt_trim"] +
        TCR_df["TRBC_nt"] +
        constants["Furin"] +
        constants["Linker_2"] +
        constants["P2A"] +
        TCR_df["TRA_nt_trim"]
    ).str.upper()

    TCR_df["complete_insert"] = (
        TCR_df["synthesis_5_arm"] +
        TCR_df["variable_synthesis"] +
        TCR_df["synthesis_3_arm"]
    ).str.upper()

    # -------------------------
    # Translation
    # -------------------------
    logger.info("Translating constructs")

    tra_c_terminal = (constants["tra_c_terminal"]).upper()
    TCR_df["mature_nt"] = (
        TCR_df["variable_synthesis"] + tra_c_terminal
    ).str.upper()

    TCR_df["mature_aa"] = TCR_df["mature_nt"].apply(lambda x: str(Seq(x).translate()))
    TCR_df["variable_synthesis_aa"] = TCR_df["variable_synthesis"].apply(lambda x: str(Seq(x).translate()))

    # -------------------------
    # Reading frame check
    # -------------------------
    logger.info("Checking reading frames")

    # - 1 since alpha trac locus starts out of frame 
    seq_len = TCR_df["variable_synthesis"].str.len() - 1
    mask = (seq_len % 3) != 0

    if mask.any():
        n_bad = mask.sum()
        logger.warning(f"{n_bad} construct(s) are not in frame")

        report = TCR_df.loc[mask, ["TCR_name"]].copy()
        report["length"] = seq_len[mask]
        report["length_mod3"] = seq_len[mask] % 3

        logger.warning("\n" + report.to_string(index=False))
    else:
        logger.info("All constructs are in frame")

    # -------------------------
    # BsmBI detection
    # -------------------------
    logger.info("Scanning for BsmBI restriction sites")

    motifs = ["CGTCTC", "GAGACG"]

    TCR_df["has_BsmBI_site"] = TCR_df["variable_synthesis"].apply(
        lambda x: has_bsmBI_site(x, motifs=motifs)
    )

    # Log if restriction site is found
    for name, has_site in zip(TCR_df["TCR_name"], TCR_df["has_BsmBI_site"]):
        if has_site:
            logger.info(f"{name}: Restriction site found")
        else:
            logger.info(f"{name}: No restriction site detected")

    # Remove sites
    logger.info("Removing BsmBI restriction sites where necessary")

    TCR_df["variable_synthesis_no_BsmBI"] = replace_restriction_sites(
        TCR_df["variable_synthesis"],
        TCR_df["variable_synthesis_aa"],
        motifs=motifs
    )

    TCR_df["complete_insert_no_BsmBI"] = (
        TCR_df["synthesis_5_arm"] +
        TCR_df["variable_synthesis_no_BsmBI"] +
        TCR_df["synthesis_3_arm"]
    ).str.upper()

    # generate final fragment with BsmBI sites 
    TCR_df["BsmBI_fragment"] = (
        constants["BsmBI_5"] + 
        TCR_df["variable_synthesis_no_BsmBI"] + 
        constants["BsmBI_3"]
    ).str.upper()

    logger.info("HDRT assembly finished successfully")

    return TCR_df, constants


# =========================
# Main pipeline
# =========================

def run_tcr_hdrt_pipeline(df: pd.DataFrame, lookup_table, logger, output_dir, use_D112K=True):
    
    logger.info("Starting TCR HDRT pipeline")

    input_path = os.path.join(output_dir, "input.tsv")
    stitchr_path = os.path.join(output_dir, "stitchr.tsv")

    df.to_csv(input_path, sep="\t", index=False)

    logger.info("Running STITCHR (thimble)")

    try:
        subprocess.run(
            ["thimble", "-in", input_path, "-o", stitchr_path, "-r", "AB", "-s", "human"],
            check=True,
            capture_output=True,
            text=True
        )

    except FileNotFoundError:
        logger.error(
            "STITCHR (thimble) was not found. "
            "Please ensure it is installed and available in your PATH."
        )
        raise SystemExit(1)

    except subprocess.CalledProcessError as e:
        logger.error("STITCHR (thimble) failed to run successfully")
        
        if e.stdout:
            logger.error("STDOUT:\n" + e.stdout.strip())
        if e.stderr:
            logger.error("STDERR:\n" + e.stderr.strip())

        raise SystemExit(1)

    TCRs = pd.read_csv(stitchr_path, sep="\t")

    # -------------------------
    # Check STITCHR warnings
    # -------------------------
    if "Warnings/Errors" in TCRs.columns:
        logger.info("Checking STITCHR warnings")

        for _, row in TCRs.iterrows():
            warning_msg = row["Warnings/Errors"]

            # Treat NaN, empty strings, and "[None]" as no warning
            if pd.notna(warning_msg):
                msg = str(warning_msg).strip()
                if msg and msg != "[None]":
                    tcr_name = row.get("TCR_name", "Unknown_TCR")

                    # Stop program if "Error:" is found
                    if "Error:" in msg:
                        logger.error(f"STITCHR error for {tcr_name}: {msg}")
                        raise SystemExit(1)

                    # Otherwise log as warning
                    logger.warning(f"STITCHR warning for {tcr_name}: {msg}")
    else:
        logger.warning("No 'Warnings/Errors' found in the STITCHR output.")

    # -------------------------
    # Cleanup data 
    # -------------------------
    columns_to_exclude = [
        "Linker", "Link_order",
        "TRA_5_prime_seq", "TRA_3_prime_seq",
        "TRB_5_prime_seq", "TRB_3_prime_seq",
        "Linked_nt", "Linked_aa"
    ]
    
    # Drop unwanted columns if they exist
    TCRs = TCRs.drop(columns=[c for c in columns_to_exclude if c in TCRs.columns])

    # -------------------------
    # Build TCR HDRTs  
    # -------------------------

    result_df, constants = assemble_tcr_hdrt(TCRs, lookup_table, logger, use_D112K=use_D112K)

    # -------------------------
    # Writing GenBank files
    # -------------------------
    gb_files = make_genbank_records(result_df, constants, logger)

    return result_df, gb_files

# =========================
# Main terminal interface
# =========================

def main():
    parser = argparse.ArgumentParser(description="TCR HDRT Pipeline")
    parser.add_argument("input_tsv", help="Path to input TSV file with TCR sequences")
    parser.add_argument("-o", "--output_dir", default="TCR_HDRT_output", help="Output folder")
    parser.add_argument("--use_D112K", action="store_true", help="Use D112K TRBC variant")
    parser.add_argument("-l","--lookup_table", default="database/TCRregion_lookup.xlsx",help="Constant TCR region used in the design")
    args = parser.parse_args()
    
    lookup_table = args.lookup_table

    os.makedirs(args.output_dir, exist_ok=True)

    log_path = os.path.join(args.output_dir, "run.log")
    logger = setup_logger(log_path)
    logger.info("Starting TCR HDRT pipeline")
    logger.info(f"Input file: {args.input_tsv}")

    df = pd.read_csv(args.input_tsv, sep="\t")

    result_df, gb_files = run_tcr_hdrt_pipeline(df, output_dir=args.output_dir, logger=logger, lookup_table=lookup_table, use_D112K=args.use_D112K)

    # Organize the order of columns found in result_df before saving it to tsv
    priority_cols = [
        "TCR_name",
        "Warnings/Errors",
        "TRAV",
        "TRAJ",
        "TRA_CDR3",
        "TRBV",
        "TRBJ",
        "TRB_CDR3",
        "TRAC",
        "TRBC",
        "TRA_leader",
        "TRB_leader",
        "complete_insert_no_BsmBI",
        "BsmBI_fragment",
    ]

    # Keep only columns that actually exist (prevents KeyError)
    priority_cols_existing = [col for col in priority_cols if col in result_df.columns]
    remaining_cols = [col for col in result_df.columns if col not in priority_cols_existing]

    # Reorder DataFrame
    result_df = result_df[priority_cols_existing + remaining_cols]

    # Save results
    tsv_path = os.path.join(args.output_dir, "result.tsv")
    result_df.to_csv(tsv_path, sep="\t", index=False)
    logger.info(f"Result saved in: {tsv_path}")

    # Save GenBank files
    gb_base_dir = os.path.join(args.output_dir, "GenBank_files")

    for subfolder, files in gb_files.items():
        subfolder_path = os.path.join(gb_base_dir, subfolder)
        os.makedirs(subfolder_path, exist_ok=True)

        for filename, content in files.items():
            path = os.path.join(subfolder_path, filename)
            with open(path, "wb") as f:
                f.write(content)

        logger.info(f"Saved {len(files)} GenBank files in {subfolder_path}")

    logger.info("TCR HDRT pipeline finished successfully")
    
    print(f"\nAll results saved in folder: {args.output_dir}")

if __name__ == "__main__":
    main()

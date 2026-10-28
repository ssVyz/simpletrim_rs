//! Built-in adapter sequences.
//!
//! The table is taken from fastp's `knownadapters.h` (MIT licensed), which in turn draws on
//! <https://github.com/stephenturner/adapters>. Entries are sorted by sequence (byte order),
//! matching the iteration order of fastp's `std::map`; adapter detection relies on that order
//! to break ties the same way fastp does.

/// `(sequence, description)` pairs, sorted by sequence.
pub const KNOWN_ADAPTERS: &[(&str, &str)] = &[
    ("AACTGTAGGCACCATCAAT", "QIASeq miRNA adapter"),
    (
        "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
        "MGI/BGI adapter (forward)",
    ),
    (
        "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG",
        "MGI/BGI adapter (reverse)",
    ),
    (
        "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA",
        "Illumina DpnII expression PCR Primer 2 | Illumina NlaIII expression PCR Primer 2 | Illumina Small RNA PCR Primer 2 | Illumina DpnII Gex PCR Primer 2 | Illumina NlaIII Gex PCR Primer 2",
    ),
    (
        "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA",
        "Illumina RNA PCR Primer",
    ),
    (
        "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        "TruSeq_Universal_Adapter | PrefixPE/1 | PCR_Primer1 | Illumina Single End PCR Primer 1 | Illumina Paried End PCR Primer 1 | Illumina Multiplexing PCR Primer 1.01 | TruSeq Universal Adapter",
    ),
    (
        "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG",
        "pcr_dimer",
    ),
    (
        "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
        "PCR_Primers",
    ),
    (
        "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        "Illumina Single End Sequencing Primer | Illumina Paired End Adapter 1 | Illumina Paried End Sequencing Primer 1 | Illumina Multiplexing Adapter 2 | Illumina Multiplexing Read1 Sequencing Primer",
    ),
    (
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        "Illumina TruSeq Adapter Read 1",
    ),
    (
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        "PE2_rc | TruSeq3_IndexedAdapter",
    ),
    (
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG",
        "Reverse_adapter",
    ),
    ("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG", "TruSeq2_PE_r"),
    (
        "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG",
        "PCR_Primer2_rc",
    ),
    (
        "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAA",
        "PhiX_read1_adapter",
    ),
    (
        "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        "Illumina TruSeq Adapter Read 2",
    ),
    (
        "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA",
        "PE1_rc | TruSeq3_UniversalAdapter",
    ),
    (
        "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "PCR_Primer1_rc",
    ),
    (
        "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA",
        "PhiX_read2_adapter",
    ),
    ("AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG", "TruSeq2_SE"),
    (
        "CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 35",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC",
        "Illumina PCR Primer Index 10",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 10",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC",
        "Illumina PCR Primer Index 2",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 2",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 38",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 27",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 25",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 31",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 44",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 37",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC",
        "Illumina PCR Primer Index 6",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 6",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC",
        "Illumina PCR Primer Index 5",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 5",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 23",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 30",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 21",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 42",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 33",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
        "PrefixPE/2 | PCR_Primer2 | Illumina Paired End PCR Primer 2",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 22",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC",
        "Illumina PCR Primer Index 1",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 1",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 17",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC",
        "Illumina PCR Primer Index 9",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 9",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 47",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 28",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 45",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC",
        "Illumina PCR Primer Index 7",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 7",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 34",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC",
        "Illumina PCR Primer Index 3",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 3",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 18",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 24",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 26",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 43",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 14",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 16",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 20",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC",
        "Illumina PCR Primer Index 11",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 11",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 39",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 41",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC",
        "Illumina PCR Primer Index 12",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 12",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 29",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC",
        "Illumina PCR Primer Index 8",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 8",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 46",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 40",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 15",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 32",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 48",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC",
        "Illumina PCR Primer Index 4",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 4",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 36",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 13",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        "RNA PCR Primer, Index 19",
    ),
    (
        "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
        "Illumina Single End Adapter 2 | Illumina Single End PCR Primer 2",
    ),
    (
        "CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT",
        "ABI Solid3 Adapter B",
    ),
    (
        "CCGACAGGTTCAGAGTTCTACAGTCCGACATG",
        "Illumina NlaIII expression Sequencing Primer | Illumina NlaIII Gex Sequencing Primer",
    ),
    (
        "CCGAGCCCACGAGACAAGAGGCAATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_and_Nextera_Enrichment_N711 | I7_Primer_Nextera_XT_Index_Kit_v2_N711",
    ),
    (
        "CCGAGCCCACGAGACACTCGCTAATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_Index_Kit_v2_N716",
    ),
    (
        "CCGAGCCCACGAGACACTGAGCGATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_Index_Kit_v2_N724",
    ),
    (
        "CCGAGCCCACGAGACAGGCAGAAATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_and_Nextera_Enrichment_N703 | I7_Primer_Nextera_XT_Index_Kit_v2_N703",
    ),
    (
        "CCGAGCCCACGAGACATCTCAGGATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_Index_Kit_v2_N715",
    ),
    (
        "CCGAGCCCACGAGACATGCGCAGATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_Index_Kit_v2_N722",
    ),
    (
        "CCGAGCCCACGAGACCAGAGAGGATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_and_Nextera_Enrichment_N708",
    ),
    (
        "CCGAGCCCACGAGACCCTAAGACATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_Index_Kit_v2_N726",
    ),
    (
        "CCGAGCCCACGAGACCGAGGCTGATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_and_Nextera_Enrichment_N710 | I7_Primer_Nextera_XT_Index_Kit_v2_N710",
    ),
    (
        "CCGAGCCCACGAGACCGATCAGTATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_Index_Kit_v2_N727",
    ),
    (
        "CCGAGCCCACGAGACCGGAGCCTATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_Index_Kit_v2_N720",
    ),
    (
        "CCGAGCCCACGAGACCGTACTAGATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_and_Nextera_Enrichment_N702 | I7_Primer_Nextera_XT_Index_Kit_v2_N702",
    ),
    (
        "CCGAGCCCACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_and_Nextera_Enrichment_N707 | I7_Primer_Nextera_XT_Index_Kit_v2_N707",
    ),
    (
        "CCGAGCCCACGAGACGCGTAGTAATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_Index_Kit_v2_N719",
    ),
    (
        "CCGAGCCCACGAGACGCTACGCTATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_and_Nextera_Enrichment_N709",
    ),
    (
        "CCGAGCCCACGAGACGCTCATGAATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_Index_Kit_v2_N714",
    ),
    (
        "CCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_and_Nextera_Enrichment_N705 | I7_Primer_Nextera_XT_Index_Kit_v2_N705",
    ),
    (
        "CCGAGCCCACGAGACGGAGCTACATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_Index_Kit_v2_N718",
    ),
    (
        "CCGAGCCCACGAGACGTAGAGGAATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_and_Nextera_Enrichment_N712 | I7_Primer_Nextera_XT_Index_Kit_v2_N712",
    ),
    (
        "CCGAGCCCACGAGACTAAGGCGAATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_and_Nextera_Enrichment_N701 | I7_Primer_Nextera_XT_Index_Kit_v2_N701",
    ),
    (
        "CCGAGCCCACGAGACTACGCTGCATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_Index_Kit_v2_N721",
    ),
    (
        "CCGAGCCCACGAGACTAGCGCTCATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_Index_Kit_v2_N723",
    ),
    (
        "CCGAGCCCACGAGACTAGGCATGATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_and_Nextera_Enrichment_N706 | I7_Primer_Nextera_XT_Index_Kit_v2_N706",
    ),
    (
        "CCGAGCCCACGAGACTCCTGAGCATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_and_Nextera_Enrichment_N704 | I7_Primer_Nextera_XT_Index_Kit_v2_N704",
    ),
    (
        "CCGAGCCCACGAGACTCGACGTCATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_Index_Kit_v2_N729",
    ),
    (
        "CCGAGCCCACGAGACTGCAGCTAATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Primer_Nextera_XT_Index_Kit_v2_N728",
    ),
    (
        "CGACAGGTTCAGAGTTCTACAGTCCGACGATC",
        "Illumina DpnII expression Sequencing Primer | Illumina Small RNA Sequencing Primer | Illumina DpnII Gex Sequencing Primer",
    ),
    (
        "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
        "Illumina Paired End Sequencing Primer 2",
    ),
    (
        "CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT",
        "Clontech Universal Primer Mix Long",
    ),
    (
        "CTGAGCGGGCTGGCAAGGCAGACCGATCTCGTATGCCGTCTTCTGCTTG",
        "I7_Adapter_Nextera_No_Barcode",
    ),
    (
        "CTGATGGCGCGAGGGAGGCGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Adapter_Nextera",
    ),
    ("CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG", "ABI Solid3 Adapter A"),
    (
        "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
        "I7_Nextera_Transposase_1 | Trans2_rc",
    ),
    (
        "CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGC",
        "I7_Nextera_Transposase_2",
    ),
    (
        "CTGTCTCTTATACACATCTCTGATGGCGCGAGGGAGGC",
        "I5_Nextera_Transposase_2",
    ),
    (
        "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
        "I5_Nextera_Transposase_1 | Trans1_rc",
    ),
    (
        "GACGCTGCCGACGAACTCTAGGGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_Index_Kit_v2_S516",
    ),
    (
        "GACGCTGCCGACGAAGAGGATAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]503 | I5_Primer_Nextera_XT_Index_Kit_v2_S503",
    ),
    (
        "GACGCTGCCGACGAAGCTAGAAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_Index_Kit_v2_S515",
    ),
    (
        "GACGCTGCCGACGAAGGCTTAGGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]508 | I5_Primer_Nextera_XT_Index_Kit_v2_S508",
    ),
    (
        "GACGCTGCCGACGAATAGAGAGGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]502 | I5_Primer_Nextera_XT_Index_Kit_v2_S502",
    ),
    (
        "GACGCTGCCGACGAATAGCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_Index_Kit_v2_S520",
    ),
    (
        "GACGCTGCCGACGAATTAGACGGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_Index_Kit_v2_S510",
    ),
    (
        "GACGCTGCCGACGACGGAGAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_Index_Kit_v2_S511",
    ),
    (
        "GACGCTGCCGACGACTAGTCGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_Index_Kit_v2_S513",
    ),
    (
        "GACGCTGCCGACGACTCCTTACGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]505 | I5_Primer_Nextera_XT_Index_Kit_v2_S505",
    ),
    (
        "GACGCTGCCGACGACTTAATAGGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_Index_Kit_v2_S518",
    ),
    (
        "GACGCTGCCGACGAGCGATCTAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]501",
    ),
    (
        "GACGCTGCCGACGATAAGGCTCGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_Index_Kit_v2_S521",
    ),
    (
        "GACGCTGCCGACGATACTCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]507 | I5_Primer_Nextera_XT_Index_Kit_v2_S507",
    ),
    (
        "GACGCTGCCGACGATATGCAGTGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]506 | I5_Primer_Nextera_XT_Index_Kit_v2_S506",
    ),
    (
        "GACGCTGCCGACGATCGCATAAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_Index_Kit_v2_S522",
    ),
    (
        "GACGCTGCCGACGATCTACTCTGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]504",
    ),
    (
        "GACGCTGCCGACGATCTTACGCGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]517 | I5_Primer_Nextera_XT_Index_Kit_v2_S517",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        "Nextera_LMP_Read1_External_Adapter | Illumina Multiplexing Index Sequencing Primer",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_5 | TruSeq Adapter, Index 5",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_25",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq Adapter, Index 25",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_8 | TruSeq Adapter, Index 8",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_13",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq Adapter, Index 13",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_14",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq Adapter, Index 14",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_1_6 | TruSeq Adapter, Index 1",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_15",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGTCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq Adapter, Index 15",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_27",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq Adapter, Index 27",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_7 | TruSeq Adapter, Index 7",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCACTCTTCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq Adapter, Index 23",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_16",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCTCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq Adapter, Index 16",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_2 | TruSeq Adapter, Index 2",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_22",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTTCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq Adapter, Index 22",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_12 | TruSeq Adapter, Index 12",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_23",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_9 | TruSeq Adapter, Index 9",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_6 | TruSeq Adapter, Index 6",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_11 | TruSeq Adapter, Index 11",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_18_7",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq Adapter, Index 18",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_19",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACTCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq Adapter, Index 19",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_20",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq Adapter, Index 20",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGAATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_21",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGTCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq Adapter, Index 21",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_10 | TruSeq Adapter, Index 10",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_4 | TruSeq Adapter, Index 4",
    ),
    (
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG",
        "TruSeq_Adapter_Index_3 | TruSeq Adapter, Index 3",
    ),
    (
        "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
        "Illumina Paired End Adapter 2",
    ),
    (
        "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        "Nextera_LMP_Read2_External_Adapter",
    ),
    (
        "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
        "Illumina Single End Adapter 1",
    ),
    (
        "GATCGTCGGACTGTAGAACTCTGAACGTGTAGA",
        "Illumina Small RNA Adapter Read 2",
    ),
    ("GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG", "Trans2"),
    (
        "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
        "PrefixPE/2 | PE2 | Illumina Multiplexing PCR Primer 2.01 | Illumina Multiplexing Read2 Sequencing Primer",
    ),
    ("TACACTCTTTCCCTACACGACGCTCTTCCGATCT", "PrefixPE/1 | PE1"),
    (
        "TCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "RNA_PCR_Primer_(RP1)_part_#_15013198",
    ),
    ("TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", "Trans1"),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_5_(RPI5)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_25_(RPI25)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_8_(RPI8)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTCAAATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_13_(RPI13)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTTCCATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_14_(RPI14)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_1_(RPI1)_2,9",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGAGCATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_26_(RPI26)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGTCAATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_15_(RPI15)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATTCCTATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_27_(RPI27)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAAAAGATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_28_(RPI28)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAACTAATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_29_(RPI29)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACCGGATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_30_(RPI30)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACGATATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_31_(RPI31)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACTCAATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_32_(RPI32)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_7_(RPI7)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGGCGATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_33_(RPI33)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATGGCATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_34_(RPI34)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATTTTATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_35_(RPI35)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCAACAATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_36_(RPI36)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCGTCCATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_16_(RPI16)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_2_(RPI2)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGGAATATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_37_(RPI37)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGTACGATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_22_(RPI22)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTAGCTATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_38_(RPI38)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTATACATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_39_(RPI39)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTCAGAATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_40_(RPI40)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_12_(RPI12)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_41_(RPI41)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGAGTGGATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_23_(RPI23)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_9_(RPI9)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_6_(RPI6)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_11_(RPI11)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGTAGCATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_24_(RPI24)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTAGAGATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_17_(RPI17)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_18_(RPI18)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGAAAATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_19_(RPI19)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGGCCATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_20_(RPI20)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTTTCGATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_21_(RPI21)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAATCGATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_42_(RPI42)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTACAGCATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_43_(RPI43)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_10_(RPI10)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTATAATATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_44_(RPI44)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCATTCATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_45_(RPI45)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCCCGAATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_46_(RPI46)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGAAGATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_47_(RPI47)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGGCAATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_48_(RPI48)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_4_(RPI4)",
    ),
    (
        "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG",
        "RNA_PCR_Primer_Index_3_(RPI3)",
    ),
    ("TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC", "FlowCell1"),
    ("TTTTTTTTTTCAAGCAGAAGACGGCATACGA", "FlowCell2"),
];

/// Named adapter presets for `--adapter-preset`: `(name, read 1 adapter, read 2 adapter, description)`.
pub const ADAPTER_PRESETS: &[(&str, &str, &str, &str)] = &[
    (
        "truseq",
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        "Illumina TruSeq",
    ),
    (
        "nextera",
        "CTGTCTCTTATACACATCT",
        "CTGTCTCTTATACACATCT",
        "Illumina Nextera / Tn5 transposase",
    ),
    (
        "smallrna",
        "TGGAATTCTCGGGTGCCAAGG",
        "GATCGTCGGACTGTAGAACTCTGAACGTGTAGA",
        "Illumina small RNA",
    ),
    (
        "mgi",
        "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
        "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG",
        "MGI / BGI DNBSEQ",
    ),
];

/// Description of a known adapter sequence, if the sequence is in the table.
pub fn describe(seq: &[u8]) -> Option<&'static str> {
    KNOWN_ADAPTERS
        .binary_search_by(|(s, _)| s.as_bytes().cmp(seq))
        .ok()
        .map(|i| KNOWN_ADAPTERS[i].1)
}

/// Looks up an adapter preset by (case-insensitive) name.
pub fn preset(name: &str) -> Option<(&'static str, &'static str)> {
    ADAPTER_PRESETS
        .iter()
        .find(|(n, ..)| n.eq_ignore_ascii_case(name))
        .map(|(_, r1, r2, _)| (*r1, *r2))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_is_sorted_and_unique() {
        assert!(KNOWN_ADAPTERS.windows(2).all(|w| w[0].0 < w[1].0));
    }

    #[test]
    fn describe_finds_truseq() {
        assert_eq!(
            describe(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"),
            Some("Illumina TruSeq Adapter Read 1")
        );
        assert_eq!(describe(b"ACGT"), None);
    }

    #[test]
    fn presets_resolve() {
        assert!(preset("TruSeq").is_some());
        assert!(preset("unknown").is_none());
    }
}

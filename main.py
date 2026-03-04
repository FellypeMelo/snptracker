"""
SNPTracker - Detector de SNPs

Propósito: Comparar uma sequência de DNA de referência com uma ou mais
amostras para detectar mutações pontuais (SNPs - Single Nucleotide
Polymorphisms) e indels (inserções e deleções).

Módulos relacionados:
    fasta_parser.py  — leitura de arquivos FASTA
    annotation.py    — anotação funcional baseada no código genético padrão

Formato do dict de SNP retornado por detect_snps():
    {
        "position":  int   — posição 1-indexed na sequência
        "reference": str   — base da referência (A/C/G/T ou '-' para inserções)
        "alternate": str   — base da amostra    (A/C/G/T ou '-' para deleções)
        "type":      str   — "TRANSITION" | "TRANSVERSION" | "INSERTION" | "DELETION"
        "annotation":str   — "SYNONYMOUS" | "NON_SYNONYMOUS" | "NONSENSE" | "NON_CODING"
        "context":   str   — contexto COSMIC ex: "T[A>G]C"  (ausente em INDELs)
    }

Uso via CLI:
    # Par único (strings ou arquivos FASTA)
    python main.py --reference "ACTG" --sample "ACTT"
    python main.py --reference ref.fasta --sample sample.fasta --output out.txt

    # Com regiões codificantes explícitas (1-indexed, inclusive)
    python main.py --reference ref.fasta --sample sample.fasta --cds "1-90"
    python main.py --reference ref.fasta --sample sample.fasta --cds "1-90,100-150"

    # Multi-amostra (primeira sequência do FASTA = referência)
    python main.py --input data/sequences.txt
"""


import argparse
import os
from fasta_parser import read_fasta
from annotation import annotate_snp, annotate_snp_with_regions


def parse_cds_regions(cds_str: str) -> list[tuple[int, int]]:
    """Parses a CDS string like '1-90,100-150' into a list of tuples.

    Args:
        cds_str: Comma-separated ranges in 'start-end' format (1-indexed).

    Returns:
        list[tuple[int, int]]: List of (start, end) tuples.

    Raises:
        ValueError: If the format is invalid or start > end.
    """
    regions = []
    for part in cds_str.split(","):
        part = part.strip()
        try:
            start_str, end_str = part.split("-")
            start, end = int(start_str), int(end_str)
        except ValueError:
            raise ValueError(
                f"Invalid CDS region '{part}'. Expected format: 'start-end' "
                f"(e.g. '1-90')."
            )
        if start > end:
            raise ValueError(
                f"Invalid CDS region '{part}': start ({start}) must be "
                f"less than or equal to end ({end})."
            )
        regions.append((start, end))
    return regions


def get_trinucleotide_context(
    reference: str,
    position: int,
    ref_base: str,
    alt_base: str,
) -> str:
    """Returns the trinucleotide context of a SNP in COSMIC format.

    Args:
        reference: Full reference sequence (uppercase).
        position: 1-indexed SNP position.
        ref_base: Reference base at this position.
        alt_base: Alternate base at this position.

    Returns:
        str: Context string in format 'X[R>A]Y', where X and Y are
             the immediate flanking bases, or '_' at sequence edges.
    """
    idx = position - 1
    prev = reference[idx - 1] if idx > 0 else "_"
    next_ = reference[idx + 1] if idx < len(reference) - 1 else "_"
    return f"{prev}[{ref_base}>{alt_base}]{next_}"


def detect_snps(
    reference: str,
    sample: str,
    cds_regions: list[tuple[int, int]] | None = None,
) -> list[dict]:
    """
    Compara duas sequências e identifica SNPs e indels.

    Realiza comparação posição a posição até o comprimento mínimo das
    sequências (SNPs e substituições). Diferenças de comprimento são
    reportadas como INSERTION ou DELETION ao final da lista.

    INDELs recebem annotation='NON_CODING' e não possuem a chave 'context'.
    SNPs recebem todas as chaves, incluindo 'context' no formato COSMIC.

    Args:
        reference: Sequência de referência (string ou uppercase).
        sample: Sequência da amostra (string ou uppercase).
        cds_regions: Lista opcional de regiões codificantes como tuplas
            (start, end) 1-indexed inclusive. Se None, toda a sequência
            é tratada como codificante (comportamento padrão). Se lista
            vazia, todos os SNPs recebem annotation='NON_CODING'.

    Returns:
        list[dict]: Lista de SNPs. Ver formato no topo deste módulo.
    """
    snps = []

    # Normaliza para maiúsculas
    ref = str(reference).upper()
    smp = str(sample).upper()

    # Determina o comprimento mínimo
    min_length = min(len(ref), len(smp))

    # Compara base a base
    for i in range(min_length):
        ref_base = ref[i]
        smp_base = smp[i]

        if ref_base != smp_base:
            snp_info = {
                "position": i + 1,  # Posição 1-indexed
                "reference": ref_base,
                "alternate": smp_base,
                "type": classify_mutation(ref_base, smp_base),
                "annotation": annotate_snp_with_regions(
                    i + 1, ref, smp, cds_regions
                ),
                "context": get_trinucleotide_context(
                    ref, i + 1, ref_base, smp_base
                ),
            }
            snps.append(snp_info)

    # Detecta diferenças de tamanho — reportadas como INDELs sem 'context'
    if len(ref) != len(smp):
        if len(ref) > len(smp):
            for i in range(min_length, len(ref)):
                snps.append(
                    {
                        "position": i + 1,
                        "reference": ref[i],
                        "alternate": "-",
                        "type": "DELETION",
                        "annotation": "NON_CODING",
                    }
                )
        else:
            for i in range(min_length, len(smp)):
                snps.append(
                    {
                        "position": i + 1,
                        "reference": "-",
                        "alternate": smp[i],
                        "type": "INSERTION",
                        "annotation": "NON_CODING",
                    }
                )

    return snps


def classify_mutation(ref_base: str, alt_base: str) -> str:
    """
    Classifica o tipo de mutação.

    Args:
        ref_base: Base de referência
        alt_base: Base alternativa

    Returns:
        str: Tipo de mutação
    """
    purines = {"A", "G"}
    pyrimidines = {"C", "T"}

    if (ref_base in purines and alt_base in purines) or (
        ref_base in pyrimidines and alt_base in pyrimidines
    ):
        return "TRANSITION"
    else:
        return "TRANSVERSION"


def print_snp_report(snps: list[dict], reference: str, sample: str) -> None:
    """
    Imprime relatório de SNPs formatado.

    Args:
        snps: Lista de SNPs detectados
        reference: Sequência de referência
        sample: Sequência da amostra
    """
    print("=" * 60)
    print("SNPTracker - Relatório de Mutações")
    print("=" * 60)
    print(f"\nReferência: {reference}")
    print(f"Amostra:    {sample}")
    print(f"\nTotal de variações encontradas: {len(snps)}")

    if snps:
        print("\n" + "-" * 70)
        print(
            f"{'Posição':<10} {'Ref':<5} {'Alt':<5} {'Tipo':<15} "
            f"{'Anotação':<20} {'Contexto'}"
        )
        print("-" * 70)

        for snp in snps:
            print(
                f"{snp['position']:<10} {snp['reference']:<5} "
                f"{snp['alternate']:<5} {snp['type']:<15} "
                f"{snp.get('annotation', ''):<20} "
                f"{snp.get('context', 'N/A')}"
            )
    else:
        print("\nNenhuma variação detectada (sequências idênticas)")


def generate_snp_file(snps: list[dict], output_file: str = "snps_report.txt") -> None:
    """
    Salva relatório em arquivo.

    Args:
        snps: Lista de SNPs
        output_file: Nome do arquivo de saída
    """
    with open(output_file, "w") as f:
        f.write("SNPTRACKER - RELATÓRIO DE SNPs\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Total de SNPs: {len(snps)}\n\n")
        f.write(f"{'Posição':<10} {'Ref':<5} {'Alt':<5} {'Tipo':<15} {'Anotação':<20} {'Contexto'}\n")
        f.write("-" * 70 + "\n")

        for snp in snps:
            f.write(
                f"{snp['position']:<10} {snp['reference']:<5} "
                f"{snp['alternate']:<5} {snp['type']:<15} "
                f"{snp.get('annotation', ''):<20} "
                f"{snp.get('context', 'N/A')}\n"
            )

    print(f"\nRelatório salvo em: {output_file}")


def run_multi_sample(
    reference: str,
    samples: list[tuple[str, str]],
) -> dict[str, list[dict]]:
    """
    Runs SNP detection for each sample against a reference sequence.

    Args:
        reference: Reference DNA sequence.
        samples: List of (name, sequence) tuples.

    Returns:
        dict[str, list[dict]]: Mapping of sample name to its SNP list.
    """
    return {name: detect_snps(reference, sequence) for name, sequence in samples}


def load_sequence(input_data: str) -> str:
    """
    Loads sequence from a file if it exists, otherwise returns the string.

    Args:
        input_data: File path or raw sequence string.

    Returns:
        str: DNA sequence.
    """
    if os.path.isfile(input_data):
        return read_fasta(input_data)
    return input_data


def parse_args(args: list[str] | None = None) -> argparse.Namespace:
    """
    Parses command line arguments.

    Supports two modes:
    - Single-pair mode: --reference and --sample (raw strings or FASTA files)
    - Multi-sample mode: --input (a multi-sequence FASTA file)

    Args:
        args: List of arguments to parse. If None, uses sys.argv.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="SNPTracker - Detector de SNPs")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--reference", help="Sequência de referência (DNA ou arquivo FASTA)"
    )
    group.add_argument(
        "--input",
        help="Arquivo FASTA com múltiplas sequências (primeira = referência)",
    )
    parser.add_argument(
        "--sample",
        default=None,
        help="Sequência da amostra (DNA ou arquivo FASTA). Obrigatório com --reference.",
    )
    parser.add_argument(
        "--output", default="snps_report.txt", help="Nome do arquivo de saída"
    )
    parser.add_argument(
        "--cds",
        default=None,
        help=(
            "Regiões codificantes (CDS) no formato 'start-end,start-end' "
            "(1-indexed, ex: '1-90,100-150'). Apenas aplicável com "
            "--reference. Se omitido, toda a sequência é tratada como "
            "codificante."
        ),
    )
    namespace = parser.parse_args(args)
    if namespace.reference is not None and namespace.sample is None:
        parser.error("--sample é obrigatório quando --reference é utilizado.")
    return namespace


def main():
    """Função principal do programa."""
    args = parse_args()

    if args.input:
        _run_multi_sample_mode(args)
    else:
        _run_single_sample_mode(args)

    print("\nAnálise concluída!")


def _run_single_sample_mode(args: argparse.Namespace) -> None:
    """Executes the original single reference vs single sample flow."""
    reference = load_sequence(args.reference)
    sample = load_sequence(args.sample)

    cds_regions = None
    if args.cds:
        cds_regions = parse_cds_regions(args.cds)

    snps = detect_snps(reference, sample, cds_regions=cds_regions)
    print_snp_report(snps, reference, sample)

    if snps:
        generate_snp_file(snps, output_file=args.output)


def _run_multi_sample_mode(args: argparse.Namespace) -> None:
    """Executes multi-sample analysis from a single multi-sequence FASTA file."""
    from fasta_parser import read_all_sequences

    sequences = read_all_sequences(args.input)
    if not sequences:
        print("Erro: nenhuma sequência encontrada no arquivo.")
        return

    ref_header, ref_seq = sequences[0]
    samples = sequences[1:]

    print("=" * 60)
    print("SNPTracker - Análise Multi-Amostra")
    print("=" * 60)
    print(f"Referência: {ref_header} ({len(ref_seq)} bp)")
    print(f"Amostras:   {len(samples)}\n")

    if not samples:
        print("Aviso: apenas uma sequência encontrada. Nenhuma amostra para comparar.")
        return

    results = run_multi_sample(ref_seq, samples)

    output_prefix = args.output.replace(".txt", "")
    for i, (name, snps) in enumerate(results.items(), start=1):
        count = len(snps)
        if count > 0:
            output_file = f"{output_prefix}_{name.split()[0]}.txt"
            generate_snp_file(snps, output_file=output_file)
            print(f"[{i}/{len(results)}] {name} → {count} SNP(s) → salvo em {output_file}")
        else:
            print(f"[{i}/{len(results)}] {name} → 0 SNPs")


if __name__ == "__main__":
    main()

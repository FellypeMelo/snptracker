"""
SNPTracker - Detector de SNPs

Propósito: Comparar uma sequência de referência com uma amostra
para detectar mutações pontuais (SNPs - Single Nucleotide Polymorphisms).

Um SNP é uma variação em uma única posição nucleotídica entre
diferentes indivíduos ou amostras.
"""


import argparse
import os
from fasta_parser import read_fasta


def detect_snps(reference: str, sample: str) -> list[dict]:
    """
    Compara duas sequências e identifica SNPs.

    Args:
        reference: Sequência de referência (string)
        sample: Sequência da amostra (string)

    Returns:
        list: Lista de dicionários com informações dos SNPs
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
            }
            snps.append(snp_info)

    # Detecta diferenças de tamanho (indels)
    if len(ref) != len(smp):
        if len(ref) > len(smp):
            for i in range(min_length, len(ref)):
                snps.append(
                    {
                        "position": i + 1,
                        "reference": ref[i],
                        "alternate": "-",
                        "type": "DELETION",
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
        print("\n" + "-" * 60)
        print(f"{'Posição':<10} {'Ref':<5} {'Alt':<5} {'Tipo':<15}")
        print("-" * 60)

        for snp in snps:
            print(
                f"{snp['position']:<10} {snp['reference']:<5} "
                f"{snp['alternate']:<5} {snp['type']:<15}"
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
        f.write(f"{'Posição':<10} {'Ref':<5} {'Alt':<5} {'Tipo':<15}\n")
        f.write("-" * 60 + "\n")

        for snp in snps:
            f.write(
                f"{snp['position']:<10} {snp['reference']:<5} "
                f"{snp['alternate']:<5} {snp['type']:<15}\n"
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

    snps = detect_snps(reference, sample)
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

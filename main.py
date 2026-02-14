"""
SNPTracker - Detector de SNPs

Propósito: Comparar uma sequência de referência com uma amostra
para detectar mutações pontuais (SNPs - Single Nucleotide Polymorphisms).

Um SNP é uma variação em uma única posição nucleotídica entre
diferentes indivíduos ou amostras.
"""


def detect_snps(reference, sample):
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


def classify_mutation(ref_base, alt_base):
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


def print_snp_report(snps, reference, sample):
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


def generate_snp_file(snps, output_file="snps_report.txt"):
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


def main():
    """Função principal do programa."""
    # Sequências de exemplo com SNPs
    reference = "ACTGCTAGCTAGCTA"
    sample = "ACTGCTGGCTAGATA"

    # Detecta SNPs
    snps = detect_snps(reference, sample)

    # Imprime relatório
    print_snp_report(snps, reference, sample)

    # Gera arquivo de saída
    if snps:
        generate_snp_file(snps)

    print("\nAnálise concluída!")


if __name__ == "__main__":
    main()

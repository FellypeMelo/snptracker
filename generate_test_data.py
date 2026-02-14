#!/usr/bin/env python3
"""
SNPTracker - Gerador de Dados de Teste

Este script gera 50+ conjuntos de dados com SNPs conhecidos.
Útil para validar se o detector está identificando mutações corretamente.

Os dados de teste são COMMITADOS no GitHub.
Para dados reais, use a pasta data/ (gitignored)
"""

import random
import os
from datetime import datetime

TEST_DATA_DIR = "test_data"
NUM_DATASETS = 55


def generate_sequence(length, gc_content=0.5):
    """Gera uma sequência de DNA aleatória."""
    seq = []
    for _ in range(length):
        if random.random() < gc_content:
            seq.append(random.choice(["G", "C"]))
        else:
            seq.append(random.choice(["A", "T"]))
    return "".join(seq)


def introduce_snps(reference, num_snps=5):
    """
    Introduz SNPs conhecidos em uma sequência.

    Returns:
        tuple: (sequência_mutada, lista_de_snps)
        lista_de_snps: [(posição, ref_base, alt_base), ...]
    """
    seq_list = list(reference)
    snps = []
    bases = ["A", "T", "G", "C"]

    positions = random.sample(range(len(reference)), min(num_snps, len(reference)))
    positions.sort()

    for pos in positions:
        ref_base = seq_list[pos]
        alt_base = random.choice([b for b in bases if b != ref_base])
        seq_list[pos] = alt_base
        snps.append((pos, ref_base, alt_base))

    return "".join(seq_list), snps


def introduce_indels(reference, num_indels=2):
    """Introduz inserções ou deleções."""
    seq_list = list(reference)
    indels = []

    for _ in range(num_indels):
        pos = random.randint(0, len(seq_list) - 1)
        if random.random() < 0.5:
            # Deleção
            deleted_base = seq_list.pop(pos)
            indels.append((pos, deleted_base, "-", "DELETION"))
        else:
            # Inserção
            new_base = random.choice(["A", "T", "G", "C"])
            seq_list.insert(pos, new_base)
            indels.append((pos, "-", new_base, "INSERTION"))

    return "".join(seq_list), indels


def generate_test_datasets():
    """Gera todos os datasets de teste."""
    datasets = []

    # 1-10: Sem SNPs (sequências idênticas)
    for i in range(10):
        length = random.randint(100, 300)
        seq = generate_sequence(length)
        datasets.append(
            (f"snptracker_test_{i + 1:02d}_identical", seq, seq, "No SNPs (identical)")
        )

    # 11-20: Poucos SNPs (1-3)
    for i in range(10, 20):
        length = random.randint(100, 300)
        ref = generate_sequence(length)
        sample, snps = introduce_snps(ref, num_snps=random.randint(1, 3))
        datasets.append(
            (f"snptracker_test_{i + 1:02d}_few_snps", ref, sample, f"{len(snps)} SNPs")
        )

    # 21-30: Médio número de SNPs (5-10)
    for i in range(20, 30):
        length = random.randint(150, 400)
        ref = generate_sequence(length)
        sample, snps = introduce_snps(ref, num_snps=random.randint(5, 10))
        datasets.append(
            (
                f"snptracker_test_{i + 1:02d}_medium_snps",
                ref,
                sample,
                f"{len(snps)} SNPs",
            )
        )

    # 31-35: Muitos SNPs (15-25)
    for i in range(30, 35):
        length = random.randint(200, 500)
        ref = generate_sequence(length)
        sample, snps = introduce_snps(ref, num_snps=random.randint(15, 25))
        datasets.append(
            (f"snptracker_test_{i + 1:02d}_many_snps", ref, sample, f"{len(snps)} SNPs")
        )

    # 36-40: Com indels
    for i in range(35, 40):
        length = random.randint(150, 300)
        ref = generate_sequence(length)
        sample, indels = introduce_indels(ref, num_indels=random.randint(1, 3))
        datasets.append(
            (
                f"snptracker_test_{i + 1:02d}_with_indels",
                ref,
                sample,
                f"{len(indels)} indels",
            )
        )

    # 41-45: SNPs + Indels
    for i in range(40, 45):
        length = random.randint(200, 400)
        ref = generate_sequence(length)
        sample, snps = introduce_snps(ref, num_snps=random.randint(3, 8))
        sample, indels = introduce_indels(sample, num_indels=random.randint(1, 2))
        total = len(snps) + len(indels)
        datasets.append(
            (
                f"snptracker_test_{i + 1:02d}_snps_indels",
                ref,
                sample,
                f"{total} variants",
            )
        )

    # 46-50: Cenários biológicos
    # Homozigoto vs Heterozigoto (simulado)
    ref = generate_sequence(200)
    sample = list(ref)
    # Muda ~50% das bases em posições aleatórias (simula heterozigoto)
    for pos in random.sample(range(200), 20):
        if random.random() < 0.5:
            bases = ["A", "T", "G", "C"]
            sample[pos] = random.choice([b for b in bases if b != ref[pos]])
    datasets.append(
        ("snptracker_test_46_heterozygous", ref, "".join(sample), "Heterozygous-like")
    )

    # Transitions predominantes
    ref = generate_sequence(200)
    sample = list(ref)
    transitions = 0
    for pos in random.sample(range(200), 15):
        ref_base = ref[pos]
        # Transitions: A<->G ou C<->T
        if ref_base == "A":
            sample[pos] = "G"
            transitions += 1
        elif ref_base == "G":
            sample[pos] = "A"
            transitions += 1
        elif ref_base == "C":
            sample[pos] = "T"
            transitions += 1
        elif ref_base == "T":
            sample[pos] = "C"
            transitions += 1
    datasets.append(
        (
            "snptracker_test_47_transitions",
            ref,
            "".join(sample),
            f"{transitions} transitions",
        )
    )

    # Transversions predominantes
    ref = generate_sequence(200)
    sample = list(ref)
    transversions = 0
    for pos in random.sample(range(200), 15):
        ref_base = ref[pos]
        purines = ["A", "G"]
        pyrimidines = ["C", "T"]
        if ref_base in purines:
            sample[pos] = random.choice(pyrimidines)
            transversions += 1
        else:
            sample[pos] = random.choice(purines)
            transversions += 1
    datasets.append(
        (
            "snptracker_test_48_transversions",
            ref,
            "".join(sample),
            f"{transversions} transversions",
        )
    )

    # CpG islands (alta mutabilidade)
    # Gera sequência rica em CpG
    cpg_rich = []
    for _ in range(100):
        if random.random() < 0.3:
            cpg_rich.append("CG")
        else:
            cpg_rich.append(generate_sequence(2))
    ref = "".join(cpg_rich)
    sample, snps = introduce_snps(ref, num_snps=10)
    datasets.append(
        ("snptracker_test_49_cpg_island", ref, sample, "CpG island (high mutation)")
    )

    # Código genético degenerado (sinônimo vs não-sinônimo)
    ref = generate_sequence(150)
    sample = list(ref)
    # Altera apenas a 3ª posição de códons (geralmente sinônimo)
    synonymous = 0
    for i in range(2, 150, 3):
        if random.random() < 0.3:
            bases = ["A", "T", "G", "C"]
            sample[i] = random.choice([b for b in bases if b != ref[i]])
            synonymous += 1
    datasets.append(
        (
            "snptracker_test_50_synonymous",
            ref,
            "".join(sample),
            f"{synonymous} synonymous changes",
        )
    )

    # 51-55: Casos extremos
    # SNP em cada posição
    ref = generate_sequence(30)
    sample = "".join([random.choice(["A", "T", "G", "C"]) for _ in range(30)])
    datasets.append(
        ("snptracker_test_51_all_mutated", ref, sample, "All positions mutated")
    )

    # Uma única mutação
    ref = generate_sequence(100)
    pos = random.randint(0, 99)
    sample = list(ref)
    bases = ["A", "T", "G", "C"]
    sample[pos] = random.choice([b for b in bases if b != ref[pos]])
    datasets.append(
        ("snptracker_test_52_single_snp", ref, "".join(sample), "Single SNP")
    )

    # Grande deleção
    ref = generate_sequence(200)
    start = 50
    end = 100
    sample = ref[:start] + ref[end:]
    datasets.append(
        ("snptracker_test_53_large_deletion", ref, sample, "Large deletion (50bp)")
    )

    # Grande inserção
    ref = generate_sequence(100)
    insertion = generate_sequence(30)
    pos = 50
    sample = ref[:pos] + insertion + ref[pos:]
    datasets.append(
        ("snptracker_test_54_large_insertion", ref, sample, "Large insertion (30bp)")
    )

    # Múltiplas amostras vs referência
    ref = generate_sequence(100)
    # Gera 3 amostras diferentes
    samples = []
    for _ in range(3):
        sample, _ = introduce_snps(ref, num_snps=random.randint(2, 5))
        samples.append(sample)
    datasets.append(
        ("snptracker_test_55_multi_sample", ref, samples[0], "Sample 1 of 3")
    )
    # Salva as outras amostras também
    datasets.append(("snptracker_test_55b_sample2", ref, samples[1], "Sample 2 of 3"))
    datasets.append(("snptracker_test_55c_sample3", ref, samples[2], "Sample 3 of 3"))

    return datasets


def save_datasets(datasets):
    """Salva os datasets em arquivos."""
    os.makedirs(TEST_DATA_DIR, exist_ok=True)

    manifest_path = os.path.join(TEST_DATA_DIR, "MANIFEST.txt")
    with open(manifest_path, "w") as manifest:
        manifest.write("=" * 70 + "\n")
        manifest.write("SNPTracker - Dados de Teste Sintéticos\n")
        manifest.write("=" * 70 + "\n\n")
        manifest.write(f"Gerado em: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        manifest.write(f"Total de datasets: {len(datasets)}\n\n")
        manifest.write("ATENÇÃO: Estes são dados FABRICADOS com SNPs conhecidos.\n")
        manifest.write("Útil para validar se o detector encontra as mutações.\n\n")
        manifest.write("Formato: Cada arquivo contém referência e amostra\n")
        manifest.write("Lista de arquivos:\n")
        manifest.write("-" * 70 + "\n")

        for filename, ref, sample, description in datasets:
            filepath = os.path.join(TEST_DATA_DIR, f"{filename}.txt")

            with open(filepath, "w") as f:
                f.write(f">reference {description}\n")
                f.write(ref + "\n")
                f.write(f">sample {description}\n")
                f.write(sample + "\n")

            manifest.write(f"{filename}.txt - {description}\n")
            manifest.write(f"  Ref: {len(ref)} bp, Sample: {len(sample)} bp\n")
            print(f"[OK] Gerado: {filename}.txt ({len(ref)}bp vs {len(sample)}bp)")

    print(f"\n[OK] Manifesto salvo em: {manifest_path}")
    print(f"[OK] Total: {len(datasets)} arquivos gerados")


def main():
    print("=" * 70)
    print("SNPTracker - Gerador de Dados de Teste")
    print("=" * 70)
    print()

    datasets = generate_test_datasets()
    save_datasets(datasets)

    print()
    print("=" * 70)
    print("Geração concluída!")
    print("=" * 70)
    print(f"\nDados em: {TEST_DATA_DIR}/")
    print("Execute: python main.py")


if __name__ == "__main__":
    main()

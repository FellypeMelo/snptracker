# SNPTracker - Detector de SNPs

## Descrição

O **SNPTracker** é uma ferramenta para detecção de **SNPs (Single Nucleotide Polymorphisms)** - variações em uma única posição nucleotídica entre diferentes amostras de DNA. SNPs são os tipos mais comuns de variação genética e são fundamentais em:

- Genética médica (predisposição a doenças)
- Farmacogenômica (resposta a medicamentos)
- Melhoramento genético
- Estudos populacionais
- Identificação forense

### O que são SNPs?

Um **SNP** é uma mudança em um único nucleotídeo na sequência de DNA:

```
Referência: A C T G C T A G C T A
Amostra:    A C T G G T A G A T A
                 ↑       ↑
               SNP 1   SNP 2
```

**Frequência**: Aproximadamente 1 SNP a cada 300 bases no genoma humano.

## Funcionalidades

- **Interface de Linha de Comando (CLI)**: Processamento via argumentos argparse.
- **Suporte a Arquivos FASTA**: Lê sequências diretamente de arquivos `.fasta` ou `.fa`.
- **Comparação Base a Base**: Análise posição por posição de alta precisão.
- **Classificação de Mutações**: Identifica transitions e transversions.
- **Detecção de Indels**: Identifica inserções e deleções.
- **Anotação Funcional**: Classifica cada SNP como SYNONYMOUS, NON_SYNONYMOUS, NONSENSE ou NON_CODING com base no código genético padrão.
- **Múltiplos Reading Frames**: Suporte a frames +1, +2, +3 (fita direta) e -1, -2, -3 (complemento reverso) via `--frame`.
- **Regiões Codificantes (CDS)**: O usuário pode definir quais regiões da sequência são codificantes. SNPs fora das regiões recebem `NON_CODING` automaticamente.
- **Contexto de Trinucleotídeo**: Reporta as bases vizinhas de cada SNP no formato COSMIC (`BASE_ANTERIOR[REF>ALT]BASE_POSTERIOR`).
- **Relatório Estruturado**: Saída formatada em tabela com colunas de anotação e contexto, incluindo o reading frame ativo.
- **Exportação**: Salva resultados em arquivo texto customizável.

## Estrutura de Dados

### 📁 `test_data/` - Dados Sintéticos (Commitados)
Contém **55+ pares referência/amostra** com SNPs conhecidos:
- ✅ **Commitados no GitHub**
- 🎯 **SNPs injetados** (quantidade e posição conhecidas)
- 📊 **Tipos variados** (transitions, transversions, indels)
- 🧪 **Casos clínicos** (CpG islands, heterozigoto, sinônimo)

**Regenerar:**
```bash
python generate_test_data.py
```

### 📁 `data/` - Dados de Pesquisa

Diretório para dados locais de trabalho. Contém `sequences.txt` como **arquivo de demonstração** do formato multi-amostra.

| Tipo de arquivo | Git | Motivo |
|---|---|---|
| `sequences.txt` | ✅ Rastreado | Demonstração do formato multi-amostra |
| `*.fasta`, `*.fa`, `*.fna` | 🚫 Ignorado | Dados reais de sequenciamento |
| `*.fastq`, `*.fq`, `*.gz` | 🚫 Ignorado | Dados brutos de sequenciador |

> Para usar seus próprios dados reais, coloque arquivos `.fasta` neste diretório — eles **não serão commitados**.

**Fontes recomendadas:**
- **dbSNP (NCBI)** - SNPs validados
- **1000 Genomes** - Variantes populacionais
- **Sequenciamento próprio** - Dados de sua pesquisa

## Instalação

### Pré-requisitos

- Python 3.7 ou superior
- Nenhuma dependência externa!

### Instalação

```bash
git clone https://github.com/FellypeMelo/snptracker.git
cd snptracker
```

Pronto! Não precisa instalar nada mais.

## Como Usar

### Execução via CLI

O SNPTracker suporta dois modos de operação: **par único** e **multi-amostra**.

---

#### Modo 1 — Par único (1 referência × 1 amostra)

Use `--reference` e `--sample`. Aceita sequência bruta ou arquivo FASTA.

```bash
# Com strings brutas
python main.py --reference "ACTGCTAGCTA" --sample "ACTGCTGGCTA"

# Com arquivos FASTA
python main.py --reference reference.fasta --sample sample.fasta

# Customizando o arquivo de saída
python main.py --reference ref.fa --sample smp.fa --output meu_relatorio.txt

# Definindo regiões codificantes (CDS) — SNPs fora recebem NON_CODING
python main.py --reference ref.fa --sample smp.fa --cds "1-90"
python main.py --reference ref.fa --sample smp.fa --cds "1-90,100-150"

# Especificando o reading frame (padrão: 1)
python main.py --reference ref.fa --sample smp.fa --frame 2
python main.py --reference ref.fa --sample smp.fa --frame -1
python main.py --reference ref.fa --sample smp.fa --frame -1 --cds "1-90"
```

> ⚠️ `--reference` e `--input` são **mutuamente exclusivos** — use um ou outro, nunca os dois.

---

#### Modo 2 — Multi-amostra (1 referência × N amostras)

Use `--input` apontando para um único arquivo FASTA com múltiplas sequências.

**Regra do arquivo:** a **primeira** sequência é sempre a referência; as demais são amostras.

```
# Formato do arquivo (ex: data/sequences.txt)
>reference
ACTGCTAGCTAGCTAGCTA
>sample1
ACTGCTGGCTAGATAGCTA
>sample2
ACTACTAGCTAGCTAGCTA
```

```bash
python main.py --input data/sequences.txt
```

Um relatório `.txt` separado é gerado para cada amostra que tiver SNPs.
Amostras idênticas à referência aparecem no terminal mas **não geram arquivo**.

---

#### Saída — Modo par único

```
============================================================
SNPTracker - Relatório de Mutações
============================================================

Referência: ACTGCTAGCTAGCTA
Amostra:    ACTGCTGGCTAGATA
Frame:      +1

Total de variações encontradas: 3

----------------------------------------------------------------------
Posição    Ref   Alt   Tipo            Anotação             Contexto
----------------------------------------------------------------------
7          A     G     TRANSVERSION    NON_SYNONYMOUS       T[A>G]C
11         C     A     TRANSVERSION    NON_SYNONYMOUS       G[C>A]T
14         T     A     TRANSVERSION    NON_CODING           C[T>A]_

Relatório salvo em: snps_report.txt

Análise concluída!
```

#### Saída — Modo multi-amostra

```
============================================================
SNPTracker - Análise Multi-Amostra
============================================================
Referência: reference (19 bp)
Amostras:   2

[1/2] sample1 → 2 SNPs → salvo em snps_report_sample1.txt
[2/2] sample2 → 1 SNPs → salvo em snps_report_sample2.txt

Análise concluída!
```

## Tipos de Mutações

### 1. Transitions (Mudança Purina→Purina ou Pirimidina→Pirimidina)

| De | Para | Tipo |
|----|------|------|
| A (Adenina) | G (Guanina) | Transition |
| G (Guanina) | A (Adenina) | Transition |
| C (Citosina) | T (Timina) | Transition |
| T (Timina) | C (Citosina) | Transition |

**Frequência**: Mais comuns (2/3 das mutações)

### 2. Transversions (Mudança Purina↔Pirimidina)

| De | Para | Tipo |
|----|------|------|
| A/G | C/T | Transversion |
| C/T | A/G | Transversion |

**Frequência**: Menos comuns (1/3 das mutações)

### 3. Indels (Inserções e Deleções)

```
Referência: ACTGCTAGCTAGCTA (15 bp)
Amostra:    ACTGCTAGCTA---- (11 bp)  # Deleção de 4 bases

Referência: ACTGCTAGCTA (11 bp)
Amostra:    ACTGCTAGCTAGCTA (15 bp)  # Inserção de 4 bases
```

## Estrutura do Projeto

```
snptracker/
├── main.py              # Código principal e CLI
├── annotation.py        # Anotação funcional de SNPs (código genético)
├── fasta_parser.py      # Utilitário de leitura FASTA
├── requirements.txt     # Sem dependências
├── README.md           # Documentação
├── tests/              # Suíte de testes unitários e integração
├── data/               # Diretório para dados reais
└── test_data/          # Dados sintéticos para validação
```

## Guia de Desenvolvimento

### Milestones do Projeto

#### Milestone 1: Detecção Básica ✅
- [x] Comparar sequências base a base
- [x] Classificar transitions e transversions
- [x] Detectar indels simples
- [x] Gerar relatório formatado
- [x] Documentação inicial

#### Milestone 2: Melhorias de Funcionalidade ✅
- [x] Ler sequências de arquivos FASTA
- [x] Implementar Argparse CLI
- [x] Criar suíte de testes (TDD)
- [x] Suporte a múltiplas amostras vs referência
- [x] Anotação de SNPs (SYNONYMOUS / NON_SYNONYMOUS / NONSENSE / NON_CODING)

#### Milestone 3: Análises Avançadas 📊
- [x] Contexto de trinucleotídeo (formato COSMIC: `X[R>A]Y`)
- [x] Regiões codificantes vs não-codificantes (`--cds "start-end,..."`)
- [x] Múltiplos reading frames (`--frame 1/2/3/-1/-2/-3`)
- [ ] Análise de qualidade (Phred scores / FASTQ)
- [ ] Predição funcional (SIFT, PolyPhen)

#### Milestone 4: Integração e Bancos de Dados 🔄
- [ ] Consulta a dbSNP (NCBI)
- [ ] Anotação com genes e transcripts
- [ ] Exportação em formato VCF
- [ ] Comparação com populações (1000 Genomes)

## Algoritmo

### Detecção Base a Base

```python
for i in range(min_length):
    if reference[i] != sample[i]:
        snp = {
            'position': i + 1,
            'reference': reference[i],
            'alternate': sample[i],
            'type': classify_mutation(reference[i], sample[i]),
            'annotation': annotate_snp(i + 1, reference, sample),
            'context': get_trinucleotide_context(reference, i + 1, reference[i], sample[i]),
        }
```

INDELs detectados pela diferença de comprimento entre as sequências **não recebem** a chave `context`.

## Conceitos Relacionados

### SNP vs Mutação
- **Mutação**: Qualquer alteração no DNA (termo geral)
- **SNP**: Mutação comum na população (>1% de frequência)

### Efeitos Funcionais

A anotação funcional é calculada automaticamente por `annotation.py` com base
no código genético padrão (64 códons). O reading frame padrão é +1, mas pode
ser configurado via `--frame` (valores aceitos: 1, 2, 3, -1, -2, -3).
Posições anteriores ao início do frame (ex.: posição 1 no frame +2) recebem
`NON_CODING` automaticamente.

| Anotação | Descrição |
|----------|-----------|
| `SYNONYMOUS` | Troca de base que **não altera** o aminoácido codificado |
| `NON_SYNONYMOUS` | Troca de base que **altera** o aminoácido codificado |
| `NONSENSE` | Troca que cria um **códon de parada** (TAA, TAG, TGA) |
| `NON_CODING` | SNP em região de trinca incompleta, ou indel (sem anotação de códon aplicável) |

### Contexto de Trinucleotídeo (formato COSMIC)

Para cada SNP (não INDEL), o relatório inclui as bases imediatamente vizinhas na sequência de referência:

```
BASE_ANTERIOR[REF>ALT]BASE_POSTERIOR
```

Exemplos:

| Posição | Contexto | Significado |
|---------|----------|-------------|
| Meio da sequência | `T[A>G]C` | Base anterior T, posterior C |
| Início (pos 1) | `_[A>G]C` | Sem base anterior |
| Final | `T[G>T]_` | Sem base posterior |

O padrão `_` é usado como sentinela nas bordas da sequência.

Exemplo:

```
GCT → GCC  →  Ala → Ala  →  SYNONYMOUS
GCT → GTT  →  Ala → Val  →  NON_SYNONYMOUS
GCT → TAA  →  Ala → STOP →  NONSENSE
```

## Limitações Atuais

- Sem informação de qualidade de leitura (suporte FASTQ planejado)
- Não identifica SNPs em repetições
- Sem exportação VCF
- `--cds` e `--frame` disponíveis apenas no modo par único (`--reference`/`--sample`); suporte multi-amostra planejado

## Próximos Passos Recomendados

1. **FASTQ / Phred**: Suporte a qualidade de leitura
2. **VCF Export**: Formato padrão da indústria
3. **SIFT / PolyPhen**: Predição funcional computacional

## Licença

MIT License - veja arquivo LICENSE

## Contato

Abra uma issue para dúvidas ou sugestões.

---

**Status**: 🟢 Funcional - Pronto para uso e expansão
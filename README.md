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
- **Relatório Estruturado**: Saída formatada em tabela.
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

### 📁 `data/` - Dados Reais (Gitignored)
Para dados reais de pesquisa:
- 🚫 **Ignorado pelo Git**
- 🧬 **Dados de sequenciamento** (WGS, WES, painéis)
- 🔬 **Amostras clínicas**

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

### Execução via CLI (Recomendado)

O SNPTracker agora suporta argumentos de linha de comando para maior flexibilidade.

#### Usando strings brutas:
```bash
python main.py --reference "ACTGCTAGCTA" --sample "ACTGCTGGCTA"
```

#### Usando arquivos FASTA:
```bash
python main.py --reference reference.fasta --sample sample.fasta
```

#### Customizando a saída:
```bash
python main.py --reference ref.fa --sample smp.fa --output meu_relatorio.txt
```

### Exemplo de Saída

```
============================================================
SNPTracker - Relatório de Mutações
============================================================

Referência: ACTGCTAGCTAGCTA
Amostra:    ACTGCTGGCTAGATA

Total de variações encontradas: 3

------------------------------------------------------------
Posição    Ref   Alt   Tipo
------------------------------------------------------------
7          A     G     TRANSVERSION
11         C     A     TRANSVERSION
14         T     A     TRANSVERSION

Relatório salvo em: snps_report.txt

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
- [ ] Suporte a múltiplas amostras vs referência
- [ ] Anotação de SNPs (sinônimo/não-sinônimo)

#### Milestone 3: Análises Avançadas 📊
- [ ] Efeito funcional predito (SIFT, Polyphen)
- [ ] Contexto de sequência (trinucleotídeos)
- [ ] Regiões codificantes vs não-codificantes
- [ ] Análise de qualidade (Phred scores)

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
            'type': classify(reference[i], sample[i])
        }
```

## Conceitos Relacionados

### SNP vs Mutação
- **Mutação**: Qualquer alteração no DNA (termo geral)
- **SNP**: Mutação comum na população (>1% de frequência)

### Efeitos Funcionais
1. **Sinônimo (Silencioso)**: Não altera o aminoácido
2. **Não-sinônimo**: Altera o aminoácido
3. **Nonsense**: Cria códon de parada prematuro

## Limitações Atuais

- Apenas duas sequências por vez
- Sem anotação funcional
- Sem informação de qualidade
- Não identifica SNPs em repetições
- Sem exportação VCF

## Próximos Passos Recomendados

1. **Múltiplas Amostras**: Comparar vários indivíduos
2. **Anotação**: Identificar efeito na proteína
3. **VCF Export**: Formato padrão da indústria
4. **Filtros**: Qualidade, profundidade, etc.

## Licença

MIT License - veja arquivo LICENSE

## Contato

Abra uma issue para dúvidas ou sugestões.

---

**Status**: 🟢 Funcional - Pronto para uso e expansão
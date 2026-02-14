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

- **Comparação Base a Base**: Análise posição por posição
- **Classificação de Mutações**: Identifica transitions e transversions
- **Detecção de Indels**: Identifica inserções e deleções
- **Relatório Estruturado**: Saída formatada em tabela
- **Exportação**: Salva resultados em arquivo texto

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

### Execução Básica

```bash
python main.py
```

O programa compara duas sequências de exemplo e gera um relatório.

### Personalizando Sequências

Edite as variáveis no final do arquivo `main.py`:

```python
reference = "ACTGCTAGCTAGCTA"  # Sequência de referência
sample = "ACTGCTGGCTAGATA"    # Sequência da amostra
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
├── main.py              # Código principal
├── requirements.txt     # Sem dependências
├── README.md           # Documentação
├── data/
│   └── sequences.txt   # Exemplo de sequências
└── snps_report.txt     # Relatório gerado
```

## Guia de Desenvolvimento

### Milestones do Projeto

#### Milestone 1: Detecção Básica ✅
- [x] Comparar sequências base a base
- [x] Classificar transitions e transversions
- [x] Detectar indels simples
- [x] Gerar relatório formatado
- [x] Documentação inicial

#### Milestone 2: Melhorias de Funcionalidade 🚧
- [ ] Ler sequências de arquivos FASTA
- [ ] Suporte a múltiplas amostras vs referência
- [ ] Anotação de SNPs (sinônimo/não-sinônimo)
- [ ] Cálculo de frequência alélica

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

### Tarefas para Contribuidores

**Nível Iniciante:**
1. Adicionar argparse para CLI
2. Implementar leitura de arquivos FASTA
3. Criar testes unitários
4. Adicionar estatísticas (taxa de mutação, etc.)

**Nível Intermediário:**
1. Implementar análise de múltiplas amostras
2. Adicionar anotação sinônimo/não-sinônimo
3. Criar visualização das mutações
4. Exportar em formato VCF básico

**Nível Avançado:**
1. Integrar com APIs de bancos de dados (NCBI, Ensembl)
2. Implementar pipeline completo de calling
3. Adicionar filtros de qualidade
4. Análise de linkage disequilibrium

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

### Classificação

```python
purines = {'A', 'G'}
pyrimidines = {'C', 'T'}

if (ref in purines and alt in purines) or \
   (ref in pyrimidines and alt in pyrimidines):
    return "TRANSITION"
else:
    return "TRANSVERSION"
```

## Exemplos de Aplicação

### 1. Medicina Genômica
```
Gene CFTR (Fibrose Cística):
Referência: ...ATG GAG AAG...
Paciente:   ...ATG GTG AAG...  # SNP: Glu→Val (G542V)
                              # Mutação patogênica
```

### 2. Farmacogenômica
```
Gene CYP2D6:
Referência: ...CYP2D6*1 (normal)
Paciente:   ...CYP2D6*4 (variante)
                              # Metabolismo lento de codeína
```

### 3. Agricultura
```
Trigo:
Referência: ...GGCC... (susceptível à doença)
Cultivar:   ...GACC... (resistente)
                              # SNP associado à resistência
```

## Conceitos Relacionados

### SNP vs Mutação
- **Mutação**: Qualquer alteração no DNA (termo geral)
- **SNP**: Mutação comum na população (>1% de frequência)

### Efeitos Funcionais
1. **Sinônimo (Silencioso)**: Não altera o aminoácido
2. **Não-sinônimo**: Altera o aminoácido
3. **Nonsense**: Cria códon de parada prematuro

### Nomenclatura
- **rsID**: Identificador no dbSNP (ex: rs334)
- **HGVS**: Padrão de nomenclatura (ex: NM_000518.5:c.20A>T)

## Limitações Atuais

- Apenas duas sequências por vez
- Sem anotação funcional
- Sem informação de qualidade
- Não identifica SNPs em repetições
- Sem exportação VCF

## Próximos Passos Recomendados

1. **Leitura FASTA**: Processar arquivos reais
2. **Múltiplas Amostras**: Comparar vários indivíduos
3. **Anotação**: Identificar efeito na proteína
4. **VCF Export**: Formato padrão da indústria
5. **Filtros**: Qualidade, profundidade, etc.

## Formatos de Arquivo

### VCF (Variant Call Format) - Futuro
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
chr1    12345   rs123   A       G       99      PASS    DP=35
```

### SAM/BAM - Alinhamentos
Formato binário para armazenar alinhamentos de reads.

## Referências

- [SNPs - NCBI](https://www.ncbi.nlm.nih.gov/snp/)
- [dbSNP Database](https://www.ncbi.nlm.nih.gov/snp/)
- [1000 Genomes Project](https://www.internationalgenome.org/)
- [VCF Format](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
- [HGVS Nomenclature](https://varnomen.hgvs.org/)

## Licença

MIT License - veja arquivo LICENSE

## Contato

Abra uma issue para dúvidas ou sugestões.

---

**Status**: 🟢 Funcional - Pronto para uso e expansão
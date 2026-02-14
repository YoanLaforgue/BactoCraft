# BactoCraft

> **BactoCraft** est une méthodologie d'analyse bioinformatique conçue pour reconstruire des génomes bactériens circulaires et complets en utilisant uniquement le séquençage Oxford Nanopore Technologies (ONT). 

Développé dans un cadre hospitalier, cet outil vise à fournir un diagnostic de précision, offrant une alternative moderne et rapide au pulsotypage (PFGE - *Pulsed-Field Gel Electrophoresis*) pour l'investigation des clusters épidémiques et la détection des facteurs de virulence/résistance.

---

## Contexte :

Dans un contexte hospitalier caractérisé par un flux patient important, l'émergence de clusters bactériens nosocomiaux est une préoccupation majeure. Les techniques traditionnelles de biologie moléculaire permettent de distinguer des espèces mais manquent souvent de la résolution pour discriminer des souches clonales dans des cas complexes.

**Les objectifs de la reconstruction génomique sont :**
1.  **Épidémiologie génomique :** Reconstruire le génome complet pour confirmer ou infirmer la présence d'un cluster hospitalier avec une résolution fine (niveau SNP/cgMLST).
2.  **Diagnostic de précision :** Détecter l'ensemble des gènes de virulence et de résistance (AMR) pour ajuster les thérapies antimicrobiennes.

> **Note :** Ce pipeline est optimisé pour des isolats bactériens (souches pures cultivées).

---

## Prérequis

*   **Système d'exploitation** : Un système basé sur Linux/Unix avec un shell `Bash`.
*   **Accès GPU requis.**
*   **Outils** 
    *   `Dorado` (v1.2.0+)
    *   `Porechop` (v0.2.3, avec SeqAn 2.1.1)
    *   `NanoPlot` (v1.42.0)
    *   `NanoFilt` (v2.7.1)
    *   `Flye` (v2.9+)
    *   `Autocycler` (v0.5.1+)
    *   `Medaka` (v2.1.1+)
    *   `Python` (v3.9)
    *   Packages : `edlib`, `biopython`
*   **Script .py**
    *   `split_fastq_coverage.py` : 

---

## Tutoriel

Ce tutoriel décrit les différentes étapes de la méthodologie `BactoCraft`. Les chemins (`$path/...`) et variables (`$nb_threads`, `$numBarcode`) doivent être adaptés à votre environnement.

Les données de séquençage utilisées pour ce tutoriel sont disponibles sur Zenodo : [ONT POD5](https://zenodo.org/records/18644616)

### Étape 1 : Basecalling & Demultiplexing

Le basecalling est réalisé avec **Dorado** en mode `duplex`.

Le séquençage `duplex` permet de lire les deux brins d'une même molécule d'ADN. Cela génère des reads avec une qualité >Q30, rivalisant avec le séquençage court.

```bash
dorado duplex --threads "${SLURM_CPUS_PER_TASK}" --device cuda:0 ${MODEL_DIR}/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 -r ${POD5_DIR} > ${FASTQ_DIR}/barcode_all.bam
```

```bash
dorado demux --emit-fastq --threads "${SLURM_CPUS_PER_TASK}" --verbose --kit-name SQK-NBD114-96 "${FASTQ_DIR}/barcode_all.bam" -o "${DEMUXED_DIR}"
```

### Étape 2 : Contrôle Qualité Initial (QC) des *Reads*

Évaluation de la qualité globale du run de séquençage.

```bash
NanoPlot -t "$nb_threads" \
        --fastq "$path/to/your/fastq/$numBarcode.fastq" \
        --title "${dateSeq}_${numBarcode}_QC" \
        --outdir "$path/to/output/NANOPLOT_REPORT" \
        --maxlength 1000000 \
        --plots dot
```

### Étape 3 : Suppression des Adaptateurs ONT

Utilisation de `porechop` pour supprimer les séquences d'adaptateurs résiduelles.

```bash
porechop -i input.fastq -o adapter_trim.fastq
```


### Étape 4 : Les populations de reads

<img width="1304" height="572" alt="Capture d’écran 2026-02-14 224856" src="https://github.com/user-attachments/assets/b45cf517-5cea-4e9d-9ce5-f5d4f22976dc" />

<img width="1318" height="640" alt="Capture d’écran 2026-02-14 2249161" src="https://github.com/user-attachments/assets/5170c1f5-7a90-4223-a3d5-5a9fc34b5fe2" />

Dans un contexte d’activité de routine, la qualité des données produites n’est pas toujours homogène. Elle peut être influencée par divers facteurs (qualité de l’extraction d’ADN, usure des flow cells) qui impactent directement le déroulement du run.

Cette variabilité est particulièrement marquée lors de la comparaison entre les bactéries Gram négatif, qui fournissent généralement des fragments d’ADN longs, et les bactéries Gram positif, dont la paroi cellulaire épaisse complique l’extraction des fragments de grande taille.

Afin de prendre en compte cette hétérogénéité et de garantir des résultats exploitables, y compris à partir de runs non optimaux, `BactoCraft` n’applique pas un filtrage global et uniforme des reads. À la place, les données sont segmentées en trois populations distinctes, chacune remplissant un rôle spécifique dans le processus de reconstruction :

- reads > 10k pb --> reads correction --> assemblage "initial"
- reads > 4k pb --> assemblage "rescue"
- reads > Q25 --> polishing
ex:
```bash
NanoFilt ${FASTQ} -q 15 --headcrop 10 --tailcrop 10 --length 10000 --maxlength 1000000 > ${TRIM > 10K}
```

### Étape 5 : Correction des erreurs (Herro)

Avant l'assemblage, les reads longs (>10kb) subissent une correction via l'algorithme Herro (intégré dans Dorado). Cela réduit drastiquement les erreurs liées au séquençage.

```bash
dorado correct --model-path ./herro-v1 input.fastq > corrected.fasta
```

### Étape 6 : Assemblage "initial"

Flye est exécuté en mode `--nano-corr`.

```bash
flye --nano-corr "$MERGED_CORRECTED" --out-dir ${ASSEMBLY_DIR}/NANO_CORR --threads ${SLURM_CPUS_PER_TASK} --deterministic --meta --genome-size ${GENOME_SIZE}
```

> **Note :** L'option `--deterministic` est cruciale en contexte clinique pour garantir que deux analyses du même jeu de données produisent exactement le même résultat (reproductibilité).

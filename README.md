# BactoCraft

<div align="center">
  <a href="[https://github.com/YoanLaforgue/BactoCraft/blob/main/imgs/logo_bactocraft.png]">
    <img src="/imgs/logo_bactocraft.png" alt="BactoCraft_logo" width="340"/>
  </a>
  <br><br>
</div>

> **BactoCraft** est une méthodologie d'analyse bioinformatique conçue pour reconstruire des génomes bactériens circulaires et complets en utilisant uniquement le séquençage Oxford Nanopore Technologies (ONT). 

Développé dans un cadre hospitalier, cet outil vise à fournir un diagnostic de précision, offrant une alternative au pulsotypage (PFGE - *Pulsed-Field Gel Electrophoresis*) pour l'investigation des clusters épidémiques.

---

## Contexte

Dans un contexte hospitalier caractérisé par un flux de patient important, l'émergence de clusters bactériens nosocomiaux est une préoccupation majeure. Les techniques traditionnelles de biologie moléculaire permettent de distinguer des espèces mais manquent souvent de la résolution pour discriminer des souches clonales dans des cas complexes.

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
    *   `split_fastq_coverage.py` : Réalise un sous-échantillonnage des fichiers FASTQ en fonction de la taille du génome et de la "depth coverage" souhaitée.

---

## Tutoriel

Ce tutoriel décrit les différentes étapes de la méthodologie `BactoCraft`. Les chemins (`$path/...`) et variables (`$nb_threads`, `$numBarcode`) doivent être adaptés à votre environnement.

Les données de séquençage utilisées pour ce tutoriel sont disponibles sur Zenodo : [ONT POD5](https://zenodo.org/records/18644616)

### Étape 1 : Basecalling & Demultiplexing

Le basecalling est réalisé avec **Dorado** en mode `duplex`.

Le séquençage `duplex` permet de lire les deux brins d'une molécule d'ADN. Cela génère des reads avec une qualité >Q30, rivalisant avec le séquençage court.

```bash
dorado duplex --threads "$nb_threads" --device cuda:0 ${MODEL_DIR}/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 -r ${POD5_DIR} > ${FASTQ_DIR}/barcode_all.bam
```

```bash
dorado demux --emit-fastq --threads "$nb_threads" --verbose --kit-name SQK-NBD114-96 "${FASTQ_DIR}/barcode_all.bam" -o "${DEMUXED_DIR}"
```

### Étape 2 : Contrôle qualité des Reads

Évaluation de la qualité globale du run de séquençage.

```bash
NanoPlot -t "$nb_threads" \
        --fastq "$path/to/your/fastq/$numBarcode.fastq" \
        --title "${dateSeq}_${numBarcode}_QC" \
        --outdir "$path/to/output/nanoplot_report" \
        --maxlength 1000000 \
        --plots dot
```

### Étape 3 : Suppression des adaptateurs

Utilisation de `porechop` pour supprimer les séquences d'adaptateurs résiduelles.

```bash
porechop -i "$path/to/your/fastq/$numBarcode.fastq" -o "$path/to/output/$numBarcode_adapter_trim.fastq"
```


### Étape 4 : Les populations de reads

<img width="1304" height="572" alt="Capture d’écran 2026-02-14 224856" src="https://github.com/user-attachments/assets/b45cf517-5cea-4e9d-9ce5-f5d4f22976dc" />

<img width="1318" height="640" alt="Capture d’écran 2026-02-14 2249161" src="https://github.com/user-attachments/assets/5170c1f5-7a90-4223-a3d5-5a9fc34b5fe2" />

Dans un contexte d’activité de routine, la qualité des données produites n’est pas toujours homogène. Elle peut être influencée par divers facteurs (qualité de l’extraction, usure des flow cells...) qui impactent directement le déroulement du run.

Afin de prendre en compte cette hétérogénéité et de garantir des résultats exploitables, y compris à partir de runs non optimaux, `BactoCraft` n’applique pas un filtrage uniforme des reads. À la place, les données sont segmentées en trois populations distinctes, chacune remplissant un rôle spécifique dans le processus de reconstruction :

- reads > **10k** pb --> reads correction --> assemblage "initial"
- reads > **4k** pb --> assemblage "rescue"
- reads > **Q25** --> polishing
ex:
```bash
NanoFilt "$path/to/your/fastq/$numBarcode_adapter_trim.fastq" -q 15 --headcrop 10 --tailcrop 10 --length 10000 --maxlength 1000000 > "$path/to/output/$numBarcode_10k_reads.fastq"
```

### Étape 5 : Correction des erreurs

Avant l'assemblage, les reads longs (**>10kb**) subissent une correction via l'algorithme `Herro` (intégré dans `Dorado`). Cela réduit drastiquement les erreurs liées au séquençage.

```bash
dorado correct --model-path "$path/to/your/model/herro-v1" "$path/to/your/fastq/$numBarcode_10k_reads.fastq" > "$path/to/output/corrected_reads.fasta"
```

### Étape 6 : Assemblage "initial"

Flye est exécuté en mode `--nano-corr`.

```bash
flye --nano-corr "$path/to/your/corrected_reads.fasta" --out-dir "$path/to/output/assembly_dir/nano-corr" --threads "$nb_threads" --deterministic --meta --genome-size ${GENOME_SIZE}
```

> **Note :** L'option `--deterministic` est cruciale en contexte clinique pour garantir que deux analyses du même jeu de données produisent exactement le même résultat (reproductibilité).

### Étape 7 : Assemblage "rescue"

Si l'assemblage "initial" est fragmenté ou non circulaire, le pipeline bascule automatiquement vers l'assemblage "rescue".

La gestion de la profondeur de séquençage est un paramètre critique en assemblage *de novo*. Il existe une fenêtre optimale de quantité de données :

> * **Insuffisante :** Le graphe d'assemblage présente des ruptures (gaps), empêchant la circularisation.

> * **Excessive :** Une surcharge de données (>500x) peut dégrader la qualité de l'assemblage. Elle introduit un bruit de fond stochastique et des artefacts de séquençage qui complexifient inutilement le graphe, conduisant à des erreurs d'assemblages.

L'objectif de cette étape est de rationaliser l'apport en données. Le pipeline effectue un **sous-échantillonnage** des reads >4kb.

Le script calcule et extrait précisément le nombre de reads nécessaire pour atteindre une couverture théorique cible de **200X**. 

1. **Sous-échantillonnage :**

```bash
python3 split_fastq_coverage.py --input_fastq "$path/to/your/fastq/$numBarcode_4k_reads.fastq" --output_folder "$path/to/output/subset_dir" --genome_size "$GENOME_SIZE_BP" --nanostats "$path/to/your/reads_4k_100k_nanoplot_report_folder/NanoStats.txt" --X 200
```   

2. **Assemblage :**

`Flye` est exécuté en mode `--nano-hq`.

```bash
for READ_SET in "$subset_dir"/*.fastq*; do
            [ -f "$READ_SET" ] || continue
            BASENAME=$(basename "$READ_SET" | cut -d. -f1)
            
            DIR_FLYE_HQ="$path/to/output/assembly_dir/nano-hq/FLYE_TEMP_${BASENAME}"

            flye --nano-hq "$READ_SET" --out-dir "$DIR_FLYE_HQ" --threads "$nb_threads" --deterministic --meta --genome-size "${GENOME_SIZE}"
            
            if ls "${DIR_FLYE_HQ}/"assembly.fasta >/dev/null 2>&1; then
                cp "${DIR_FLYE_HQ}/"assembly.fasta "$path/to/output/assembly_dir/nano-hq/assembly_flye_${BASENAME}.fasta"
            fi

            rm -rf "$DIR_FLYE_HQ" "$READ_SET"
done
```
3. **Consensus calling :**

Afin de regrouper les différents assemblages et de générer une séquence consensus, nous utilisons l’outil `Autocycler`.

```bash
autocycler compress -i "${ASSEMBLY_DIR}/NANO_HQ" -a "$AUTOCYCLER_OUT"

autocycler cluster -a "$AUTOCYCLER_OUT"
for c in "${AUTOCYCLER_OUT}/clustering/qc_pass/cluster_"*; do
            autocycler trim -c "$c" 
            autocycler resolve -c "$c"
done

autocycler combine -a "$AUTOCYCLER_OUT" -i "${AUTOCYCLER_OUT}/clustering/qc_pass/cluster_"*"/5_final.gfa"
```

### Étape 8 : Polissage

L'assemblage final est poli avec `Medaka` avec les reads > Q25.

Le `--model` r1041_e82_400bps_bacterial_methylation prend en compte la méthylation bactérienne.

```bash
medaka_consensus -i "$path/to/your/fastq/$numBarcode_1k_reads.fastq" -d "$CONSENSUS_RAW" -o "${CONSENSUS_DIR}" -m r1041_e82_400bps_bacterial_methylation -t "$nb_threads"  --bacteria    
```

<div align="center">
  <a href="[https://github.com/YoanLaforgue/BactoCraft/blob/main/imgs/result_atcc_700603.png]">
    <img src="/imgs/result_atcc_700603.png" alt="ATCC" width="850"/>
  </a>
  <br><br>
</div>

---

## Perspectives

Les avancées récentes dans le domaine de la reconstruction génomique ont permis une amélioration significative de la qualité des génomes bactériens obtenus. 

À titre personnel, il semble que nous ayons atteint un plafond technologique. Le principal défi ne réside désormais plus dans l’analyse elle-même, mais plutôt au niveau organisationnel et clinique : faire accepter le séquençage du génome bactérien complet comme un outil de routine en pratique clinique. Un grand nombre d’espèces bactériennes n’ont pas encore été séquencées et assemblées selon des méthodologies rigoureuses. L’étude approfondie de ces organismes représente une opportunité majeure pour améliorer notre compréhension.

Enfin, des progrès restent nécessaires concernant l’interprétation des données. Bien que des outils performants, tels que chewBBACA, permettent l’analyse du cgMLST, le nombre de schémas actuellement disponibles demeure insuffisant. L’élargissement de ces référentiels constitue donc un enjeu important pour une exploitation des données de séquençage.

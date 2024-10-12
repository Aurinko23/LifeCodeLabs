# LifeCodeLabs

## Project Description

LifeCodeLabs is a comprehensive software package for biological data analysis, specializing in working with nucleic acids: DNA and RNA. The project is intended for scientists and researchers who require a reliable tool for genetic data analysis.

The main goal of the project is to simplify complex sequence processing by providing users with a flexible set of tools. These tools allow for quality data filtering and perform basic operations with nucleic acids, such as transcription and translation of DNA, obtaining complementary and reverse complementary sequences for further research and experiments.

LifeCodeLabs is implemented in the Python programming language. The project is geared toward use in educational and scientific research where reliability and reproducibility of results are required.

## Project Structure

- `LifeCodeLabs.py` — the main module coordinating the operation of the tools.
- `bio_files_processor.py` — a tool for preliminary data viewing and analysis.
- `Moduls/`
  - `moduls_dna_rna_tools.py` — tools for working with DNA and RNA.
  - `moduls_filter_fastq.py` — tools for filtering FASTQ data.

## Functionality

### FASTQ Filtering

The function filters sequences based on quality, GC content, and length.
The program accepts 5 arguments: `input_fastq`, `output_fastq`, `gc_bounds`, `length_bounds`, `quality_threshold`:
- `input_fastq` - a file containing FASTQ sequences with the following structure: sequence identifier, sequence, sequence identifier, quality.
- `gc_bounds` - GC content range (in percentages) for filtering (default is (0, 100)). You can pass a single number as an argument; in that case, the number will be used as the upper bound, and the program will filter out reads with quality less than the given number.
- `length_bounds` - length range for filtering (default is (0, 2**32)). You can pass a single number as an argument; in that case, the number will be used as the upper bound, and the program will filter out reads with length less than the given number.
- `quality_threshold` - threshold value of average read quality for filtering (default is 0, phred33 scale). Reads with an average quality across all nucleotides below the threshold are discarded.

As a result, the program returns a similar file consisting only of those sequences that meet all the conditions.

### Operations with DNA and RNA

- `transcribe` — returns the transcribed sequence.
- `reverse` — returns the reversed sequence.
- `complement` — returns the complementary sequence.
- `reverse_complement` — returns the reverse complementary sequence.

The program requires an arbitrary number of arguments with DNA or RNA sequences and the name of the procedure to perform.
- If one sequence is provided, it returns a string with the result. If multiple sequences are provided, it returns a list of strings.
- The program checks that **only** nucleic acid sequences (RNA or DNA) are provided. Otherwise, it returns an input error.

### Reading Bioinformatics Files

- `convert_multiline_fasta_to_oneline` - accepts 2 arguments (`input_fasta` and `output_fasta`). Reads the input FASTA file containing a sequence (DNA/RNA/protein), converts the given sequence into a single line, and writes it to `output_fasta`.
- `output_fasta` - accepts 2 arguments (`input_file`, `output_file`), reads the specified txt file, finds the best matches with the reference database, and saves the protein name in `output_file`.

### Usage Examples

### FASTQ Filtering

```python
seqs = {
    'seq1': ('ATCG', '!!!!'),
    'seq2': ('GCTA', '++++')
}
filtered_seqs = filter_fastq(seqs, gc_bounds=(40, 60), quality_threshold=20)
```
### Operations with DNA and RNA
run_dna_rna_tools('ATG', 'transcribe')    # 'AUG'
run_dna_rna_tools('ATG', 'reverse')       # 'GTA'
run_dna_rna_tools('AtG', 'complement')    # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']

### Authors
Evgeniya Tsymbalova, Saint Petersburg, 2024.

### Software Requirements
The project requires Python 3.6 or higher.

<!--
# LifeCodeLabs

## Описание проекта

LifeCodeLabs — это комплексный программный пакет для анализа биологических данных, специализирующийся на работе с нуклеиновыми кислотами: ДНК и РНК. Проект предназначен для ученых и исследователей, которым требуется надежный инструмент для анализа генетических данных.

Основная цель проекта — упростить сложные процессы обработки последовательностей, предоставляя пользователю гибкий набор инструментов. Эти инструменты позволяют проводить качественную фильтрацию данных, осуществлять основные операций с нуклеиновыми кислотами, такие как транскрипция и трансляция ДНК, получение комплеменатарной и обратно комплементарной последовательностей для дальнейших исследований и экспериментов.

LifeCodeLabs реализован на языке программирования Python. Проект ориентирован на использование в учебных и научных исследованиях, где требуется надежность и воспроизводимость результатов.

## Структура проекта

- `LifeCodeLabs.py` — основной модуль, координирующий работу инструментов.
- `bio_files_processor.py` — инструмент, для предварительного просмотра и анализа данных 
- `Moduls/`
  - `moduls_dna_rna_tools.py` — инструменты для работы с ДНК и РНК.
  - `moduls_filter_fastq.py` — инструменты для фильтрации данных FASTQ.
  

## Функционал

### Фильтрация FASTQ

Функция фильтрует последовательности на основе качества, содержания GC и длины.
Программа принимает на вход  5 аргументов: `input_fastq`, `output_fastq`, `gc_bounds`, `length_bounds`, `quality_threshold`:
    - `input_fastq` - файл, содержащий fastq-сиквенсы со следующей струкурой: идентификатор последовательности, последовательность, идентификатор последовательности, качество. 
    - `gc_bounds` - интервал GC состава (в процентах) для фильтрации (по-умолчанию равен (0, 100)). В аргумент можно передать одно число, тогда по умолчанию данное число будет использоваться в качестве верхней границы и программа отфильтрует риды с качеством меньше введенного числа.
    - `length_bounds` - интервал длины для фильтрации, по-умолчанию равен (0, 2**32).В аргумент можно передать одно число, тогда по умолчанию данное число будет использоваться в качестве верхней границы и программа отфильтрует риды с длиной меньше введенного числа.
    - `quality_threshold` - пороговое значение среднего качества рида для фильтрации, по-умолчанию равно 0 (шкала phred33). Риды со средним качеством по всем нуклеотидам ниже порогового отбрасываются. </br></br>
   
  По итогам работы программа возвращает аналогичный файл, состоящий только из тех сиквенсов, которые удовлетвопмли всем условиям. 


### Операции с ДНК и РНК

- `transcribe` — возвращает транскрибированную последовательность
- `reverse` — возвращает развёрнутую последовательность
- `complement` — возвращает комплементарную последовательность
- `reverse_complement` — возвращает обратную комплементарную последовательность

 Программа требует произвольное количество аргументов с последовательностями ДНК или РНК, и название процедуры которую необходимо выполнить. 
 - Если подана одна последовательность - возвращается строка с результатом. Если подано несколько - возвращается список из строк. 
- Программа проверяет, что поданы **только** последовательности нуклеиновых кислот (РНК или ДНК). В противном случае возвращает ошибку введения.

### Чтение биоинформатических файлов

-`convert_multiline_fasta_to_oneline` - принимает на вход 2 аргумента (input_fasta и output_fasta). Читает поданный на вход fasta-файл, содержащий последовательность (ДНК/РНК/белка) преобразуе данную последоваельность в единую строку и записывае в output_fasta. 
-`output_fasta` - принимает на вход 2 аргумента (input_file, output_file), читает заданный txt файл, находит наилучшее совпадения с референсной базой и сохранить название белка в output_file.
  
## Примеры использования

### Фильтрация FASTQ

```python
seqs = {
    'seq1': ('ATCG', '!!!!'),
    'seq2': ('GCTA', '++++')
}
filtered_seqs = filter_fastq(seqs, gc_bounds=(40, 60), quality_threshold=20)
```

### Операции с ДНК и РНК


```python
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
```

## Авторы

Цымбалова Евгения, Санкт-Петербург, 2024г.

## Требования к программному обеспечению

Проект требует установки Python 3.6 или выше.
-->

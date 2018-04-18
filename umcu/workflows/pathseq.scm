(define-module (umcu workflows pathseq)
  #:use-module (guix workflows)
  #:use-module (guix processes)
  #:use-module (guix gexp)
  #:use-module (gnu packages bioinformatics)
  #:use-module (gnu packages java)
  #:use-module (umcu packages picard)
  #:use-module (umcu packages pathseq))

;;
;; This file contains a workflow that can be used with the Guix Workflow
;; Language.  See: https://www.guixwl.org.
;;

;;
;; DATA LOCATIONS
;; ----------------------------------------------------------------------------
;;
;; Adjust the following variables to your paths.
;;

(define %project-root "/hpc/cog_bioinf/cuppen/project_data/Roel_Jens_SU2C")
(define %bacteria-reference
  (string-append %project-root "/reference/Bacteria_Sequences_removedPlasmid"
                 "_Contigs_Scaffolds_Oct232015_Fusonec_Fusoval.fasta"))

;;
;; For each sample, the relevant processes will be created.  It expects
;; to find the files for the samples in:
;; %project-root/bams/<sample>_dedup.realigned.bam
;;
(define %samples '("FR11123472" "FR11123159"))

;;
;; INDEX REFERENCE GENOME
;; ----------------------------------------------------------------------------
;;
;; The reference genome needs to be Burrow-Wheeler Transformed so that
;; "bwa aln" and "bwa sampe" can work (faster).
;;

(define-public index-bacteria-reference
  (process
     (name "index-bacteria-reference")
     (package-inputs (list bwa))
     (run-time (complexity
                 (space (gigabytes 50))
                 (time (hours 36))))
     (procedure
      `(begin
         ;; Create a Burrows-Wheeler index of the reference genome.
         (system (string-append "bwa index -a bwtsw " ,%bacteria-reference))
         ;; Create an additional index for 'bwa sampe'.
         (system (string-append "bwa bwt2sa "
                                ,%bacteria-reference ".bwt "
                                ,%bacteria-reference ".sa"))))))

;;
;; Additionally, to BLAST against this reference genome, we must build
;; a BLAST database.  This process does exactly that.
;;

(define-public blastdb-index-reference
  (process
   (name "blastdb-index-reference")
   (package-inputs (list rmblast))
   (run-time (complexity
              (space (gigabytes 40))
              (time (hours 8))))
   (procedure
    `(system
      ,(string-append "makeblastdb -in " %bacteria-reference " -dbtype nucl")))))

;;
;; GET UNMAPPED READS FROM HUMAN MAPPING
;; ----------------------------------------------------------------------------
;;
;; This process extracts unmapped reads after mapping reads to GRCh37.
;;

(define (extract-unmapped-reads input-file output-file)
  (process
   (name (string-append "extract-unmapped-reads-for-" (basename input-file ".bam")))
   (package-inputs (list samtools))
   (data-inputs input-file)
   (outputs output-file)
   (run-time (complexity
              (space (gigabytes 10))
              (time (hours 16))))
   (procedure
    #~(system (string-append "samtools view -b -f 4 "
                             #$data-inputs " > " #$outputs)))))

(for-each (lambda (sample)
            (define-dynamically
              (symbol-append 'extract-unmapped-for- (string->symbol sample))
              (extract-unmapped-reads
               (string-append %project-root "/bams/" sample "_dedup.realigned.bam")
               (string-append %project-root "/bams/" sample "_unmapped.bam"))))
          %samples)

;;
;; RUN BWA ALN
;; ----------------------------------------------------------------------------
;;
;; This process maps reads using the bwa-aln algorithm.
;;

(define (run-bwa-aln input-file output-file)
  (process
   (name (string-append "run-bwa-aln-for-" (basename input-file)))
   (package-inputs (list bwa))
   (run-time (complexity
              (space (gigabytes 50))
              (time (hours 1)) ; Unmapped reads aren't so many.
              (threads 8)))
   (procedure 
    `(system ,(string-append "bwa aln -t 8 " %bacteria-reference " "
                             input-file " > " output-file)))))

(for-each (lambda (sample)
            (define-dynamically
              (symbol-append 'run-bwa-aln-for- (string->symbol sample))
              (run-bwa-aln (string-append %project-root "/bams/" sample
                                          "_unmapped.fastq1")
                           (string-append %project-root "/bams/" sample
                                          "_unmapped.fastq1.sai")))
            (define-dynamically
              (symbol-append 'run-bwa-aln-for- (string->symbol sample))
              (run-bwa-aln (string-append %project-root "/bams/" sample
                                          "_unmapped.fastq2")
                           (string-append %project-root "/bams/" sample
                                          "_unmapped.fastq2.sai"))))
          %samples)

;;
;; MAPPING USING BWA SAMPE
;; ------------------------------------------------------------------------------
;;

(define (bwa-sampe sample)
  (process
   (name (string-append "bwa-sampe-" (basename sample)))
   (package-inputs (list bwa samtools))
   (run-time (complexity
              (space (gigabytes 50))
              (time (hours 48))))
   (procedure
    #~(system 
       (string-append
        "bwa sampe " #$%bacteria-reference " "
        #$%project-root "/bams/" #$sample "_unmapped.fastq1.sai "
        #$%project-root "/bams/" #$sample "_unmapped.fastq2.sai "
        #$%project-root "/bams/" #$sample "_unmapped.fastq1 "
        #$%project-root "/bams/" #$sample "_unmapped.fastq2 "
        " > " #$%project-root "/bams/" #$sample "_unmapped.aln.sam")))))

(for-each (lambda (sample)
            (define-dynamically
              (symbol-append 'run-sampe-for- (string->symbol sample))
              (bwa-sampe sample)))
          %samples)

;;
;; CONVERT TO FASTA
;; ----------------------------------------------------------------------------
;;

(define (fastq-to-fasta input-file output-file)
  (process
   (name (string-append "fastq-to-fasta-for-"
                        (basename (basename input-file ".fastq") ".fastq1")))
   (package-inputs (list seqtk))
   (run-time (complexity
              (space (gigabytes 2))
              (time (minutes 30))))
   (procedure
    `(system ,(string-append "seqtk seq -a " input-file " > " output-file)))))

(for-each (lambda (sample)
            (define-dynamically
              (symbol-append 'fastq-to-fasta-for- (string->symbol sample))
              (fastq-to-fasta
               (string-append %project-root "/bams/" sample "_unmapped.fastq1")
               (string-append %project-root "/bams/" sample "_unmapped.1.fasta"))))
          %samples)

(define (blastn sample)
  (process
   (name (string-append "blastn-for-" (basename sample)))
   (package-inputs (list rmblast))
   (run-time (complexity
              (space (gigabytes 1))
              (time (hours 24))))
   (procedure
    (let ((output-dir (string-append %project-root "/blasts"))
          (bam-dir    (string-append %project-root "/bams")))
      `(system
        ,(string-append
          "blastn -task megablast"
          ;; Only blast on forward reads.
          " -query " bam-dir "/" sample "_unmapped.1.fasta"
          " -db " %bacteria-reference
          " -outfmt 5"
          " -evalue 0.0000001"
          " -word_size 16"
          " -max_target_seqs 5"
          " -dust no"
          " -num_threads " (number->string (complexity-threads run-time))
          " -out " output-dir "/" sample ".blast.out"))))))

(for-each (lambda (sample)
            (define-dynamically
              (symbol-append 'blastn-for- (string->symbol sample))
              (blastn sample)))
          %samples)

;; (define (pathseq-bwaunmapped-june2016 sample)
;;   (process
;;    (name (string-append "pathseq-bwaunmapped-june2016-for-" sample))
;;    (package-inputs (list icedtea pathseq-pipeline-tools))
;;    (run-time (complexity
;;               (space (gigabytes 27))
;;               (time (hours 4))))
;;    (procedure
;;     #~(system (string-append
;;                "java -Xmx25G -cp "
;;                #$pathseq-pipeline-tools "/share/java/user-classes/PathSeq.jar"
;;                " BWAunmapped_June2016 "
;;                #$%project-root "/bams/" #$sample "_unmapped.fastq1 "
;;                #$%project-root "/bams/" #$sample "_unmapped.aln.sam "
;;                #$%project-root "/bams/" #$sample "_unmapped.psp.1.out")))))

;; (for-each (lambda (sample)
;;             (define-dynamically
;;               (symbol-append 'psp-bwaunmapped- (string->symbol sample))
;;               (pathseq-bwaunmapped-june2016 sample)))
;;           %samples)

(define (pathseq-blastxml sample)
  (process
   (name (string-append "pathseq-blastxml-for-" sample))
   (package-inputs (list icedtea pathseq-pipeline-tools))
   (run-time (complexity
              (space (gigabytes 28))
              (time (minutes 30))))
   (procedure
    #~(system (string-append
               "java -Xmx24G -cp "
               #$pathseq-pipeline-tools "/share/java/user-classes/PathSeq.jar"
               " blastxml "
               #$%project-root "/blasts/" #$sample ".blast.out "
               #$%project-root "/blasts/" #$sample ".blast.hit")))))

(for-each (lambda (sample)
            (define-dynamically
              (symbol-append 'psp-blastxml- (string->symbol sample))
              (pathseq-blastxml sample)))
          %samples)

(define (pathseq-extract-hits sample)
  (process
   (name (string-append "pathseq-extract-hits-for-" sample))
   (package-inputs (list icedtea pathseq-pipeline-tools))
   (run-time (complexity
              (space (gigabytes 28))
              (time (minutes 30))))
   (procedure
    #~(system (string-append
               "java -Xmx24G -cp "
               #$pathseq-pipeline-tools "/share/java/user-classes/PathSeq.jar"
               " extractFullQuert4BHitTable "
               #$%project-root "/bams/" #$sample "_unmapped.fastq1 "
               #$%project-root "/blasts/" #$sample ".blast.hit "
               #$%project-root "/blasts/" #$sample ".blast.hittable")))))

(for-each (lambda (sample)
            (define-dynamically
              (symbol-append 'psp-extract-hits- (string->symbol sample))
              (pathseq-extract-hits sample)))
          %samples)

;;
;; REMOVE ALIGNMENTS FROM THE ALIGNMENT MAP
;; ----------------------------------------------------------------------------
;;
;; This operation shuffles the reads as well, because it sorts the reads based
;; on their name.
;;

(define (convert-to-unaligned-bam input-file output-file)
  (process
   (name (string-append "convert-to-unaligned-bam-for-"
                        (basename input-file ".bam")))
   (package-inputs (list picard-bin-1.141))
   (data-inputs input-file)
   (outputs output-file)
   (run-time (complexity
              (space (gigabytes 10))
              (time (hours 16))))
   (procedure
    #~(system (string-append
               "java -Xmx8G "
               "-Djava.io.tmpdir=" (getcwd) "/tmp "
               "-jar " #$picard-bin-1.141 "/share/java/picard/picard.jar "
               "RevertSam "
               "I=" #$data-inputs " "
               "O=" #$outputs " "
               "SANITIZE=true "
               "MAX_DISCARD_FRACTION=0.005 "
               "ATTRIBUTE_TO_CLEAR=XT "
               "ATTRIBUTE_TO_CLEAR=XN "
               "ATTRIBUTE_TO_CLEAR=AS "
               "ATTRIBUTE_TO_CLEAR=OC "
               "ATTRIBUTE_TO_CLEAR=OP "
               "SORT_ORDER=queryname "
               "RESTORE_ORIGINAL_QUALITIES=true "
               "REMOVE_DUPLICATE_INFORMATION=true "
               "REMOVE_ALIGNMENT_INFORMATION=true")))))

;;
;; SORT BAM FILES
;; ----------------------------------------------------------------------------
;;

(define (sort-bam-file input-file output-file)
  (process
   (name (string-append "sort-bam-" (basename input-file ".bam")))
   (package-inputs (list samtools))
   (run-time (complexity
              (space (gigabytes 30))
              (time (hours 4))))
   (procedure
    #~(system (string-append "samtools sort -n " #$input-file " -o " #$output-file)))
   (synopsis "Sort reads by name in a BAM file.")
   (description "After unmapping reads, the pairs may not be next to each other
in the BAM file.  Before remapping, we need to sort the reads so that the pairs
are lined up properly.")))

;;
;; CONVERT BAM TO FASTQ
;; ----------------------------------------------------------------------------
;;

(define (bam-to-fastq input-file output-file)
  (process
   (name (string-append "bam-to-fastq-" (basename input-file ".bam")))
   (package-inputs (list bedtools))
   (run-time (complexity
              (space (megabytes 100))
              (time (hours 6))))
   (procedure
    #~(system 
       (string-append 
        "bedtools bamtofastq -i " #$input-file " -fq " #$output-file "1 "
        "-fq2 " #$output-file "2")))))

;;
;; MAPPING USING BWA MEM.
;; ------------------------------------------------------------------------------
;;

(define (bwa-mem sample reference output-file)
  (process
   (name (string-append
          "bwa-mem-" (basename sample) "-" (basename reference)))
   (package-inputs (list bwa samtools))
   (run-time (complexity
              (space (gigabytes 20))
              (time (hours 16))
              (threads 16)))
   (procedure
    #~(system 
       (string-append
        "bwa mem -M -t 16 " #$reference " " #$sample "1 " #$sample "2 "
        "| samtools view -bo " #$output-file)))))

;; (define-public bwa-mem-CP022124
;;   (bwa-mem
;;    %unaligned-fastq
;;    (string-append %project-root "/CP022124.1.fasta")
;;    (string-append %project-root "/" %sample-name "_CP022124.1.bam")))

;; (define-public bwa-mem-AE009951
;;   (bwa-mem
;;    %unaligned-fastq
;;    (string-append %project-root "/AE009951.2.fasta")
;;    (string-append %project-root "/" %sample-name "_AE009951.2.bam")))

;;
;; INDEXING USING BWA INDEX.
;; ------------------------------------------------------------------------------
;;

(define (bwa-index filename)
  (process
   (name (string-append "bwa-index-" (basename filename)))
   (package-inputs (list bwa))
   (run-time (complexity
              (space (gigabytes 8))
              (time (hours 4))))
   (procedure
    #~(system (string-append "bwa index " #$filename)))))

;; (define-public bwa-index-AE009951
;;   (bwa-index (string-append %project-root "/AE009951.2.fasta")))

;; (define-public bwa-index-CP022124
;;   (bwa-index (string-append %project-root "/CP022124.1.fasta")))

;;
;; STRIP UNMAPPED DATA
;; ------------------------------------------------------------------------------
;;

(define (strip-unmapped-reads input-file output-file)
  (process
   (name (string-append "strip-unmapped-reads-" (basename input-file)))
   (package-inputs (list samtools))
   (data-inputs input-file)
   (outputs output-file)
   (run-time (complexity
              (space (gigabytes 6))
              (time (hours 6))))
   (procedure
    #~(system (string-append "samtools view -b -F 4 " #$data-inputs  " > " #$outputs)))))

;; (define-public strip-unmapped-reads-CP022124
;;   (strip-unmapped-reads
;;    (string-append %project-root "/" %sample-name "_CP022124.1.bam")
;;    (string-append %project-root "/" %sample-name "_CP022124.1_mapped.bam")))

;; (define-public strip-unmapped-reads-AE009951
;;   (strip-unmapped-reads
;;    (string-append %project-root "/" %sample-name "_AE009951.2.bam")
;;    (string-append %project-root "/" %sample-name "_AE009951.2_mapped.bam")))

;;
;; SORT THE BAM OUTPUT
;; ----------------------------------------------------------------------------
;;

(define (samtools-sort input-file output-file)
  (process
   (name (string-append "samtools-sort-" (basename input-file)))
   (package-inputs (list samtools))
   (run-time (complexity
              (space (gigabytes 2))
              (time (minutes 10))))
   (procedure
    #~(system (string-append "samtools sort " #$input-file
                             " > " #$output-file)))))

;; (define-public samtools-sort-AE009951
;;   (samtools-sort
;;    (string-append %project-root "/" %sample-name "_AE009951.2_mapped.bam")
;;    (string-append %project-root "/" %sample-name "_AE009951.2_mapped_sorted.bam")))

;; (define-public samtools-sort-CP022124
;;   (samtools-sort
;;    (string-append %project-root "/" %sample-name "_CP022124.1_mapped.bam")
;;    (string-append %project-root "/" %sample-name "_CP022124.1_mapped_sorted.bam")))

;;
;; DETERMINE READ DEPTH
;; ----------------------------------------------------------------------------
;;

(define (samtools-depth input-file output-file)
  (process
   (name (string-append "samtools-depth-" (basename input-file)))
   (package-inputs (list samtools))
   (run-time (complexity
              (space (gigabytes 2))
              (time (minutes 10))))
   (procedure
    #~(system (string-append "samtools depth " #$input-file
                             " > " #$output-file)))))


(define-public pathseq-workflow
  (let ((bwa-sampe-procs (map bwa-sampe %samples))
        )
    (workflow
     (name "pathseq-workflow")
     (processes (append (list index-bacteria-reference) bwa-sampe-procs))
     (restrictions
      (map (lambda (proc)
             `(,proc ,index-bacteria-reference))
           bwa-sampe-procs)))))

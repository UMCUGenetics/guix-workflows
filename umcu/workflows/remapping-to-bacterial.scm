;; Copyright © 2017  Roel Janssen <roel@gnu.org>

;; This program is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.

;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see <http://www.gnu.org/licenses/>.

;; ----------------------------------------------------------------------------
;;
;; This workflow can be used to re-map any BAM file to bacterial genomes
;; CP022124.1 and AE009951.2.  This code has been developed for a single
;; experiment, but may contain useful processes for other purposes as well.
;;
;; ----------------------------------------------------------------------------

(define-module (umcu workflows remapping-to-bacterial)
  #:use-module (guix workflows)
  #:use-module (guix processes)
  #:use-module (guix gexp)
  #:use-module (gnu packages bioinformatics)
  #:use-module (umcu packages picard))

;; CONSTANTS / DATA PATHS
;; ----------------------------------------------------------------------------
;; This is the folder that will be used to save the output to.  It is set
;; to the current working directory, because that's where we need to run
;; the "guixr workflow -r <sample-name>-analysis from.
(define %project-root (getcwd))

;; This will be used to make the %sample-name "-Fusobacterium" files.
;; Adjust this to the relevant sample name.
(define %sample-name "")

;; This file is the already-mapped-to-human BAM file.
;; Adjust to the location of the new BAM file.
(define %sample-bam "")

;; These files will be created after the first step in this pipeline.
(define %unaligned-bam (string-append %project-root "/" %sample-name "_unaligned.bam"))

;; This will automatically replace .bam with .fastq, so no need to adjust this.
(define %unaligned-fastq
  (string-append
   (substring %unaligned-bam 0 (string-index-right %unaligned-bam #\.))
   ".fastq"))

;;
;; REMOVE ALIGNMENTS FROM THE ALIGNMENT MAP
;; ----------------------------------------------------------------------------
;;
;; This operation shuffles the reads as well, because it sorts the reads based
;; on their name.
;;

(define-public convert-to-unaligned-bam
  (process
   (name "convert-to-unaligned-bam")
   (package-inputs (list picard-bin-1.141))
   (data-inputs %sample-bam)
   (outputs %unaligned-bam)
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
;; CONVERT UNALIGNED BAM TO FASTQ
;; ----------------------------------------------------------------------------
;;

(define (bam-to-fastq input-file output-file)
  (process
   (name (string-append "bam-to-fastq-" (basename input-file)))
   (package-inputs (list bedtools))
   (run-time (complexity
              (space (megabytes 100))
              (time (hours 6))))
   (procedure
    #~(system 
       (string-append 
        "bedtools bamtofastq -i " #$input-file " -fq " #$output-file "1 "
        "-fq2 " #$output-file "2")))))

(define-public bam-to-fastq-sample
  (bam-to-fastq %unaligned-bam %unaligned-fastq))

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

(define-public bwa-mem-CP022124
  (bwa-mem
   %unaligned-fastq
   (string-append %project-root "/CP022124.1.fasta")
   (string-append %project-root "/" %sample-name "_CP022124.1.bam")))

(define-public bwa-mem-AE009951
  (bwa-mem
   %unaligned-fastq
   (string-append %project-root "/AE009951.2.fasta")
   (string-append %project-root "/" %sample-name "_AE009951.2.bam")))

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

(define-public bwa-index-AE009951
  (bwa-index (string-append %project-root "/AE009951.2.fasta")))

(define-public bwa-index-CP022124
  (bwa-index (string-append %project-root "/CP022124.1.fasta")))

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

(define-public strip-unmapped-reads-CP022124
  (strip-unmapped-reads
   (string-append %project-root "/" %sample-name "_CP022124.1.bam")
   (string-append %project-root "/" %sample-name "_CP022124.1_mapped.bam")))

(define-public strip-unmapped-reads-AE009951
  (strip-unmapped-reads
   (string-append %project-root "/" %sample-name "_AE009951.2.bam")
   (string-append %project-root "/" %sample-name "_AE009951.2_mapped.bam")))

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

(define-public samtools-sort-AE009951
  (samtools-sort
   (string-append %project-root "/" %sample-name "_AE009951.2_mapped.bam")
   (string-append %project-root "/" %sample-name "_AE009951.2_mapped_sorted.bam")))

(define-public samtools-sort-CP022124
  (samtools-sort
   (string-append %project-root "/" %sample-name "_CP022124.1_mapped.bam")
   (string-append %project-root "/" %sample-name "_CP022124.1_mapped_sorted.bam")))

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

(define-public samtools-depth-AE009951
  (samtools-depth
   (string-append %project-root "/" %sample-name "_AE009951.2_mapped_sorted.bam")
   (string-append %project-root "/" %sample-name "_AE009951.2_mapped_sorted.coverage")))

(define-public samtools-depth-CP022124
  (samtools-depth
   (string-append %project-root "/" %sample-name "_CP022124.1_mapped_sorted.bam")
   (string-append %project-root "/" %sample-name "_CP022124.1_mapped_sorted.coverage")))

;;
;; COMBINE THE PROCESSES ABOVE IN A PIPELINE
;; ----------------------------------------------------------------------------

(define-public remap-to-bacterial
  (workflow
   (name (string-append %sample-name "-remap-to-bacterial"))
   (processes
    (list convert-to-unaligned-bam
          bam-to-fastq-sample
          bwa-index-CP022124
          bwa-index-AE009951
          bwa-mem-CP022124
          bwa-mem-AE009951
          strip-unmapped-reads-CP022124
          strip-unmapped-reads-AE009951
          samtools-depth-CP022124
          samtools-depth-AE009951
          samtools-sort-CP022124
          samtools-sort-AE009951))
   (restrictions
    `((,bam-to-fastq-sample ,convert-to-unaligned-bam)
      (,bwa-mem-AE009951 ,bwa-index-AE009951)
      (,bwa-mem-CP022124 ,bwa-index-CP022124)
      (,bwa-mem-CP022124 ,bam-to-fastq-sample)
      (,bwa-mem-AE009951 ,bam-to-fastq-sample)
      (,strip-unmapped-reads-AE009951 ,bwa-mem-AE009951)
      (,strip-unmapped-reads-CP022124 ,bwa-mem-CP022124)
      (,samtools-depth-AE009951 ,samtools-sort-AE009951)
      (,samtools-depth-CP022124 ,samtools-sort-CP022124)
      (,samtools-sort-AE009951 ,strip-unmapped-reads-AE009951)
      (,samtools-sort-CP022124 ,strip-unmapped-reads-CP022124)))))

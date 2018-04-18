(define-module (umcu workflows histone-mark-analysis)
  #:use-module (guix processes)
  #:use-module (guix workflows)
  #:use-module (guix gexp)
  #:use-module (gnu packages bioinformatics)
  #:use-module (gnu packages python)
  #:use-module (gnu packages statistics)
  #:use-module (umcu packages picard)
  #:use-module (umcu packages python)
  #:use-module (umcu packages bioconductor)
  #:use-module (srfi srfi-1)
  #:use-module (srfi srfi-34)
  #:use-module (ice-9 ftw)     ; 'scandir' functie.
  #:use-module (ice-9 rdelim)  ; read file
  #:use-module (ice-9 regex))

;; Author: Enrico Schmitz
;; Project: Histone mark analysis

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;; This pipeline will convert your fastq to peakfiles and then calculate phenotype cell specificity. 
 ;; To start you need 1 Directory containing at least 1 treatment file and 1 control file both in fastq(can be gzipped), bam, bed/tagAlign or sra
 ;; The name needs to be the same except 1 part that is used for all control files e.g. treatment file Human_Lung.fastq.gz and the control file Human_Lung_control.fastq.gz
;; REF needs to be the path to a indexed reference genome file 
(define REF "")
;; Genome is the path to a genome file, only needed if input files are .sra
(define Genome "")
;; ControlString is the part of a the controle file names that is unique for control files.
(define ControlString "")
;; The Directory containing all the fastq files that need to be converted
(define Directory "")
;; Output_Directory where all output can be stored, should differ from Directory
(define Output_Directory "")
;; Directory with beagle files
(define dir1000G "")
;; (tab delimited) file listing individuals present in the Beagle files along with population information (popName) in the third column
(define BeaglePopList "")
;; popName is the origin of the population
(define popName "")
;;  (tab delimited) SNP name, chromosome (chrN format), base pair position. These are the SNPs with location (location needs to be the same as in BeaglePopList)
(define SNPmappings "")
;;  computes LD between the index SNP and SNPs within this window for LD 
(define LDwindow "")
;; cut­off for Ld calculation
(define r2 "")
;; Name of the histonemark
(define histoneName "")
;; How many times does a permutation needs to be done.
(define permutations "")


;; This directory can be downloaded from http://archive.broadinstitute.org/mpg/epigwas/
;; IMPORTANT: tissue indexing must be the same for background and phenotype SNPs.    
(define bgDir "")
;;Nothing else needs to be personalised
;;Make sure enough time and memory is provided for the jobs  (run-time ( complexity (space (gigabytes ?))(time (minutes ?))))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; (run-time ( complexity (space (gigabytes ?))(time    (minutes ?)))) this line is for setting the jobs memory and time 
;; List with marks where -broad were used
(define broad-list '(
    "H3K27me3"
    "H3K36me3"
    "H3K9me3"
    "H3K4me1"
))
;; Adds the broad mark if the given histonename is inside the broad-list
(define broad (if ( null? ( delete "" (filter (lambda (item) (string-ci=? item "--broad")  )(map (lambda (item) (if (string-ci=? item histoneName) "--broad" "" ))broad-list)))) "" "--broad"))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Home-made functions
;; Get extension from a file
(define (file-extension file-name)
	(last (string-split file-name #\.)))
  
;; Remove extensions		 
(define (remove-extension item)
	(basename (basename (basename (basename (basename (basename item ".gz") ".fastq") ".tagAlign") ".sra") ".bam") ".bed"))

;; Create Output filename
(define (Create-Name dir Name n)
		(string-append Output_Directory dir Name  "_output" n ))
		
;; Get frame-size from PPQT
(define (get-frame file)
 (number->string (string->number (first (string-split (third (string-split (read-line (open-input-file file)) #\tab))#\,)))))

;; =======================================================================================================================
;; 	Files and Directories	
;; =======================================================================================================================
;; Make Output_Directory if it doesn't exist
(if (not(file-exists? Output_Directory)) (mkdir Output_Directory))

;; Define all directory names
(define dirBeagleFiles (string-append Output_Directory "dir1000G/" ))
(define TMP (string-append Output_Directory "tmp/" ))
(define PEAK (string-append Output_Directory "peak/" ))
(define LD (string-append Output_Directory "LD/" ))
(define MACS2 (string-append Output_Directory "macs2/" ))

;; Check if all the directories exist
(if (not(file-exists? TMP)) (mkdir TMP))
(if (not(file-exists? PEAK)) (mkdir PEAK))
(if (not(file-exists? dirBeagleFiles)) (mkdir dirBeagleFiles))
(if (not(file-exists? LD)) (mkdir LD))
(if (not(file-exists? MACS2)) (mkdir MACS2))

;; Gets all the file names from the directory in a list
(define file-list (scandir Directory (lambda (file)
                                      ;; No hidden files
                             (not (or (eq? (string-ref file 0) #\.))))))
							 
;; filter out all control file names and make other list with only treatment file names
(define Name-list  (reverse (remove (lambda (file)
                                      ;; No files containing the ControlString
                                      (string-contains file ControlString)  )file-list)))
 ;; Read in the dir1000G directory and get all names in a list					  
(define dir1000GFiles (scandir dir1000G (lambda (file)
                                      ;; No hidden files
                             (not (or (eq? (string-ref file 0) #\.))))))
;; rename all beaglefiles in a directory so computeLD can use them.
 (map (lambda (name) (if (not(file-exists? (string-append dirBeagleFiles (regexp-substitute/global #f (third (string-split (basename name) #\.))   (basename name) 'pre "phase1_release_v2" 'post))))
 (copy-file (string-append dir1000G name) (string-append dirBeagleFiles (regexp-substitute/global #f (third (string-split (basename name) #\.))   (basename name) 'pre "phase1_release_v2" 'post)))))dir1000GFiles)
(define dir1000G dirBeagleFiles)
     
;; Make tab delimited tissue list with "complete name of file" and name of file without extension 
(define tissueList  (string-append Output_Directory "tissueList.txt"))
(let ((output-port (open-file tissueList "w")))
  (map(lambda (file)
    (display (string-append (basename (Create-Name "peak/" (remove-extension file ) ""))"_peaks.xls" "\t\t\t"  (remove-extension file )) output-port)
    (newline output-port))Name-list)
  (close output-port))


;;==================================================================================================================================
;; Process
;;==================================================================================================================================
;;Convert to FASTQ
(define-public  (run-FASTQ File)
  (process
   (name (string-append "FASTQ-Convert-" (remove-extension File )))
   (package-inputs ( list bedtools sra-tools))
   (run-time ( complexity (space (gigabytes 15))(time    (minutes 60))))
    (procedure #~(begin
	(use-modules (srfi srfi-1))
	;; Convert sra file to fastq file
	(if (string-ci=?  #$(file-extension File) "sra" ) 
	(system (string-append #$sra-tools "/bin/fastq-dump.2.8.2 -O " #$Output_Directory " -I "  #$Directory  #$File " --gzip")))
	;; convert bed and tagAlign to fastq
	(if (or (string-ci=?  #$(file-extension File) "tagAlign" ) (string-ci=?  #$(file-extension File) "bed" )) (begin 
	(system (string-append #$bedtools "/bin/bedToBam -i "  #$Directory  #$File " -g "  #$Genome " > " #$Output_Directory #$(basename (basename (basename File ".gz" ) ".tagAlign") ".bed") ".bam"))
	(system (string-append #$bedtools "/bin/bamToFastq -i " #$Output_Directory #$(basename (basename (basename File ".gz" ) ".tagAlign") ".bed") ".bam -fq "  #$Output_Directory  #$(basename (basename (basename File ".gz") ".tagAlign") ".bed")".fastq"))))
	;; convert bam to fastq
	(if (string-ci=?  #$(file-extension File) "bam" ) 
	(system (string-append #$bedtools "/bin/bamToFastq -i " #$Directory  #$File " -fq "  #$Output_Directory  #$(basename (basename File ".gz") ".bam") ".fastq")))))))
	
;; Do BWA
;; Mapping
(define-public  (run-BWA File)
  (process
   (name (string-append "bwa-" (remove-extension File )))
   (package-inputs ( list bwa sambamba))
   (run-time ( complexity (space (gigabytes 70))(time    (minutes 300))))
	;;BWA and sambamba and check where the file is because if it wasn't in fastq format it is inside the output directory
	(procedure #~(system (string-append "bwa mem -t 5 -c 100 -M "  #$REF " "   (if (string-contains #$File ".fastq")  #$(string-append Directory File) #$(string-append Output_Directory (basename (basename (basename (basename (basename File ".gz") ".tagAlign") ".bed") ".bam") ".sra") ".fastq")) " | sambamba view --format=bam -S -o "  #$(Create-Name "" (remove-extension File ) "1") ".bam /dev/stdin")))))


;; Do SortSam
;; Mapping
(define-public (run-SortSam File)
  (process
   (name (string-append "sortsam-" (remove-extension File )))
   (package-inputs ( list picard-bin-1.141))
   (run-time ( complexity (space (gigabytes 50))(time    (minutes 120))))
    (procedure
     #~(system (string-append 
                            "java -Xmx10g -Djava.io.tmpdir=" #$TMP " -jar " #$picard-bin-1.141
                            "/share/java/picard/picard.jar SortSam INPUT=" #$(Create-Name "" (remove-extension File ) "1") ".bam"
							" OUTPUT="  #$(Create-Name "" (remove-extension File ) "2") ".bam"
							" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=" #$TMP "")))))

 
;; do MarkDup
;; Mapping
(define-public (run-MarkDup File)
  (process
   (name (string-append "MarkDup-" (remove-extension File )))
   (package-inputs ( list picard-bin-1.141))
   (run-time ( complexity (space (gigabytes 40))(time    (minutes 120))))
    (procedure
     #~(system (string-append 
                            "java -Xmx10g -Djava.io.tmpdir=" #$TMP " -jar " #$picard-bin-1.141
                            "/share/java/picard/picard.jar MarkDuplicates INPUT="  #$(Create-Name "" (remove-extension File ) "2")  ".bam"
							" OUTPUT="  #$(Create-Name "" (remove-extension File ) "3") ".bam"
							" METRICS_FILE="  #$(Create-Name "" (remove-extension File ) "3") ".dup_metrics ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT TMP_DIR=" #$TMP )))))

 
;; Do SamTools
;; Filtering
(define-public (run-SamTools File)
  (process
   (name (string-append "SamTools-" (remove-extension File )))
   (package-inputs ( list samtools))
   (run-time ( complexity (space (gigabytes 40))(time    (minutes 120))))
    (procedure `(system (string-append 
                                     "samtools view -b -F 4 -q 5 "  ,(Create-Name "" (remove-extension File ) "3")".bam" " | samtools view -b -F 1024 > "  ,(Create-Name "" (remove-extension File ) "4") ".bam")))))

;; Bam to tagAlign
;; Convert to tagAlign so PhantomPeakQualTools can use the files.
(define-public (run-BedTools File)
  (process
   (name (string-append "BedTools-" (remove-extension File )))
   (package-inputs ( list bedtools))
   (run-time ( complexity (space (gigabytes 40))(time    (minutes 180))))
    (procedure `(system (string-append 
                                     "bamToBed -i " ,(Create-Name "" (remove-extension File ) "4") ".bam" " > " ,(Create-Name "" (remove-extension File ) "4") ".tagAlign")))))
;; Do PPQT
;; Get the framesize for macs2
(define-public (run-PhantomPeakQualTools File)
  (process
   (name (string-append "PPQT-" (remove-extension File)))
   (package-inputs ( list r r-phantompeakqualtools))
   (run-time ( complexity (space (gigabytes 40))(time    (minutes 120))))
 (procedure #~(system (string-append "Rscript " #$r-phantompeakqualtools "/share/scripts/run_spp.R -c="  #$(Create-Name "" (remove-extension File ) "4")".tagAlign" " -rf -out=" #$(Create-Name "" (remove-extension File ) "4") "_params.out")))))

;; Do macs2
;; Calculate peaks
(define-public (run-MACS2 Co In) 
  (process
   (name (string-append "MACS2-" (remove-extension In) ))
   (package-inputs ( list python-macs2))
   (run-time ( complexity (space (gigabytes 40))(time    (minutes 120))))
   (procedure `(begin
		(use-modules (ice-9 rdelim)) 
		(use-modules (srfi srfi-1))
		(define (get-frame file)
			(number->string (string->number (first (string-split (third (string-split (read-line (open-input-file file)) #\tab))#\,)))))
		(system (string-append "macs2 callpeak --tempdir " ,TMP " -t "  ,(Create-Name "" (remove-extension In ) "4") ".bam"  " -n "  ,(Create-Name "macs2/" (remove-extension In ) "") " --gsize hs -c " ,(Create-Name "" (remove-extension Co) "4")".bam" " --nomodel --extsize=" (get-frame ,(string-append(Create-Name "" (remove-extension In ) "4") "_params.out")) " " ,broad))))))
    
;; prepare files
;; Changes the format form the macs files so epigwas can use them.
(define-public (run-prepFiles File)
  (process
   (name (string-append "run-prepFiles-"(remove-extension File )))
   (run-time ( complexity (space (gigabytes 10))(time    (minutes 30))))
   (procedure

    `(begin
	(use-modules (ice-9 rdelim))
	(use-modules (ice-9 regex))
	(use-modules(ice-9 format)) 
	(use-modules (srfi srfi-1))
	;; rewrite the macs2 output so makeFiles and makeFilesBG can use them. Change a tab to ;
	(define peakFile (open-file ,(string-append (Create-Name "peak/" (remove-extension File ) "") "_peaks.xls") "w"))
	(define openfile (open-input-file ,(string-append (Create-Name "macs2/" (remove-extension File ) "") "_peaks.xls")))
	(define (my-reader port line-number)
	(let ((line (read-line port)))
		(if (not (eof-object? line))
			(begin
			;; Start after line 26 to remove all the input information
			(if (not(< line-number 27))
				(format peakFile "~a~&" (string-append (first (string-split (regexp-substitute/global #f "[\t]+"  line 'pre ";" 'post) #\;)) ";" (second (string-split (regexp-substitute/global #f "[\t]+"  line 'pre ";" 'post) #\;)) ";" (third (string-split (regexp-substitute/global #f "[\t]+"  line 'pre ";" 'post) #\;)) ";" (number->string(- (string->number(fifth (string-split (regexp-substitute/global #f "[\t]+"  line 'pre ";" 'post) #\;))) (string->number(second (string-split (regexp-substitute/global #f "[\t]+"  line 'pre ";" 'post) #\;))))) ";-;" (eighth (string-split (regexp-substitute/global #f "[\t]+"  line 'pre ";" 'post) #\;)) ";-")))
			(my-reader port (+ 1 line-number))))))
	(my-reader openfile 0)))))
	; )

;; Do computeLD
;; Calculate ld for all snps
(define-public run-computeLD
  (process
   (name "run-computeLD" )
   (package-inputs ( list python-3.5 epigwas))
   (run-time ( complexity (space (gigabytes 45))(time  (minutes 300))))
   (procedure #~(system (string-append "python " #$epigwas "/share/epigwas/scripts/computeLd/computeLd.py "  #$dir1000G " " #$BeaglePopList " " #$popName " " #$SNPmappings " " #$LDwindow " " #$r2 " " #$LD " " #$Output_Directory "computeLD" )))))
	
;; Do makeFiles
;; Make the files needed for phenoCellSpecif
(define-public run-makeFiles
  (process
   (name  "run-makeFiles")
   (package-inputs ( list python-3.5 epigwas))
   (run-time ( complexity (space (gigabytes 40))(time    (minutes 120))))
   (procedure #~(system (string-append "python " #$epigwas "/share/epigwas/scripts/prepFiles_v01/makeFiles.py " #$Output_Directory "computeLD.1kg.map.txt " #$LD " " #$PEAK " " #$tissueList " " #$histoneName " " #$r2 " " #$Output_Directory "treatment")))))

;; Makefiles for the BG the files made can be used in phenoCellSpecif   
(define-public run-makeFilesBG
  (process
   (name  "run-makeFilesBG")
   (package-inputs ( list python-3.5 epigwas))
   (run-time ( complexity (space (gigabytes 40))(time  (minutes 240))))
   (procedure #~(system (string-append "python " #$epigwas "/share/epigwas/scripts/prepFiles_v01/makeFiles.py " #$bgDir "bg.1kg.map.txt " #$bgDir "ld/ " #$PEAK " " #$tissueList " " #$histoneName " " #$r2 " " #$Output_Directory "BG")))))
   

;; Do phenoCellSpecif
;; Calculate specificity for all tissues
(define-public run-phenoCellSpecif
  (process
   (name "run-phenoCellSpecif")
   (package-inputs ( list python-3.5 epigwas))
   (run-time ( complexity (space (gigabytes 30))(time    (minutes 720))))
   (procedure #~(system (string-append "python " #$epigwas "/share/epigwas/scripts/phenoCellSpec_v01/phenoCellSpecif.py " #$Output_Directory "treatment." #$histoneName ".snpPeak.txt " #$Output_Directory "treatment.ld.txt " #$Output_Directory "treatment." #$histoneName ".tisIndex.txt " #$Output_Directory "BG." #$histoneName ".snpPeak.txt " #$Output_Directory "BG.ld.txt "  #$permutations " " #$Output_Directory "phenoCell")))))
    
;; This process makes the paired processes for macs2   
(define* (make-macs2-processes items #:optional (output '()))
  (if (null? items)
      output 
      (if (even? (length items))
          (make-macs2-processes (cdr (cdr items))
                                (cons (run-MACS2 (first items) (second items)) output))
          (make-macs2-processes (cdr (cdr items)) output))))

;; This process pairs the treatment file to the control file.		  
(define* (get-pairs items #:optional (output '()))
  (if (null? items)
      output
      (if (even? (length items))
          (get-pairs (cdr (cdr items))
                     (cons (list (first items) (second items)) output))
          (get-pairs (cdr (cdr items)) output))))
	
;;==================================================================================================================================
;; Workflow
;; First all processes are made
;; And then the restrictions make sure everything goes in order.
;;==================================================================================================================================
   
(define-public pipeline-workflow
(let (
;; Make all the jobs
	   (FASTQ-procs     (map run-FASTQ file-list))
	   (bwa-procs     (map run-BWA file-list))
       (sortsam-procs (map run-SortSam file-list))
	   (markdup-procs (map run-MarkDup file-list))
	   (samtools-procs (map run-SamTools file-list))
	   (bedtools-procs (map run-BedTools Name-list))
	   (ppqt-procs (map run-PhantomPeakQualTools Name-list)) 
	   (macs2-procs (make-macs2-processes file-list))  
	   (prepFile-procs (map run-prepFiles Name-list))
	  )
  (workflow 
   (name "Pipeline-workflow")
   ;; define all processes
   (processes (append FASTQ-procs bwa-procs sortsam-procs markdup-procs samtools-procs bedtools-procs ppqt-procs macs2-procs prepFile-procs `(,run-computeLD ,run-makeFiles ,run-makeFilesBG ,run-phenoCellSpecif))) 
    ;; Make sure all jobs go in order
	(restrictions (append
					(zip bwa-procs FASTQ-procs)
					(zip sortsam-procs bwa-procs )
					(zip markdup-procs sortsam-procs )
					(zip samtools-procs markdup-procs)
					(map cons bedtools-procs (get-pairs samtools-procs))	
					(zip ppqt-procs bedtools-procs )					
					(zip macs2-procs ppqt-procs)
					(zip prepFile-procs macs2-procs )
					`((,run-makeFiles ,run-computeLD . ,(append prepFile-procs ))) ;;Zip all prepare file processes to one make-file process. First all prepFile processes finish and then make-files can start
					`((,run-makeFilesBG . ,(append prepFile-procs)))
					(list `(,run-phenoCellSpecif ,run-makeFilesBG ,run-makeFiles))
					))
					)))

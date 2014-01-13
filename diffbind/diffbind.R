## Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
# we need that to not crash galaxy with an UTF8 error on German LC settings.
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library('getopt');
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'verbose', 'v', 2, "integer",
	'help' , 'h', 0, "logical",
	'outfile' , 'o', 1, "character",
	'plots' , 'p', 2, "character",
	'infile' , 'i', 1, "character",
    'format', 'f', 1, 'character',
), byrow=TRUE, ncol=4);
opt = getopt(spec);

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
	cat(getopt(spec, usage=TRUE));
	q(status=1);
}


library(DiffBind)
# used to save to BED, GFF or WIG format
library(rtracklayer)

if ( !is.null(opt$plots) ) {
	pdf(opt$plots)
}


sample = dba(sampleSheet=opt$infile)
sample_count = dba.count(sample)
sample_contrast = dba.contrast(sample_count, categories=DBA_CONDITION)
sample_analyze = dba.analyze(sample_contrast)
diff_bind = dba.report(sample_analyze)


export(diff_bind, opt$outfile, format=opt$format)

dev.off()
sessionInfo()

ALL_TESTS=archaea-default archaea-probabilistic haca.snoRNA-default
ALL_TESTS_STDOUT=$(ALL_TESTS:%=%.stdout)

COMMON_DEFAULTS=--tau 50 --plfold-span 150 --plfold-winsize 300 -p	\
0.0005 --max-diff-am 30 --max-diff 60 --max-diff-at-am -1		\
--max-bps-length-ratio 0.0 --noLP --ribofit true

all: $(ALL_TESTS_STDOUT)

## goals for single tests
## convention: 
## -- produce stdout text output in <testname>.stdout
## -- produce result directory <testname>.out
archaea-default.stdout: archaea.fa
	mlocarna ${COMMON_DEFAULTS} $< > $@ --tgtdir $(@:.stdout=.out)

archaea-probabilistic.stdout: archaea.fa
	mlocarna  --probabilistic --consistency-transformation ${COMMON_DEFAULTS} \
	 $< > $@ --tgtdir $(@:.stdout=.out)

haca.snoRNA-default.stdout: haca.snoRNA.fa
	mlocarna ${COMMON_DEFAULTS} $< > $@ --tgtdir $(@:.stdout=.out)

clean:
	$(RM) -r $(ALL_TESTS_STDOUT:.stdout=.out)

realclean: clean
	$(RM) $(ALL_TESTS_STDOUT)

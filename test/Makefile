PY = python

PYTESTS = $(wildcard *.py)
PYNAMES  = $(patsubst %.py,.%,$(PYTESTS))
PYCHECK  = $(patsubst %.py,.%.check,$(PYTESTS))



all: $(PYCHECK)


.%.check: %.py
	$(PY) $<
	#touch $@

%: .%.check
	echo $<


clean:
	rm -f *.dat *.txt *.png .*.check


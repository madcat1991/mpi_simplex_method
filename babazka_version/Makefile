simplex: simplex.c
	mpicc-TV -lm -fopenmp -o simplex simplex.c
sub: submit
submit: simplex
	rm -f *.std*
	rm -f trace.ttf
	rm -f trace
	llsubmit small.job

clean:
	rm -f *.std*
	rm -f trace.ttf
	rm -f trace
	rm simplex

ttf: trace
	tracesort trace

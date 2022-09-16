READS_URL=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR666/004/ERR6668574/ERR6668574.fastq.gz
READS_TEST=reads.fq
READS_LINES=1000000 # 250,000 reads
BTLLIB_VERSION=1.4.7

.PHONY: clean scalability.png

all:
	$(MAKE) input_plot.png
	$(MAKE) cbf_plot.png

btllib:
	wget https://github.com/bcgsc/btllib/releases/download/v$(BTLLIB_VERSION)/btllib-$(BTLLIB_VERSION).tar.gz
	tar zxf btllib-$(BTLLIB_VERSION).tar.gz
	mv btllib-$(BTLLIB_VERSION) btllib
	export AR=gcc-ar
	btllib/compile

$(READS_TEST):
	wget -O- $(READS_URL) | gzip -dc | head -n$(READS_LINES) >$@

test_cbf: test_cbf.cpp btllib
	g++ -std=c++17 -fopenmp -pthread -march=native -O3 -flto $< -o $@ -Ibtllib/install/include -Lbtllib/install/lib -lbtllib

test_seqreader: test_seqreader.cpp btllib
	g++ -std=c++17 -fopenmp -pthread -march=native -O3 -flto $< -o $@ -Ibtllib/install/include -Lbtllib/install/lib -lbtllib

test_kseq: test_kseq.cpp
	g++ -std=c++17 -fopenmp -pthread -march=native -O3 -flto -lz -lm $< -o $@ -I/projects/btl/vnikolic/miniconda/include -L/projects/btl/vnikolic/miniconda/lib

input_data.csv: $(READS_TEST) test_kseq test_seqreader
	./run-test-input $<

input_plot.png: input_data.csv
	./plot-input-data

cbf_data.csv: $(READS_TEST) test_cbf
	./run-test-cbf $<
	./join-cbf-logs test_cbf.log
	./gather-cbf-data test_cbf.log $@

cbf_plot.png: cbf_data.csv
	./plot-cbf-data

scalability.png:
	convert cbf_plot.png input_plot.png +append $@

clean:
	$(RM) *.png *.csv *.log $(READS_TEST) test_kseq test_seqreader test_cbf

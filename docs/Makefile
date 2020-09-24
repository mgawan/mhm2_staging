title="MHM User Guide"
revision="Version 2.0 (09/30/2020)"

TARGET = mhm_guide
INPUT = copy-ack-legal.md mhm_guide.md

pdfprops=-title=$(title)", "$(revision) \
	-pdf:subject=$(subtitle) \
	-createdate="2020:09:30 00:00:00" -modifydate="2020:09:30 00:00:00" \
	-author="S. Hofmeyr, R. Egan" \
	-keywords="Metagenome Assembly, UPC++, Terabase-scale Coassembly, Exascale Computing"

exiftool:=$(shell command -v exiftool 2> /dev/null)
pdflatex:=$(shell command -v pdflatex 2> /dev/null)

#pandocopts= --toc --number-sections --highlight-style=haddock
pandocopts= --highlight-style=haddock

all: $(TARGET).pdf $(TARGET).html

    # embed-image replaces all images with inline base64 encodings so this is truly standalone
$(TARGET).html: $(INPUT)
#	printf "<center><h1>%s<br>%s</h1><h3>%s</h3></center>\n<h2>Contents</h2>\n" $(title) $(revision) $(subtitle) > header.html
	printf "<center><h1>%s</h1><h3>%s</h3></center>\n" $(title) $(revision) >> header.html
	pandoc --metadata pagetitle=$(title)", "$(revision) $(pandocopts) \
               pandoc.css -H header.html \
               -o $@ $^
	./embed-image.sh mhm_guide.html

$(TARGET).pdf: $(INPUT)
ifndef pdflatex
	$(warning WARNING: pdflatex is not available - cannot generate pdf)
else
	perl -pe 's/\*Figure [0-9]+\. ?/*/' -- $(TARGET).md > $(TARGET)_pdf.md
	pandoc -V title=$(title)'\\[15pt]\large{}'$(revision) \
               -V fontsize=10pt -V geometry:margin=1in -V linkcolor=blue \
	       $(pandocopts) \
               -o $@ $(^:$(TARGET).md=$(TARGET)_pdf.md)
	rm -f $(TARGET)_pdf.md
ifndef exiftool
	$(warning WARNING: exiftool is not available - cannot add pdf metadata)
else
	exiftool $(pdfprops) mhm_guide.pdf
endif
endif

clean:
	rm -f $(TARGET).html $(TARGET).pdf *_original header.html *.base64

.PHONY: clean force

#    Continuous Wavelet Transform Library
#    Copyright (C) 2004-2009 Stepan V. Karpenko <carp@mail.ru>
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the
#    Free Software Foundation, Inc., 59 Temple Place, Suite 330,
#    Boston, MA  02111-1307  USA

all: cwtlib_demo.jar


# Intermediate output folder
OUTPUT_DIR := output
# Documentation output folder
DOC_OUTPUT_DIR := docs

# cwtlib .class files
CWTLIB_CLASS_FILES := \
	$(OUTPUT_DIR)/cwtlib/RangeFunctor/RangeFunctor.class \
	$(OUTPUT_DIR)/cwtlib/RangeFunctor/LinearRangeFunctor.class \
	$(OUTPUT_DIR)/cwtlib/Wavelet/Wavelet.class \
	$(OUTPUT_DIR)/cwtlib/Wavelet/MexicanHat.class \
	$(OUTPUT_DIR)/cwtlib/Wavelet/ComplexMorlet.class \
	$(OUTPUT_DIR)/cwtlib/Signal.class \
	$(OUTPUT_DIR)/cwtlib/WTransform.class \
	$(OUTPUT_DIR)/cwtlib/CWTalgorithm.class \
	$(OUTPUT_DIR)/cwtlib_demo/Main.class

# cwtlib .java files
CWTLIB_JAVA_FILES := \
	cwtlib/RangeFunctor/RangeFunctor.java \
	cwtlib/RangeFunctor/LinearRangeFunctor.java \
	cwtlib/Wavelet/Wavelet.java \
	cwtlib/Wavelet/MexicanHat.java \
	cwtlib/Wavelet/ComplexMorlet.java \
	cwtlib/Signal.java \
	cwtlib/WTransform.java \
	cwtlib/CWTalgorithm.java \
	cwtlib_demo/Main.java


# Build demo program binary jar file
cwtlib_demo.jar:  $(OUTPUT_DIR) $(CWTLIB_CLASS_FILES)
	jar -cfe cwtlib_demo.jar cwtlib_demo.Main -C output .


# cwtlib compilation rules
$(OUTPUT_DIR)/cwtlib/RangeFunctor/%.class: cwtlib/RangeFunctor/%.java
	javac -cp output -d output $<

$(OUTPUT_DIR)/cwtlib/Wavelet/%.class: cwtlib/Wavelet/%.java
	javac -cp output -d output $<

$(OUTPUT_DIR)/cwtlib/%.class: cwtlib/%.java
	javac -cp output -d output $<

$(OUTPUT_DIR)/cwtlib_demo/%.class: cwtlib_demo/%.java
	javac -cp output -d output $<

# Create intermediate output folder
$(OUTPUT_DIR):
	mkdir -p $(OUTPUT_DIR)


# Run self-tests
test: cwtlib_demo.jar testcase/sine.txt testcase/delta.txt testcase/freqmod.txt
	java -jar cwtlib_demo.jar testcase/sine.txt    wt_sine.txt
	java -jar cwtlib_demo.jar testcase/delta.txt   wt_delta.txt
	java -jar cwtlib_demo.jar testcase/freqmod.txt wt_freqmod.txt


# Build cwtlib documentation
doc: $(DOC_OUTPUT_DIR) $(CWTLIB_JAVA_FILES)
	javadoc -d $(DOC_OUTPUT_DIR) $(CWTLIB_JAVA_FILES)

$(DOC_OUTPUT_DIR):
	mkdir -p $(DOC_OUTPUT_DIR)


# Remove results of compilation and self-tests execution
clean:
	rm -f $(OUTPUT_DIR)/cwtlib/RangeFunctor/*.class
	rm -f $(OUTPUT_DIR)/cwtlib/Wavelet/*.class
	rm -f $(OUTPUT_DIR)/cwtlib/*.class
	rm -f $(OUTPUT_DIR)/cwtlib_demo/*.class
	rm -f cwtlib_demo.jar
	rm -f wt_*.txt
	rmdir $(OUTPUT_DIR)/cwtlib/RangeFunctor
	rmdir $(OUTPUT_DIR)/cwtlib/Wavelet
	rmdir $(OUTPUT_DIR)/cwtlib
	rmdir $(OUTPUT_DIR)/cwtlib_demo
	rmdir $(OUTPUT_DIR)
